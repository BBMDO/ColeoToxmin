#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
03_clade_enrichment.py

Fecha Lacuna C: padrões por clado + estatística.

Inputs:
- --in_candidates_tsv: TSV principal (ex.: 04_results/tables/all_candidates_features.tsv
  ou ColeoTox_min_catalog.tsv), precisa ter 'species' e idealmente:
    - is_secreted (0/1) ou signalp_call/is_secreted
    - n_cys ou cys_density
    - hmm_toxinlike / pfam_is_toxinlike / is_toxinlike (qualquer um)
    - priority_score_simple (opcional, mas ótimo para comparar distribuições)

- --taxonomy_tsv (opcional mas recomendado): um TSV com colunas:
    species, family, order, suborder
  Isso torna o pipeline 100% reprodutível.

Outputs:
- <out_prefix>.taxonomy_mapped.tsv
- <out_prefix>.family_enrichment.tsv
- <out_prefix>.clade_score_tests.tsv
- <out_prefix>.plots.pdf
- <out_prefix>.summary.txt

Modo taxonomia:
- Se taxonomy_tsv for fornecido: usa ele
- Senão: tenta inferir via ete3 NCBITaxa (se disponível)
"""

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from scipy.stats import fisher_exact, kruskal


def normalize_species(s: str) -> str:
    # aceita "Genus_species" ou "Genus species"
    s = str(s).strip()
    return s.replace("_", " ").strip()


def pick_first_existing(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


def fdr_bh(pvals):
    """Benjamini-Hochberg FDR."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        prev = min(prev, val)
        q[i] = prev
    out = np.empty(n, dtype=float)
    out[order] = q
    return out


def load_taxonomy_table(path):
    tx = pd.read_csv(path, sep="\t")
    needed = {"species", "family", "order", "suborder"}
    missing = needed - set(tx.columns)
    if missing:
        raise ValueError(f"taxonomy_tsv falta colunas: {missing}")
    tx["species_norm"] = tx["species"].apply(normalize_species)
    return tx[["species_norm", "family", "order", "suborder"]].drop_duplicates()


def try_ete3_taxonomy(species_list):
    """
    Tenta mapear via ete3.NCBITaxa.
    Requer ete3 instalado e banco NCBI local (pode precisar download prévio).
    """
    try:
        from ete3 import NCBITaxa
    except Exception as e:
        return None, f"ete3 não disponível: {e}"

    ncbi = NCBITaxa()

    rows = []
    failures = 0

    for sp in species_list:
        spn = normalize_species(sp)
        # tenta resolver taxid por nome
        try:
            name2tax = ncbi.get_name_translator([spn])
            if spn not in name2tax:
                failures += 1
                continue
            taxid = name2tax[spn][0]
            lineage = ncbi.get_lineage(taxid)
            ranks = ncbi.get_rank(lineage)
            names = ncbi.get_taxid_translator(lineage)

            fam = ord_ = subord = "NA"
            for tid in lineage:
                r = ranks.get(tid, "")
                if r == "family":
                    fam = names.get(tid, "NA")
                elif r == "order":
                    ord_ = names.get(tid, "NA")
                elif r == "suborder":
                    subord = names.get(tid, "NA")

            rows.append({"species_norm": spn, "family": fam, "order": ord_, "suborder": subord})
        except Exception:
            failures += 1

    if not rows:
        return None, f"ete3 disponível, mas falhou em mapear (falhas={failures})."

    tx = pd.DataFrame(rows).drop_duplicates()
    msg = f"ete3 ok: mapeou {tx.shape[0]} espécies; falhas={failures}"
    return tx, msg


def define_interest_group(df, cys_rule="n_cys>=6"):
    """
    Define grupo (0/1) de candidatos “toxins-like secretados”.
    Regra padrão:
      is_secreted==1 AND Cys-rich AND toxinlike==1 (por HMM/Pfam/flag)
    """
    # secreted flag
    col_secreted = pick_first_existing(df, ["is_secreted", "signalp_call"])
    if col_secreted is None:
        raise ValueError("Não achei coluna de secreção (is_secreted ou signalp_call).")

    # toxinlike flag
    col_tox = pick_first_existing(df, ["hmm_toxinlike", "pfam_is_toxinlike", "is_toxinlike", "toxinlike", "toxinlike_label"])
    if col_tox is None:
        # se não existir, assume 1 para todos (para ainda permitir enriquecimento “secreted + cys-rich”)
        df["_toxinlike_tmp"] = 1
        col_tox = "_toxinlike_tmp"

    # Cys rule
    # suporta: "n_cys>=6" ou "cys_density>=0.08"
    if "n_cys" in df.columns and "n_cys" in cys_rule:
        pass
    elif "cys_density" in df.columns and "cys_density" in cys_rule:
        pass
    elif "n_cys" in df.columns:
        cys_rule = "n_cys>=6"
    elif "cys_density" in df.columns:
        cys_rule = "cys_density>=0.08"
    else:
        raise ValueError("Não achei n_cys nem cys_density para definir Cys-rich.")

    # coerce numeric
    df[col_secreted] = pd.to_numeric(df[col_secreted], errors="coerce").fillna(0).astype(int)
    df[col_tox] = pd.to_numeric(df[col_tox], errors="coerce").fillna(0).astype(int)

    # cys mask
    try:
        cys_mask = df.eval(cys_rule)
    except Exception as e:
        raise ValueError(f"Falha ao aplicar cys_rule='{cys_rule}': {e}")

    df["interest_group"] = ((df[col_secreted] == 1) & (df[col_tox] == 1) & (cys_mask)).astype(int)
    return df, col_secreted, col_tox, cys_rule


def family_enrichment(df):
    """
    Fisher por família:
      interest_group=1 vs 0 (em cada família) comparado ao resto.
    """
    if "family" not in df.columns:
        raise ValueError("Sem coluna family após mapear taxonomia.")

    fams = [f for f in df["family"].fillna("NA").unique() if f != "NA"]
    rows = []

    total_interest = int(df["interest_group"].sum())
    total_non = int((df["interest_group"] == 0).sum())

    for fam in fams:
        sub = df[df["family"] == fam]
        a = int(sub["interest_group"].sum())                  # fam & interest
        b = int((sub["interest_group"] == 0).sum())           # fam & not
        c = total_interest - a                                # not fam & interest
        d = total_non - b                                     # not fam & not

        # Fisher
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")

        rows.append({
            "family": fam,
            "a_interest_in_family": a,
            "b_non_interest_in_family": b,
            "odds_ratio": odds,
            "p_value": p,
            "family_size": a + b,
            "interest_fraction_in_family": a / max(1, (a + b))
        })

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    out["fdr_bh"] = fdr_bh(out["p_value"].values)
    out = out.sort_values(["fdr_bh", "p_value", "odds_ratio"], ascending=[True, True, False])
    return out


def clade_score_tests(df):
    """
    Compara priority_score_simple entre famílias (ou ordem/subordem se quiser).
    Usa Kruskal-Wallis (não paramétrico).
    """
    score_col = pick_first_existing(df, ["priority_score_simple", "simple_rank_score", "priority_score"])
    if score_col is None:
        return None, "Sem coluna de score (priority_score_simple)."

    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df2 = df.dropna(subset=[score_col]).copy()

    # só famílias com n>=10 para estabilidade
    fam_counts = df2["family"].value_counts()
    keep_fams = fam_counts[fam_counts >= 10].index.tolist()
    df2 = df2[df2["family"].isin(keep_fams)].copy()
    if df2["family"].nunique() < 2:
        return None, "Poucas famílias com n>=10 para teste de score."

    groups = [df2.loc[df2["family"] == f, score_col].values for f in df2["family"].unique()]
    stat, p = kruskal(*groups)

    out = pd.DataFrame({
        "test": ["Kruskal-Wallis"],
        "by": ["family"],
        "score_col": [score_col],
        "n_families": [df2["family"].nunique()],
        "n_total": [df2.shape[0]],
        "statistic": [stat],
        "p_value": [p],
    })
    return out, None


def make_plots(df, enr, out_pdf):
    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(out_pdf) as pdf:
        # 1) Top enriched families (bar odds or fraction)
        if enr is not None and not enr.empty:
            top = enr.head(15).copy()
            fig = plt.figure()
            plt.barh(top["family"][::-1], top["odds_ratio"][::-1])
            plt.xlabel("Odds ratio (Fisher, greater)")
            plt.title("Top enriched families (interest_group)")
            plt.tight_layout()
            pdf.savefig(fig); plt.close(fig)

        # 2) Boxplot score by family (top N families)
        score_col = pick_first_existing(df, ["priority_score_simple", "simple_rank_score", "priority_score"])
        if score_col is not None and "family" in df.columns:
            dfc = df.dropna(subset=[score_col]).copy()
            dfc[score_col] = pd.to_numeric(dfc[score_col], errors="coerce")
            fam_counts = dfc["family"].value_counts()
            fams = fam_counts.head(10).index.tolist()
            dfc = dfc[dfc["family"].isin(fams)].copy()
            if dfc["family"].nunique() >= 2:
                fig = plt.figure()
                data = [dfc.loc[dfc["family"] == f, score_col].values for f in fams]
                plt.boxplot(data, labels=fams, vert=True)
                plt.xticks(rotation=60, ha="right")
                plt.ylabel(score_col)
                plt.title("Score distribution across top families")
                plt.tight_layout()
                pdf.savefig(fig); plt.close(fig)

        # 3) Interest fraction by order (se existir)
        if "order" in df.columns:
            tmp = df.copy()
            tmp["order"] = tmp["order"].fillna("NA")
            ord_tab = tmp.groupby("order")["interest_group"].agg(["sum", "count"]).reset_index()
            ord_tab["fraction"] = ord_tab["sum"] / ord_tab["count"].clip(lower=1)
            ord_tab = ord_tab[ord_tab["order"] != "NA"].sort_values("fraction", ascending=False).head(15)
            if not ord_tab.empty:
                fig = plt.figure()
                plt.barh(ord_tab["order"][::-1], ord_tab["fraction"][::-1])
                plt.xlabel("interest fraction")
                plt.title("Interest-group fraction by order (top 15)")
                plt.tight_layout()
                pdf.savefig(fig); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_candidates_tsv", required=True)
    ap.add_argument("--out_prefix", required=True)
    ap.add_argument("--taxonomy_tsv", default="", help="TSV com species,family,order,suborder (recomendado)")
    ap.add_argument("--cys_rule", default="n_cys>=6", help="Ex.: n_cys>=6 ou cys_density>=0.08")
    args = ap.parse_args()

    df = pd.read_csv(args.in_candidates_tsv, sep="\t")
    if "species" not in df.columns:
        print("[ERROR] TSV precisa ter coluna 'species'.", file=sys.stderr)
        sys.exit(1)

    df["species_norm"] = df["species"].apply(normalize_species)

    # Taxonomy mapping
    tx_msg = ""
    if args.taxonomy_tsv.strip():
        tx = load_taxonomy_table(args.taxonomy_tsv.strip())
        tx_msg = f"taxonomy_tsv usado: {args.taxonomy_tsv.strip()}"
    else:
        unique_species = sorted(df["species_norm"].dropna().unique().tolist())
        tx, tx_msg = try_ete3_taxonomy(unique_species)
        if tx is None:
            print("[ERROR] Sem taxonomy_tsv e não consegui mapear via ete3.", file=sys.stderr)
            print("        Motivo:", tx_msg, file=sys.stderr)
            print("        Solução: crie um TSV species,family,order,suborder e rode com --taxonomy_tsv.", file=sys.stderr)
            sys.exit(1)

    df = df.merge(tx, on="species_norm", how="left")

    # Interest group
    df, col_secreted, col_tox, cys_rule_used = define_interest_group(df, cys_rule=args.cys_rule)

    # Save taxonomy mapped
    out_tax = f"{args.out_prefix}.taxonomy_mapped.tsv"
    df.to_csv(out_tax, sep="\t", index=False)

    # Enrichment by family
    enr = family_enrichment(df)
    out_enr = f"{args.out_prefix}.family_enrichment.tsv"
    enr.to_csv(out_enr, sep="\t", index=False)

    # Score test by clade
    score_test, msg = clade_score_tests(df)
    out_score = f"{args.out_prefix}.clade_score_tests.tsv"
    if score_test is not None:
        score_test.to_csv(out_score, sep="\t", index=False)
    else:
        pd.DataFrame([{"note": msg}]).to_csv(out_score, sep="\t", index=False)

    # Plots
    out_pdf = f"{args.out_prefix}.plots.pdf"
    make_plots(df, enr, out_pdf)

    # Summary
    out_sum = f"{args.out_prefix}.summary.txt"
    with open(out_sum, "w", encoding="utf-8") as f:
        f.write(f"taxonomy_source\t{tx_msg}\n")
        f.write(f"n_rows\t{df.shape[0]}\n")
        f.write(f"n_species\t{df['species_norm'].nunique()}\n")
        f.write(f"secreted_col\t{col_secreted}\n")
        f.write(f"toxinlike_col\t{col_tox}\n")
        f.write(f"cys_rule\t{cys_rule_used}\n")
        f.write(f"interest_count\t{int(df['interest_group'].sum())}\n")
        f.write(f"interest_fraction\t{df['interest_group'].mean():.4f}\n")
        if enr is not None and not enr.empty:
            top = enr.head(5)
            f.write("top5_families\t" + ";".join([f"{r.family}(OR={r.odds_ratio:.2g},FDR={r.fdr_bh:.2g})" for _, r in top.iterrows()]) + "\n")

    print("[OK] Wrote:")
    print(" -", out_tax)
    print(" -", out_enr)
    print(" -", out_score)
    print(" -", out_pdf)
    print(" -", out_sum)


if __name__ == "__main__":
    main()
