#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

from utils_config import load_config


def norm_id(x: str) -> str:
    return str(x).strip().lower()


def safe_read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {path}")
    return pd.read_csv(path, sep="\t")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)

    base = Path(".").resolve()
    feat_dir = Path(cfg["paths"]["features_dir"]).resolve()
    tmp_dir = Path(cfg["paths"].get("tmp_dir", str(base / "01_intermediate"))).resolve()

    # -----------------------------
    # Inputs esperados
    # -----------------------------
    seq_pq    = feat_dir / "seq" / "candidates_seq_features.parquet"
    struct_pq = feat_dir / "struct" / "candidates_struct_features.parquet"
    pfam_pq   = feat_dir / "hmm" / "candidates_pfam_summary.parquet"

    # Fold clustering (E self-clustering)
    fold_assign_annot = feat_dir / "fold" / "candidates_fold_cluster_annot.tsv"
    fold_assign_plain = feat_dir / "fold" / "candidates_fold_cluster.tsv"
    fold_cluster_sum  = feat_dir / "fold" / "fold_cluster_summary.tsv"

    missing = []
    for p in [seq_pq, struct_pq, pfam_pq, fold_cluster_sum]:
        if not p.exists():
            missing.append(str(p))
    if missing:
        raise FileNotFoundError("Faltando:\n" + "\n".join(missing))

    # assignments: tenta annot primeiro, senão plain
    if fold_assign_annot.exists():
        assign_path = fold_assign_annot
    elif fold_assign_plain.exists():
        assign_path = fold_assign_plain
    else:
        raise FileNotFoundError(f"Faltando assignments de fold cluster: {fold_assign_annot} ou {fold_assign_plain}")

    # -----------------------------
    # Load
    # -----------------------------
    seq = pd.read_parquet(seq_pq)
    st  = pd.read_parquet(struct_pq)
    pf  = pd.read_parquet(pfam_pq)

    asg = safe_read_tsv(assign_path)
    clu = safe_read_tsv(fold_cluster_sum)

    # -----------------------------
    # Normaliza IDs para merge robusto
    # -----------------------------
    for df, col in [(seq, "candidate_id"), (st, "candidate_id"), (pf, "candidate_id")]:
        if col not in df.columns:
            raise KeyError(f"Coluna '{col}' não existe em {df}")
        df["candidate_id_norm"] = df[col].map(norm_id)

    # assignments: detectar colunas automaticamente
    # Esperado: candidate_id + cluster_id (ou nomes equivalentes)
    cand_col = None
    for c in ["candidate_id", "cand_id", "id", "query", "member"]:
        if c in asg.columns:
            cand_col = c
            break
    if cand_col is None:
        raise KeyError(f"Não achei coluna de candidato em {assign_path}. Colunas: {list(asg.columns)}")

    cluster_col = None
    for c in ["cluster_id", "cluster", "rep", "representative"]:
        if c in asg.columns:
            cluster_col = c
            break
    if cluster_col is None:
        raise KeyError(f"Não achei coluna cluster_id em {assign_path}. Colunas: {list(asg.columns)}")

    asg["candidate_id_norm"] = asg[cand_col].map(norm_id)
    asg["cluster_id"] = asg[cluster_col].astype(str)

    # cluster summary: garantir colunas mínimas
    if "cluster_id" not in clu.columns:
        # às vezes vem como rep/representative
        for c in ["rep", "representative", "cluster"]:
            if c in clu.columns:
                clu = clu.rename(columns={c: "cluster_id"})
                break
    if "cluster_id" not in clu.columns:
        raise KeyError(f"Não achei cluster_id em {fold_cluster_sum}. Colunas: {list(clu.columns)}")

    # colunas opcionais no summary
    if "n_members" not in clu.columns:
        # tenta alternativa
        for c in ["n", "size", "cluster_size", "members"]:
            if c in clu.columns:
                clu = clu.rename(columns={c: "n_members"})
                break
    if "n_members" not in clu.columns:
        clu["n_members"] = np.nan

    if "pfam_mode" not in clu.columns:
        clu["pfam_mode"] = "no_pfam_record"
    if "pfam_mode_frac" not in clu.columns:
        clu["pfam_mode_frac"] = 0.0

    clu["pfam_mode"] = clu["pfam_mode"].astype(str).str.lower()

    # -----------------------------
    # Merge master table (candidatos)
    # -----------------------------
    df = (
        seq.merge(st.drop(columns=["candidate_id"], errors="ignore"), on="candidate_id_norm", how="left")
           .merge(pf.drop(columns=["candidate_id"], errors="ignore"), on="candidate_id_norm", how="left")
           .merge(asg[["candidate_id_norm", "cluster_id"]].drop_duplicates(), on="candidate_id_norm", how="left")
           .merge(clu[["cluster_id","n_members","pfam_mode","pfam_mode_frac"]], on="cluster_id", how="left")
    )

    # restaura candidate_id original (prefer seq)
    if "candidate_id" not in df.columns:
        df["candidate_id"] = df["candidate_id_norm"]

    # -----------------------------
    # Regras PU: cria labels heurísticos (paper-friendly)
    # -----------------------------
    # Positivo heurístico = secretado + rico em Cys + PFAM toxin-ish OU cluster com PFAM toxin-ish
    toxin_like_pfams = set([
        "kunitz_bpti", "kazal_2", "defensin_2", "pacifastin_i", "til",
        "trypsin", "peptidase_s9", "astacin", "cap", "phospholip_a2_2"
    ])

    if "signalp_call" not in df.columns:
        df["signalp_call"] = 0
    if "n_cys" not in df.columns:
        df["n_cys"] = 0
    if "len_aa" not in df.columns:
        df["len_aa"] = 0

    # melhor PFAM individual (do seu summary)
    best_pfam = None
    for c in ["best_pfam_model", "best_pfam", "best_pfam_name"]:
        if c in df.columns:
            best_pfam = c
            break
    if best_pfam is None:
        df["best_pfam_model"] = ""
        best_pfam = "best_pfam_model"

    df[best_pfam] = df[best_pfam].astype(str).str.lower()

    df["pfam_is_toxinlike"] = (
        df[best_pfam].isin(toxin_like_pfams) | df["pfam_mode"].isin(toxin_like_pfams)
    ).astype(int)

    # heurística básica
    df["is_secreted"] = (df["signalp_call"].astype(int) == 1).astype(int)
    df["is_cys_rich"] = ((df["n_cys"].astype(int) >= 6) & (df["len_aa"].astype(int) <= 250)).astype(int)

    # label PU: 1 = positivo conhecido (heurístico), 0 = unlabeled
    df["pu_label"] = ((df["is_secreted"] == 1) & (df["is_cys_rich"] == 1) & (df["pfam_is_toxinlike"] == 1)).astype(int)

    # -----------------------------
    # Splits sem leakage: por cluster_id (fold clusters)
    # -----------------------------
    pu_cfg = cfg.get("pu", {})
    seed = int(pu_cfg.get("random_seed", 13))
    test_frac = float(pu_cfg.get("test_frac", 0.2))
    val_frac  = float(pu_cfg.get("val_frac", 0.1))

    rng = np.random.default_rng(seed)

    # clusters válidos (descarta NaN)
    clusters = sorted([c for c in df["cluster_id"].dropna().unique().tolist()])
    rng.shuffle(clusters)

    n = len(clusters)
    n_test = int(round(n * test_frac))
    n_val  = int(round(n * val_frac))

    test_set = set(clusters[:n_test])
    val_set  = set(clusters[n_test:n_test+n_val])
    train_set= set(clusters[n_test+n_val:])

    def split_of(cid):
        if cid in test_set:
            return "test"
        if cid in val_set:
            return "val"
        return "train"

    df["split"] = df["cluster_id"].map(split_of)
    df["split"] = df["split"].fillna("train")

    # -----------------------------
    # Save outputs
    # -----------------------------
    out_dir = tmp_dir / "pu"
    out_dir.mkdir(parents=True, exist_ok=True)

    out_table = out_dir / "pu_master_table.parquet"
    df.sort_values("candidate_id_norm").to_parquet(out_table, index=False)

    out_tsv = out_dir / "pu_master_table.tsv"
    df.sort_values("candidate_id_norm").to_csv(out_tsv, sep="\t", index=False)

    # split files
    for sp in ["train", "val", "test"]:
        ids = df.loc[df["split"] == sp, "candidate_id"].astype(str).tolist()
        (out_dir / f"{sp}.tsv").write_text("\n".join(ids) + "\n", encoding="utf-8")

    print(f"[OK] PU table: {out_table}")
    print(f"[OK] PU table TSV: {out_tsv}")
    print("[OK] Splits:", out_dir / "train.tsv", out_dir / "val.tsv", out_dir / "test.tsv")
    print("[INFO] n candidates:", len(df))
    print("[INFO] positives (heuristic):", int(df["pu_label"].sum()))
    print("[INFO] clusters:", n, "| train/val/test:", len(train_set), len(val_set), len(test_set))
    print("[INFO] assign file used:", assign_path)


if __name__ == "__main__":
    main()
