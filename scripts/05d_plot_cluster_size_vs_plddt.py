#!/usr/bin/env python3
"""
Figura 2 — Cluster size × pLDDT_mean

Lê:
- 02_features/fold/candidates_fold_cluster_annot.tsv   (ou candidates_fold_cluster.tsv)
- 02_features/struct/candidates_struct_features.parquet

Gera:
- 04_results/figures/Fig2_cluster_size_vs_plddt.png
- 04_results/figures/Fig2_cluster_size_vs_plddt.pdf

Uso:
python scripts/05d_plot_cluster_size_vs_plddt.py --config 09_reproducibility/configs/config.yaml
"""
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

from utils_config import load_config


def pick_cluster_file(fold_dir: Path) -> Path:
    # Preferir o arquivo anotado (com pfam_mode etc), mas aceitar o básico
    p1 = fold_dir / "candidates_fold_cluster_annot.tsv"
    p2 = fold_dir / "candidates_fold_cluster.tsv"
    if p1.exists():
        return p1
    if p2.exists():
        return p2
    raise FileNotFoundError(f"Não achei {p1} nem {p2}")


def normalize_id(s: str) -> str:
    # evita o seu problema de interseção 0 (maiúsc/minúsc)
    return str(s).strip().lower()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--topk", type=int, default=0, help="0 = todos os clusters; senão plota só top-K por tamanho")
    args = ap.parse_args()

    cfg = load_config(args.config)

    feat_dir = Path(cfg["paths"]["features_dir"])
    fold_dir = feat_dir / "fold"
    struct_path = feat_dir / "struct" / "candidates_struct_features.parquet"

    cluster_path = pick_cluster_file(fold_dir)
    if not struct_path.exists():
        raise FileNotFoundError(f"Não achei struct features: {struct_path}")

    # --- Load ---
    cl = pd.read_csv(cluster_path, sep="\t")
    st = pd.read_parquet(struct_path)

    # Esperado no cluster file: candidate_id, cluster_id (ou cluster_rep)
    if "cluster_id" not in cl.columns:
        # fallback comum: cluster_rep / rep / representative
        for alt in ["cluster_rep", "rep", "representative"]:
            if alt in cl.columns:
                cl = cl.rename(columns={alt: "cluster_id"})
                break
    if "cluster_id" not in cl.columns:
        raise ValueError(f"{cluster_path} precisa ter coluna cluster_id (ou rep/cluster_rep). Colunas: {list(cl.columns)}")
    if "candidate_id" not in cl.columns:
        raise ValueError(f"{cluster_path} precisa ter coluna candidate_id. Colunas: {list(cl.columns)}")
    if "candidate_id" not in st.columns or "plddt_mean" not in st.columns:
        raise ValueError(f"{struct_path} precisa ter candidate_id e plddt_mean. Colunas: {list(st.columns)[:20]}...")

    # Normaliza IDs (case-insensitive)
    cl["candidate_id_norm"] = cl["candidate_id"].map(normalize_id)
    st["candidate_id_norm"] = st["candidate_id"].map(normalize_id)

    # Merge
    m = cl.merge(st[["candidate_id_norm", "plddt_mean", "struct_pass"]], on="candidate_id_norm", how="left")

    # Sanity
    miss = m["plddt_mean"].isna().mean()
    if miss > 0.05:
        print(f"[WARN] {miss:.1%} dos candidatos no cluster não acharam pLDDT no struct parquet. Verifique IDs/arquivos.")

    # Agrega por cluster
    agg = (
        m.groupby("cluster_id", as_index=False)
         .agg(
             n_members=("candidate_id_norm", "nunique"),
             plddt_mean_cluster=("plddt_mean", "mean"),
             plddt_median_cluster=("plddt_mean", "median"),
             frac_struct_pass=("struct_pass", lambda x: float((x.fillna(0) >= 1).mean())),
         )
    )

    # opcional: topK clusters por tamanho
    if args.topk and args.topk > 0:
        agg = agg.sort_values("n_members", ascending=False).head(args.topk)

    # --- Plot ---
    out_dir = Path("04_results/figures")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "Fig2_cluster_size_vs_plddt.png"
    out_pdf = out_dir / "Fig2_cluster_size_vs_plddt.pdf"

    # tamanho do ponto proporcional ao tamanho do cluster (com clamp)
    sizes = (agg["n_members"].astype(float) ** 1.2) * 6.0
    sizes = sizes.clip(lower=30, upper=1200)

    plt.figure(figsize=(7.2, 5.2))
    plt.scatter(
        agg["n_members"],
        agg["plddt_mean_cluster"],
        s=sizes,
        alpha=0.75,
        edgecolors="black",
        linewidths=0.5,
    )

    # linhas-guia úteis (paper-friendly)
    plt.axhline(70, linestyle="--", linewidth=1.0)
    plt.xlabel("Cluster size (n members)")
    plt.ylabel("Mean pLDDT (cluster)")
    plt.title("Fold self-clusters: size vs. structural confidence (pLDDT)")

    # rótulos só para clusters muito grandes (evita poluição)
    big = agg.sort_values("n_members", ascending=False).head(5)
    for _, r in big.iterrows():
        plt.text(
            r["n_members"] + 0.4,
            r["plddt_mean_cluster"] + 0.2,
            str(r["cluster_id"])[:28],
            fontsize=8,
        )

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()

    # salva tabela para suplementar / debug
    out_tsv = out_dir / "Fig2_cluster_size_vs_plddt_data.tsv"
    agg.sort_values(["n_members", "plddt_mean_cluster"], ascending=[False, False]).to_csv(out_tsv, sep="\t", index=False)

    print(f"[OK] Figura salva: {out_png}")
    print(f"[OK] Figura salva: {out_pdf}")
    print(f"[OK] Dados da figura: {out_tsv}")
    print(f"[INFO] clusters plotados: {len(agg)}")


if __name__ == "__main__":
    main()
