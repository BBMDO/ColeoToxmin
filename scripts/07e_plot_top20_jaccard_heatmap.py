#!/usr/bin/env python3
import pandas as pd
from glob import glob
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np


def main():
    K = 20
    files = sorted(glob("04_results/tables/seeds/candidates_pu_scored_seed*.tsv"))
    if not files:
        raise SystemExit("[ERR] Não achei 04_results/tables/seeds/candidates_pu_scored_seed*.tsv")

    # extrai seeds e conjuntos topK
    top_sets = {}
    seeds = []
    for f in files:
        seed = f.split("seed")[-1].split(".tsv")[0]
        seeds.append(seed)
        df = pd.read_csv(f, sep="\t")
        if "pu_prob" not in df.columns:
            raise SystemExit(f"[ERR] Coluna pu_prob não encontrada em {f}")

        df = df.sort_values("pu_prob", ascending=False)
        top = df.head(K)["candidate_id"].astype(str).str.strip().str.lower().tolist()
        top_sets[seed] = set(top)

    # matriz Jaccard
    n = len(seeds)
    J = np.zeros((n, n), dtype=float)
    for i, si in enumerate(seeds):
        for j, sj in enumerate(seeds):
            a, b = top_sets[si], top_sets[sj]
            J[i, j] = len(a & b) / max(1, len(a | b))

    outdir = Path("04_results/figures")
    outdir.mkdir(parents=True, exist_ok=True)

    # salvar TSV
    jdf = pd.DataFrame(J, index=[f"seed{z}" for z in seeds], columns=[f"seed{z}" for z in seeds])
    tsv_path = outdir / "top20_jaccard_heatmap.tsv"
    jdf.to_csv(tsv_path, sep="\t")

    # plot heatmap (sem seaborn, só matplotlib)
    plt.figure(figsize=(6, 5))
    im = plt.imshow(J, vmin=0, vmax=1)
    plt.colorbar(im, label="Jaccard overlap (TOP20)")

    plt.xticks(range(n), [f"seed{z}" for z in seeds], rotation=45, ha="right")
    plt.yticks(range(n), [f"seed{z}" for z in seeds])

    plt.title("TOP20 overlap across seeds (Jaccard index)")
    plt.tight_layout()

    png_path = outdir / "top20_jaccard_heatmap.png"
    pdf_path = outdir / "top20_jaccard_heatmap.pdf"
    plt.savefig(png_path, dpi=200)
    plt.savefig(pdf_path)
    plt.close()

    print("[OK] Saved:")
    print(f" - {png_path}")
    print(f" - {pdf_path}")
    print(f" - {tsv_path}")


if __name__ == "__main__":
    main()
