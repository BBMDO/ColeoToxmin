#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt


def main():
    tab = Path("04_results/tables/top20_stability_across_seeds.tsv")
    ens = Path("04_results/tables/candidates_pu_scored_ensemble.tsv")
    outdir = Path("04_results/figures")
    outdir.mkdir(parents=True, exist_ok=True)

    if not tab.exists():
        raise SystemExit(f"[ERR] Missing: {tab}. Run scripts/07a_topk_stability.py first.")

    df = pd.read_csv(tab, sep="\t")

    # opcional: juntar com ensemble mean para desempate (fica mais bonito)
    if ens.exists():
        e = pd.read_csv(ens, sep="\t")
        keep = ["candidate_id", "pu_prob_ens_mean", "pu_prob_ens_std"]
        e = e[[c for c in keep if c in e.columns]].copy()
        df = df.merge(e, on="candidate_id", how="left")

    # ------------- Figura 1: Top estáveis (barplot) -------------
    # mostre, por padrão, candidatos que aparecem em >=2 seeds
    df_plot = df[df["n_seeds_in_top"] >= 2].copy()

    # ordena: mais estáveis primeiro; se tiver ensemble, usa como desempate
    if "pu_prob_ens_mean" in df_plot.columns:
        df_plot = df_plot.sort_values(["n_seeds_in_top", "pu_prob_ens_mean", "candidate_id"],
                                      ascending=[False, False, True])
    else:
        df_plot = df_plot.sort_values(["n_seeds_in_top", "candidate_id"],
                                      ascending=[False, True])

    # limite para não ficar enorme
    max_bars = 25
    df_plot = df_plot.head(max_bars)

    # plot
    plt.figure(figsize=(10, max(4, 0.35 * len(df_plot))))
    y = range(len(df_plot))[::-1]  # inverte para ficar do mais alto em cima
    plt.barh(list(y), df_plot["n_seeds_in_top"].values)

    labels = df_plot["candidate_id"].astype(str).tolist()
    plt.yticks(list(y), labels)
    plt.xlabel("Number of seeds where candidate appears in TOP20")
    plt.title("TOP20 stability across seeds (most stable candidates)")
    plt.xlim(0, 5)
    plt.tight_layout()

    plt.savefig(outdir / "top20_stability_bar.png", dpi=200)
    plt.savefig(outdir / "top20_stability_bar.pdf")
    plt.close()

    # ------------- Figura 2: Distribuição de estabilidade -------------
    counts = df["n_seeds_in_top"].value_counts().sort_index()

    plt.figure(figsize=(7, 4))
    plt.bar(counts.index.astype(int), counts.values)
    plt.xlabel("n_seeds_in_top (out of 5)")
    plt.ylabel("Number of candidates")
    plt.title("Distribution of TOP20 stability across seeds")
    plt.xticks([1, 2, 3, 4, 5])
    plt.tight_layout()

    plt.savefig(outdir / "top20_stability_distribution.png", dpi=200)
    plt.savefig(outdir / "top20_stability_distribution.pdf")
    plt.close()

    # salvar um TSV “do plot” (o top mostrado)
    df_plot.to_csv(outdir / "top20_stability_bar_data.tsv", sep="\t", index=False)

    print("[OK] Saved:")
    print(" - 04_results/figures/top20_stability_bar.(png/pdf)")
    print(" - 04_results/figures/top20_stability_distribution.(png/pdf)")
    print(" - 04_results/figures/top20_stability_bar_data.tsv")


if __name__ == "__main__":
    main()
