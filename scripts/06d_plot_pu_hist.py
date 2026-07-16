#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("04_results/tables/candidates_pu_scored.tsv", sep="\t")
x = df["pu_prob"].astype(float).values

plt.figure(figsize=(6,4))
plt.hist(x, bins=30)
plt.xlabel("PU probability")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig("04_results/figures/pu_prob_hist.png", dpi=300)
plt.savefig("04_results/figures/pu_prob_hist.pdf")
print("[OK] 04_results/figures/pu_prob_hist.(png/pdf)")
