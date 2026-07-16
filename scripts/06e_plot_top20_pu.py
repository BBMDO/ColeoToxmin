#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("04_results/tables/candidates_pu_scored.tsv", sep="\t").sort_values("pu_prob", ascending=False).head(20)
plt.figure(figsize=(7,4))
plt.barh(df["candidate_id"][::-1], df["pu_prob"][::-1])
plt.xlabel("PU probability")
plt.ylabel("Candidate")
plt.tight_layout()
plt.savefig("04_results/figures/top20_pu.png", dpi=300)
plt.savefig("04_results/figures/top20_pu.pdf")
print("[OK] 04_results/figures/top20_pu.(png/pdf)")
