#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Caminho do arquivo gerado na Etapa E
file_path = Path("02_features/fold/fold_cluster_summary.tsv")

if not file_path.exists():
    raise FileNotFoundError(
        f"{file_path} not found. Run Etapa E first."
    )

# Load data
df = pd.read_csv(file_path, sep="\t")

# Select top 15 clusters by size
top15 = df.sort_values("n_members", ascending=False).head(15)

# Plot (single plot, no specific colors)
plt.figure()
plt.bar(top15["cluster_id"], top15["n_members"])
plt.xticks(rotation=90)
plt.xlabel("Cluster ID")
plt.ylabel("Number of Members")
plt.title("Top 15 Structural Fold Clusters (Self-Clustering)")
plt.tight_layout()

# Save
out_path = Path("02_features/fold/fold_cluster_top15_barplot.png")
plt.savefig(out_path, dpi=300)

print(f"[OK] Saved figure: {out_path}")
