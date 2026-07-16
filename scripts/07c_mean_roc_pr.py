import numpy as np
import pandas as pd
from glob import glob
from sklearn.metrics import roc_curve, precision_recall_curve, auc
import matplotlib.pyplot as plt

files = sorted(glob("04_results/tables/seeds/candidates_pu_scored_seed*.tsv"))
assert files, "Não achei candidates_pu_scored_seed*.tsv"

# grid comum
fpr_grid = np.linspace(0, 1, 500)
rec_grid = np.linspace(0, 1, 500)

tprs = []
rocs = []
precs = []
prs  = []
aucs_roc = []
aucs_pr  = []

for f in files:
    df = pd.read_csv(f, sep="\t")
    # foco no TEST
    df = df[df["split"] == "test"].copy()
    y = df["pu_label"].astype(int).values
    s = df["pu_prob"].astype(float).values

    fpr, tpr, _ = roc_curve(y, s)
    prec, rec, _ = precision_recall_curve(y, s)

    # interp em grid
    tpr_i = np.interp(fpr_grid, fpr, tpr)
    # PR: precision como função do recall (rec é crescente)
    # precision_recall_curve retorna recall crescente, então interp funciona
    prec_i = np.interp(rec_grid, rec[::-1], prec[::-1])  # garante monotonicidade

    tprs.append(tpr_i)
    precs.append(prec_i)

    aucs_roc.append(auc(fpr, tpr))
    aucs_pr.append(auc(rec, prec))

tprs = np.vstack(tprs)
precs = np.vstack(precs)

tpr_mean = tprs.mean(axis=0)
tpr_std  = tprs.std(axis=0, ddof=1)

prec_mean = precs.mean(axis=0)
prec_std  = precs.std(axis=0, ddof=1)

print(f"ROC AUC mean±sd: {np.mean(aucs_roc):.4f} ± {np.std(aucs_roc, ddof=1):.4f}")
print(f"PR  AUC mean±sd: {np.mean(aucs_pr ):.4f} ± {np.std(aucs_pr , ddof=1):.4f}")

# --- plot ROC
plt.figure()
plt.plot(fpr_grid, tpr_mean)
plt.fill_between(fpr_grid, np.clip(tpr_mean - tpr_std, 0, 1), np.clip(tpr_mean + tpr_std, 0, 1), alpha=0.2)
plt.plot([0,1],[0,1], linestyle="--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Mean ROC (test) across seeds")
plt.tight_layout()
plt.savefig("04_results/figures/roc_mean_across_seeds.png", dpi=200)
plt.savefig("04_results/figures/roc_mean_across_seeds.pdf")

# --- plot PR
plt.figure()
plt.plot(rec_grid, prec_mean)
plt.fill_between(rec_grid, np.clip(prec_mean - prec_std, 0, 1), np.clip(prec_mean + prec_std, 0, 1), alpha=0.2)
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("Mean PR (test) across seeds")
plt.tight_layout()
plt.savefig("04_results/figures/pr_mean_across_seeds.png", dpi=200)
plt.savefig("04_results/figures/pr_mean_across_seeds.pdf")

# save curves
out = pd.DataFrame({
    "fpr": fpr_grid, "tpr_mean": tpr_mean, "tpr_sd": tpr_std
})
out.to_csv("04_results/figures/roc_mean_across_seeds.tsv", sep="\t", index=False)

out2 = pd.DataFrame({
    "recall": rec_grid, "precision_mean": prec_mean, "precision_sd": prec_std
})
out2.to_csv("04_results/figures/pr_mean_across_seeds.tsv", sep="\t", index=False)

print("[OK] Saved mean curves + figures in 04_results/figures/")
