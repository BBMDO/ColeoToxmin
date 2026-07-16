import pandas as pd
from glob import glob

files = sorted(glob("04_results/tables/seeds/candidates_pu_scored_seed*.tsv"))
assert files, "Não achei candidates_pu_scored_seed*.tsv"

dfs = []
for f in files:
    seed = f.split("seed")[-1].split(".tsv")[0]
    d = pd.read_csv(f, sep="\t")[["candidate_id","pu_prob","split","pu_label"]].copy()
    d["candidate_id"] = d["candidate_id"].astype(str)
    d = d.rename(columns={"pu_prob": f"pu_prob_seed{seed}"})
    dfs.append(d)

# merge por candidate_id (split/pu_label devem ser idênticos entre seeds? podem variar se split muda)
# então vamos manter split/pu_label do próprio arquivo "features" ou do primeiro seed, e salvar os demais como cols
base = dfs[0].copy()
base_cols = ["candidate_id","split","pu_label", base.columns[1]]
base = base[base_cols]

for d in dfs[1:]:
    base = base.merge(d[["candidate_id", d.columns[1]]], on="candidate_id", how="inner")

prob_cols = [c for c in base.columns if c.startswith("pu_prob_seed")]
base["pu_prob_ens_mean"] = base[prob_cols].mean(axis=1)
base["pu_prob_ens_std"]  = base[prob_cols].std(axis=1, ddof=1)

# Ranking ensemble
ens = base.sort_values("pu_prob_ens_mean", ascending=False)

out = "04_results/tables/candidates_pu_scored_ensemble.tsv"
ens.to_csv(out, sep="\t", index=False)
print(f"[OK] Saved: {out}")
print("\nTop20 ensemble (candidate_id, mean, std):")
print(ens[["candidate_id","pu_prob_ens_mean","pu_prob_ens_std"]].head(20).to_string(index=False))
