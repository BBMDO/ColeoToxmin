import pandas as pd
from glob import glob

K = 20
files = sorted(glob("04_results/tables/seeds/candidates_pu_scored_seed*.tsv"))
assert files, "Não achei candidates_pu_scored_seed*.tsv"

top_sets = {}
top_lists = {}

for f in files:
    seed = f.split("seed")[-1].split(".tsv")[0]
    df = pd.read_csv(f, sep="\t")
    # coluna padrão: pu_prob (se no seu for diferente, troque aqui)
    df = df.sort_values("pu_prob", ascending=False)
    top = df.head(K)["candidate_id"].astype(str).tolist()
    top_lists[seed] = top
    top_sets[seed] = set(top)

# Interseção e união
all_union = set().union(*top_sets.values())
all_inter  = set.intersection(*top_sets.values())

print(f"\nSeeds: {list(top_sets.keys())}")
print(f"TOP{K} union size: {len(all_union)}")
print(f"TOP{K} intersection size (present in all seeds): {len(all_inter)}")

# Frequência de aparição no TOPK
freq = {}
for seed, sset in top_sets.items():
    for cid in sset:
        freq[cid] = freq.get(cid, 0) + 1

freq_df = pd.DataFrame({"candidate_id": list(freq.keys()), "n_seeds_in_top": list(freq.values())})
freq_df = freq_df.sort_values(["n_seeds_in_top","candidate_id"], ascending=[False, True])

print("\nCandidatos mais estáveis (top por frequência):")
print(freq_df.head(30).to_string(index=False))

# Jaccard entre seeds
seeds = list(top_sets.keys())
print("\nJaccard entre TOPs (quanto maior, mais estável):")
for i in range(len(seeds)):
    for j in range(i+1, len(seeds)):
        a, b = top_sets[seeds[i]], top_sets[seeds[j]]
        jac = len(a & b) / len(a | b)
        print(f" seed {seeds[i]} vs {seeds[j]}: {jac:.3f}")

# salvar tabela
out = "04_results/tables/top20_stability_across_seeds.tsv"
freq_df.to_csv(out, sep="\t", index=False)
print(f"\n[OK] Saved: {out}")
