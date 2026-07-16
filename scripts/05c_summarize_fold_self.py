#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
from utils_config import load_config

def norm_id(x: str) -> str:
    p = Path(str(x))
    return p.stem

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)

    feat_dir = Path(cfg["paths"]["features_dir"])
    hits_path = feat_dir / "fold" / "foldseek_allvsall.tsv"
    clu_path  = feat_dir / "fold" / "fold_clusters.tsv"

    hits = pd.read_csv(hits_path, sep="\t")
    hits["candidate_id"] = hits["query"].map(norm_id)
    hits["target_id"] = hits["target"].map(norm_id)

    clu = pd.read_csv(clu_path, sep="\t", header=None, names=["rep","member"])
    clu["rep_id"] = clu["rep"].map(norm_id)
    clu["member_id"] = clu["member"].map(norm_id)

    # cluster_id = rep_id
    clu["cluster_id"] = clu["rep_id"]

    # size por cluster
    clusize = clu.groupby("cluster_id", as_index=False).agg(cluster_size=("member_id","nunique"))

    # mapeia candidato -> cluster_id (pode haver 1)
    cand2clu = clu.groupby("member_id", as_index=False).first()[["member_id","cluster_id"]]
    cand2clu = cand2clu.rename(columns={"member_id":"candidate_id"})

    # hit do candidato contra seu representante (rep) para ter um score "intra-cluster"
    # pega best score candidate->rep (ou rep->candidate), o que existir
    reps = set(clu["rep_id"].unique())
    hits_rep = hits[hits["target_id"].isin(reps)].copy()
    hits_rep = hits_rep.sort_values(["candidate_id","evalue","score"], ascending=[True, True, False])
    best_to_rep = hits_rep.groupby("candidate_id", as_index=False).first()
    best_to_rep = best_to_rep.rename(columns={
        "target_id":"best_rep_hit",
        "evalue":"best_rep_evalue",
        "score":"best_rep_score",
        "qcov":"best_rep_qcov",
        "tcov":"best_rep_tcov"
    })[["candidate_id","best_rep_hit","best_rep_evalue","best_rep_score","best_rep_qcov","best_rep_tcov"]]

    summary = cand2clu.merge(clusize, on="cluster_id", how="left").merge(best_to_rep, on="candidate_id", how="left")
    summary["cluster_size"] = summary["cluster_size"].fillna(1).astype(int)

    out = feat_dir / "fold" / "candidates_fold_self_summary.tsv"
    summary.sort_values(["cluster_size","candidate_id"], ascending=[False, True]).to_csv(out, sep="\t", index=False)
    print(f"[OK] saved: {out} (n={len(summary)})")

if __name__ == "__main__":
    main()
