#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import yaml


def load_cfg(p):
    with open(p, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def norm_id(s):
    return (
        s.astype(str)
        .str.replace("\r", "", regex=False)
        .str.replace("\t", "", regex=False)
        .str.strip()
        .str.lower()
    )


def ensure_unique(df, name):
    if "candidate_id" not in df.columns:
        raise SystemExit(f"[ERR] {name} missing candidate_id")

    df = df.copy()
    df["candidate_id"] = norm_id(df["candidate_id"])

    if df["candidate_id"].is_unique:
        return df

    print(f"[WARN] {name} duplicated candidate_id — collapsing deterministically")

    agg = {}
    for c in df.columns:
        if c == "candidate_id":
            continue
        if pd.api.types.is_numeric_dtype(df[c]):
            agg[c] = "max"
        else:
            agg[c] = "first"

    return df.groupby("candidate_id", as_index=False).agg(agg)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--top", type=int, default=50)
    args = ap.parse_args()

    cfg = load_cfg(args.config)

    results_dir = Path(cfg["paths"]["results_dir"])
    tables_dir = results_dir / "tables"

    features_path = tables_dir / "all_candidates_features.tsv"
    if not features_path.exists():
        raise SystemExit("[ERR] Run 04_merge_final.py first.")

    ensemble_path = tables_dir / "candidates_pu_scored_ensemble.tsv"
    if not ensemble_path.exists():
        raise SystemExit(
            "[ERR] Ensemble file not found.\n"
            "Run: python scripts/07b_pu_ensemble.py"
        )

    print("[OK] Using ENSEMBLE probabilities only")

    feats = pd.read_csv(features_path, sep="\t")
    feats = ensure_unique(feats, "features")

    pu = pd.read_csv(ensemble_path, sep="\t")
    pu = ensure_unique(pu, "ensemble")

    required_cols = ["candidate_id", "pu_prob_ens_mean"]
    for c in required_cols:
        if c not in pu.columns:
            raise SystemExit(f"[ERR] Missing required ensemble column: {c}")

    pu_keep = pu[["candidate_id", "pu_prob_ens_mean", "pu_prob_ens_std"]].copy()

    df = feats.merge(pu_keep, on="candidate_id", how="left", validate="one_to_one")

    df["pu_prob"] = df["pu_prob_ens_mean"]
    df["pu_prob_sd"] = df["pu_prob_ens_std"]

    df = df.sort_values("pu_prob", ascending=False)

    out_all = tables_dir / "all_candidates_final.tsv"
    out_top = tables_dir / f"top{args.top}_candidates_final.tsv"

    df.to_csv(out_all, sep="\t", index=False)
    df.head(args.top).to_csv(out_top, sep="\t", index=False)

    print("[OK] saved:")
    print(" ", out_all)
    print(" ", out_top)
    print("[INFO] Ensemble ranking applied.")


if __name__ == "__main__":
    main()
