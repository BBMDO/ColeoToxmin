#!/usr/bin/env python3
"""06c_make_final_table.py

Builds the *final* publishable tables by merging:
1) Features table produced by 04_merge_final.py (canonical features per candidate)
2) PU scores produced by 06b_train_pu.py

Important:
- Does NOT overwrite the canonical features file.
- Keeps structural columns as NaN when missing (no silent fill->0), except pu_prob.
"""

import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from utils_config import load_config

# Columns to keep in final table (if present).
# Feel free to add/remove here for the paper.
KEEP_COLS = [
    # identifiers / metadata
    "candidate_id","species",

    # PU outputs
    "pu_prob","pu_label","split",

    # sequence features
    "len_aa","n_cys","cys_density","net_charge_pH7","gravy",
    "is_secreted","is_cys_rich",

    # structural summary (keep NaNs if missing)
    "plddt_source","plddt_mean","plddt_std","plddt_frac_ge_70","plddt_frac_ge_80","plddt_frac_ge_90",
    "struct_pass",
    "dssp_n_res","sec_helix_frac","sec_sheet_frac","sec_coil_frac",
    "sasa_total","rsa_mean","rsa_std","core_fraction_rsa_le_thr","hydrophobic_exposed_frac",
    "n_ca","rg","contact_density","n_contact_pairs",
    "n_cys_sg","n_disulfide_pairs","disulfide_pairs",

    # pfam summary
    "best_pfam_model","best_pfam_dom_evalue","n_pfam_hits","pfam_is_toxinlike",

    # fold cluster summary (if present in features table)
    "cluster_id","n_members","pfam_mode","pfam_mode_frac",
]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--top", type=int, default=50)
    args=ap.parse_args()
    cfg=load_config(args.config)

    outd = Path(cfg["paths"]["results_dir"])
    tables_dir = outd / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    # Canonical features (from 04_merge_final.py)
    features_path = tables_dir / "all_candidates_features.tsv"
    if not features_path.exists():
        # backward compat: allow old name
        legacy = tables_dir / "all_candidates_final.tsv"
        raise SystemExit(
            f"[ERR] Missing canonical features table: {features_path}\n"
            f"(Tip: run 04_merge_final.py first; do not use 06c before that.)\n"
            f"Legacy found? {legacy.exists()}")

    df = pd.read_csv(features_path, sep="\t")

    # PU scores (from 06b_train_pu.py)
    scored = tables_dir / "candidates_pu_scored.tsv"
    if scored.exists():
        sc = pd.read_csv(scored, sep="\t")
        # keep only needed pu cols
        keep = [c for c in ["candidate_id","pu_prob","pu_label","split"] if c in sc.columns]
        sc = sc[keep].copy()
        out = df.merge(sc, on="candidate_id", how="left")
    else:
        print(f"[WARN] Missing {scored}. Writing final tables WITHOUT pu_prob.")
        out = df.copy()

    # numeric coercion for pu_prob only (keep NaNs elsewhere)
    if "pu_prob" in out.columns:
        out["pu_prob"] = pd.to_numeric(out["pu_prob"], errors="coerce").fillna(0.0)
        out = out.sort_values("pu_prob", ascending=False)
    else:
        out = out.sort_values("candidate_id")

    # choose cols
    cols = [c for c in KEEP_COLS if c in out.columns]
    out = out[cols].copy()

    # write outputs (final publishable tables)
    out.to_csv(tables_dir / "all_candidates_final.tsv", sep="\t", index=False)
    out.head(args.top).to_csv(tables_dir / f"top{args.top}_candidates_final.tsv", sep="\t", index=False)

    print("[OK] saved:")
    print(f"  {tables_dir / 'all_candidates_final.tsv'}")
    print(f"  {tables_dir / f'top{args.top}_candidates_final.tsv'}")
    print(f"[INFO] canonical features kept at: {features_path}")

if __name__ == "__main__":
    main()
