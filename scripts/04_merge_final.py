#!/usr/bin/env python3
"""
04_merge_final_ready_v2.py

Merge final feature tables into a single, 1-row-per-candidate table.

Fixes included in this version:
1) Robust selection of structural features parquet (prefers candidates_struct_features.parquet).
2) Normalizes candidate_id across all tables (strip whitespace/control chars) BEFORE merging,
   preventing silent merge-miss that turns DSSP into all-zeros later.
3) Avoids "masking" missing structural values: we do NOT globally fill numeric NaNs with 0.
   We only coerce/fill the small set of columns needed to compute priority_score_simple.
4) PFAM hits are collapsed to 1 row/candidate_id before merge (prevents cartesian duplication).
"""

import argparse
import re
from pathlib import Path

import pandas as pd
import yaml


# -----------------------------
# IO helpers
# -----------------------------
def load_cfg(p: str) -> dict:
    with open(p, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def read_parquet_if_exists(p: Path) -> pd.DataFrame:
    return pd.read_parquet(p) if p.exists() else pd.DataFrame()


def _pick_first_existing(paths):
    for p in paths:
        if p is not None and Path(p).exists():
            return Path(p)
    return None


def _norm_id(df: pd.DataFrame, col: str = "candidate_id") -> pd.DataFrame:
    """Normalize candidate_id to maximize merge matches."""
    if df is None or df.empty or col not in df.columns:
        return df
    out = df.copy()
    out[col] = (
        out[col]
        .astype(str)
        .str.replace("\r", "", regex=False)
        .str.replace("\t", "", regex=False)
        .str.strip()
        .str.lower()
    )
    return out


# -----------------------------
# PFAM collapsing
# -----------------------------
def collapse_pfam_hits(pfam_df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse a 1->N PFAM hits table to 1 row per candidate_id.

    Produces numeric summaries and simple booleans such as pfam_is_toxinlike if present.
    """
    df = pfam_df.copy()
    if "candidate_id" not in df.columns:
        return pd.DataFrame(columns=["candidate_id"])

    df = _norm_id(df)

    # Common columns across hmmscan/pfam exports
    # We'll keep "best" hit by e-value if possible and compute simple aggregates.
    # If your PFAM table already has 1 row/candidate, this still works.
    ev_cols = [c for c in df.columns if re.search(r"(evalue|e_value|e\-value)", c, flags=re.I)]
    evalue_col = ev_cols[0] if ev_cols else None

    # Make a stable ordering so "first" is deterministic
    if evalue_col:
        df[evalue_col] = pd.to_numeric(df[evalue_col], errors="coerce")
        df = df.sort_values(["candidate_id", evalue_col], ascending=[True, True])
    else:
        df = df.sort_values(["candidate_id"])

    agg = {"candidate_id": "first"}

    # Numeric columns: take min/max/mean depending on meaning; default to max (signal presence)
    for c in df.columns:
        if c == "candidate_id":
            continue
        if pd.api.types.is_numeric_dtype(df[c]):
            # For evalue-like: min is best. For scores/counts: max keeps signal.
            if evalue_col and c == evalue_col:
                agg[c] = "min"
            else:
                agg[c] = "max"

    # Frequently used boolean flag in this project
    if "pfam_is_toxinlike" in df.columns:
        agg["pfam_is_toxinlike"] = lambda x: int(pd.to_numeric(x, errors="coerce").fillna(0).astype(int).max())

    out = df.groupby("candidate_id", as_index=False).agg(agg)
    return out


# -----------------------------
# Dedup helper (deterministic)
# -----------------------------
def ensure_unique_candidate_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse to 1 row/candidate_id deterministically (keeps max for numeric, first for text)."""
    df = df.copy()
    df = _norm_id(df)
    agg = {}
    for c in df.columns:
        if c == "candidate_id":
            continue
        if pd.api.types.is_numeric_dtype(df[c]):
            agg[c] = "max"
        else:
            agg[c] = "first"
    return df.groupby("candidate_id", as_index=False).agg(agg)


def gravy_fallback(df: pd.DataFrame) -> pd.DataFrame:
    """
    Backward-compat: if gravy_1_20 is missing but gravy exists, map it.
    """
    if "gravy_1_20" not in df.columns and "gravy" in df.columns:
        df = df.copy()
        df["gravy_1_20"] = pd.to_numeric(df["gravy"], errors="coerce")
    return df


# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_cfg(args.config)

    feat = Path(cfg["paths"]["features_dir"])
    outd = Path(cfg["paths"]["results_dir"])
    tables_dir = outd / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    # ---- SEQ (required)
    seq_path = _pick_first_existing([
        feat / "seq" / "seq_features.parquet",
        feat / "seq" / "candidates_seq_features.parquet",
    ])
    if seq_path is None:
        raise SystemExit(f"[ERR] missing seq features parquet under: {feat/'seq'}")

    seq = pd.read_parquet(seq_path)
    if "candidate_id" not in seq.columns:
        raise SystemExit(f"[ERR] seq features missing 'candidate_id': {seq_path}")

    seq = _norm_id(seq)
    print(f"[OK] Loaded seq: {seq_path} (n={len(seq)}) unique={seq['candidate_id'].nunique()}")

    if not seq["candidate_id"].is_unique:
        print("[WARN] seq_features has duplicate candidate_id. Collapsing to 1 row/candidate (deterministic).")
        seq = ensure_unique_candidate_rows(seq)
        print(f"[OK] seq collapsed: rows={len(seq)} unique={seq['candidate_id'].nunique()}")

    # ---- STRUCT (optional) — prefer merged candidates_struct_features (03h output)
    struct_path = _pick_first_existing([
        feat / "struct" / "candidates_struct_features.parquet",
        feat / "struct" / "struct_features.parquet",
    ])
    struct = pd.read_parquet(struct_path) if struct_path else pd.DataFrame()
    struct = _norm_id(struct)

    if not struct.empty:
        if "candidate_id" not in struct.columns:
            raise SystemExit(f"[ERR] struct features missing 'candidate_id': {struct_path}")
        if not struct["candidate_id"].is_unique:
            raise SystemExit("[ERR] struct_features has duplicate candidate_id (should not happen).")
        print(f"[INFO] Structural table used: {struct_path}")
        print(f"[OK] Loaded struct: {struct_path} (n={len(struct)}) unique={struct['candidate_id'].nunique()}")
        # Debug overlap
        overlap = len(set(seq["candidate_id"]) & set(struct["candidate_id"]))
        print(f"[DEBUG] overlap seq vs struct candidate_id: {overlap}")
        if "dssp_n_res" in struct.columns:
            print(f"[DEBUG] struct dssp_n_res>0: {(struct['dssp_n_res'].fillna(0)>0).sum()}/{len(struct)}")
    else:
        print("[INFO] No struct table found (skipping).")

    # ---- PFAM (optional)
    pfam_path = _pick_first_existing([
        feat / "pfam" / "pfam_hits.parquet",
        feat / "pfam" / "pfam_hits.tsv",
        feat / "pfam" / "hmmscan_hits.parquet",
        feat / "pfam" / "hmmscan_hits.tsv",
        feat / "pfam" / "pfam_features.parquet",
        feat / "pfam" / "pfam_best.parquet",
    ])

    pfam_1row = pd.DataFrame(columns=["candidate_id"])
    if pfam_path is not None:
        pfam = pd.read_parquet(pfam_path) if pfam_path.suffix == ".parquet" else pd.read_csv(pfam_path, sep="\t")
        pfam = _norm_id(pfam)
        if "candidate_id" in pfam.columns:
            print(f"[OK] Loaded PFAM: {pfam_path} (rows={len(pfam)} unique={pfam['candidate_id'].nunique()})")
            pfam_1row = collapse_pfam_hits(pfam)
            print(f"[OK] Collapsed PFAM -> 1 row/candidate (rows={len(pfam_1row)} unique={pfam_1row['candidate_id'].nunique()})")
        else:
            print(f"[WARN] PFAM file lacks candidate_id: {pfam_path} (skipping)")
    else:
        print("[INFO] No PFAM table found (skipping).")

    # ---- MERGE
    df = seq.copy()
    if not struct.empty:
        df = df.merge(struct, on="candidate_id", how="left", validate="one_to_one")
        if "plddt_mean" in df.columns:
            print(f"[DEBUG] rows with struct matched (plddt_mean notna): {df['plddt_mean'].notna().sum()}/{len(df)}")
        if "dssp_n_res" in df.columns:
            print(f"[DEBUG] rows with dssp_n_res>0 after merge: {(df['dssp_n_res'].fillna(0)>0).sum()}/{len(df)}")

    if not pfam_1row.empty and pfam_1row.shape[1] > 1:
        df = df.merge(pfam_1row, on="candidate_id", how="left", validate="one_to_one")

    # ---- Minimal default filling:
    # Keep NaNs for structural analysis columns (DSSP/SASA/etc.) to avoid masking missing data.
    # For text columns, fill empty string.
    for c in df.columns:
        if c == "candidate_id":
            continue
        if pd.api.types.is_numeric_dtype(df[c]):
            df[c] = pd.to_numeric(df[c], errors="coerce")
        else:
            df[c] = df[c].fillna("")

    df = gravy_fallback(df)

    # ---- Compute priority score (coerce only needed columns to numeric; fill missing with 0 in the score)
    def col(name, default=0.0):
        if name in df.columns:
            return pd.to_numeric(df[name], errors="coerce").fillna(default)
        return pd.Series(default, index=df.index)

    df["priority_score_simple"] = (
        0.25 * col("gravy_1_20") +
        0.25 * col("cys_density") +
        0.20 * (col("plddt_mean") / 100.0) +
        0.15 * col("n_disulfide_pairs") +
        0.10 * col("is_secreted") +
        0.05 * col("pfam_is_toxinlike")
    )

    # ---- Hard safety: uniqueness
    if not df["candidate_id"].is_unique:
        print("[WARN] Duplicates detected after merges. Collapsing deterministically...")
        df = ensure_unique_candidate_rows(df)

    if not df["candidate_id"].is_unique:
        dup = df.loc[df["candidate_id"].duplicated(keep=False), "candidate_id"].value_counts().head(20)
        raise SystemExit(f"[ERR] Still duplicated after collapse. Top dups:\n{dup}")

    df = df.sort_values("priority_score_simple", ascending=False)

    out_all = tables_dir / "all_candidates_features.tsv"
    out_top50 = tables_dir / "top50_candidates_features.tsv"
    df.to_csv(out_all, sep="\t", index=False)
    df.head(50).to_csv(out_top50, sep="\t", index=False)

    print(f"[OK] Wrote: {out_all} (rows={len(df)} unique={df['candidate_id'].nunique()})")
    print(f"[OK] Wrote: {out_top50} (rows={min(50, len(df))})")


if __name__ == "__main__":
    main()
