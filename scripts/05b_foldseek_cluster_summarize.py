#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
from collections import Counter

import pandas as pd
import yaml


def load_config(cfg_path: str) -> dict:
    cfg_path = Path(cfg_path).expanduser().resolve()
    with open(cfg_path, "r") as f:
        cfg = yaml.safe_load(f)
    if cfg is None:
        cfg = {}
    if "paths" not in cfg:
        cfg["paths"] = {}
    return cfg


def canon_str(x) -> str:
    """Canonical string for joining IDs (lowercase + strip)."""
    if pd.isna(x):
        return ""
    return str(x).strip().lower()


def mode_with_count(values):
    """Return (mode_value, count). Empty -> ("no_hit", 0)."""
    values = [v for v in values if v not in ("", None) and not pd.isna(v)]
    if not values:
        return "no_hit", 0
    c = Counter(values)
    val, cnt = c.most_common(1)[0]
    return val, cnt


def main():
    ap = argparse.ArgumentParser(
        description="Summarize Foldseek self-clusters and annotate with Pfam (case-insensitive safe)."
    )
    ap.add_argument("--config", required=True, help="YAML config path")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = cfg.get("paths", {})

    # Defaults (match your repo layout)
    fold_assign_path = Path(paths.get("fold_cluster_assignments", "02_features/fold/candidates_fold_cluster.tsv")).expanduser()
    pfam_summary_path = Path(paths.get("pfam_summary", "02_features/hmm/candidates_pfam_summary.parquet")).expanduser()
    out_summary_path = Path(paths.get("fold_cluster_summary", "02_features/fold/fold_cluster_summary.tsv")).expanduser()
    out_annot_path = Path(paths.get("fold_cluster_assignments_annot", "02_features/fold/candidates_fold_cluster_annot.tsv")).expanduser()

    fold_assign_path = fold_assign_path.resolve()
    pfam_summary_path = pfam_summary_path.resolve()
    out_summary_path = out_summary_path.resolve()
    out_annot_path = out_annot_path.resolve()

    if not fold_assign_path.exists():
        raise FileNotFoundError(f"[ERR] Missing fold assignments TSV: {fold_assign_path}")
    if not pfam_summary_path.exists():
        raise FileNotFoundError(f"[ERR] Missing Pfam summary parquet: {pfam_summary_path}")

    out_summary_path.parent.mkdir(parents=True, exist_ok=True)
    out_annot_path.parent.mkdir(parents=True, exist_ok=True)

    # Load fold cluster assignments
    cl = pd.read_csv(fold_assign_path, sep="\t")
    # Expect columns: candidate_id, cluster_id (anything else is preserved)
    if "candidate_id" not in cl.columns or "cluster_id" not in cl.columns:
        raise ValueError(
            f"[ERR] {fold_assign_path} must contain candidate_id and cluster_id columns. Found: {list(cl.columns)}"
        )

    # Canonicalize IDs for join + stable summarization
    cl["candidate_id_raw"] = cl["candidate_id"].astype(str)
    cl["cluster_id_raw"] = cl["cluster_id"].astype(str)
    cl["candidate_id"] = cl["candidate_id"].map(canon_str)
    cl["cluster_id"] = cl["cluster_id"].map(canon_str)

    # Load Pfam summary
    pf = pd.read_parquet(pfam_summary_path)

    if "candidate_id" not in pf.columns:
        raise ValueError(f"[ERR] Pfam summary missing candidate_id column: {pfam_summary_path}")

    # Canonicalize Pfam IDs
    pf["candidate_id_raw"] = pf["candidate_id"].astype(str)
    pf["candidate_id"] = pf["candidate_id"].map(canon_str)

    # Columns we will try to use
    # (matches what your pfam summary showed: n_pfam_hits, best_pfam_model, best_pfam_dom_evalue)
    best_model_col = "best_pfam_model" if "best_pfam_model" in pf.columns else None
    n_hits_col = "n_pfam_hits" if "n_pfam_hits" in pf.columns else None
    best_eval_col = "best_pfam_dom_evalue" if "best_pfam_dom_evalue" in pf.columns else None

    # If pfam summary doesn't have these exact columns, we still merge candidate_id only.
    keep_cols = ["candidate_id"]
    for c in (best_model_col, n_hits_col, best_eval_col):
        if c and c in pf.columns:
            keep_cols.append(c)

    pf_small = pf[keep_cols].copy()

    # Merge (LEFT: clusters) so we keep all clustered candidates
    m = cl.merge(pf_small, on="candidate_id", how="left")

    # Make a safe "pfam_label" per candidate for cluster mode:
    # - if no pfam record (merge NaN) -> no_pfam_record
    # - if pfam record exists but n_pfam_hits == 0 -> no_hit
    # - else use best_pfam_model (or "hit" fallback)
    def pfam_label(row):
        # no record in Pfam summary
        if best_model_col is None and n_hits_col is None:
            # we cannot infer hits reliably; only know record existence by merge:
            # if all merged pfam cols absent, then everything is "no_pfam_record"
            return "no_pfam_record"

        # record existence test: if n_hits_col exists, NaN => missing record
        if n_hits_col is not None:
            if pd.isna(row.get(n_hits_col)):
                return "no_pfam_record"
            try:
                nh = int(row.get(n_hits_col))
            except Exception:
                nh = 0
            if nh <= 0:
                return "no_hit"
            # has hits
            if best_model_col is not None:
                v = row.get(best_model_col)
                return canon_str(v) if canon_str(v) else "hit"
            return "hit"

        # if no n_hits_col but best_model exists:
        if best_model_col is not None:
            v = row.get(best_model_col)
            if pd.isna(v) or canon_str(v) == "":
                # ambiguous: could be no-hit or missing; treat as no_pfam_record only if merge had no record:
                # We check candidate_id presence in pf_small
                return "no_pfam_record"
            return canon_str(v)

        return "no_pfam_record"

    m["pfam_label"] = m.apply(pfam_label, axis=1)

    # Save annotated assignments
    # Restore original-looking ids for readability (optional):
    # We'll keep both raw and canonical
    m.to_csv(out_annot_path, sep="\t", index=False)

    # Build cluster summary
    rows = []
    for cid, g in m.groupby("cluster_id", dropna=False):
        n = len(g)
        # mode computed on pfam_label
        mode_val, mode_cnt = mode_with_count(g["pfam_label"].tolist())
        mode_frac = (mode_cnt / n) if n > 0 else 0.0

        rows.append(
            {
                "cluster_id": cid,
                "n_members": n,
                "pfam_mode": mode_val,
                "pfam_mode_n": mode_cnt,
                "pfam_mode_frac": round(mode_frac, 6),
            }
        )

    summary = pd.DataFrame(rows).sort_values(["n_members", "pfam_mode_frac"], ascending=[False, False])
    summary.to_csv(out_summary_path, sep="\t", index=False)

    print(f"[OK] Annotated assignments: {out_annot_path}")
    print(f"[OK] Cluster summary saved: {out_summary_path}")
    # quick stats
    if len(summary) > 0:
        n_nohit = (summary["pfam_mode"] == "no_hit").sum()
        n_norec = (summary["pfam_mode"] == "no_pfam_record").sum()
        print(f"[INFO] clusters: {len(summary)} | no_hit: {n_nohit} | no_pfam_record: {n_norec}")


if __name__ == "__main__":
    main()
