#!/usr/bin/env python3
"""
06b_train_pu_fixed.py (v3)

Fixes:
- seed label creation robust when PFAM/struct columns are missing
- grouped split to reduce leakage by near-duplicates
- excludes seed-derived columns from training features
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression


def load_cfg(p: str) -> dict:
    with open(p, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _pick_first_existing(paths):
    for p in paths:
        if p.exists():
            return p
    return None


_rx_suffix = re.compile(r"^(.*)_\d+$")


def group_id_from_candidate_id(cid: str) -> str:
    m = _rx_suffix.match(str(cid))
    return m.group(1) if m else str(cid)


def col_series(df: pd.DataFrame, name: str, default=0):
    """Return a Series aligned to df.index even if column is missing."""
    if name in df.columns:
        return df[name]
    return pd.Series(default, index=df.index)


def make_seed_label(df: pd.DataFrame, cfg: dict) -> pd.Series:
    """
    Creates a PU seed label (1 = positive seed, 0 = unlabeled).
    Robust to missing columns.
    """
    pu = (cfg or {}).get("pu", {})
    rule = pu.get("seed_rules", "pfam_or_struct")

    y = pd.Series(0, index=df.index, dtype=int)

    if rule == "custom":
        q = pu.get("custom_positive_query", "").strip()
        if not q:
            raise SystemExit("[ERR] pu.seed_rules=custom but no query provided")
        pos = df.query(q, engine="python").index
        y.loc[pos] = 1
        return y

    if rule == "secreted_cysrich":
        is_sec = pd.to_numeric(col_series(df, "is_secreted", 0), errors="coerce").fillna(0).astype(int)
        cys = pd.to_numeric(col_series(df, "cys_density", 0.0), errors="coerce").fillna(0.0)
        length = pd.to_numeric(col_series(df, "len_aa", 0), errors="coerce").fillna(0).astype(int)
        y[(is_sec == 1) & (cys >= 0.08) & (length.between(30, 140))] = 1
        return y

    # Default: PFAM/struct-driven seed (works even if cols absent -> no positives from that rule)
    tox = pd.to_numeric(col_series(df, "pfam_is_toxinlike", 0), errors="coerce").fillna(0).astype(int)
    ev = pd.to_numeric(col_series(df, "best_pfam_dom_evalue", np.nan), errors="coerce")
    has_good_pfam = (ev.notna()) & (ev <= 1e-10)

    plddt = pd.to_numeric(col_series(df, "plddt_mean", np.nan), errors="coerce")
    has_good_struct = (plddt.notna()) & (plddt >= 80)

    y[(tox == 1) | (has_good_pfam) | (has_good_struct)] = 1
    return y


def leaky_columns_used_for_seed(cfg: dict):
    base = [
        "pu_label", "pu_prob", "split", "seed_positive",
        "pfam_is_toxinlike", "best_pfam_model", "best_pfam_dom_evalue",
        "pfam_mode", "pfam_mode_frac", "n_pfam_hits",
        "plddt_mean", "sec_helix_frac", "sec_sheet_frac", "rg", "n_disulfide_pairs",
        "priority_score_simple",
    ]
    pu = (cfg or {}).get("pu", {})
    rule = pu.get("seed_rules", "pfam_or_struct")
    if rule == "secreted_cysrich":
        base += ["is_secreted", "cys_density", "len_aa"]
    return sorted(set(base))


def build_feature_matrix(df: pd.DataFrame, exclude):
    exclude = set(exclude or [])
    exclude |= {"candidate_id", "species", "tox_fasta_path", "group_id"}

    num_cols = [c for c in df.columns if c not in exclude and pd.api.types.is_numeric_dtype(df[c])]

    coerce_candidates = [
        "len_aa", "n_cys", "cys_density", "net_charge_pH7", "gravy",
        "aromaticity", "instability_index", "signalp_call", "signalp_score",
        "is_secreted",
    ]
    for c in coerce_candidates:
        if c in df.columns and c not in exclude and c not in num_cols:
            df[c] = pd.to_numeric(df[c], errors="coerce")
            if pd.api.types.is_numeric_dtype(df[c]):
                num_cols.append(c)

    num_cols = sorted(set(num_cols))
    X = df[num_cols].copy()
    return X, num_cols


def grouped_train_val_test_split(df, groups,
                                 test_size=0.15, val_size=0.15,
                                 min_test=25, min_val=25, random_state=0):
    n = len(df)
    test_size = max(test_size, min_test / max(n, 1))
    val_size = max(val_size, min_val / max(n, 1))
    test_size = min(test_size, 0.3)
    val_size = min(val_size, 0.3)

    gss1 = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=random_state)
    train_idx, test_idx = next(gss1.split(df, groups=groups))

    df_train = df.iloc[train_idx].copy()
    df_test = df.iloc[test_idx].copy()

    groups_train = groups.iloc[train_idx]
    gss2 = GroupShuffleSplit(n_splits=1, test_size=val_size, random_state=random_state + 1337)
    tr_idx2, val_idx2 = next(gss2.split(df_train, groups=groups_train))

    df_tr = df_train.iloc[tr_idx2].copy()
    df_val = df_train.iloc[val_idx2].copy()
    return df_tr, df_val, df_test


def proxy_metrics(y_true_seed, y_prob):
    out = {}
    if len(np.unique(y_true_seed)) >= 2:
        out["auprc_proxy"] = float(average_precision_score(y_true_seed, y_prob))
        try:
            out["auroc_proxy"] = float(roc_auc_score(y_true_seed, y_prob))
        except Exception:
            out["auroc_proxy"] = np.nan
    else:
        out["auprc_proxy"] = np.nan
        out["auroc_proxy"] = np.nan
    out["n_pos_seed"] = int(np.sum(y_true_seed == 1))
    out["n_total"] = int(len(y_true_seed))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--test_size", type=float, default=0.15)
    ap.add_argument("--val_size", type=float, default=0.15)
    args = ap.parse_args()

    cfg = load_cfg(args.config)
    outd = Path(cfg["paths"]["results_dir"])
    tables_dir = outd / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    in_path = _pick_first_existing([
        tables_dir / "all_candidates_features.tsv",
        tables_dir / "all_candidates_final.tsv",  # fallback (older)
        outd / "ColeoTox_min_catalog.tsv",
    ])
    if in_path is None:
        raise SystemExit("[ERR] Could not find input table")

    df = pd.read_csv(in_path, sep="\t")
    if "candidate_id" not in df.columns:
        raise SystemExit("[ERR] input table missing candidate_id")

    # Ensure 1 row per candidate_id
    if not df["candidate_id"].is_unique:
        df = df.sort_values("candidate_id").drop_duplicates("candidate_id", keep="first")

    df["pu_label"] = make_seed_label(df, cfg).astype(int)
    df["group_id"] = df["candidate_id"].apply(group_id_from_candidate_id)

    tr, va, te = grouped_train_val_test_split(
        df,
        groups=df["group_id"],
        test_size=args.test_size,
        val_size=args.val_size,
        min_test=25,
        min_val=25,
        random_state=args.seed
    )

    tr["split"] = "train"
    va["split"] = "val"
    te["split"] = "test"
    df_split = pd.concat([tr, va, te], axis=0).reset_index(drop=True)

    train_pos = int(df_split[df_split["split"] == "train"]["pu_label"].sum())
    val_pos = int(df_split[df_split["split"] == "val"]["pu_label"].sum())
    test_pos = int(df_split[df_split["split"] == "test"]["pu_label"].sum())

    print("[OK] Splits: train={}, val={}, test={}".format(
        sum(df_split["split"] == "train"),
        sum(df_split["split"] == "val"),
        sum(df_split["split"] == "test")
    ))
    print("[OK] Seed positives: train={}, val={}, test={}".format(train_pos, val_pos, test_pos))

    if train_pos < 5:
        raise SystemExit("[ERR] Too few seed positives in TRAIN. "
                         "Either your seed rule needs PFAM/struct cols present, "
                         "or set pu.seed_rules=secreted_cysrich/custom in config.yaml.")

    leaky = leaky_columns_used_for_seed(cfg)
    X, feat_cols = build_feature_matrix(df_split, exclude=leaky)

    if X.shape[1] == 0:
        raise SystemExit("[ERR] No numeric features left after leakage exclusion")

    clf = Pipeline(steps=[
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("lr", LogisticRegression(
            solver="lbfgs",
            max_iter=5000,
            C=0.5,
            class_weight="balanced",
        ))
    ])

    train_mask = (df_split["split"] == "train").values
    X_train = X.loc[train_mask].values
    y_train = df_split.loc[train_mask, "pu_label"].astype(int).values

    clf.fit(X_train, y_train)

    df_split["pu_prob"] = clf.predict_proba(X.values)[:, 1]

    metrics = []
    for split in ["train", "val", "test"]:
        mask = (df_split["split"] == split).values
        y = df_split.loc[mask, "pu_label"].values
        p = df_split.loc[mask, "pu_prob"].values
        m = {"split": split}
        m.update(proxy_metrics(y, p))
        metrics.append(m)

    met_df = pd.DataFrame(metrics)

    out_scored = tables_dir / "candidates_pu_scored.tsv"
    out_metrics = tables_dir / "pu_metrics.tsv"

    df_split.to_csv(out_scored, sep="\t", index=False)
    met_df.to_csv(out_metrics, sep="\t", index=False)

    print("[OK] Saved:", out_scored)
    print("[OK] Saved:", out_metrics)


if __name__ == "__main__":
    main()
