#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
from utils_config import load_config

def rp(p): 
    return pd.read_parquet(p) if p.exists() else pd.DataFrame()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)

    out_dir = Path(cfg["paths"]["features_dir"]) / "struct"
    out_dir.mkdir(parents=True, exist_ok=True)

    plddt = rp(out_dir / "plddt_features.parquet")
    dssp  = rp(out_dir / "dssp_features.parquet")
    sasa  = rp(out_dir / "sasa_rsa_features.parquet")
    rgc   = rp(out_dir / "contacts_rg_features.parquet")
    ss    = rp(out_dir / "disulfide_features.parquet")

    dfs = [plddt, dssp, sasa, rgc, ss]
    df = None
    for d in dfs:
        if d.empty:
            continue
        df = d if df is None else df.merge(d, on="candidate_id", how="outer")

    if df is None:
        raise SystemExit("[ERR] no structural tables to merge")

    pmin = float(cfg["gating"]["plddt_min_mean"])
    df["struct_pass"] = (df["plddt_mean"].fillna(0) >= pmin).astype(int)

    # sanity checks
    if "sec_helix_frac" in df.columns and "sec_sheet_frac" in df.columns:
        frac_bad = ((df["sec_helix_frac"].fillna(0) + df["sec_sheet_frac"].fillna(0)) < 0.01).mean()
        print(f"[CHECK] (helix+sheet)<0.01 fraction: {frac_bad:.2f}")
        #if frac_bad > 0.8:
            #print("[WARN] DSSP suspicious. Check mkdssp install / input PDB validity.")
        if "dssp_n_res" in df.columns:
            bad = (df["dssp_n_res"].fillna(0) == 0).mean()
            if bad > 0.8:
                raise SystemExit("[ERR] DSSP praticamente ausente (dssp_n_res==0 em >80%). "
                         "Verifique mkdssp e o input (.cif/.pdb).")

    nan_rg = df["rg"].isna().mean() if "rg" in df.columns else 1.0
    print(f"[CHECK] rg NaN fraction: {nan_rg:.2f}")

    out1 = out_dir / "candidates_struct_features.parquet"
    out2 = out_dir / "struct_features.parquet"  # <- sobrescreve o “antigo” e vira o padrão

    df.sort_values("candidate_id").to_parquet(out1, index=False)
    df.sort_values("candidate_id").to_parquet(out2, index=False)

    print(f"[OK] Saved: {out1} and {out2} (n={len(df)})")
    #out = out_dir / "candidates_struct_features.parquet"
    #df.sort_values("candidate_id").to_parquet(out, index=False)
    #print(f"[OK] Saved FINAL: {out} (n={len(df)})")

if __name__ == "__main__":
    main()
