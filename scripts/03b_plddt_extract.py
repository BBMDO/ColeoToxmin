#!/usr/bin/env python3
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import gemmi
from utils_config import load_config

def find_json(cand_dir: Path):
    j = list(cand_dir.glob("fold_*_full_data_0.json"))
    return j[0] if j else None

def plddt_from_json(json_path: Path):
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    arr = data.get("atom_plddts", None)
    if arr is None:
        return None
    a = np.asarray(arr, dtype=float)
    return a if a.size else None

def plddt_from_pdb_bfactor(pdb_path: Path):
    st = gemmi.read_structure(str(pdb_path))
    vals = []
    for model in st:
        for chain in model:
            for res in chain:
                for atom in res:
                    if atom.element.name != "H":
                        vals.append(atom.b_iso)
        break
    return np.asarray(vals, dtype=float) if vals else None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)

    alphafold_dir = Path(cfg["paths"]["alphafold_dir"])
    pdb_dir = Path(cfg["paths"]["tmp_dir"]) / "pdb"
    out_dir = Path(cfg["paths"]["features_dir"]) / "struct"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for cand_dir in sorted(alphafold_dir.iterdir()):
        if not cand_dir.is_dir():
            continue
        cand_id = cand_dir.name
        pdb_path = pdb_dir / f"{cand_id}.pdb"
        if not pdb_path.exists():
            continue

        plddt = None
        source = None

        json_path = find_json(cand_dir)
        if json_path and json_path.exists():
            plddt = plddt_from_json(json_path)
            if plddt is not None:
                source = "json"

        if plddt is None:
            plddt = plddt_from_pdb_bfactor(pdb_path)
            if plddt is not None:
                source = "pdb_bfactor"

        if plddt is None:
            print(f"[WARN] No pLDDT for {cand_id}")
            continue

        rows.append({
            "candidate_id": cand_id,
            "plddt_source": source,
            "plddt_mean": float(np.mean(plddt)),
            "plddt_std": float(np.std(plddt)),
            "plddt_frac_ge_70": float(np.mean(plddt >= 70.0)),
            "plddt_frac_ge_80": float(np.mean(plddt >= 80.0)),
            "plddt_frac_ge_90": float(np.mean(plddt >= 90.0)),
        })

    df = pd.DataFrame(rows).sort_values("candidate_id")
    out = out_dir / "plddt_features.parquet"
    df.to_parquet(out, index=False)
    print(f"[OK] Saved {out} (n={len(df)})")

if __name__ == "__main__":
    main()
