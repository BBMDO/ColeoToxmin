#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import gemmi
from utils_config import load_config

def ca_coords(pdb: Path):
    st = gemmi.read_structure(str(pdb))
    coords=[]
    for model in st:
        for chain in model:
            for res in chain:
                if res.is_water():
                    continue
                ca = res.find_atom("CA", '\0')
                if ca is not None:
                    coords.append([ca.pos.x, ca.pos.y, ca.pos.z])
        break
    return np.asarray(coords, dtype=float) if coords else None

def rg(coords):
    c = coords.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum((coords-c)**2, axis=1))))

def contact_density(coords, cutoff):
    n = coords.shape[0]
    if n < 2:
        return 0.0, 0
    cutoff2 = cutoff*cutoff
    pairs=0
    for i in range(n):
        d = coords[i+1:] - coords[i]
        dist2 = np.sum(d*d, axis=1)
        pairs += int(np.sum(dist2 <= cutoff2))
    return float((2.0*pairs)/n), int(pairs)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)

    cutoff = float(cfg["structure"]["contact_cutoff_A"])
    pdb_dir = Path(cfg["paths"]["tmp_dir"]) / "pdb"
    out_dir = Path(cfg["paths"]["features_dir"]) / "struct"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows=[]
    for pdb in sorted(pdb_dir.glob("*.pdb")):
        cid = pdb.stem
        coords = ca_coords(pdb)
        if coords is None:
            print(f"[WARN] no CA: {cid}")
            continue
        rgv = rg(coords)
        cd, npairs = contact_density(coords, cutoff)
        rows.append({
            "candidate_id": cid,
            "n_ca": int(coords.shape[0]),
            "rg": rgv,
            "contact_density": cd,
            "n_contact_pairs": npairs
        })

    df = pd.DataFrame(rows).sort_values("candidate_id")
    out = out_dir / "contacts_rg_features.parquet"
    df.to_parquet(out, index=False)
    print(f"[OK] Saved {out} (n={len(df)})")

if __name__ == "__main__":
    main()
