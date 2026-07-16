#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import gemmi
from utils_config import load_config

def cys_sg(pdb: Path):
    st = gemmi.read_structure(str(pdb))
    pts=[]
    for model in st:
        for chain in model:
            for res in chain:
                if res.name not in ("CYS","CYX"):
                    continue
                sg = res.find_atom("SG", '\0')
                if sg is None:
                    continue
                pts.append((chain.name, int(res.seqid.num), np.array([sg.pos.x, sg.pos.y, sg.pos.z], float)))
        break
    return pts

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)

    cutoff = float(cfg["structure"]["disulfide_cutoff_A"])
    cutoff2 = cutoff*cutoff

    pdb_dir = Path(cfg["paths"]["tmp_dir"]) / "pdb"
    out_dir = Path(cfg["paths"]["features_dir"]) / "struct"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows=[]
    for pdb in sorted(pdb_dir.glob("*.pdb")):
        cid = pdb.stem
        pts = cys_sg(pdb)
        pairs=[]
        for i in range(len(pts)):
            for j in range(i+1, len(pts)):
                d = pts[i][2] - pts[j][2]
                if float(np.dot(d,d)) <= cutoff2:
                    a = pts[i]; b = pts[j]
                    pairs.append(f"{a[0]}{a[1]}-{b[0]}{b[1]}")
        rows.append({
            "candidate_id": cid,
            "n_cys_sg": int(len(pts)),
            "n_disulfide_pairs": int(len(pairs)),
            "disulfide_pairs": ";".join(pairs)
        })

    df = pd.DataFrame(rows).sort_values("candidate_id")
    out = out_dir / "disulfide_features.parquet"
    df.to_parquet(out, index=False)
    print(f"[OK] Saved {out} (n={len(df)})")

if __name__ == "__main__":
    main()
