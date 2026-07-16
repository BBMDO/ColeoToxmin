#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import freesasa
import gemmi
from utils_config import load_config

MAX_ASA = {
    "A":121.0,"R":265.0,"N":187.0,"D":187.0,"C":148.0,"Q":214.0,"E":214.0,"G":97.0,"H":216.0,
    "I":195.0,"L":191.0,"K":230.0,"M":203.0,"F":228.0,"P":154.0,"S":143.0,"T":163.0,"W":264.0,
    "Y":255.0,"V":165.0
}

def seq_from_pdb(pdb: Path):
    st = gemmi.read_structure(str(pdb))
    seq=[]
    for model in st:
        for chain in model:
            for res in chain:
                if res.is_water(): 
                    continue
                aa = gemmi.find_tabulated_residue(res.name).one_letter_code
                if aa != "X":
                    seq.append(aa)
        break
    return "".join(seq)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)

    hydrophobic = set(cfg["structure"]["hydrophobic_residues"])
    rsa_thr = float(cfg["structure"]["rsa_core_threshold"])

    pdb_dir = Path(cfg["paths"]["tmp_dir"]) / "pdb"
    out_dir = Path(cfg["paths"]["features_dir"]) / "struct"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows=[]
    for pdb in sorted(pdb_dir.glob("*.pdb")):
        cid = pdb.stem
        try:
            fs = freesasa.Structure(str(pdb))
            res = freesasa.calc(fs)
            sasa_total = float(res.totalArea())
            areas = res.residueAreas()

            per_res=[]
            for chain_id, residues in areas.items():
                for resnum, a in residues.items():
                    per_res.append(float(a.total))

            seq = seq_from_pdb(pdb)
            n = min(len(seq), len(per_res))
            if n == 0:
                continue

            rsa=[]
            for i in range(n):
                aa = seq[i]
                denom = MAX_ASA.get(aa, None)
                if denom:
                    rsa.append(per_res[i]/denom)
            rsa = np.asarray(rsa, dtype=float)
            if rsa.size == 0:
                continue

            hydro_mask = np.array([1 if seq[i] in hydrophobic else 0 for i in range(min(n, len(rsa)))], dtype=int)
            exposed_mask = (rsa[:hydro_mask.size] > rsa_thr).astype(int)
            hydro_exposed_frac = float((hydro_mask*exposed_mask).sum() / hydro_mask.sum()) if hydro_mask.sum() else 0.0

            rows.append({
                "candidate_id": cid,
                "sasa_total": sasa_total,
                "rsa_mean": float(rsa.mean()),
                "rsa_std": float(rsa.std()),
                "core_fraction_rsa_le_thr": float((rsa <= rsa_thr).mean()),
                "hydrophobic_exposed_frac": hydro_exposed_frac,
            })
        except Exception as e:
            print(f"[WARN] SASA/RSA failed: {cid} :: {e}")

    df = pd.DataFrame(rows).sort_values("candidate_id")
    out = out_dir / "sasa_rsa_features.parquet"
    df.to_parquet(out, index=False)
    print(f"[OK] Saved {out} (n={len(df)})")

if __name__ == "__main__":
    main()
