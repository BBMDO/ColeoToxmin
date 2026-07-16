#!/usr/bin/env python3
import argparse
from pathlib import Path
import gemmi
from utils_config import load_config

def find_cif_files(alphafold_dir: Path):
    for cand_dir in sorted(alphafold_dir.iterdir()):
        if not cand_dir.is_dir():
            continue
        cifs = list(cand_dir.glob("fold_*_model_0.cif"))
        if cifs:
            yield cand_dir.name, cifs[0]

def cif_to_pdb(cif_path: Path, pdb_path: Path):
    st = gemmi.read_structure(str(cif_path))
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    st.write_pdb(str(pdb_path))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)

    alphafold_dir = Path(cfg["paths"]["alphafold_dir"])
    out_pdb_dir = Path(cfg["paths"]["tmp_dir"]) / "pdb"

    n_ok = 0
    for cand_id, cif_path in find_cif_files(alphafold_dir):
        pdb_path = out_pdb_dir / f"{cand_id}.pdb"
        try:
            cif_to_pdb(cif_path, pdb_path)
            n_ok += 1
        except Exception as e:
            print(f"[WARN] CIF->PDB failed: {cand_id} :: {e}")

    print(f"[OK] Converted {n_ok} CIFs to PDB in {out_pdb_dir}")

if __name__ == "__main__":
    main()
