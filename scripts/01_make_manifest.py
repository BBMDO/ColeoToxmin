#!/usr/bin/env python3
import argparse
from pathlib import Path
import yaml
import pandas as pd

def load_cfg(p):
    with open(p, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def parse_fasta_headers(fasta_path: Path):
    ids = []
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids

def find_alphafold_files(af_root: Path, cand_id: str):
    # tenta: exato / lower / upper
    candidates = [cand_id, cand_id.lower(), cand_id.upper()]
    for name in candidates:
        d = af_root / name
        if d.is_dir():
            cif = next(iter(d.glob("fold_*_model_0.cif")), None)
            js  = next(iter(d.glob("fold_*_full_data_0.json")), None)
            return str(d), str(cif) if cif else "", str(js) if js else ""
    return "", "", ""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_cfg(args.config)

    toxin_dir = Path(cfg["paths"]["toxin_dir"])
    af_dir    = Path(cfg["paths"]["alphafold_dir"])
    work_dir  = Path(cfg["paths"]["work_dir"])
    work_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for fasta in sorted(toxin_dir.rglob("*.toxin_renomeado.fasta")):
        # tenta inferir species como a pasta logo após toxin/
        try:
            species = fasta.relative_to(toxin_dir).parts[0]
        except Exception:
            species = "unknown"

        for cid in parse_fasta_headers(fasta):
            af_folder, cif_path, json_path = find_alphafold_files(af_dir, cid)
            rows.append({
                "candidate_id": cid,
                "candidate_id_lower": cid.lower(),
                "species": species,
                "tox_fasta_path": str(fasta),
                "alphafold_folder": af_folder,
                "alphafold_cif_path": cif_path,
                "alphafold_json_path": json_path,
                "has_structure": 1 if cif_path else 0
            })

    df = pd.DataFrame(rows).drop_duplicates(subset=["candidate_id","tox_fasta_path"])
    out = work_dir / "manifest.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"[OK] manifest: {out}  (n={len(df)})")
    print(df[["candidate_id","species","has_structure"]].head(10).to_string(index=False))

if __name__ == "__main__":
    main()
