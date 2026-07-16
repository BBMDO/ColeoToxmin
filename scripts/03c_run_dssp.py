#!/usr/bin/env python3
import argparse, subprocess
from pathlib import Path
import pandas as pd
from utils_config import load_config

def run_dssp(dssp_exec: str, cif_in: Path, dssp_out: Path):
    dssp_out.parent.mkdir(parents=True, exist_ok=True)
    # mkdssp [options] input-file [output-file]
    cmd = [dssp_exec, "--output-format", "dssp", str(cif_in), str(dssp_out)]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        raise RuntimeError((res.stderr or res.stdout).strip()[:800])

def parse_dssp_ss(dssp_path: Path):
    helix = set(["H","G","I"])
    sheet = set(["E","B"])
    total = 0
    h = e = c = 0

    lines = dssp_path.read_text(encoding="utf-8", errors="ignore").splitlines()

    start = None
    for i, line in enumerate(lines):
        if line.startswith("  #  RESIDUE"):
            start = i + 1
            break
    if start is None:
        return None

    for line in lines[start:]:
        if len(line) < 18:
            continue
        if line[13] == "!":   # chain break
            continue
        ss = line[16]
        total += 1
        if ss in helix:
            h += 1
        elif ss in sheet:
            e += 1
        else:
            c += 1

    if total == 0:
        return None
    return {
        "dssp_n_res": total,
        "sec_helix_frac": h/total,
        "sec_sheet_frac": e/total,
        "sec_coil_frac": c/total,
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)

    dssp_exec = cfg["structure"]["dssp_exec"]
    alphafold_dir = Path(cfg["paths"]["alphafold_dir"])
    tmp_dir = Path(cfg["paths"]["tmp_dir"])
    dssp_dir = tmp_dir / "dssp"
    out_dir = Path(cfg["paths"]["features_dir"]) / "struct"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for cand_dir in sorted(alphafold_dir.iterdir()):
        if not cand_dir.is_dir():
            continue
        cand_id = cand_dir.name
        cifs = list(cand_dir.glob("fold_*_model_0.cif"))
        if not cifs:
            continue
        cif_path = cifs[0]
        dssp_path = dssp_dir / f"{cand_id}.dssp"
        try:
            run_dssp(dssp_exec, cif_path, dssp_path)
            parsed = parse_dssp_ss(dssp_path)
            if parsed is None:
                print(f"[WARN] DSSP parse failed for {cand_id}")
                continue
            rows.append({"candidate_id": cand_id, **parsed})
        except Exception as e:
            print(f"[WARN] DSSP failed for {cand_id}: {e}")

    df = pd.DataFrame(rows).sort_values("candidate_id")
    out_path = out_dir / "dssp_features.parquet"
    df.to_parquet(out_path, index=False)
    print(f"[OK] DSSP features saved: {out_path} ({len(df)} candidates)")

if __name__ == "__main__":
    main()
