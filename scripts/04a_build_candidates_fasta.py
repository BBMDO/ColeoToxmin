#!/usr/bin/env python3
import argparse
from pathlib import Path
from utils_config import load_config

def iter_candidate_fastas(toxin_dir: Path):
    for p in toxin_dir.rglob("*.toxin_renomeado.fasta"):
        yield p

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)

    toxin_dir = Path(cfg["paths"]["toxin_dir"])
    tmp_dir = Path(cfg["paths"]["tmp_dir"])
    out_dir = tmp_dir / "hmm"
    out_dir.mkdir(parents=True, exist_ok=True)

    out_fasta = out_dir / "candidates_all.faa"

    n_seq = 0
    with open(out_fasta, "w", encoding="utf-8") as out:
        for fasta_path in sorted(iter_candidate_fastas(toxin_dir)):
            txt = fasta_path.read_text(encoding="utf-8", errors="ignore").strip()
            if not txt:
                continue
            if not txt.endswith("\n"):
                txt += "\n"
            out.write(txt)
            n_seq += sum(1 for line in txt.splitlines() if line.startswith(">"))

    print(f"[OK] candidates fasta: {out_fasta} (n={n_seq} sequences)")

if __name__ == "__main__":
    main()
