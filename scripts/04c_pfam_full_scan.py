#!/usr/bin/env python3
import argparse, subprocess
from pathlib import Path
import pandas as pd
from utils_config import load_config

def run(cmd):
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        raise RuntimeError(res.stderr.strip()[:2000])
    return res

def parse_domtblout(domtbl_path: Path):
    rows = []
    with open(domtbl_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            # domtblout stable columns:
            # 0 target_name, 3 query_name, 6 full_Evalue, 12 dom_Evalue, 13 dom_score, 17 ali_from, 18 ali_to
            target = parts[0]          # sequence id in candidates fasta = candidate_id
            query  = parts[3]          # Pfam model name
            full_e = float(parts[6])
            dom_e  = float(parts[12])
            dom_score = float(parts[13])
            ali_from = int(parts[17])
            ali_to   = int(parts[18])

            rows.append({
                "candidate_id": target,
                "pfam_model": query,
                "full_evalue": full_e,
                "dom_evalue": dom_e,
                "dom_score": dom_score,
                "ali_from": ali_from,
                "ali_to": ali_to,
            })
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)

    tmp_dir = Path(cfg["paths"]["tmp_dir"])
    feat_dir = Path(cfg["paths"]["features_dir"]) / "hmm"
    feat_dir.mkdir(parents=True, exist_ok=True)

    candidates_faa = tmp_dir / "hmm" / "candidates_all.faa"
    if not candidates_faa.exists():
        raise SystemExit(f"[ERR] candidates fasta not found: {candidates_faa}. Run 04a_build_candidates_fasta.py first.")

    hmm_cfg = cfg.get("hmm", {})
    pfam_db = Path(hmm_cfg["pfam_db_hmm"])
    if not pfam_db.exists():
        raise SystemExit(f"[ERR] pfam_db_hmm not found: {pfam_db}")

    cpu = str(hmm_cfg.get("cpu", 8))
    evalue = str(hmm_cfg.get("evalue", 1e-3))

    work_dir = tmp_dir / "hmm"
    work_dir.mkdir(parents=True, exist_ok=True)

    domtbl = work_dir / "pfam_full.domtblout"

    # hmmsearch: Pfam HMMs (queries) vs sequences (targets)
    cmd = [
        "hmmsearch",
        "--cpu", cpu,
        "--domtblout", str(domtbl),
        "-E", evalue,
        str(pfam_db),
        str(candidates_faa),
    ]
    print("[RUN]", " ".join(cmd))
    run(cmd)

    rows = parse_domtblout(domtbl)
    hits_path = feat_dir / "pfam_hits.tsv"
    if rows:
        hits = pd.DataFrame(rows)
        hits.sort_values(["candidate_id", "dom_evalue", "dom_score"], ascending=[True, True, False], inplace=True)
        hits.to_csv(hits_path, sep="\t", index=False)
        print(f"[OK] pfam hits saved: {hits_path} ({len(hits)} domain hits)")
    else:
        pd.DataFrame(columns=["candidate_id","pfam_model","full_evalue","dom_evalue","dom_score","ali_from","ali_to"])\
          .to_csv(hits_path, sep="\t", index=False)
        print(f"[OK] pfam hits saved (empty): {hits_path}")

    # Summary per candidate
    if rows:
        best = hits.groupby("candidate_id", as_index=False).first()
        summary = (hits.groupby("candidate_id", as_index=False)
                     .agg(n_pfam_hits=("dom_evalue","count"),
                          best_dom_evalue=("dom_evalue","min"),
                          best_dom_score=("dom_score","max")))
        summary = summary.merge(best[["candidate_id","pfam_model","dom_evalue","dom_score"]], on="candidate_id", how="left")
        summary = summary.rename(columns={
            "pfam_model":"best_pfam_model",
            "dom_evalue":"best_pfam_dom_evalue",
            "dom_score":"best_pfam_dom_score",
        })
    else:
        summary = pd.DataFrame(columns=[
            "candidate_id","n_pfam_hits","best_dom_evalue","best_dom_score",
            "best_pfam_model","best_pfam_dom_evalue","best_pfam_dom_score"
        ])

    summary_path = feat_dir / "candidates_pfam_summary.parquet"
    summary.sort_values("candidate_id").to_parquet(summary_path, index=False)
    print(f"[OK] pfam summary saved: {summary_path} ({len(summary)} candidates)")

if __name__ == "__main__":
    main()
