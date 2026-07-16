#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
from utils_config import load_config

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)

    feat_dir = Path(cfg["paths"]["features_dir"])
    fold_sum = pd.read_csv(feat_dir / "fold" / "candidates_fold_self_summary.tsv", sep="\t")
    pfam = pd.read_csv(feat_dir / "hmm" / "pfam_hits.tsv", sep="\t")

    # Ajuste aqui os nomes das colunas se seu pfam_hits.tsv tiver cabeçalho diferente:
    # esperamos algo como:
    # candidate_id, pfam_model (ou hmm_query), dom_evalue (ou similar)
    # vou inferir:
    col_model = "pfam_model" if "pfam_model" in pfam.columns else ("hmm_query" if "hmm_query" in pfam.columns else None)
    col_eval  = "dom_evalue" if "dom_evalue" in pfam.columns else ("best_dom_evalue" if "best_dom_evalue" in pfam.columns else None)
    if col_model is None:
        raise SystemExit("[ERR] Não achei coluna de modelo Pfam (pfam_model/hmm_query).")

    # pega “melhor pfam por candidato” = menor evalue
    if "dom_evalue" in pfam.columns:
        best = pfam.sort_values(["candidate_id","dom_evalue"]).groupby("candidate_id", as_index=False).first()
        best = best[["candidate_id", col_model, "dom_evalue"]].rename(columns={col_model:"best_pfam", "dom_evalue":"best_pfam_evalue"})
    else:
        best = pfam.groupby("candidate_id", as_index=False).first()[["candidate_id", col_model]].rename(columns={col_model:"best_pfam"})
        best["best_pfam_evalue"] = None

    x = fold_sum.merge(best, on="candidate_id", how="left").fillna({"best_pfam":"NO_PFAM"})

    # por cluster, Pfam dominante
    dom = (x.groupby(["cluster_id","best_pfam"], as_index=False)
             .agg(n=("candidate_id","count"))
             .sort_values(["cluster_id","n"], ascending=[True, False])
             .groupby("cluster_id", as_index=False).first()
             .rename(columns={"best_pfam":"cluster_pfam_label","n":"cluster_pfam_count"}))

    out = feat_dir / "fold" / "fold_cluster_pfam_enrichment.tsv"
    dom.to_csv(out, sep="\t", index=False)
    print(f"[OK] saved: {out} (n={len(dom)})")

if __name__ == "__main__":
    main()
