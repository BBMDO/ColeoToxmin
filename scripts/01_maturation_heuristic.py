#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
01_maturation_heuristic.py

Heurística simples e documentável para estimar arquitetura:
Signal peptide -> (propeptídeo) -> peptídeo maduro

Inputs (do seu pacote):
- all_candidates_features.tsv (features por candidato)
- candidates_all.faa (sequências)

Outputs:
- <out_prefix>.maturation.tsv
- <out_prefix>.summary.txt
- <out_prefix>.plots.pdf

Obs:
- Como o TSV não contém 'signalp_end', o script estima o corte do SP por:
  (i) motivo AXA (Ala-X-Ala) na janela 15..35 (típico de SPI)
  (ii) fallback default (20 aa)
"""

import argparse
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Optional, Tuple, List

import pandas as pd
import matplotlib.pyplot as plt


DIBASIC_RE = re.compile(r"(KR|RR|RK|KK)")
AXA_RE = re.compile(r"A.A")  # Ala-any-Ala


def norm_id(s: str) -> str:
    """Normaliza IDs para cruzamento candidato_id <-> FASTA header."""
    return s.strip().lstrip(">").strip().lower()


def read_fasta_as_dict(fasta_path: str) -> Dict[str, str]:
    """Lê FASTA e retorna dict {norm_header: sequence_without_stop}."""
    seqs: Dict[str, str] = {}
    cur_id = None
    cur_seq = []

    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[norm_id(cur_id)] = "".join(cur_seq).replace("*", "")
                cur_id = line[1:].strip()
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id is not None:
            seqs[norm_id(cur_id)] = "".join(cur_seq).replace("*", "")
    return seqs


def estimate_signal_peptide_end(
    seq: str,
    default_sp_len: int = 20,
    axa_window: Tuple[int, int] = (15, 35),
) -> int:
    """
    Estima SP_end (1-based na biologia, mas aqui retornamos 0-based index do corte):
    retorno = número de aa do SP (ex.: 20 significa que SP é seq[0:20] e pós-SP começa em 20)
    """
    n = len(seq)
    w0, w1 = axa_window
    w0 = max(0, min(w0, n))
    w1 = max(0, min(w1, n))

    # Busca AXA na janela típica de clivagem SPI
    window = seq[w0:w1]
    m = AXA_RE.search(window)
    if m:
        # corte geralmente após o motivo; colocamos o corte no final do motivo (posição +3)
        sp_end = w0 + m.end()
        # limitações de sanidade
        if 12 <= sp_end <= 35:
            return sp_end

    # fallback
    return min(default_sp_len, n)


def find_dibasic_site(
    seq_after_sp: str,
    search_window: int = 60
) -> Optional[Tuple[str, int]]:
    """
    Procura KR/RR/RK/KK na região pós-signal (primeiros `search_window` aa).
    Retorna (motif, pos_end) onde pos_end é o índice (0-based) do final do motivo na string seq_after_sp.
    Ex.: se motif em 10..12, pos_end=12.
    """
    region = seq_after_sp[:search_window]
    m = DIBASIC_RE.search(region)
    if not m:
        return None
    motif = m.group(1)
    return motif, m.end()


@dataclass
class MaturationResult:
    candidate_id: str
    species: str
    seq: str
    sp_end: int
    sp_method: str
    dibasic_motif: str
    dibasic_pos_after_sp: int
    mature_start: int
    propeptide_len: int
    mature_len: int
    mature_seq: str
    architecture: str


def predict_maturation(
    candidate_id: str,
    species: str,
    seq: str,
    signalp_call: int,
    default_sp_len: int,
    axa_window: Tuple[int, int],
    dibasic_window: int,
    fallback_pro_len: int,
    min_mature_len: int,
    max_mature_len: int
) -> MaturationResult:
    """
    Define:
    - sp_end (aa): fim do signal peptide
    - dibasic (se tiver)
    - mature_start (0-based index no seq)
    - propeptide_len
    - mature_len/seq
    - architecture label
    """
    seq = seq.replace("*", "")
    n = len(seq)

    if int(signalp_call) == 1 and n > 0:
        sp_end = estimate_signal_peptide_end(seq, default_sp_len=default_sp_len, axa_window=axa_window)
        sp_method = f"AXA_or_default({default_sp_len})"
    else:
        sp_end = 0
        sp_method = "no_signalp"

    after_sp = seq[sp_end:]

    dibasic = find_dibasic_site(after_sp, search_window=dibasic_window)
    if dibasic:
        motif, endpos = dibasic
        dibasic_motif = motif
        dibasic_pos_after_sp = endpos  # pos na string after_sp (fim do motivo)
        mature_start = sp_end + endpos
        propeptide_len = max(0, mature_start - sp_end)
        architecture = "signal+pro+mat (dibasic)"
    else:
        dibasic_motif = ""
        dibasic_pos_after_sp = -1
        mature_start = min(n, sp_end + fallback_pro_len)
        propeptide_len = max(0, mature_start - sp_end)
        architecture = "signal+pro?+mat (fallback)"

    mature_seq = seq[mature_start:]
    mature_len = len(mature_seq)

    # marca como fora de faixa (mas ainda reporta)
    if mature_len < min_mature_len or mature_len > max_mature_len:
        architecture = architecture + " [mature_len_out_of_range]"

    return MaturationResult(
        candidate_id=candidate_id,
        species=species,
        seq=seq,
        sp_end=sp_end,
        sp_method=sp_method,
        dibasic_motif=dibasic_motif,
        dibasic_pos_after_sp=dibasic_pos_after_sp,
        mature_start=mature_start,
        propeptide_len=propeptide_len,
        mature_len=mature_len,
        mature_seq=mature_seq,
        architecture=architecture
    )


def save_summary(path: str, stats: Dict[str, str]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for k, v in stats.items():
            f.write(f"{k}\t{v}\n")


def make_plots(df_out: pd.DataFrame, pdf_path: str) -> None:
    """
    Gera PDF com:
    - hist mature_len
    - hist propeptide_len
    - bar dibasic motifs
    - scatter mature_len vs n_cys (se houver coluna)
    """
    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(pdf_path) as pdf:
        # 1) mature_len
        fig = plt.figure()
        plt.hist(df_out["mature_len"].astype(int), bins=30)
        plt.xlabel("mature_len (aa)")
        plt.ylabel("count")
        plt.title("Distribution of predicted mature peptide length")
        pdf.savefig(fig); plt.close(fig)

        # 2) propeptide_len
        fig = plt.figure()
        plt.hist(df_out["propeptide_len"].astype(int), bins=30)
        plt.xlabel("propeptide_len (aa)")
        plt.ylabel("count")
        plt.title("Distribution of predicted propeptide length")
        pdf.savefig(fig); plt.close(fig)

        # 3) dibasic motifs
        motifs = df_out["dibasic_motif"].fillna("").astype(str)
        c = Counter([m for m in motifs if m])
        if c:
            fig = plt.figure()
            xs = list(c.keys())
            ys = [c[k] for k in xs]
            plt.bar(xs, ys)
            plt.xlabel("dibasic motif")
            plt.ylabel("count")
            plt.title("Dibasic cleavage motifs (post-signal)")
            pdf.savefig(fig); plt.close(fig)

        # 4) scatter mature_len vs n_cys (se existir n_cys)
        if "n_cys" in df_out.columns:
            fig = plt.figure()
            x = df_out["mature_len"].astype(int)
            y = df_out["n_cys"].astype(int)
            plt.scatter(x, y, s=15)
            plt.xlabel("mature_len (aa)")
            plt.ylabel("n_cys (full precursor)")
            plt.title("Mature length vs cysteine count (precursor)")
            pdf.savefig(fig); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_features_tsv", required=True, help="TSV com features (ex.: all_candidates_features.tsv)")
    ap.add_argument("--in_fasta", required=True, help="FASTA com sequências (ex.: candidates_all.faa)")
    ap.add_argument("--out_prefix", required=True, help="Prefixo de saída (sem extensão)")

    ap.add_argument("--default_sp_len", type=int, default=20, help="Fallback para tamanho do signal peptide")
    ap.add_argument("--axa_window", default="15,35", help="Janela (ini,fim) para buscar motivo AXA")
    ap.add_argument("--dibasic_window", type=int, default=60, help="Janela pós-signal para buscar KR/RR/RK/KK")
    ap.add_argument("--fallback_pro_len", type=int, default=20, help="Se não achar dibásico, assume propeptídeo ~20 aa")

    ap.add_argument("--min_mature_len", type=int, default=8, help="Mínimo do maduro para considerar plausível")
    ap.add_argument("--max_mature_len", type=int, default=80, help="Máximo do maduro para considerar plausível")

    args = ap.parse_args()

    axa0, axa1 = [int(x) for x in args.axa_window.split(",")]

    df = pd.read_csv(args.in_features_tsv, sep="\t")
    seqs = read_fasta_as_dict(args.in_fasta)

    rows: List[dict] = []
    missing = 0

    for _, r in df.iterrows():
        cid = str(r["candidate_id"])
        key = norm_id(cid)
        seq = seqs.get(key)

        if not seq:
            missing += 1
            continue

        res = predict_maturation(
            candidate_id=cid,
            species=str(r.get("species", "")),
            seq=seq,
            signalp_call=int(r.get("signalp_call", 0)),
            default_sp_len=args.default_sp_len,
            axa_window=(axa0, axa1),
            dibasic_window=args.dibasic_window,
            fallback_pro_len=args.fallback_pro_len,
            min_mature_len=args.min_mature_len,
            max_mature_len=args.max_mature_len,
        )

        out = {
            "candidate_id": res.candidate_id,
            "species": res.species,
            "len_aa": len(res.seq),
            "signalp_call": int(r.get("signalp_call", 0)),
            "signalp_prediction": r.get("signalp_prediction", ""),
            "sp_end": res.sp_end,
            "sp_method": res.sp_method,
            "dibasic_motif": res.dibasic_motif,
            "dibasic_pos_after_sp": res.dibasic_pos_after_sp,
            "mature_start": res.mature_start,
            "propeptide_len": res.propeptide_len,
            "mature_len": res.mature_len,
            "architecture": res.architecture,
            "mature_seq": res.mature_seq,
        }

        # anexa colunas úteis do features TSV (se existirem)
        for col in ["n_cys", "cys_density", "cys_grammar", "net_charge_pH7", "gravy", "priority_score_simple", "is_secreted"]:
            if col in df.columns:
                out[col] = r.get(col)

        rows.append(out)

    df_out = pd.DataFrame(rows)

    out_tsv = f"{args.out_prefix}.maturation.tsv"
    out_txt = f"{args.out_prefix}.summary.txt"
    out_pdf = f"{args.out_prefix}.plots.pdf"

    df_out.to_csv(out_tsv, sep="\t", index=False)

    # Summary
    n_total = df.shape[0]
    n_reported = df_out.shape[0]
    n_dibasic = (df_out["dibasic_motif"].fillna("") != "").sum()
    in_range = df_out["architecture"].astype(str).str.contains("out_of_range").sum()
    stats = {
        "N_total_features": str(n_total),
        "N_with_sequence": str(n_reported),
        "N_missing_sequence": str(missing),
        "dibasic_fraction": f"{(n_dibasic / max(1,n_reported)):.3f}",
        "mature_len_out_of_range_count": str(in_range),
        "default_sp_len": str(args.default_sp_len),
        "axa_window": f"{axa0},{axa1}",
        "dibasic_window": str(args.dibasic_window),
        "fallback_pro_len": str(args.fallback_pro_len),
        "min_mature_len": str(args.min_mature_len),
        "max_mature_len": str(args.max_mature_len),
    }
    save_summary(out_txt, stats)

    # Plots
    make_plots(df_out, out_pdf)

    print(f"[OK] Wrote: {out_tsv}")
    print(f"[OK] Wrote: {out_txt}")
    print(f"[OK] Wrote: {out_pdf}")
    print(f"[SUMMARY] N_total={n_total} | with_seq={n_reported} | missing_seq={missing} | dibasic={n_dibasic}")


if __name__ == "__main__":
    main()
