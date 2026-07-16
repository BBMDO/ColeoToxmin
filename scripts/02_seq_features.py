#!/usr/bin/env python3
import argparse
from pathlib import Path
import re
import pandas as pd
import hashlib
import numpy as np

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from utils_config import load_config

from pathlib import Path
import pandas as pd

import re

def fasta_headers(fasta_path: Path):
    """Retorna lista de headers (sem '>') na ordem do FASTA."""
    headers = []
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                headers.append(line[1:].strip().split()[0])  # pega só o primeiro token
    return headers

def read_signalp_summary(signalp_file: Path):
    """Parse SignalP5 summary (.signalp5) robustly.

    The file format often looks like:
      #_SignalP-5.0   Organism:_euk   Timestamp:_...
      #_ID    Prediction      SP(Sec/SPI)     OTHER   CS_Position
      <ID>    SP(Sec/SPI)     0.999902        0.000098        CS_pos:_19-20._ALA-AA._Pr:_0.4021

    Pandas sometimes raises ParserError because the last column can contain extra
    whitespace (or the file mixes tabs/spaces). We parse line-by-line and split
    with maxsplit=4 so any remainder stays in the last field.

    Returns:
      dict[id] = (prediction, sp_score)
    """
    signalp_file = Path(signalp_file)
    out: dict[str, tuple[str, float]] = {}

    with open(signalp_file, "r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            # ID, Prediction, SP_score, OTHER_score, CS_Position(remainder)
            parts = line.split(None, 4)
            if len(parts) < 3:
                continue

            sid = parts[0].strip()
            pred = parts[1].strip()

            try:
                sp_score = float(parts[2])
            except Exception:
                sp_score = 0.0

            out[sid] = (pred, sp_score)

    return out


def read_fasta_seqs(fasta_path: Path):
    """Return dict[id] = sequence (uppercase, no spaces) from a FASTA."""
    fasta_path = Path(fasta_path)
    seqs = {}
    cur_id = None
    cur = []
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(cur).upper().replace(" ", "")
                cur_id = line[1:].strip().split()[0]
                cur = []
            else:
                cur.append(line)
        if cur_id is not None:
            seqs[cur_id] = "".join(cur).upper().replace(" ", "")
    return seqs

def seq_hash(seq: str):
    """Stable hash for mapping sequences."""
    return hashlib.sha1(seq.encode("utf-8")).hexdigest()

def find_signalp_file_for_fasta(fasta_path: Path):
    """
    Heurística: procura o .signalp5 no mesmo diretório do FASTA,
    ou no diretório da espécie em 00_inputs/toxin/<species>/.
    """
    fasta_path = Path(fasta_path)
    # 1) mesmo diretório
    cand = list(fasta_path.parent.glob("*.toxin_signalp_summary.signalp5"))
    if cand:
        return cand[0]
    # 2) diretório da espécie (um nível acima do fasta pode variar)
    cand = list(fasta_path.parent.rglob("*.toxin_signalp_summary.signalp5"))
    if cand:
        return cand[0]
    return None

def load_signalp_map(signalp_root: Path):
    """
    Lê arquivos SignalP5 no formato:
      #_SignalP-5.0 ...
      #_ID Prediction SP(Sec/SPI) OTHER CS_Position
      <ID> <Prediction> <SP_score> <OTHER_score> <CS_Position...>

    Retorna:
      dict[id] = {"signalp_prediction": str, "signalp_score": float, "signalp_call": int}
    """
    signalp_root = Path(signalp_root)
    out = {}

    files = list(signalp_root.rglob("*.toxin_signalp_summary.signalp5"))
    if not files:
        print(f"[WARN] No SignalP summary files found under: {signalp_root}")
        return out

    for fp in files:
        try:
            # IGNORA linhas que começam com '#'
            # e usa separação por whitespace (tabs ou espaços)
            df = pd.read_csv(
                fp,
                sep=r"\s+",
                comment="#",
                engine="python",
                dtype=str
            )
            if df.empty:
                continue

            # colunas esperadas
            # primeiro campo geralmente é o ID
            id_col = df.columns[0]

            # achar Prediction
            pred_col = None
            for c in df.columns:
                if c.lower() == "prediction":
                    pred_col = c
                    break
            if pred_col is None:
                # fallback: procurar substring
                for c in df.columns:
                    if "prediction" in c.lower():
                        pred_col = c
                        break

            # score de SP: pode ser "SP(Sec/SPI)" ou parecido
            sp_col = None
            for c in df.columns:
                if "sp(" in c.lower() or c.lower().startswith("sp"):
                    sp_col = c
                    break

            for _, row in df.iterrows():
                cid = str(row.get(id_col, "")).strip()
                if not cid:
                    continue

                pred = str(row.get(pred_col, "")).strip() if pred_col else ""
                pred_upper = pred.upper()

                # chama secretado se predição for SP (Sec/SPI)
                call = 1 if "SP" in pred_upper else 0

                sp_score = 0.0
                if sp_col is not None:
                    try:
                        sp_score = float(str(row.get(sp_col, "0")).strip())
                    except Exception:
                        sp_score = 0.0

                out[cid] = {
                    "signalp_prediction": pred,
                    "signalp_score": sp_score,
                    "signalp_call": call
                }

        except Exception as e:
            print(f"[WARN] Could not parse SignalP file: {fp} ({e})")
            continue

    return out

def parse_fasta(path: Path):
    sid = None
    buf = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sid is not None:
                    yield sid, "".join(buf)
                sid = line[1:].strip()
                buf = []
            else:
                buf.append(line)
        if sid is not None:
            yield sid, "".join(buf)

def guess_species_from_path(toxin_dir: Path, fasta_path: Path) -> str:
    # tenta pegar o primeiro nível depois de toxin_dir
    try:
        rel = fasta_path.relative_to(toxin_dir)
        parts = rel.parts
        if len(parts) >= 2:
            return parts[0]
    except Exception:
        pass
    # fallback
    return "unknown"

def cys_grammar(seq: str):
    # retorna: string tipo C-x(3)-C-x(6)-C...
    pos = [i+1 for i,a in enumerate(seq) if a == "C"]
    if len(pos) < 2:
        return ("", "")
    gaps = [pos[i+1] - pos[i] - 1 for i in range(len(pos)-1)]
    grammar = "C" + "".join([f"-x({g})-C" for g in gaps])
    gaps_str = ",".join(map(str, gaps))
    return grammar, gaps_str

def net_charge_pH7(seq: str):
    # aproximação padrão
    seq = seq.upper()
    return (seq.count("K") + seq.count("R") + 0.1*seq.count("H")) - (seq.count("D") + seq.count("E"))

def compute_features(seq_raw: str):
    seq = seq_raw.replace("*","").upper()
    seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq)  # remove X etc
    L = len(seq)
    if L == 0:
        return None

    n_cys = seq.count("C")
    cys_density = n_cys / L

    pa = ProteinAnalysis(seq)
    try:
        gravy = float(pa.gravy())
    except Exception:
        gravy = np.nan
    try:
        arom = float(pa.aromaticity())
    except Exception:
        arom = np.nan
    try:
        instab = float(pa.instability_index())
    except Exception:
        instab = np.nan

    grammar, gaps_str = cys_grammar(seq)

    return {
        "len_aa": int(L),
        "n_cys": int(n_cys),
        "cys_density": float(cys_density),
        "cys_grammar": grammar,
        "cys_gaps": gaps_str,
        "net_charge_pH7": float(net_charge_pH7(seq)),
        "gravy": gravy,
        "aromaticity": arom,
        "instability_index": instab,
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)

    toxin_dir = Path(cfg["paths"]["toxin_dir"])
    feat_dir = Path(cfg["paths"]["features_dir"]) / "seq"
    feat_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    fasta_files = sorted(toxin_dir.rglob("*.toxin_renomeado.fasta"))
    if not fasta_files:
        raise SystemExit(f"[ERR] Nenhum *.toxin_renomeado.fasta encontrado em: {toxin_dir}")

    for fp in fasta_files:
        species = guess_species_from_path(toxin_dir, fp)
        for cand_id, seq in parse_fasta(fp):
            feats = compute_features(seq)
            if feats is None:
                continue
            rows.append({
                "candidate_id": cand_id,
                "species": species,
                "tox_fasta_path": str(fp),
                # SignalP (se você não tiver por enquanto, fica 0)
                "signalp_call": 0,
                "signalp_score": 0.0,
                **feats
            })

    df = pd.DataFrame(rows)
    if df.empty:
        raise SystemExit("[ERR] dataframe vazio. Verifique seus fastas.")
    # --- SignalP mapping por sequência: original fasta (Gene...) -> renomeado fasta (candidate_id) ---
    df["signalp_call"] = 0
    df["signalp_score"] = 0.0
    df["signalp_prediction"] = "OTHER"

    for fasta_ren, sub in df.groupby("tox_fasta_path", dropna=False):
        if not fasta_ren:
            continue
        fasta_ren = Path(fasta_ren)
        if not fasta_ren.exists():
            continue

        # tenta achar o fasta original correspondente
        fasta_orig = Path(str(fasta_ren).replace(".toxin_renomeado.fasta", ".toxin.fasta"))
        if not fasta_orig.exists():
            fasta_orig = Path(str(fasta_ren).replace("_renomeado", ""))
        if not fasta_orig.exists():
            continue

        # SignalP geralmente foi rodado no fasta original
        sig_file = find_signalp_file_for_fasta(fasta_orig) or find_signalp_file_for_fasta(fasta_ren)
        if sig_file is None:
            continue

        sig_map = read_signalp_summary(sig_file)  # dict[orig_id] = (pred, sp_score)
        orig = read_fasta_seqs(fasta_orig)        # Gene... -> seq
        ren  = read_fasta_seqs(fasta_ren)         # candidate_id -> seq

        ren_hash_to_id = {}
        for rid, s in ren.items():
            ren_hash_to_id[seq_hash(s)] = rid

        for orig_id, (pred, score) in sig_map.items():
            seq = orig.get(orig_id)
            if seq is None:
                continue
            h = seq_hash(seq)
            ren_id = ren_hash_to_id.get(h)
            if ren_id is None:
                continue

            call = 1 if "SP" in str(pred).upper() else 0
            mask = (df["candidate_id"] == ren_id)
            df.loc[mask, "signalp_call"] = call
            df.loc[mask, "signalp_score"] = float(score)
            df.loc[mask, "signalp_prediction"] = pred


    df["is_secreted"] = (pd.to_numeric(df["signalp_call"], errors="coerce")
                         .fillna(0).astype(int).eq(1).astype(int))
    df.loc[df["signalp_prediction"].astype(str).str.len() == 0, "signalp_prediction"] = df["signalp_call"].apply(lambda x: "SP" if int(x) == 1 else "OTHER")
    out = feat_dir / "candidates_seq_features.parquet"
    df.sort_values("candidate_id").to_parquet(out, index=False)
    print(f"[OK] Saved: {out} (n={len(df)})")
    print("[INFO] columns:", ", ".join(df.columns))

if __name__ == "__main__":
    main()
