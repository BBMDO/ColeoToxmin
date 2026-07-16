#!/usr/bin/env python3
import argparse, subprocess
from pathlib import Path
import yaml
import pandas as pd
import numpy as np
import gemmi

def load_cfg(p):
    with open(p, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def cif_to_pdb(cif_path: Path, pdb_path: Path):
    st = gemmi.read_structure(str(cif_path))
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    st.write_pdb(str(pdb_path))

def plddt_from_pdb_bfactor(pdb_path: Path):
    st = gemmi.read_structure(str(pdb_path))
    vals = []
    for model in st:
        for chain in model:
            for res in chain:
                for atom in res:
                    if atom.element.name != "H":
                        vals.append(atom.b_iso)
        break
    if not vals:
        return None
    a = np.asarray(vals, dtype=float)
    return float(a.mean()), float(a.std())

def get_ca_coords(pdb_path: Path):
    st = gemmi.read_structure(str(pdb_path))
    coords = []
    for model in st:
        for chain in model:
            for res in chain:
                ca = res.find_atom("CA", '\0')
                if ca is not None:
                    coords.append([ca.pos.x, ca.pos.y, ca.pos.z])
        break
    if not coords:
        return None
    return np.asarray(coords, dtype=float)

def rg(coords: np.ndarray):
    c = coords.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum((coords-c)**2, axis=1))))

def disulfides_sg(pdb_path: Path, cutoff=2.2):
    st = gemmi.read_structure(str(pdb_path))
    pts=[]
    for model in st:
        for chain in model:
            for res in chain:
                if res.name not in ("CYS","CYX"):
                    continue
                sg = res.find_atom("SG", '\0')
                if sg:
                    pts.append((chain.name, int(res.seqid.num), np.array([sg.pos.x,sg.pos.y,sg.pos.z], float)))
        break
    pairs=[]
    c2=cutoff*cutoff
    for i in range(len(pts)):
        for j in range(i+1,len(pts)):
            d = pts[i][2]-pts[j][2]
            if float(np.dot(d,d)) <= c2:
                pairs.append((pts[i][0],pts[i][1],pts[j][0],pts[j][1]))
    return len(pts), len(pairs), ";".join([f"{a}{ra}-{b}{rb}" for a,ra,b,rb in pairs])

def run_dssp(dssp_exec: str, pdb_in: Path, dssp_out: Path):
    dssp_out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [dssp_exec, "-i", str(pdb_in), "-o", str(dssp_out)]
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[:300])

def parse_dssp_fractions(dssp_path: Path):
    helix=set("HGI")
    sheet=set("EB")
    with open(dssp_path, "r", encoding="utf-8", errors="ignore") as f:
        lines=f.readlines()
    start=None
    for i,l in enumerate(lines):
        if l.startswith("  #  RESIDUE"):
            start=i+1; break
    if start is None:
        return None
    h=e=c=tot=0
    for l in lines[start:]:
        if len(l)<18: 
            continue
        if l[13]=="!":
            continue
        ss=l[16]
        tot += 1
        if ss in helix: h += 1
        elif ss in sheet: e += 1
        else: c += 1
    if tot==0:
        return None
    return tot, h/tot, e/tot, c/tot

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_cfg(args.config)

    work = Path(cfg["paths"]["work_dir"])
    feat_dir = Path(cfg["paths"]["features_dir"]) / "struct"
    feat_dir.mkdir(parents=True, exist_ok=True)

    dssp_exec = cfg["structure"]["dssp_exec"]
    plddt_min = float(cfg["gating"]["plddt_min_mean"])

    man = pd.read_csv(work / "manifest.tsv", sep="\t")
    pdb_dir = work / "pdb"
    dssp_dir = work / "dssp"

    rows=[]
    for _, r in man[man["has_structure"]==1].iterrows():
        cid = r["candidate_id"]
        cif = Path(r["alphafold_cif_path"])
        if not cif.exists():
            continue

        pdb = pdb_dir / f"{cid}.pdb"
        try:
            cif_to_pdb(cif, pdb)
        except Exception:
            # tenta com lowercase (caso cid tenha caixa diferente)
            pdb = pdb_dir / f"{cid.lower()}.pdb"
            cif_to_pdb(cif, pdb)

        pld = plddt_from_pdb_bfactor(pdb)
        if pld is None:
            continue
        pld_mean, pld_std = pld

        coords = get_ca_coords(pdb)
        rg_val = rg(coords) if coords is not None else np.nan

        n_cys_sg, n_ss, ss_pairs = disulfides_sg(pdb, cutoff=float(cfg["structure"]["disulfide_cutoff_A"]))

        dssp_out = dssp_dir / f"{cid}.dssp"
        try:
            run_dssp(dssp_exec, pdb, dssp_out)
            parsed = parse_dssp_fractions(dssp_out)
        except Exception:
            parsed = None

        if parsed is None:
            nres = np.nan; h=e=c = np.nan
        else:
            nres, h, e, c = parsed

        rows.append({
            "candidate_id": cid,
            "plddt_mean": float(pld_mean),
            "plddt_std": float(pld_std),
            "struct_pass": 1 if pld_mean >= plddt_min else 0,
            "rg": float(rg_val) if not np.isnan(rg_val) else np.nan,
            "dssp_n_res": nres,
            "sec_helix_frac": h,
            "sec_sheet_frac": e,
            "sec_coil_frac": c,
            "n_cys_sg": int(n_cys_sg),
            "n_disulfide_pairs": int(n_ss),
            "disulfide_pairs": ss_pairs
        })

    df = pd.DataFrame(rows).drop_duplicates(subset=["candidate_id"])
    out = feat_dir / "struct_features.parquet"
    df.to_parquet(out, index=False)

    # sanity check simples (pega o seu erro “H/E=0 em tudo”)
    if len(df) > 0 and "sec_helix_frac" in df.columns:
        bad = ((df["sec_helix_frac"].fillna(0) + df["sec_sheet_frac"].fillna(0)) < 0.01).mean()
        print(f"[CHECK] frac(helix+sheet < 0.01) = {bad:.2f}")
        if bad > 0.8:
            print("[WARN] DSSP suspeito (cheque mkdssp / PDB).")

    print(f"[OK] struct features: {out} (n={len(df)})")

if __name__ == "__main__":
    main()
