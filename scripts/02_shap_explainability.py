#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
02_shap_explainability.py

Fecha a Lacuna B: gera outputs de SHAP (global + local) para explicar
quais features dirigem a classificação toxin-like.

Input:
- all_candidates_features.tsv (precisa ter candidate_id e alguma coluna-alvo: pu_label ou is_toxinlike ou similar)

Outputs:
- <out_prefix>.shap_global_bar.pdf
- <out_prefix>.shap_summary_beeswarm.pdf
- <out_prefix>.shap_top5_waterfalls.pdf
- <out_prefix>.shap_importance_meanabs.tsv
- <out_prefix>.top5_cases.tsv
- <out_prefix>.model_metrics.txt

Estratégia:
- Treina modelo proxy (XGBoost) com as features do TSV
- Calcula SHAP (TreeExplainer)
- Salva figuras e tabelas
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, average_precision_score

# XGBoost + SHAP
from xgboost import XGBClassifier
import shap


DEFAULT_FEATURES = [
    # sequence / biophysics (ajuste se quiser)
    "len_aa",
    "n_cys",
    "cys_density",
    "net_charge_pH7",
    "gravy",
    "aromaticity",
    "instability_index",
    "boman_index",
    # secretion / signal
    "signalp_score",
    "is_secreted",
    "signalp_call",
    # structure (se tiver)
    "plddt_mean",
    "pLDDT_mean",
    "helix_frac",
    "sheet_frac",
    "coil_frac",
    "contact_density",
    # ranking
    "priority_score_simple",
]

# nomes alternativos comuns
ALT_COLS = {
    "plddt_mean": ["plddt_mean", "pLDDT_mean", "plddt_avg", "mean_plddt"],
    "helix_frac": ["helix_frac", "sec_helix_frac", "helix_fraction"],
    "sheet_frac": ["sheet_frac", "sec_sheet_frac", "sheet_fraction"],
    "coil_frac": ["coil_frac", "sec_coil_frac", "coil_fraction"],
}


def pick_existing_columns(df, wanted):
    """Retorna apenas colunas existentes, com fallback para aliases em ALT_COLS."""
    cols = []
    for c in wanted:
        if c in df.columns:
            cols.append(c)
            continue
        if c in ALT_COLS:
            for alt in ALT_COLS[c]:
                if alt in df.columns:
                    cols.append(alt)
                    break
    # remove duplicadas preservando ordem
    seen = set()
    out = []
    for c in cols:
        if c not in seen:
            out.append(c)
            seen.add(c)
    return out


def resolve_label_column(df, label_candidates):
    for c in label_candidates:
        if c in df.columns:
            return c
    return None


def coerce_numeric(df, cols):
    out = df.copy()
    for c in cols:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_tsv", required=True, help="TSV com features (ex.: all_candidates_features.tsv)")
    ap.add_argument("--out_prefix", required=True, help="Prefixo de saída (sem extensão)")
    ap.add_argument("--label_col", default="", help="Coluna alvo (ex.: pu_label). Se vazio, tenta detectar.")
    ap.add_argument("--label_candidates", default="pu_label,is_toxinlike,toxinlike_label,label",
                    help="Se label_col vazio, tenta uma dessas (separadas por vírgula).")
    ap.add_argument("--features", default="", help="Lista de features separadas por vírgula. Se vazio, usa default.")
    ap.add_argument("--top_k", type=int, default=5, help="Quantos casos top para waterfall")
    ap.add_argument("--test_size", type=float, default=0.2)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max_rows", type=int, default=0, help="Opcional: limitar linhas (0 = sem limite)")
    args = ap.parse_args()

    df = pd.read_csv(args.in_tsv, sep="\t")

    if args.max_rows and args.max_rows > 0:
        df = df.head(args.max_rows).copy()

    if "candidate_id" not in df.columns:
        print("[ERROR] Precisa ter coluna candidate_id no TSV.", file=sys.stderr)
        sys.exit(1)

    # label
    if args.label_col.strip():
        ycol = args.label_col.strip()
        if ycol not in df.columns:
            print(f"[ERROR] label_col='{ycol}' não existe no TSV.", file=sys.stderr)
            sys.exit(1)
    else:
        candidates = [x.strip() for x in args.label_candidates.split(",") if x.strip()]
        ycol = resolve_label_column(df, candidates)
        if ycol is None:
            print("[ERROR] Não achei coluna-alvo. Use --label_col pu_label (ou equivalente).", file=sys.stderr)
            print(f"       Candidatas tentadas: {candidates}", file=sys.stderr)
            sys.exit(1)

    # features
    if args.features.strip():
        wanted = [x.strip() for x in args.features.split(",") if x.strip()]
    else:
        wanted = DEFAULT_FEATURES

    feat_cols = pick_existing_columns(df, wanted)
    if len(feat_cols) < 4:
        print("[ERROR] Poucas features disponíveis após checar existência/aliases.", file=sys.stderr)
        print("        Colunas detectadas:", feat_cols, file=sys.stderr)
        sys.exit(1)

    # preparar X/y
    # y precisa ser 0/1
    y = pd.to_numeric(df[ycol], errors="coerce").fillna(0).astype(int).values

    Xraw = df[feat_cols].copy()
    Xraw = coerce_numeric(Xraw, feat_cols)

    # imputação simples (mediana)
    X = Xraw.fillna(Xraw.median(numeric_only=True))

    # split
    X_train, X_test, y_train, y_test, df_train, df_test = train_test_split(
        X, y, df[["candidate_id", "species"]].copy() if "species" in df.columns else df[["candidate_id"]].copy(),
        test_size=args.test_size, random_state=args.seed, stratify=y if len(np.unique(y)) > 1 else None
    )

    # modelo proxy (XGBoost)
    # parâmetros seguros e estáveis
    model = XGBClassifier(
        n_estimators=500,
        max_depth=4,
        learning_rate=0.05,
        subsample=0.9,
        colsample_bytree=0.9,
        reg_lambda=1.0,
        random_state=args.seed,
        eval_metric="logloss",
        n_jobs=4
    )
    model.fit(X_train, y_train)

    # métricas (para você reportar como sanity check)
    proba_test = model.predict_proba(X_test)[:, 1]
    try:
        auc = roc_auc_score(y_test, proba_test) if len(np.unique(y_test)) > 1 else float("nan")
    except Exception:
        auc = float("nan")
    try:
        ap_score = average_precision_score(y_test, proba_test) if len(np.unique(y_test)) > 1 else float("nan")
    except Exception:
        ap_score = float("nan")

    # SHAP
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_test)

    # 1) importância global (mean |SHAP|)
    mean_abs = np.abs(shap_values).mean(axis=0)
    imp = pd.DataFrame({"feature": feat_cols, "mean_abs_shap": mean_abs}).sort_values("mean_abs_shap", ascending=False)
    imp_path = f"{args.out_prefix}.shap_importance_meanabs.tsv"
    imp.to_csv(imp_path, sep="\t", index=False)

    # 2) Figura bar global
    bar_pdf = f"{args.out_prefix}.shap_global_bar.pdf"
    plt.figure()
    plt.bar(imp["feature"].head(20), imp["mean_abs_shap"].head(20))
    plt.xticks(rotation=75, ha="right")
    plt.ylabel("mean(|SHAP|)")
    plt.title("Global feature importance (SHAP, proxy XGBoost)")
    plt.tight_layout()
    plt.savefig(bar_pdf)
    plt.close()

    # 3) Beeswarm (summary plot) — precisa do shap.summary_plot
    beeswarm_pdf = f"{args.out_prefix}.shap_summary_beeswarm.pdf"
    plt.figure()
    shap.summary_plot(shap_values, X_test, show=False)
    plt.tight_layout()
    plt.savefig(beeswarm_pdf)
    plt.close()

    # 4) Top-K casos por score (probabilidade) + waterfall
    # pega maiores proba como “casos top”
    order = np.argsort(-proba_test)
    topk = order[: args.top_k]
    top_cases = df_test.iloc[topk].copy()
    top_cases["proxy_proba"] = proba_test[topk]
    top_cases_path = f"{args.out_prefix}.top{args.top_k}_cases.tsv"
    top_cases.to_csv(top_cases_path, sep="\t", index=False)

    # Waterfalls em um PDF
    wf_pdf = f"{args.out_prefix}.shap_top{args.top_k}_waterfalls.pdf"
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(wf_pdf) as pdf:
        for i in topk:
            # base values / expected value
            # shap.Explanation facilita waterfall
            ex = shap.Explanation(
                values=shap_values[i],
                base_values=explainer.expected_value,
                data=X_test.iloc[i].values,
                feature_names=feat_cols
            )
            plt.figure()
            shap.plots.waterfall(ex, max_display=20, show=False)
            cid = df_test.iloc[i]["candidate_id"]
            plt.title(f"Waterfall (candidate_id={cid})")
            plt.tight_layout()
            pdf.savefig()
            plt.close()

    # 5) métricas
    metrics_path = f"{args.out_prefix}.model_metrics.txt"
    with open(metrics_path, "w", encoding="utf-8") as f:
        f.write(f"label_col\t{ycol}\n")
        f.write(f"n_rows\t{df.shape[0]}\n")
        f.write(f"n_features\t{len(feat_cols)}\n")
        f.write(f"features_used\t{','.join(feat_cols)}\n")
        f.write(f"test_size\t{args.test_size}\n")
        f.write(f"seed\t{args.seed}\n")
        f.write(f"ROC_AUC\t{auc}\n")
        f.write(f"AveragePrecision\t{ap_score}\n")

    print("[OK] Wrote:")
    print(" -", imp_path)
    print(" -", bar_pdf)
    print(" -", beeswarm_pdf)
    print(" -", top_cases_path)
    print(" -", wf_pdf)
    print(" -", metrics_path)
    print(f"[INFO] label={ycol} | AUC={auc} | AP={ap_score}")


if __name__ == "__main__":
    main()
