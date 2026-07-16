#!/usr/bin/env bash
set -euo pipefail

CFG="${1:-09_reproducibility/configs/config.yaml}"

pyget () {
python - <<PY
import yaml
cfg=yaml.safe_load(open("$CFG","r"))
x=cfg
for k in "$1".split("."):
    x=x[k]
print(x)
PY
}

TMP_DIR="$(pyget paths.tmp_dir)"
FEAT_DIR="$(pyget paths.features_dir)"

WORK_DIR="${TMP_DIR}/foldseek_self"
DB="${WORK_DIR}/db"
CLU="${WORK_DIR}/clu"
TMPFS="${WORK_DIR}/tmp_clu"

mkdir -p "${WORK_DIR}"

# easy-cluster faz clustering estrutural direto no DB
# Ajuste de sensibilidade:
#   -c (coverage) e --cov-mode impactam clusterização
#   --min-seq-id não é o parâmetro aqui; Foldseek usa score/align internamente
foldseek easy-cluster "${DB}" "${CLU}" "${TMPFS}" -e 1e-3

# Exporta cluster membership (rep<TAB>member)
OUT="${FEAT_DIR}/fold/fold_clusters.tsv"
foldseek createtsv "${DB}" "${DB}" "${CLU}_cluster.tsv" "${OUT}"

echo "[OK] clusters: ${OUT}"
