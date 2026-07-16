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

PDB_DIR="${TMP_DIR}/pdb"
WORK_DIR="${TMP_DIR}/foldseek_self"
OUT_DIR="${WORK_DIR}/out"
TMP_FS="${WORK_DIR}/tmp"

MIN_SEQ_ID="$(pyget fold.min_seq_id)"
COV="$(pyget fold.cov)"
COV_MODE="$(pyget fold.cov_mode)"
CLUSTER_MODE="$(pyget fold.cluster_mode)"
THREADS="$(pyget fold.threads)"

mkdir -p "${WORK_DIR}" "${OUT_DIR}" "${TMP_FS}" "$(dirname "${FEAT_DIR}/fold/")" "${FEAT_DIR}/fold"

# Sanity: tem pdb?
N_PDB="$(ls -1 "${PDB_DIR}"/*.pdb 2>/dev/null | wc -l || true)"
if [[ "${N_PDB}" -eq 0 ]]; then
  echo "[ERR] No PDBs found in: ${PDB_DIR}"
  exit 1
fi
echo "[OK] Found ${N_PDB} PDBs in ${PDB_DIR}"

# Jeito recomendado: passa a pasta diretamente (evita erro "No structures found")
echo "[RUN] foldseek easy-cluster ${PDB_DIR} ${OUT_DIR}/clusters ${TMP_FS} ..."
foldseek easy-cluster \
  "${PDB_DIR}" \
  "${OUT_DIR}/clusters" \
  "${TMP_FS}" \
  --min-seq-id "${MIN_SEQ_ID}" \
  -c "${COV}" \
  --cov-mode "${COV_MODE}" \
  --cluster-mode "${CLUSTER_MODE}" \
  --threads "${THREADS}"

echo "[OK] Foldseek clustering done."
echo "[INFO] Output prefix: ${OUT_DIR}/clusters"
echo "[INFO] Next: python scripts/05b_foldseek_cluster_summarize.py --config ${CFG}"
