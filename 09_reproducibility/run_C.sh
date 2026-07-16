#!/usr/bin/env bash
set -euo pipefail
CFG="${1:-09_reproducibility/configs/config.yaml}"

python scripts/03a_cif_to_pdb.py --config "$CFG"
python scripts/03b_plddt_extract.py --config "$CFG"
python scripts/03c_run_dssp.py --config "$CFG"
python scripts/03d_sasa_rsa.py --config "$CFG"
python scripts/03e_contacts_rg.py --config "$CFG"
python scripts/03f_disulfides.py --config "$CFG"
python scripts/03h_merge_struct_features.py --config "$CFG"
