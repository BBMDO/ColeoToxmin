# ColeoTox_min pipeline notes (fixed)

Canonical outputs:
- `04_merge_final.py` writes **feature-only** tables:
  - `04_results/tables/all_candidates_features.tsv`
  - `04_results/tables/top50_candidates_features.tsv`

PU stage:
- `06b_train_pu.py` reads `all_candidates_features.tsv` and writes:
  - `04_results/tables/candidates_pu_scored.tsv`
  - `04_results/tables/pu_metrics.tsv`

Final publishable tables:
- `06c_make_final_table.py` merges features + PU scores and writes:
  - `04_results/tables/all_candidates_final.tsv`
  - `04_results/tables/top50_candidates_final.tsv`

Why:
- Previously, `06c_make_final_table.py` overwrote `all_candidates_final.tsv`, deleting DSSP columns.
- This fix makes the file naming unambiguous and preserves structural features.
