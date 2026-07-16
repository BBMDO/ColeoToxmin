# ColeoTox-StructAtlas

**Structure-guided discovery and prioritization of toxin-like peptide architectures across Coleoptera**

This repository contains the computational workflow, processed data, predicted structures, machine-learning outputs, and figure source files associated with the study:

> **Structure-guided discovery reveals recurrent bioactive peptide architectures across Coleoptera**  
> Thaís Caroline Gonçalves and Danilo T. Amaral

The study integrates transcriptome mining, physicochemical characterization, AlphaFold-based structural modeling, Foldseek structural clustering, Pfam annotation, and positive–unlabeled (PU) learning to identify and prioritize candidate toxin-like peptides across Coleoptera.

---

## Repository contents

```text
.
├── 01_intermediate/        # Intermediate tables and structural-analysis outputs
│   ├── dssp/               # DSSP secondary-structure assignments
│   ├── foldseek_self/      # Foldseek self-search and clustering outputs
│   ├── hmm/                # HMM/Pfam search outputs
│   └── pu/                 # PU-learning input and split tables
├── 01_work/                # Working structural files used during analysis
├── 02_features/            # Sequence, structural, HMM, and Foldseek-derived features
│   ├── fold/
│   ├── hmm/
│   └── structure/
├── 03_models/              # Trained PU-learning model objects
├── 04_results/             # Final tables, figures, and candidate catalogues
│   ├── figures/
│   └── tables/
├── 09_reproducibility/     # Example configuration files and pipeline launcher
└── scripts/                # Analysis and visualization scripts
```

---

## Main outputs

The principal manuscript-ready outputs are located in `04_results/`.

### Candidate catalogues

- `04_results/ColeoTox_min_catalog.tsv`  
  Compact catalogue of candidate peptides.

- `04_results/tables/all_candidates_final.tsv`  
  Complete candidate table containing sequence-derived, structural, Foldseek, Pfam, and PU-learning information.

- `04_results/tables/top50_candidates_final.tsv`  
  Highest-ranked candidates according to the final prioritization workflow.

- `04_results/tables/top20_stability_across_seeds.tsv`  
  Stability of top-ranked candidates across independent random seeds.

### Structural clustering

- `02_features/fold/candidates_fold_cluster.tsv`  
  Assignment of each candidate to a Foldseek-derived structural cluster.

- `02_features/fold/fold_cluster_summary.tsv`  
  Summary statistics for Foldseek-derived structural clusters.

- `04_results/figures/Figure_S1_Foldseek_clusters.pdf`  
  Supplementary summary of recurrent structural clusters.

### Machine learning

- `04_results/tables/candidates_pu_scored_ensemble.tsv`  
  Ensemble PU-learning predictions.

- `04_results/tables/pu_metrics.tsv`  
  Predictive-performance metrics.

- `04_results/tables/colex_shap_proxy.shap_importance_meanabs.tsv`  
  Mean absolute SHAP feature importance.

- `03_models/pu/pu_model.joblib`  
  Serialized trained model. Compatibility depends on the Python and scikit-learn versions used during training.

### Predicted structures

Predicted structural models are provided as PDB files in the working and structural-output directories. These models were generated computationally and should be interpreted as structural hypotheses rather than experimentally determined structures.

---

## Analysis overview

The workflow comprises the following major steps:

1. Candidate peptide identification and sequence filtering.
2. Prediction of secretion signals and peptide maturation features.
3. Physicochemical characterization.
4. AlphaFold-based structural modeling.
5. Structural feature extraction using DSSP, solvent-accessibility calculations, contact analysis, radius of gyration, and disulfide-bond prediction.
6. Foldseek all-versus-all structural comparison and clustering.
7. Pfam/HMM-based annotation.
8. Positive–unlabeled learning and ensemble candidate prioritization.
9. Stability analysis across independent random seeds.
10. Generation of final catalogues, figures, and structural representatives.

---

## Reproducibility

### Software requirements

The workflow uses the following main software packages and tools:

- Python 3
- pandas
- NumPy
- scikit-learn
- joblib
- matplotlib
- PyYAML
- Biopython
- Foldseek
- HMMER
- DSSP / `mkdssp`
- FreeSASA
- MDTraj

Exact package versions should be recorded in a release-specific environment file before archival publication.

### Configuration

Example configuration files are provided in:

```text
09_reproducibility/configs/
```

The archived configuration files may contain machine-specific absolute paths. Before running the workflow, replace them with local paths or create a portable configuration file such as:

```yaml
paths:
  alphafold_dir: "00_inputs/AlphaFold"
  toxin_dir: "00_inputs/toxin"
  tmp_dir: "01_intermediate"
  features_dir: "02_features"
  results_dir: "04_results"
```

### Example execution

From the repository root:

```bash
bash 09_reproducibility/run_C.sh 09_reproducibility/configs/config.yaml
```

Individual analysis stages can also be executed directly from the `scripts/` directory. See `scripts/PIPELINE_NOTES.md` for notes on canonical outputs and table generation.

---

## Input data

The original transcriptomic datasets were obtained from public repositories, primarily the NCBI Transcriptome Shotgun Assembly and Sequence Read Archive databases. Accession information and dataset provenance are reported in the manuscript supplementary material.

Raw public sequencing reads and third-party reference databases are not redistributed in this repository. Users should obtain those resources from their original repositories according to the relevant database terms.

---

## Data interpretation

The reported candidates are computational predictions. Structural similarity, secretion signals, cysteine patterns, PU-learning scores, and structural-cluster membership support candidate prioritization but do not constitute experimental evidence of toxicity or biological activity.

Functional assignments should therefore be treated as hypotheses requiring proteomic, biochemical, and biological validation.

---

## Zenodo archive

A versioned archival snapshot containing the complete processed dataset, predicted structures, final figures, tables, and reproducibility files should be deposited on Zenodo.

After deposition, replace the placeholder below with the final DOI:

> **Zenodo DOI:** `10.5281/zenodo.XXXXXXXX`

For the manuscript, the recommended data-availability statement is:

> The processed datasets, predicted structural models, analysis scripts, and source data supporting the findings of this study are available in the ColeoTox-StructAtlas GitHub repository and in the associated Zenodo record at [DOI].

---

## Citation

Please cite the associated article when using this repository or its data:

> Gonçalves, T. C.; Amaral, D. T. *Structure-guided discovery reveals recurrent bioactive peptide architectures across Coleoptera*. Manuscript in preparation.

A machine-readable `CITATION.cff` file should be added after the article DOI and bibliographic details become available.

---

## Authors

- **Thaís Caroline Gonçalves**  
  ORCID: 0009-0001-3039-7366

- **Danilo T. Amaral**  
  ORCID: 0000-0002-8940-6546  
  Universidade Federal do ABC, Brazil

---

## License

Before public release, add explicit licenses for both code and data.

Recommended choices:

- **Code:** MIT License
- **Processed data, tables, figures, and predicted structures:** Creative Commons Attribution 4.0 International (CC BY 4.0)

Third-party software, public database records, and external reference resources remain subject to their original licenses and terms of use.

---

## Contact

For questions regarding the workflow or dataset, contact:

**Danilo T. Amaral**  
Centro de Ciências Naturais e Humanas  
Universidade Federal do ABC  
Santo André, São Paulo, Brazil  
Email: danilo.trabuco@ufabc.edu.br
