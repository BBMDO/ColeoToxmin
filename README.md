# ColeoTox-StructAtlas

**Structure-guided discovery and prioritization of toxin-like peptide architectures across Coleoptera**

This repository contains the computational workflow, processed data, predicted structures, machine-learning outputs, and figure source files associated with the study:

> **Structure-guided discovery reveals recurrent bioactive peptide architectures across Coleoptera**  
> Thaís Caroline Gonçalves, João Alfredo Teodoro, and Danilo T. Amaral

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
├── 04_results/             # Final tables, figures, candidate catalogues, and revised mature-peptide analyses
│   ├── figures/
│   ├── tables/
│   └── analysis_mature/    # Revision-specific reanalysis of 76 high-confidence mature peptides
│       ├── 01_models_cif/
│       ├── 01_models_pdb/
│       ├── 02_foldseek/
│       ├── 03_blast/
│       ├── 04_GO/
│       ├── 05_foldseek/
│       ├── 06_functional_clusters/
│       ├── 07_novelty_ranking/
│       ├── 08_external_structure_search/
│       ├── 09_external_structure_interpretation/
│       ├── 10_confirmatory_structure_search/
│       ├── 11_final_structural_atlas/
│       ├── 12_final_figures/
│       ├── 13_structure_gallery/
│       ├── 14_structure_gallery_publication/
│       ├── 15_supplementary_tables/
│       ├── figures/
│       ├── inputs/
│       ├── logs/
│       └── tables/
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

## Revised mature-peptide analysis

The complete revision-specific reanalysis requested during peer review is located in:

```text
04_results/analysis_mature/
```

The initial discovery workflow identified **291 candidate peptide precursors**. Following maturation assessment, a high-confidence subset of **76 inferred mature peptides** was retained for revised physicochemical and structural analyses.

For this reanalysis:

1. predicted signal peptide regions were removed;
2. propeptide regions were removed when applicable;
3. mature peptide sequences were inferred;
4. sequence-derived and physicochemical descriptors were recalculated from mature sequences;
5. AlphaFold 3 models were independently generated from mature peptide sequences;
6. structural descriptors were recalculated;
7. Foldseek clustering was repeated using mature-peptide structures;
8. BLASTP searches were performed against UniProtKB/Swiss-Prot;
9. Gene Ontology annotations were evaluated for representative structural clusters;
10. representative structures were compared against PDB and AlphaFold/Swiss-Prot, including sensitive confirmatory searches for initially unresolved clusters.

### Mature peptide sequences

The FASTA file containing the 76 mature peptide sequences used in the revised analyses is:

```text
04_results/analysis_mature/inputs/mature_peptides_76_for_blast.fasta
```

Sequence-derived features are available in:

```text
04_results/analysis_mature/tables/mature_sequence_features.tsv
04_results/analysis_mature/tables/mature_all_features.tsv
```

All descriptors in these files refer to the inferred mature peptide rather than the full precursor.

### AlphaFold 3 mature-peptide models

Models reconstructed from the 76 mature peptide sequences are available in:

```text
04_results/analysis_mature/01_models_cif/
04_results/analysis_mature/01_models_pdb/
```

The repository contains **76 CIF models** and **76 PDB models**.

Model provenance and confidence information are summarized in:

```text
04_results/analysis_mature/tables/model_manifest.tsv
04_results/analysis_mature/tables/mature_structural_features.tsv
```

These are computationally predicted structures and should not be interpreted as experimentally determined structures.

### Foldseek clustering

Foldseek structural comparisons and clustering of the mature-peptide models are available in:

```text
04_results/analysis_mature/05_foldseek/
```

Key files include:

```text
cluster_membership.tsv
cluster_summary.tsv
all_vs_all.tsv
clusters_rep_seq.fasta
clusters_all_seqs.fasta
```

The 76 mature peptides were organized into **26 Foldseek-derived structural clusters**.

### BLASTP against UniProtKB/Swiss-Prot

Sequence-similarity analyses are available in:

```text
04_results/analysis_mature/03_blast/
```

Key files include:

```text
mature_vs_swissprot_all_hits.tsv
mature_vs_swissprot_all_hits_parsed.tsv
mature_vs_swissprot_best_hits.tsv
blast_similarity_summary.tsv
```

These outputs provide, when available, UniProt accession numbers, percentage identity, query coverage, E-values, bit scores, and protein descriptions.

BLASTP similarity is interpreted as supporting sequence evidence rather than experimental functional validation.

### Gene Ontology and functional interpretation

GO-supported annotations are available in:

```text
04_results/analysis_mature/04_GO/
04_results/analysis_mature/06_functional_clusters/
```

Important files include:

```text
04_results/analysis_mature/04_GO/candidate_GO.tsv
04_results/analysis_mature/04_GO/uniprot_GO.tsv
04_results/analysis_mature/06_functional_clusters/cluster_representatives_GO.tsv
04_results/analysis_mature/06_functional_clusters/cluster_representatives_GO_long.tsv
04_results/analysis_mature/06_functional_clusters/cluster_representatives_annotated.tsv
04_results/analysis_mature/06_functional_clusters/cluster_functional_summary.tsv
```

GO terms provide functional context for representative structural clusters when supported by sequence-based annotation. They should not be interpreted as direct evidence of biological activity.

### External structural searches

Initial structural searches against PDB and AlphaFold/Swiss-Prot are available in:

```text
04_results/analysis_mature/08_external_structure_search/
```

Main outputs:

```text
representatives_vs_PDB.tsv
representatives_vs_AlphaFold_SwissProt.tsv
```

Interpretation of these searches is provided in:

```text
04_results/analysis_mature/09_external_structure_interpretation/
```

Sensitive confirmatory searches for prioritized unresolved candidates are available in:

```text
04_results/analysis_mature/10_confirmatory_structure_search/
```

Important final outputs include:

```text
priority3_final_structural_assessment.tsv
priority3_final_structural_conclusions.tsv
priority3_sensitive_all_hits_classified.tsv
priority3_sensitive_best_hits_final.tsv
```

Some initially unresolved candidates gained close structural counterparts in the sensitive searches. **MC08 remained without a detected close structural counterpart under the search conditions used against PDB and AlphaFold/Swiss-Prot.**

This result is reported conservatively as **“no close structural counterpart detected”** and not as proof of a definitively novel fold.

### Final structural atlas and revised figures

The integrated mature-peptide structural atlas is available in:

```text
04_results/analysis_mature/11_final_structural_atlas/
```

Publication-oriented figures and galleries are available in:

```text
04_results/analysis_mature/12_final_figures/
04_results/analysis_mature/13_structure_gallery/
04_results/analysis_mature/14_structure_gallery_publication/
```

Publication-ready supplementary tables generated during the reanalysis are available in:

```text
04_results/analysis_mature/15_supplementary_tables/
```

These tables connect candidate identifiers, species, mature peptide sequences, physicochemical descriptors, AlphaFold 3 confidence metrics, Foldseek clusters, BLAST/Swiss-Prot evidence, Gene Ontology annotations, external structural matches, and final structural interpretation.

### Interpretation of structural novelty

Structural-cluster membership alone was not treated as evidence of novelty. Low-representation Foldseek clusters were further compared against external structural databases to distinguish close, moderate, and distant structural relationships from cases in which no close counterpart was detected.

Accordingly, the revised manuscript uses conservative terminology such as:

> candidate structural architecture lacking a close structural counterpart under the search conditions employed.


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

The original transcriptomic datasets were obtained from public repositories, primarily the NCBI Transcriptome Shotgun Assembly and Sequence Read Archive databases. Species, TSA/SRA accession numbers, and tissue, organ, body region, or developmental-stage metadata are reported in Supplementary Table S1 whenever such information was available from the original TSA/SRA/BioSample records. When the original public metadata were unavailable or insufficiently specific, this information is reported as unavailable rather than inferred.

Raw public sequencing reads and third-party reference databases are not redistributed in this repository. Users should obtain those resources from their original repositories according to the relevant database terms.

---

## Data interpretation

The reported candidates are computational predictions. Structural similarity, secretion signals, cysteine patterns, PU-learning scores, and structural-cluster membership support candidate prioritization but do not constitute experimental evidence of toxicity or biological activity.

Functional assignments should therefore be treated as hypotheses requiring proteomic, biochemical, and biological validation.

---

## Important interpretation note

The sequences described in this study are **computationally inferred peptide candidates derived from publicly available transcriptomes**. They should not be interpreted as experimentally confirmed proteins or experimentally validated toxins.

Likewise:

- AlphaFold 3 structures are computational predictions;
- structural similarity does not establish biological activity;
- BLAST similarity does not establish functional identity;
- GO annotations represent computational support rather than experimental validation;
- PU-learning scores represent candidate prioritization rather than proof of toxicity.

Experimental proteomic, biochemical, and functional validation will be required to establish peptide production, maturation, structure, and biological activity.

---

## Data availability and Zenodo archive

Processed datasets, mature peptide sequences, predicted structural models, analysis scripts, structural-clustering results, BLAST similarity results, functional annotations, machine-learning inputs and outputs, and figure source data supporting this study are available through this GitHub repository.

A permanent archived version is available through Zenodo:

> **Zenodo DOI:** `10.5281/zenodo.21386392`

Raw transcriptomic datasets were obtained from public repositories and are identified by accession numbers in Supplementary Table S1.

---

## Citation

Please cite the associated article when using this repository or its data:

> Gonçalves, T. C.; Teodoro, J. A.; Amaral, D. T. *Structure-guided discovery reveals recurrent bioactive peptide architectures across Coleoptera*.

A machine-readable `CITATION.cff` file should be added after the article DOI and bibliographic details become available.

---

## Authors

- **Thaís Caroline Gonçalves**  
  ORCID: 0009-0001-3039-7366

- **João Alfredo Teodoro**  
  ORCID: 0009-0003-5601-0814

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
