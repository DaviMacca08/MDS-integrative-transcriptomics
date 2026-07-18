Cross-Platform Transcriptomic Analysis Framework for Myelodysplastic
Syndromes: Comparative RNA-seq and Microarray Workflows
================

![R](https://img.shields.io/badge/R-4.6.0-blue)
![Transcriptomics](https://img.shields.io/badge/Transcriptomics-RNA--seq%20%7C%20Microarray-purple)
![DESeq2](https://img.shields.io/badge/RNA--seq-DESeq2-darkgreen)
![limma](https://img.shields.io/badge/Microarray-limma-darkred)
![Status](https://img.shields.io/badge/status-finished-brightgreen)

------------------------------------------------------------------------

# 💼 Bioinformatics Service Demonstration

This repository presents a reproducible transcriptomic analysis
framework developed for the comparative investigation of
**Myelodysplastic Syndromes (MDS)** across independent public cohorts
generated using different transcriptomic technologies.

Rather than focusing on a single dataset, the project integrates
multiple independent studies to evaluate the reproducibility of
transcriptional alterations and biological pathway perturbations across
both **RNA-seq** and **microarray** platforms.

The workflow provides a standardized computational framework for
differential expression analysis and downstream functional
interpretation while employing statistical methodologies specifically
optimized for each experimental technology.

The framework is suitable for:

- systematic transcriptomic characterization of disease-associated
  molecular alterations
- comparative analyses across independent public cohorts
- cross-platform validation of differential expression signatures
- functional interpretation through complementary enrichment approaches
- reproducible computational analyses with transparent methodological
  reporting
- generation of publication-quality figures and analytical summaries

**Core analytical functionalities supported by the framework include:**

- comprehensive quality control and exploratory data analysis
- platform-specific preprocessing and normalization
- differential gene expression analysis using state-of-the-art
  statistical models
- functional enrichment through Over-Representation Analysis (ORA)
- Gene Set Enrichment Analysis (GSEA) using ranked transcriptome-wide
  statistics
- comparative interpretation of biological pathways across independent
  cohorts
- standardized visualization and reproducible reporting

**📤 Deliverables**

- Differential expression results
- Functional enrichment analyses
- Publication-quality visualizations
- Reproducible analytical reports
- Cross-dataset biological comparisons

------------------------------------------------------------------------

# 📊 Case Study

| Dataset   | Platform   | Healthy Controls | MDS Patients | Cell Population |
|-----------|------------|-----------------:|-------------:|-----------------|
| GSE19429  | Microarray |               17 |          183 | CD34+ cells     |
| GSE58831  | Microarray |               17 |          159 | CD34+ cells     |
| GSE114922 | RNA-seq    |                8 |           82 | CD34+ cells     |

------------------------------------------------------------------------

# 🧬 Biological Context

Myelodysplastic Syndromes comprise a heterogeneous group of
hematological malignancies characterized by ineffective hematopoiesis,
bone marrow dysplasia and an increased risk of progression to acute
myeloid leukemia.

------------------------------------------------------------------------

# 🎯 Project Objectives

- identify differentially expressed genes associated with MDS
- characterize dysregulated biological pathways through complementary
  enrichment approaches
- evaluate reproducibility across independent cohorts
- compare RNA-seq and microarray findings
- establish a standardized and reproducible workflow

------------------------------------------------------------------------

# ⚙️ Analytical Workflow

## Quality Control

- sample-level quality assessment
- exploratory data analysis
- PCA
- outlier detection

## Differential Expression

### Microarray

- limma
- empirical Bayes
- Benjamini-Hochberg correction

### RNA-seq

- DESeq2
- size-factor normalization
- negative binomial GLM
- Wald test
- Benjamini-Hochberg correction

## Functional Analysis

- Over-Representation Analysis (ORA)
- Gene Set Enrichment Analysis (GSEA)
- KEGG pathway visualization using Pathview

------------------------------------------------------------------------

# 🔬 Statistical Framework

| Analysis                 | Microarray         | RNA-seq            |
|--------------------------|--------------------|--------------------|
| Differential Expression  | limma              | DESeq2             |
| Multiple Testing         | Benjamini-Hochberg | Benjamini-Hochberg |
| ORA                      | ✓                  | ✓                  |
| GSEA                     | ✓                  | ✓                  |
| Cross-Dataset Comparison | ✓                  | ✓                  |

------------------------------------------------------------------------

# 🧪 Batch Effect Assessment (SVA)

Surrogate Variable Analysis was performed diagnostically on all three
cohorts to evaluate the presence of latent, unmodeled technical
variation prior to differential expression testing. Two estimation
strategies were compared, `be` and `leek`, to assess the stability of
the surrogate variable estimates independently of the differential
expression model.

Across datasets, the two methods did not converge on consistent
estimates. For GSE114922, the `leek` method estimated a large number of
surrogate variables relative to sample size, which destabilized the
downstream model. For GSE19429, correction with the `be` method visibly
reduced the proportion of variance attributable to biological signal,
along with a marked reduction in the number of differentially expressed
genes, consistent with overcorrection rather than genuine noise removal.

Given this instability, **SVA-adjusted surrogate variables were not
incorporated into the final differential expression models**.
`removeBatchEffect()` outputs were used exclusively for visualization
purposes (PCA, heatmaps) and never as input to statistical testing.
Baseline (unadjusted) models were retained as the definitive results for
each dataset.

As an indirect line of support for this decision, the strength of
biological replication observed across the three independent cohorts —
generated on different platforms, in different laboratories, at
different times — was taken as evidence that the unadjusted models
capture reproducible disease-associated signal rather than being
dominated by technical variation. This is presented as a qualitative,
cross-dataset argument rather than a formal quantitative diagnostic of
batch effect magnitude.

------------------------------------------------------------------------

# 🔗 Cross-Dataset Comparative Analysis

Differential expression and enrichment results were compared across the
three independent cohorts to assess reproducibility of the MDS
transcriptional signature across platforms. Datasets were analyzed
independently throughout the pipeline; integration was performed
exclusively at the level of summary statistics and enrichment results,
without merging raw or normalized expression matrices across platforms.

Two complementary comparison strategies were used:

- **Differential expression correlation**: pairwise correlation of
  gene-level statistics (e.g., log-fold-change or moderated
  t-statistics) across shared genes between datasets.
- **Enrichment concordance**: overlap and correlation of Normalized
  Enrichment Scores (NES) across datasets for shared gene sets (KEGG,
  MSigDB Hallmark), visualized as comparative heatmaps.

------------------------------------------------------------------------

# 📁 Project Structure

``` text
MDS-integrative-transcriptomics/
│
├── README.md
├── LICENSE
│
├── SetupEnvironment/
│   ├── 00_paths.R
│   ├── 01_environment.R
│   ├── 02_seed.R
│   └── 03_helper_functions.R
│
├── scripts/
│   │
│   ├── GSE114922_RNAseq/
│   │
│   ├── GSE58831_microarray/
│   │
│   └── GSE19429_microarray/
│
├── scripts_sva/
│   ├── Comparison_DEG_sva.R
│   ├── DEG_sva_GSE19429.R
│   ├── DEG_sva_GSE58831.R
│   ├── DEG_sva_GSE114922.R
│
├── results/
│   │
│   ├── GSE114922_RNAseq/
│   │
│   ├── GSE58831_microarray/
│   │
│   └── GSE19429_microarray/
│   │
│   └── Comparison/
│   │
│   └── logs/
│   │
│   └── tables/
│
├── results_sva/
│   │
│   ├── GSE114922_RNAseq/
│   │
│   ├── GSE58831_microarray/
│   │
│   └── GSE19429_microarray/
│   │
│   └── logs/
│   │
│   └── tables/
│
├── Report/
│   ├── Transcriptomic_analysis_report.md
│   ├── SVA_assessment_report.md
│
└── 
```

------------------------------------------------------------------------

# 📈 Analysis Outputs

The workflow generates standardized outputs for each dataset:

- Quality control and exploratory analysis reports
- Differential expression analysis results
- Functional enrichment analyses (ORA and GSEA)
- KEGG pathway visualization using Pathview
- Cross-platform comparison of transcriptional signatures
- Publication-quality figures and analytical summaries

------------------------------------------------------------------------

# 🔁 Reproducibility

The project was developed following reproducible bioinformatics workflow
principles:

- modular and dataset-specific R scripts
- centralized environment setup and dependency management
- standardized directory structure for inputs, scripts, and outputs
- reproducible randomization through controlled seed setting
- automated session information reporting
- transparent analytical workflow from preprocessing to biological
  interpretation

------------------------------------------------------------------------

# 🚀 Why This Project?

Unlike conventional analyses focused on a single dataset, this project
emphasizes biological reproducibility through independent validation
across multiple cohorts and transcriptomic platforms.

The framework demonstrates the ability to:

- integrate heterogeneous transcriptomic datasets
- apply technology-specific statistical methodologies
- compare biological findings across independent studies
- generate reproducible transcriptomic analyses
- translate computational results into biologically meaningful
  conclusions

------------------------------------------------------------------------

# 📚 References

## Public Transcriptomic Datasets

### GSE19429 — Microarray

- **Dataset:** Gene Expression Omnibus accession GSE19429
- **Platform:** Affymetrix Human Genome U133 Plus 2.0 Array
- **Biological material:** CD34+ bone marrow cells
- **Study reference:** Pellagatti A, Cazzola M, Giagounidis A, Perry J
  et al. Deregulated gene expression pathways in myelodysplastic
  syndrome hematopoietic stem cells. Leukemia 2010 Apr;24(4):756-64.
  PMID: 20220779

------------------------------------------------------------------------

### GSE58831 — Microarray

- **Dataset:** Gene Expression Omnibus accession GSE58831
- **Platform:** Affymetrix Human Genome U133 Plus 2.0 Array
- **Biological material:** CD34+ bone marrow cells
- **Study reference:** Gerstung M, Pellagatti A, Malcovati L,
  Giagounidis A et al. Combining gene mutation with gene expression data
  improves outcome prediction in myelodysplastic syndromes. Nat Commun
  2015 Jan 9;6:5901. PMID: 25574665

------------------------------------------------------------------------

### GSE114922 — RNA-seq

- **Dataset:** Gene Expression Omnibus accession GSE114922
- **Platform:** Illumina HiSeq 4000 (Homo sapiens)
- **Biological material:** CD34+ bone marrow cells
- **Study reference:** Pellagatti A, Armstrong RN, Steeples V, Sharma E
  et al. Impact of spliceosome mutations on RNA splicing in
  myelodysplasia: dysregulated genes/pathways and clinical associations.
  Blood 2018 Sep 20;132(12):1225-1240. PMID: 29930011

------------------------------------------------------------------------

# Tools

- R (v.4.6.0)
- Bioconductor / CRAN packages: GEOquery (v2.80.0), hgu133plus2.db
  (v3.13.0), limma (v3.68.4), DESeq2 (v1.52.0), clusterProfiler
  (v4.20.0), ReactomePA (v1.56.0), org.Hs.eg.db (v3.23.1), msigdbr
  (v26.1.0), fgsea (v1.38.0), pathview (v1.52.0), ComplexHeatmap
  (v2.28.0), sva (v3.60.0)

------------------------------------------------------------------------

# 📬 Contact

For questions, collaborations, or bioinformatics consulting inquiries:

**Davide Maccarrone**

- GitHub: <https://github.com/DaviMacca08>
- LinkedIn: www.linkedin.com/in/davidemaccarrone
- Email: <davide_maccarrone@icloud.com>
