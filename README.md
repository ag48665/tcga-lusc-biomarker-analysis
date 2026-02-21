# TCGA Lung Squamous Cell Carcinoma (LUSC) RNA-seq Biomarker Discovery

A fully reproducible bioinformatics pipeline for identifying candidate biomarkers in lung cancer using TCGA bulk RNA-seq data (tumor vs normal).

---

## Analysis Workflow

![Workflow](figures/workflow.png)

This repository is designed so the entire analysis can be reproduced from raw TCGA data using a single R command.

---

## Project Overview

Lung cancer is one of the leading causes of cancer-related mortality worldwide.  
**Lung Squamous Cell Carcinoma (LUSC)** is a major histological subtype characterized by strong transcriptional dysregulation and tumor–immune interactions.

This project performs a **reproducible transcriptomic analysis of TCGA-LUSC RNA-seq data** to:

- identify genes differentially expressed between tumor and normal lung tissue
- interpret affected biological pathways (Gene Ontology)
- evaluate clinical relevance using survival analysis
- propose candidate diagnostic / prognostic biomarkers

The analysis compares **primary tumor samples (TP)** vs **normal tissue samples (NT)**.

---

## Dataset

Source: **The Cancer Genome Atlas (TCGA)**

| Attribute | Value |
|---|---|
| Cohort | TCGA-LUSC |
| Data type | Bulk RNA-seq gene expression counts |
| Tumor samples | 511 |
| Normal samples | 51 |
| Genes analyzed | ~60,000 |

Data are downloaded programmatically using the `TCGAbiolinks` R package.

---

## Key Results

### Tumor vs Normal Separation (PCA)

![PCA](figures/pca_lusc_tp_vs_nt.png)

Principal Component Analysis shows clear separation between tumor and normal lung samples, indicating widespread transcriptomic reprogramming in LUSC.

---

### Differential Gene Expression

![Volcano](figures/volcano_lusc.png)

Thousands of genes are significantly differentially expressed (FDR-corrected), consistent with major molecular alterations in tumor tissue.

---

### Survival Analysis (example gene)

![Survival](figures/survival_KRT6A.png)

High expression of **KRT6A** stratifies TCGA-LUSC patients into different survival groups (Kaplan–Meier; log-rank p-value shown), suggesting potential prognostic relevance.

---

## Tumor Gene Signature

The most variable significantly differentially expressed genes were selected and visualized across tumor (TP) and normal (NT) samples (z-score scaled per gene).

![Tumor Gene Signature Heatmap](figures/heatmap_lusc_signature.png)

---

## Functional Interpretation

Gene Ontology enrichment highlights biological programs characteristic of squamous tumors:

- keratinization
- epidermis development
- epithelial cell differentiation
- immune activation pathways

These results align with known LUSC biology and tumor microenvironment involvement.

---

## Quality Control

Detailed QC report: `reports/qc_report.md`

QC plots:

**QC PCA**  
![](figures/qc_pca.png)

**Sample correlation**  
![](figures/qc_sample_correlation.png)

**Library sizes**  
![](figures/qc_library_sizes.png)

---

## Repository Structure

scripts/   -> analysis scripts (R)
figures/   -> generated plots
results/   -> result tables (e.g., DE genes)
reports/   -> QC report (markdown)
run_analysis.R -> one-command full pipeline


Main result file (if saved): `results/DE_genes_LUSC.csv` (FDR-corrected DE results)

---

## Reproducibility

Run the entire pipeline from the project root:
source("run_analysis.R")


The pipeline automatically:

downloads TCGA-LUSC RNA-seq data

performs quality control

computes PCA

identifies differentially expressed genes

generates figures

performs Gene Ontology enrichment

runs survival analysis

Software versions are recorded in session_info.txt.

---

## Tumor Gene Signature

The most variable significantly differentially expressed genes were selected and visualized across tumor (TP) and normal (NT) samples (z-score scaled per gene).

![Tumor Gene Signature](figures/heatmap_lusc_signature.png) 

This heatmap highlights a clear transcriptomic signature distinguishing LUSC tumors from normal lung tissue and supports the biological validity of the candidate biomarkers.
---

## Tumor Gene Signature

The most variable significantly differentially expressed genes were selected and visualized across tumor (TP) and normal (NT) samples (z-score scaled per gene).

![Tumor Gene Signature](figures/heatmap_lusc_signature.png)

This heatmap highlights a clear transcriptomic signature distinguishing LUSC tumors from normal lung tissue and supports the biological validity of the candidate biomarkers.


## Skills Demonstrated

RNA-seq data analysis

statistical modeling (DESeq2)

cancer transcriptomics

functional enrichment analysis

Kaplan–Meier survival analysis

data visualization in R

reproducible research workflow design

automated bioinformatics pipelines


## Author

Agata Gabara
Bioinformatics / Computational Biology
