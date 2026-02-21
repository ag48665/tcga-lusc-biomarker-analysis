# TCGA Lung Squamous Cell Carcinoma (LUSC) RNA-seq Biomarker Discovery

A fully reproducible bioinformatics pipeline for identifying candidate biomarkers in lung cancer using TCGA bulk RNA-seq data.

---

## Analysis Workflow

![Workflow](figures/workflow.png)

The entire study can be reproduced from raw TCGA data using a single R command.  
The pipeline automatically performs data download, preprocessing, statistical analysis, visualization, functional interpretation, and clinical survival analysis.

---

## Project Overview

Lung cancer is one of the leading causes of cancer mortality worldwide.  
**Lung Squamous Cell Carcinoma (LUSC)** is a major histological subtype characterized by strong transcriptional dysregulation, epithelial differentiation, and tumor–immune interactions.

The goal of this project is to perform a **reproducible transcriptomic analysis of TCGA-LUSC RNA-seq data** in order to:

- identify genes differentially expressed between tumor and normal lung tissue
- determine affected biological pathways
- evaluate clinical relevance of selected genes
- propose candidate diagnostic and prognostic biomarkers

The analysis compares **primary tumor samples (TP)** with **normal tissue samples (NT)** and links gene expression to biological function and patient survival.

---

## Dataset

Source: **The Cancer Genome Atlas (TCGA)**

| Attribute      | Value                               |
|-------------- |------------------------------------|
| Cohort         | TCGA-LUSC                           |
| Data type      | Bulk RNA-seq gene expression counts |
| Tumor samples  | 511                                 |
| Normal samples | 51                                  |
| Genes analyzed | ~60,000                             |

Data were downloaded programmatically using the `TCGAbiolinks` R package.

---

## Key Results

### Tumor vs Normal Separation (PCA)

![PCA](figures/pca_lusc_tp_vs_nt.png)

Principal Component Analysis shows clear separation between tumor and normal lung samples, demonstrating widespread transcriptomic reprogramming in LUSC.

---

### Differential Gene Expression

![Volcano](figures/volcano_lusc.png)

Thousands of genes are significantly differentially expressed (FDR-corrected), indicating major molecular alterations associated with tumor development.

---

### Survival Analysis

![Survival](figures/survival_KRT6A.png)

High expression of **KRT6A** stratifies TCGA-LUSC patients into different survival groups (Kaplan–Meier analysis), suggesting clinical relevance and potential prognostic value.

---

## Functional Interpretation

Gene Ontology enrichment revealed biological programs characteristic of squamous tumors:

- keratinization
- epidermis development
- epithelial cell differentiation
- immune activation pathways

These findings are consistent with the known histopathology of squamous lung carcinoma and indicate strong epithelial and immune involvement in tumor biology.

---

## Quality Control

Detailed QC report:  
`reports/qc_report.md`

Quality checks performed:

- library size distribution
- PCA clustering
- sample–sample correlation

**QC PCA**
![](figures/qc_pca.png)

**Sample correlation**
![](figures/qc_sample_correlation.png)

**Library sizes**
![](figures/qc_library_sizes.png)

---

## Repository Structure
scripts/ → analysis scripts (R)
figures/ → generated plots
results/ → differential expression tables
reports/ → QC report

---
Main result file:
results/DE_genes_LUSC.csv
---

(FDR-corrected differential expression results)

---
## Reproducibility

---


## Skills Demonstrated

This project demonstrates:

- RNA-seq data analysis
- statistical modeling (DESeq2)
- cancer transcriptomics
- functional enrichment analysis
- Kaplan–Meier survival analysis
- data visualization in R
- reproducible research workflow design
- automated bioinformatics pipelines

---

## Author

**Agata Gabara**  
Bioinformatics / Computational Biology

