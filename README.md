# TCGA Lung Squamous Cell Carcinoma (LUSC) RNA-seq Biomarker Discovery

A fully reproducible bioinformatics pipeline for identifying candidate biomarkers in lung cancer using TCGA bulk RNA-seq data.

---

## Analysis Workflow

![Workflow](figures/workflow.png)

This study can be reproduced from raw TCGA data using a single R command (`source("run_analysis.R")`).
The pipeline performs data download, preprocessing, statistical analysis, visualization, functional interpretation, and clinical survival analysis.

---

## Project Overview

Lung cancer is one of the leading causes of cancer mortality worldwide.  
**Lung Squamous Cell Carcinoma (LUSC)** is a major histological subtype characterized by strong transcriptional dysregulation, epithelial differentiation, and tumor–immune interactions.

This project performs a **reproducible transcriptomic analysis of TCGA-LUSC RNA-seq data** to:

- identify genes differentially expressed between tumor and normal lung tissue
- determine affected biological pathways (Gene Ontology)
- evaluate clinical relevance of selected genes (Kaplan–Meier survival)
- propose candidate diagnostic and prognostic biomarkers

The analysis compares **primary tumor samples (TP)** with **normal tissue samples (NT)**.

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

Principal Component Analysis shows clear separation between tumor and normal lung samples, indicating widespread transcriptomic reprogramming in LUSC.

---

### Differential Gene Expression

![Volcano](figures/volcano_lusc.png)

Thousands of genes are significantly differentially expressed (FDR-corrected), consistent with major molecular alterations associated with tumor development.

---

### Survival Analysis (Kaplan–Meier)

![Survival](figures/survival_KRT6A.png)

High expression of **KRT6A** stratifies TCGA-LUSC patients into different survival groups (log-rank p-value shown), suggesting potential prognostic relevance.

---

## Functional Interpretation

Gene Ontology enrichment revealed biological programs characteristic of squamous tumors:

- keratinization
- epidermis development
- epithelial cell differentiation
- immune activation pathways

These findings align with the known histopathology of squamous lung carcinoma and indicate strong epithelial and immune involvement in tumor biology.

---

## Quality Control

Detailed QC report: `reports/qc_report.md`

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
results/ → differential expression tables (optional / can be ignored if generated)
reports/ → QC report
run_analysis.R → one-command pipeline runner


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

To identify candidate biomarkers, the most significantly differentially expressed genes between tumor and normal lung tissue were selected and visualized as a heatmap.

![Tumor Gene Signature](figures/heatmap_lusc_signature.png)

The gene expression signature clearly separates tumor and normal samples, indicating a robust cancer-specific transcriptional program.

Many of the top genes are associated with epithelial differentiation and keratinization — hallmark biological processes of squamous cell carcinoma.

This gene set represents a potential **diagnostic and prognostic biomarker panel** for LUSC and may serve as a starting point for future experimental validation.


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
