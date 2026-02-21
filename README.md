# TCGA Lung Squamous Cell Carcinoma (LUSC) RNA-seq Biomarker Discovery

This project implements a fully reproducible bioinformatics pipeline for analysis of TCGA lung squamous cell carcinoma (LUSC) RNA-seq data.

The workflow identifies tumor-specific transcriptional programs, performs functional interpretation, and evaluates the prognostic relevance of a multi-gene expression signature using Kaplan–Meier survival analysis and multivariate Cox regression.

The analysis demonstrates that a transcriptome-derived gene signature is associated with overall survival in TCGA patients.

---

## Analysis Workflow

![Workflow](figures/workflow.png)

This repository is designed so the entire analysis can be reproduced from raw TCGA data using a single R command.

---

## Project Overview

Lung cancer is one of the leading causes of cancer-related mortality worldwide.  
**Lung Squamous Cell Carcinoma (LUSC)** is a major histological subtype characterized by strong transcriptional dysregulation and tumor–immune interactions.

This project performs a **reproducible transcriptomic analysis of TCGA-LUSC RNA-seq data** to:

- identify differentially expressed genes between tumor and normal tissue
- characterize biological pathways in LUSC
- construct a transcriptome-derived prognostic gene signature
- evaluate clinical relevance using multivariate survival modeling

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
## Repository Structure

scripts/   -> analysis scripts (R)
figures/   -> generated plots
results/   -> result tables (e.g., DE genes)
reports/   -> QC report (markdown)
run_analysis.R -> one-command full pipeline


Main result file (if saved): `results/DE_genes_LUSC.csv` (FDR-corrected DE results)

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
## Pathway Analysis (KEGG)

![GO enrichment](figures/GO_BP_simplified.png)

KEGG pathway enrichment analysis highlighted immune- and microenvironment-related pathways (e.g., cytokine–cytokine receptor interaction, complement/coagulation, hematopoietic cell lineage, neutrophil extracellular trap formation) as well as extracellular matrix and adhesion pathways (ECM–receptor interaction, CAMs). Enrichment of cornified envelope formation supports the squamous epithelial differentiation program characteristic of LUSC.

Full KEGG results: `results/KEGG_enrichment_results.csv`

---

## Functional Pathway Analysis (Gene Ontology)

![GO enrichment](figures/GO_BP_simplified.png)

To interpret the biological meaning of the differentially expressed genes, Gene Ontology (GO) enrichment analysis was performed using the `clusterProfiler` package.

The tumor gene set showed strong enrichment for epithelial and squamous differentiation processes, including:

* epidermis development
* keratinization
* keratinocyte differentiation
* epithelial cell differentiation

These pathways are characteristic of lung squamous cell carcinoma histopathology and confirm that the transcriptomic analysis recovers the expected lineage program of airway basal epithelial cells.

In addition, immune-related pathways were significantly enriched:

* leukocyte-mediated immunity
* chemotaxis
* humoral immune response

This suggests active tumor–immune microenvironment interactions in TCGA-LUSC tumors.

Developmental programs such as cell fate commitment and pattern specification were also observed, consistent with aberrant differentiation and tumor progression.

Together, these results demonstrate that the identified gene signature captures both tumor cell identity and immune microenvironment activity.

The full enrichment results are provided in:

`results/GO_enrichment_results.csv`



---

### Survival Analysis (example gene)

![Survival](figures/survival_KRT6A.png)

High expression of **KRT6A** stratifies TCGA-LUSC patients into different survival groups (Kaplan–Meier; log-rank p-value shown), suggesting potential prognostic relevance.

---

### Prognostic Gene Signature

![Signature KM](https://github.com/ag48665/tcga-lusc-biomarker-analysis/blob/main/figures/signature_KM.png?raw=true)


A multi-gene transcriptional signature derived from tumor vs normal differential expression significantly stratifies TCGA-LUSC patients into survival groups (Cox HR ≈ 0.75, p < 0.05).

---

## Tumor Gene Signature
Differential expression analysis of TCGA-LUSC RNA-seq data identified a reproducible transcriptional program distinguishing tumor from normal lung tissue.

A large proportion of the top upregulated genes belonged to the Cancer-Testis Antigen (CTA) family, including MAGEA3, MAGEA4, MAGEA9, MAGEA10, MAGEA11, MAGEB1, and MAGEB2.
CTAs are normally restricted to germline tissues but become aberrantly activated in tumors due to epigenetic dysregulation and are well-known immunotherapy targets. Their recovery from an unbiased bulk RNA-seq analysis indicates that the computational pipeline captures biologically meaningful tumor signals.

In parallel, multiple epithelial differentiation markers (keratinization and squamous lineage genes) were detected, consistent with the known histopathology of lung squamous cell carcinoma.

The convergence of: Cancer-Testis Antigen activationepithelial differentiation programs, tumor–normal transcriptomic separation, demonstrates that transcriptome-level analysis alone can recover canonical molecular features of LUSC.

![Tumor Gene Signature Heatmap](https://raw.githubusercontent.com/ag48665/tcga-lusc-biomarker-analysis/main/figures/heatmap.png)

The heatmap shows z-score scaled expression of the most variable significantly differentially expressed genes across selected tumor (TP) and normal (NT) samples.

Genes were selected using:
adjusted p-value (FDR) < 0.05  ranked by log2 fold-change

The resulting candidate biomarker panel is provided here: results/top_LUSC_signature_genes.csv

This gene set represents a candidate transcriptional signature capable of distinguishing LUSC tumors from normal lung tissue and provides potential targets for diagnostic and immunotherapeutic investigation.

---

## Prognostic Gene Signature (Multivariate Survival Model)

To evaluate whether the transcriptomic signature provides clinical value beyond standard pathology, we performed a multivariate Cox proportional hazards analysis including:

- signature score
- patient age
- AJCC pathological stage

![Multivariate Cox forest plot](https://raw.githubusercontent.com/ag48665/tcga-lusc-biomarker-analysis/main/figures/cox_forest_plot.png)

The gene signature remained a significant predictor of overall survival after adjusting for clinical covariates.

Model performance:

- Concordance index (C-index): ~0.60
- Likelihood ratio test: significant
- Log-rank test: significant

This indicates that the transcriptomic signal captured by RNA-seq contains prognostic information not fully explained by tumor staging alone.

In other words, two patients with the same clinical stage can still have different risk profiles based on tumor gene expression.

This suggests potential application of the signature in **risk stratification** and patient prognosis assessment.


### Prognostic Transcriptomic Signature

![Signature survival](figures/signature_KM.png)

Patients stratified by the TCGA-derived gene expression signature show significantly different overall survival.

Importantly, the effect persists after adjustment for clinical covariates (age and pathological stage) in a multivariate Cox model, indicating that the signature captures independent biological risk.

---

## External validation (GEO – GSE33479)

To assess whether the TCGA-derived signature generalizes beyond TCGA, we validated it in an independent GEO microarray cohort (**GSE33479**).

**PCA using signature genes**  
![External validation PCA](figures/GEO_PCA_signature.png)

**Signature score / expression differences**  
![External validation boxplot](figures/GEO_signature_boxplot.png)

Reproducible script: reports/08_external_validation_GEO.R

---
## Prognostic Gene Signature (TCGA-LUSC)

To determine whether tumor transcriptional programs have clinical relevance, we constructed a multi-gene expression signature derived from significantly differentially expressed genes and evaluated its prognostic value in the TCGA-LUSC cohort.

For each patient, a signature score was computed as the mean z-score of normalized expression across signature genes. Patients were stratified into HIGH and LOW groups using the cohort median.

Kaplan–Meier survival analysis
Patients with high signature scores showed significantly improved overall survival compared to the low-score group.

Cox proportional hazards model:

Hazard Ratio (HR) ≈ 0.75

95% CI: 0.57–0.98

log-rank p ≈ 0.036

This indicates the gene signature functions as a prognostic biomarker in lung squamous cell carcinoma.

Importantly, the signature is derived solely from tumor vs normal transcriptomic differences, demonstrating that tumor-specific expression programs can predict patient outcomes.

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

## Reproducibility

Run the entire pipeline from the project root:
source("run_analysis.R")


The pipeline automatically:

downloads TCGA-LUSC RNA-seq data

performs quality control

computes PCA

identifies differentially expressed genes

generates figures

performs Gene Ontology enrichment and functional pathway interpretation

runs survival analysis

Software versions are recorded in session_info.txt.

---

## Methods and Techniques

RNA-seq differential expression analysis (DESeq2)
TCGA data acquisition (TCGAbiolinks)
Gene Ontology enrichment analysis
Kaplan–Meier survival modeling
Cox proportional hazards regression
Prognostic biomarker construction
External dataset validation (GEO)
Transcriptomic data visualization (ggplot2, pheatmap)
Reproducible bioinformatics workflows

---
## Use the tumor signature on your own dataset

You can compute the TCGA-LUSC signature score on any RNA-seq dataset:

source("scripts/predict_signature_score.R")

score <- predict_signature_score(expression_matrix, signature_genes)

---

## Author

Agata Gabara
Bioinformatics / Computational Biology
