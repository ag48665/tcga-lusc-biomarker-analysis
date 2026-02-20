# TCGA Lung Squamous Cell Carcinoma (LUSC) Biomarker Analysis

This project analyzes transcriptomic alterations in **TCGA Lung Squamous Cell Carcinoma (TCGA-LUSC)** using bulk RNA-seq data from The Cancer Genome Atlas (TCGA).  
The goal is to identify genes differentially expressed between **tumor (TP)** and **normal lung (NT)** tissue and highlight candidate biomarker signatures.

---

## Dataset
- Cohort: **TCGA-LUSC**
- Samples: **TP = 511**, **NT = 51** (after filtering)
- Data type: bulk RNA-seq gene expression (counts)

---

## Methods Overview
1. Data acquisition using **TCGAbiolinks**
2. Filtering **tumor (TP)** vs **normal (NT)** samples
3. Variance stabilizing transformation (**DESeq2 VST**)
4. Principal Component Analysis (PCA)
5. Differential gene expression analysis (**DESeq2**)
6. Visualization: volcano plot + heatmaps
7. Functional interpretation: **Gene Ontology enrichment**

---

## Results

### PCA: Tumor vs Normal
![PCA](figures/pca_lusc_tp_vs_nt.png)

PCA shows separation between tumor and normal lung tissue, indicating strong global transcriptomic differences.

### Differential Expression (Volcano Plot)
![Volcano](figures/volcano_lusc.png)

Thousands of genes are significantly up- or down-regulated in tumor samples, consistent with large-scale transcriptional reprogramming in LUSC.

### Heatmap: Tumor vs Normal Gene Signature
Heatmap of a selected gene set that best separates tumor and normal samples based on VST expression.

![Tumor vs Normal gene signature](figures/heatmap_lusc_signature.png)

---

## Biological Interpretation

Gene Ontology enrichment of differentially expressed genes shows strong activation of **epithelial differentiation** and **immune-related** processes in TCGA-LUSC tumors.

Enriched terms such as **keratinization**, **epidermis development**, and **skin development** are consistent with hallmark features of squamous differentiation, aligning with the histopathology of lung squamous cell carcinoma.

Enrichment of immune-associated pathways (e.g., **leukocyte-mediated immunity**, **humoral immune response**, **chemotaxis**) suggests immune infiltration and tumor–immune interactions in the tumor microenvironment, consistent with the known immunogenicity of LUSC.

---

## Quality Control
QC summary and plots are in `reports/qc_report.md`.

Key QC visuals:
- PCA (VST, top variable genes): `figures/qc_pca.png`
- Sample–sample correlation heatmap: `figures/qc_sample_correlation.png`
- Library sizes: `figures/qc_library_sizes.png`

---

## Reproducibility
Run the full pipeline from the project root:

```r
source("run_analysis.R")
