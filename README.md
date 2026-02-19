# TCGA Lung Squamous Cell Carcinoma Biomarker Analysis

This project analyzes transcriptomic alterations in Lung Squamous Cell Carcinoma (TCGA-LUSC) using bulk RNA-seq data from The Cancer Genome Atlas (TCGA).

The aim is to identify genes differentially expressed between tumor and normal lung tissue and explore candidate biomarkers associated with lung cancer.

---

## Methods Overview
1. Data acquisition using TCGAbiolinks
2. Filtering tumor (TP) vs normal (NT) samples
3. Variance stabilizing transformation (DESeq2)
4. Principal Component Analysis (PCA)
5. Differential gene expression analysis (DESeq2)
6. Visualization using volcano plot

---

## Results

### PCA: Tumor vs Normal
![PCA](figures/pca_lusc_tp_vs_nt.png)

The PCA shows clear separation between tumor and normal lung tissue, indicating global transcriptomic differences.

### Differential Expression (Volcano Plot)
![Volcano](figures/volcano_lusc.png)

Thousands of genes are significantly up- or down-regulated in tumor samples, consistent with large-scale transcriptional reprogramming in lung cancer.

---

## Tools
R, Bioconductor, TCGAbiolinks, DESeq2, ggplot2

---

## Author
Agata Gabara
