# TCGA Lung Squamous Cell Carcinoma Biomarker Analysis

This project investigates transcriptomic alterations in Lung Squamous Cell Carcinoma (TCGA-LUSC) using bulk RNA-seq data from The Cancer Genome Atlas (TCGA).

The objective is to identify differentially expressed genes and evaluate their prognostic significance using survival modeling. The analysis includes data processing, normalization, differential expression (DESeq2), functional enrichment (GO/KEGG), and Kaplan–Meier survival analysis in R/Bioconductor.

## Methods Overview
1. Data acquisition using TCGAbiolinks
2. Filtering tumor vs normal samples
3. Variance stabilizing transformation (DESeq2)
4. Principal Component Analysis (PCA)
5. Differential gene expression analysis
6. Functional enrichment analysis (GO/KEGG)
7. Survival analysis using clinical data
