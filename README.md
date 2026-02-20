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
## Biological interpretation

Differential expression analysis followed by Gene Ontology enrichment revealed strong activation of epithelial differentiation and immune-related processes in TCGA-LUSC tumors.

The most significant enriched terms included keratinization, epidermis development, and skin development, which are hallmark features of squamous cell carcinoma. This indicates that tumor cells undergo squamous differentiation and acquire transcriptional programs characteristic of stratified epithelium. Importantly, this molecular signal matches histopathological definitions of lung squamous cell carcinoma.

Additionally, enrichment of immune-related pathways such as leukocyte mediated immunity, humoral immune response, chemotaxis, and cell recognition suggests substantial immune cell infiltration within the tumor microenvironment. This finding is consistent with the known immunogenicity of LUSC and helps explain the clinical responsiveness of these tumors to immune checkpoint inhibitors.

Overall, RNA-seq data alone was sufficient to reconstruct both the histological identity of the tumor and the presence of tumor-immune interactions.

## Results

### PCA: Tumor vs Normal
![PCA](figures/pca_lusc_tp_vs_nt.png)

The PCA shows clear separation between tumor and normal lung tissue, indicating global transcriptomic differences.

### Differential Expression (Volcano Plot)
![Volcano](figures/volcano_lusc.png)

Thousands of genes are significantly up- or down-regulated in tumor samples, consistent with large-scale transcriptional reprogramming in lung cancer.

---
## Quality control

QC summary and plots are in `reports/qc_report.md`.

Key QC visuals:
- PCA (VST, top variable genes): `figures/qc_pca.png`
- Sample–sample correlation heatmap: `figures/qc_sample_correlation.png`


## Tools
R, Bioconductor, TCGAbiolinks, DESeq2, ggplot2

---

## Author
Agata Gabara
