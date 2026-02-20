cat("==== TCGA LUSC RNA-seq pipeline ====\n")

# 1. Download + preprocessing
source("scripts/01_download_tcga.R")
cat("Download finished\n")

# 2. Quality control
source("scripts/01_qc_report.R")
cat("QC finished\n")

# 3. PCA
source("scripts/02_pca.R")
cat("PCA finished\n")

# 4. Differential expression
source("scripts/03_differential_expression.R")
cat("Differential expression finished\n")

# 5. Volcano plot
source("scripts/04_volcano_plot.R")
cat("Volcano plot finished\n")

# 6. Survival
source("scripts/05_survival_analysis.R")
cat("Survival analysis finished\n")

# 7. GO enrichment
source("scripts/06_enrichment_GO.R")
cat("GO enrichment finished\n")

# 8. Heatmap
source("scripts/07_heatmap_top_genes.R")
cat("Heatmap finished\n")

cat("==== PIPELINE COMPLETED SUCCESSFULLY ====\n")