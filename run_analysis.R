cat("==== TCGA LUSC RNA-seq pipeline ====\n")

source("scripts/01_download_tcga.R")
cat("Download finished\n")

source("scripts/02_pca.R")
cat("PCA finished\n")

source("scripts/03_differential_expression.R")
cat("Differential expression finished\n")

source("scripts/04_volcano_plot.R")
cat("Volcano plot finished\n")

source("scripts/05_survival_analysis.R")
cat("Survival analysis finished\n")

source("scripts/06_enrichment_GO.R")
cat("GO enrichment finished\n")

cat("==== PIPELINE COMPLETED SUCCESSFULLY ====\n")