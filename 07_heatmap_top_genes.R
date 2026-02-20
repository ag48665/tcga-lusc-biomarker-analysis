# =========================================================

# 07_heatmap_top_genes.R

# Heatmap of top differentially expressed genes (TCGA-LUSC)

# =========================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(dplyr)
})

dir.create("figures", showWarnings = FALSE)

cat("=== HEATMAP: TOP DE GENES ===\n")

# -------------------------------

# 1) Load objects

# -------------------------------

vsd <- readRDS("results/vsd.rds")
res <- read.csv("results/DE_genes_LUSC.csv")

# remove NA p-values

res <- res[!is.na(res$padj), ]

# -------------------------------

# 2) Select significant genes

# -------------------------------

sig <- res %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

cat("Significant genes:", nrow(sig), "\n")

# take TOP 30 most significant

top_genes <- sig %>%
  arrange(padj) %>%
  slice(1:30)

genes <- top_genes$ensembl_id

# -------------------------------

# 3) Expression matrix

# -------------------------------

mat <- assay(vsd)[genes, ]

# Z-score normalization (CRUCIAL for heatmap)

mat <- t(scale(t(mat)))

# -------------------------------

# 4) Sample annotation

# -------------------------------

annotation <- as.data.frame(colData(vsd)[, "shortLetterCode", drop = FALSE])
colnames(annotation) <- "Type"

# -------------------------------

# 5) Plot heatmap

# -------------------------------

pheatmap(
  mat,
  annotation_col = annotation,
  show_rownames = FALSE,
  fontsize_col = 8,
  clustering_method = "complete",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  filename = "figures/heatmap_top30_DE_genes.png",
  width = 8,
  height = 10
)

cat("Saved: figures/heatmap_top30_DE_genes.png\n")
cat("=== HEATMAP DONE ===\n")
