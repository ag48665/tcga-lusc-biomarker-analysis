# =========================================================
# 03_differential_expression.R (TCGA-LUSC, DESeq2)
# =========================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(dplyr)
  library(tibble)
})

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Expect `lusc` in GlobalEnv from 01_download_tcga.R
stopifnot(exists("lusc"))

cat("\n=== DIFFERENTIAL EXPRESSION (DESeq2) ===\n")

# -------------------------------
# 1) counts + metadata
# -------------------------------
counts <- assay(lusc, "unstranded")   # raw counts (CRITICAL)
meta <- as.data.frame(colData(lusc))

# keep only Primary Tumor (TP) and Solid Tissue Normal (NT)
meta$sample_type <- as.character(meta$sample_type)
keep_samples <- meta$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")

counts <- counts[, keep_samples, drop = FALSE]
meta   <- meta[keep_samples, , drop = FALSE]

# align order
stopifnot(all(colnames(counts) == rownames(meta)))

# set group factor (NT as reference)
meta$group <- factor(
  ifelse(meta$sample_type == "Primary Tumor", "TP", "NT"),
  levels = c("NT", "TP")
)

cat("Samples kept:\n")
print(table(meta$group))

# -------------------------------
# 2) build DESeq2 dataset
# -------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ group
)

# -------------------------------
# 3) low-count filtering
# -------------------------------
keep_genes <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep_genes, ]

cat("Genes after filtering:", nrow(dds), "\n")

# -------------------------------
# 4) run DESeq2
# -------------------------------
dds <- DESeq(dds)

# -------------------------------
# 5) results + shrinkage
# -------------------------------
res <- results(dds, contrast = c("group", "TP", "NT"))

# shrink log2FC (apeglm if available, otherwise fallback)
if (requireNamespace("apeglm", quietly = TRUE)) {
  suppressPackageStartupMessages(library(apeglm))
  res_shr <- lfcShrink(dds, coef = "group_TP_vs_NT", type = "apeglm")
  cat("✓ LFC shrinkage: apeglm\n")
} else {
  cat("⚠ apeglm not installed — using unshrunken results\n")
  res_shr <- res
}

# -------------------------------
# 6) export
# -------------------------------
res_df <- as.data.frame(res_shr)

# ⭐⭐⭐ NAJWAŻNIEJSZE — DLA GO ANALYSIS ⭐⭐⭐
# add gene IDs as a column
res_df$ensembl_id <- rownames(res_shr)

# move ensembl_id to first column
res_df <- res_df[, c("ensembl_id", setdiff(colnames(res_df), "ensembl_id"))]

# remove version numbers from ENSEMBL (ENSG000001234.5 -> ENSG000001234)
res_df$ensembl_id <- sub("\\..*$", "", res_df$ensembl_id)

write.csv(res_df, "results/DE_genes_LUSC.csv", row.names = FALSE)

cat("✓ Saved DE results with gene IDs: results/DE_genes_LUSC.csv\n")

# save objects for later scripts
saveRDS(dds, "results/dds.rds")
saveRDS(res_df, "results/de_results.rds")

cat("=== DESeq2 DONE ===\n")