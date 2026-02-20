# scripts/01_qc_report.R
suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
})

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("reports", showWarnings = FALSE)

cat("=== QC REPORT ===\n")

# ----------------------------
# Helpers
# ----------------------------
get_type_col <- function(cd) {
  if ("shortLetterCode" %in% colnames(cd)) return("shortLetterCode")
  if ("group" %in% colnames(cd)) return("group")
  return(NULL)
}

# Bulletproof PNG saver (no width/height args to avoid list->integer issues)
save_png <- function(filename, expr) {
  png(filename = filename, width = 1200, height = 800, res = 150)
  on.exit(dev.off(), add = TRUE)
  eval.parent(substitute(expr))
}

# ----------------------------
# 1) Load VST object
# ----------------------------
if (!file.exists("results/vsd.rds")) {
  stop("Missing results/vsd.rds. Run your VST/PCA step first (scripts/02_pca.R).")
}
vsd <- readRDS("results/vsd.rds")
mat_vst <- assay(vsd)
cd <- as.data.frame(colData(vsd))

type_col <- get_type_col(cd)
if (is.null(type_col)) {
  stop("No 'shortLetterCode' or 'group' in colData(vsd). Add it in your preprocessing step.")
}
cd$Type <- as.character(cd[[type_col]])
cd$Sample <- rownames(cd)

cat("VST matrix:", nrow(mat_vst), "genes x", ncol(mat_vst), "samples\n")
cat("Type column:", type_col, "\n")

type_table <- sort(table(cd$Type), decreasing = TRUE)
print(type_table)

# ----------------------------
# 2) Try to load raw counts (optional but recommended)
# ----------------------------
has_counts <- FALSE
counts_mat <- NULL

if (file.exists("results/dds.rds")) {
  dds <- readRDS("results/dds.rds")
  counts_mat <- counts(dds)
  has_counts <- TRUE
  cat("Loaded raw counts from results/dds.rds\n")
} else if (file.exists("results/dds_object.rds")) {
  dds <- readRDS("results/dds_object.rds")
  counts_mat <- counts(dds)
  has_counts <- TRUE
  cat("Loaded raw counts from results/dds_object.rds\n")
} else {
  cat("⚠ Raw counts object (dds.rds) not found. Library-size QC will be skipped.\n")
}

# ----------------------------
# 3) QC: Library sizes (if raw counts available)
# ----------------------------
libsize_file <- "figures/qc_library_sizes.png"
if (has_counts) {
  libsize <- colSums(counts_mat)
  lib_df <- data.frame(
    Sample = names(libsize),
    LibrarySize = as.numeric(libsize),
    Type = cd[names(libsize), "Type"]
  )
  
  lib_df <- lib_df[order(lib_df$LibrarySize, decreasing = TRUE), ]
  
  save_png(libsize_file, {
    par(mar = c(10, 4, 3, 1))
    barplot(
      height = lib_df$LibrarySize,
      names.arg = rep("", nrow(lib_df)),
      las = 2,
      main = "Library sizes (raw counts)",
      ylab = "Total counts per sample"
    )
  })
  
  write.csv(lib_df, "results/qc_library_sizes.csv", row.names = FALSE)
} else {
  libsize_file <- NULL
}

# ----------------------------
# 4) QC: VST density plot
# ----------------------------
dens_file <- "figures/qc_vst_density.png"
set.seed(1)
genes_sub <- sample(rownames(mat_vst), min(5000, nrow(mat_vst)))
x <- as.vector(mat_vst[genes_sub, ])

save_png(dens_file, {
  hist(
    x, breaks = 80,
    main = "VST expression distribution (subset of genes)",
    xlab = "VST value",
    border = NA
  )
})

# ----------------------------
# 5) QC: PCA (top variable genes)
# ----------------------------
pca_file <- "figures/qc_pca.png"
vars <- apply(mat_vst, 1, var)
top <- names(sort(vars, decreasing = TRUE))[1:min(2000, length(vars))]
pca <- prcomp(t(mat_vst[top, , drop = FALSE]), scale. = FALSE)

pc1 <- pca$x[, 1]
pc2 <- pca$x[, 2]
pve <- (pca$sdev^2) / sum(pca$sdev^2)

save_png(pca_file, {
  par(mar = c(5, 5, 4, 1))
  plot(
    pc1, pc2,
    pch = 19,
    xlab = paste0("PC1 (", round(100 * pve[1], 1), "%)"),
    ylab = paste0("PC2 (", round(100 * pve[2], 1), "%)"),
    main = "PCA (top variable genes)"
  )
  types_f <- factor(cd[colnames(mat_vst), "Type"])
  cols <- as.numeric(types_f)
  points(pc1, pc2, pch = 19, col = cols)
  legend("topright", legend = levels(types_f), pch = 19,
         col = seq_along(levels(types_f)), bty = "n")
})

# ----------------------------
# 6) QC: Sample correlation heatmap
# ----------------------------
corr_file <- "figures/qc_sample_correlation.png"
cor_mat <- cor(mat_vst, method = "pearson")

ann <- data.frame(Type = cd[colnames(mat_vst), "Type"])
rownames(ann) <- colnames(mat_vst)

pheatmap(
  cor_mat,
  annotation_col = ann,
  annotation_row = ann,
  show_colnames = FALSE,
  show_rownames = FALSE,
  filename = corr_file,
  width = 10,
  height = 9
)

# ----------------------------
# 7) Outlier detection
# ----------------------------
dist_mat <- as.matrix(1 - cor_mat)
avg_dist <- rowMeans(dist_mat)
outliers <- sort(avg_dist, decreasing = TRUE)

outlier_df <- data.frame(
  Sample = names(outliers),
  AvgDistance = as.numeric(outliers),
  Type = cd[names(outliers), "Type"]
)
write.csv(outlier_df, "results/qc_outlier_ranking.csv", row.names = FALSE)
top_outliers <- head(outlier_df, 5)

# ----------------------------
# 8) Markdown report
# ----------------------------
report_path <- "reports/qc_report.md"

md <- c(
  "# QC Report (TCGA-LUSC RNA-seq)",
  "",
  "## Summary",
  paste0("- Samples (total): **", ncol(mat_vst), "**"),
  paste0("- Genes (VST matrix): **", nrow(mat_vst), "**"),
  "",
  "### Sample counts by Type",
  "",
  paste0("- ", paste(names(type_table), as.integer(type_table), sep=": ", collapse = " | ")),
  "",
  "## Plots",
  ""
)

if (!is.null(libsize_file)) {
  md <- c(md,
          "### Library sizes (raw counts)",
          "",
          "![](../figures/qc_library_sizes.png)",
          "")
} else {
  md <- c(md,
          "### Library sizes (raw counts)",
          "",
          "_Skipped: raw counts object (dds.rds) not found._",
          "")
}

md <- c(md,
        "### VST distribution",
        "",
        "![](../figures/qc_vst_density.png)",
        "",
        "### PCA (top variable genes)",
        "",
        "![](../figures/qc_pca.png)",
        "",
        "### Sample–sample correlation (VST)",
        "",
        "![](../figures/qc_sample_correlation.png)",
        "",
        "## Potential outliers (highest average distance)",
        "",
        paste0("Top 5: ", paste(top_outliers$Sample, "(", top_outliers$Type, ")", sep = "", collapse = ", ")),
        "",
        "Outputs saved to:",
        "- `figures/` (plots)",
        "- `results/qc_outlier_ranking.csv`",
        if (has_counts) "- `results/qc_library_sizes.csv`" else NULL
)

writeLines(md, report_path)

cat("Saved plots:\n")
if (exists("libsize_file") && !is.null(libsize_file)) cat("-", libsize_file, "\n")
if (exists("dens_file")) cat("-", dens_file, "\n")
if (exists("pca_file")) cat("-", pca_file, "\n")
if (exists("corr_file")) cat("-", corr_file, "\n")
if (exists("report_path")) cat("Saved report:", report_path, "\n")
cat("=== QC REPORT DONE ===\n")