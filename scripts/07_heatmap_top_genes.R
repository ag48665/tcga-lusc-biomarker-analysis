# scripts/07_heatmap_top_genes.R
suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(matrixStats)
})

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

cat("=== HEATMAP: TOP GENES (PORTFOLIO VERSION) ===\n")

# ----------------------------
# 1) Load VST object
# ----------------------------
if (!file.exists("results/vsd.rds")) {
  stop("Missing results/vsd.rds. Run scripts/02_pca.R first.")
}
vsd <- readRDS("results/vsd.rds")
mat_all <- assay(vsd)
cat("VSD matrix:", nrow(mat_all), "genes x", ncol(mat_all), "samples\n")

# ----------------------------
# 2) Load DE results
# ----------------------------
if (!file.exists("results/DE_genes_LUSC.csv")) {
  stop("Missing results/DE_genes_LUSC.csv. Run DE step first.")
}
res <- read.csv("results/DE_genes_LUSC.csv", stringsAsFactors = FALSE)

required <- c("ensembl_id", "padj", "log2FoldChange")
missing <- setdiff(required, colnames(res))
if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]

# ----------------------------
# 3) Clean IDs same way (DE + VSD)
# ----------------------------
clean_id <- function(x) {
  x <- as.character(x)
  x <- sub("\\|.*$", "", x)
  x <- sub("\\..*$", "", x)
  trimws(x)
}

res$ensembl_clean <- clean_id(res$ensembl_id)
rownames(mat_all) <- clean_id(rownames(mat_all))

# remove duplicate gene IDs after cleaning
dup <- duplicated(rownames(mat_all))
if (any(dup)) mat_all <- mat_all[!dup, , drop = FALSE]

# ----------------------------
# 4) Balanced subsample (NT vs TP)
# ----------------------------
cd <- as.data.frame(colData(vsd))

type_col <- if ("shortLetterCode" %in% colnames(cd)) "shortLetterCode" else
  if ("group" %in% colnames(cd)) "group" else
    stop("No shortLetterCode/group in colData(vsd).")

types <- cd[[type_col]]
tumor_samples  <- rownames(cd)[types == "TP"]
normal_samples <- rownames(cd)[types == "NT"]

cat("Tumor samples:", length(tumor_samples), "\n")
cat("Normal samples:", length(normal_samples), "\n")

set.seed(42)

# For a clean README figure, fewer samples looks MUCH better in TCGA
n_show <- 10   # <- change to 20 if you want a wider heatmap

tumor_pick  <- sample(tumor_samples,  min(n_show, length(tumor_samples)))
normal_pick <- sample(normal_samples, min(n_show, length(normal_samples)))

# keep only samples that exist in the matrix
tumor_pick  <- intersect(tumor_pick,  colnames(mat_all))
normal_pick <- intersect(normal_pick, colnames(mat_all))

selected_samples <- c(normal_pick, tumor_pick)  # NT first, then TP
mat_sub <- mat_all[, selected_samples, drop = FALSE]

# annotation for selected samples
annotation <- data.frame(Type = cd[selected_samples, type_col], stringsAsFactors = FALSE)
rownames(annotation) <- selected_samples

# enforce order NT then TP
ord <- order(annotation$Type)   # NT then TP
mat_sub <- mat_sub[, ord, drop = FALSE]
annotation <- annotation[ord, , drop = FALSE]

# position for a visible gap between groups
gap_pos <- sum(annotation$Type == "NT")

# ----------------------------
# 5) Choose genes for heatmap
#    Significant genes -> take top by variance across selected samples
# ----------------------------
sig <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
if (nrow(sig) < 50) {
  cat("⚠ Few sig genes; relaxing filter to padj<0.1\n")
  sig <- res[res$padj < 0.1, ]
}

sig_ids <- unique(sig$ensembl_clean)
sig_ids <- sig_ids[sig_ids %in% rownames(mat_sub)]

if (length(sig_ids) < 10) stop("Too few significant genes matched VST matrix after cleaning.")

N <- 30
vars <- rowVars(mat_sub[sig_ids, , drop = FALSE])
top_ids <- names(sort(vars, decreasing = TRUE))[1:min(N, length(vars))]

mat <- mat_sub[top_ids, , drop = FALSE]

# ----------------------------
# 6) Z-score by gene + clip
# ----------------------------
mat_z <- t(scale(t(mat)))
mat_z[is.na(mat_z)] <- 0
mat_z[mat_z > 2] <- 2
mat_z[mat_z < -2] <- -2

# ----------------------------
# 7) Plot (portfolio-friendly)
# ----------------------------
# Save also as PNG for GitHub README
png_file <- "figures/heatmap_lusc_signature.png"

png(png_file, width = 1600, height = 900, res = 150)

pheatmap(
  mat_z,
  annotation_col = annotation,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_method = "complete",
  breaks = seq(-2, 2, length.out = 101),
  fontsize_row = 8
)

dev.off()

cat("Saved:", png_file, "\n")

pheatmap(
  mat_z,
  annotation_col = annotation,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_cols = FALSE,       # no mixing of groups
  cluster_rows = TRUE,
  clustering_method = "complete",
  breaks = seq(-2, 2, length.out = 101),
  gaps_col = gap_pos,         # visible split between NT and TP
  filename = out_file,
  width = 10,
  height = 8,
  fontsize_row = 8
)

cat("Saved:", out_file, "\n")
cat("=== HEATMAP DONE ===\n")