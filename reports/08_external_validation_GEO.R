# ================================
# External validation of TCGA-LUSC signature (GEO)
# ================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

dir.create("external_validation", showWarnings = FALSE)
dir.create("external_validation/figures", showWarnings = FALSE)

cat("Downloading GEO dataset...\n")

gset_list <- getGEO("GSE33479", GSEMatrix = TRUE)
gset <- gset_list[[1]]

expr  <- exprs(gset)
pheno <- pData(gset)

cat("Raw matrix:", nrow(expr), "probes x", ncol(expr), "samples\n")

# ---------------------------------
# Load your signature genes (ENSEMBL -> SYMBOL)
# ---------------------------------
sig <- read.csv("results/top_LUSC_signature_genes.csv", stringsAsFactors = FALSE)

# Clean ENSG (remove version if exists)
sig$ensembl_clean <- sub("\\..*$", "", sig$ensembl_id)

sym <- suppressWarnings(mapIds(
  org.Hs.eg.db,
  keys      = sig$ensembl_clean,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
))

sig$symbol <- unname(sym)
sig <- sig[!is.na(sig$symbol) & sig$symbol != "", ]
genes <- unique(sig$symbol)

cat("Signature genes (SYMBOL) loaded:", length(genes), "\n")

# ---------------------------------
# Map microarray probes -> Gene SYMBOL (GPL annotation)
# ---------------------------------
cat("Downloading platform annotation...\n")
platform <- annotation(gset)  # e.g. "GPL6480"
gpl <- getGEO(platform, AnnotGPL = TRUE)
annot <- Table(gpl)

sym_candidates <- grep("symbol|gene.?symbol|genesymbol", colnames(annot),
                       ignore.case = TRUE, value = TRUE)

if (length(sym_candidates) == 0) {
  stop("Could not find a SYMBOL column in GPL table. Available columns:\n",
       paste(colnames(annot), collapse = ", "))
}

symbol_col <- sym_candidates[1]
cat("Using SYMBOL column:", symbol_col, "\n")

annot2 <- annot[, c("ID", symbol_col)]
colnames(annot2) <- c("probe", "symbol")
annot2$symbol <- as.character(annot2$symbol)

# Some platforms have multiple symbols separated by " /// " or ";"
annot2$symbol <- sub(" ///.*$", "", annot2$symbol)
annot2$symbol <- sub(";.*$", "", annot2$symbol)
annot2$symbol <- trimws(annot2$symbol)
annot2 <- annot2[annot2$symbol != "" & !is.na(annot2$symbol), ]

# Join expression with annotation
expr_df <- as.data.frame(expr, stringsAsFactors = FALSE)
expr_df$probe <- rownames(expr_df)

expr_annot <- merge(expr_df, annot2, by = "probe")
expr_annot <- expr_annot[expr_annot$symbol != "" & !is.na(expr_annot$symbol), ]

cat("Annotated rows:", nrow(expr_annot), "\n")

# ---------------------------------
# Collapse multiple probes -> 1 gene (mean) -- FAST
# ---------------------------------
cat("Collapsing probes to genes (fast)...\n")
flush.console()

mat_expr <- as.matrix(expr_annot[, setdiff(colnames(expr_annot), c("probe", "symbol"))])
storage.mode(mat_expr) <- "numeric"

sum_by_gene <- rowsum(mat_expr, group = expr_annot$symbol, reorder = TRUE)
n_by_gene   <- as.vector(table(expr_annot$symbol)[rownames(sum_by_gene)])

expr_gene <- sweep(sum_by_gene, 1, n_by_gene, "/")

cat("Mapped genes:", nrow(expr_gene), "genes x", ncol(expr_gene), "samples\n")

# ---------------------------------
# Keep only signature genes
# ---------------------------------
common <- intersect(genes, rownames(expr_gene))
cat("Genes found in GEO:", length(common), "\n")

if (length(common) < 5) {
  stop("Too few signature genes found in GEO (<5). Check gene symbols and GPL mapping.")
}

mat <- expr_gene[common, , drop = FALSE]

# ---------------------------------
# Create Tumor vs Normal labels (safer + keep Unknown)
# ---------------------------------
cat("\n===== DEBUG PHENO (what GEO says) =====\n")
cat("pheno columns:\n")
print(colnames(pheno))

cols_show <- intersect(c("title", "source_name_ch1", "characteristics_ch1"), colnames(pheno))
if (length(cols_show) > 0) {
  cat("\nFirst 10 rows of key columns:\n")
  print(head(pheno[, cols_show, drop = FALSE], 10))
}

get_text_col <- function(colname) {
  if (!colname %in% colnames(pheno)) return(rep("", nrow(pheno)))
  x <- pheno[[colname]]
  if (is.list(x)) x <- sapply(x, function(v) paste(v, collapse = " "))
  as.character(x)
}

txt <- paste(
  get_text_col("title"),
  get_text_col("source_name_ch1"),
  get_text_col("characteristics_ch1"),
  sep = " | "
)

is_normal <- grepl(
  "normal|non-?tumou?r|non-?cancer|control|healthy|adjacent|histology\\s*:\\s*normal",
  txt, ignore.case = TRUE
)
is_tumor <- grepl(
  "tumou?r|cancer|carcinoma|malignan|neoplasm|squamous|scc|lusc|nsclc|histology\\s*:\\s*tumou?r",
  txt, ignore.case = TRUE
)

group <- ifelse(is_normal & !is_tumor, "Normal",
                ifelse(is_tumor & !is_normal, "Tumor",
                       ifelse(is_tumor & is_normal, "Tumor", "Unknown")))

group <- factor(group, levels = c("Normal", "Tumor", "Unknown"))

cat("\nGroup counts (after rules):\n")
print(table(group, useNA = "ifany"))
cat(">>> FINISHED GROUP LABELING SECTION <<<\n")
flush.console()

# Safety: must have both Normal and Tumor for comparison
if (sum(group == "Normal") == 0 || sum(group == "Tumor") == 0) {
  stop("\nCould not detect BOTH Normal and Tumor groups in GSE33479. ",
       "Check DEBUG PHENO output above and adjust keyword rules.\n")
}

# ---------------------------------
# PCA (signature only) - can include Unknown
# ---------------------------------
pca <- prcomp(t(mat), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = group
)

p1 <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_classic(base_size = 14) +
  ggtitle("External validation (GSE33479): PCA using TCGA-LUSC signature")

ggsave("external_validation/figures/GEO_PCA_signature.png", p1,
       width = 7, height = 5, dpi = 200)

# ---------------------------------
# Signature score (mean z-score across genes)
# ---------------------------------
z <- t(scale(t(mat)))
z[is.na(z)] <- 0
z[z > 2] <- 2
z[z < -2] <- -2

score <- colMeans(z)

# Boxplot only for Normal vs Tumor (exclude Unknown)
keep_nt <- group %in% c("Normal", "Tumor")
score_df <- data.frame(
  score = score[keep_nt],
  Group = droplevels(group[keep_nt])
)
# statistical test
test <- wilcox.test(score ~ Group, data = score_df)
cat("\nWilcoxon p-value (Tumor vs Normal):", test$p.value, "\n")

p2 <- ggplot(score_df, aes(Group, score, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.6) +
  theme_classic(base_size = 14) +
  ggtitle("External validation: TCGA-LUSC signature score")

ggsave("external_validation/figures/GEO_signature_boxplot.png", p2,
       width = 5.5, height = 5, dpi = 200)

cat("External validation DONE\n")
cat("Saved:\n")
cat(" - external_validation/figures/GEO_PCA_signature.png\n")
cat(" - external_validation/figures/GEO_signature_boxplot.png\n")