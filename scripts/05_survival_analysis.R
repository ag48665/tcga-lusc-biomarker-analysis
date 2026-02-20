# scripts/05_survival_analysis.R
# Survival (Kaplan–Meier) for TCGA-LUSC using expression (VST) + clinical data

library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)

dir.create("figures", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

# -----------------------------
# Load expression (VST)
# -----------------------------
vsd <- readRDS("data/vsd_lusc_tp_nt.rds")
expr <- assay(vsd)

# -----------------------------
# Clinical cache (IMPORTANT: do this before merge)
# -----------------------------
CLIN_PATH <- "data/clinical_lusc.rds"

if (file.exists(CLIN_PATH)) {
  message("✓ Using cached clinical data: ", CLIN_PATH)
  clinical <- readRDS(CLIN_PATH)
} else {
  message("Downloading clinical data from GDC (first run only)...")
  clinical <- TCGAbiolinks::GDCquery_clinic(project = "TCGA-LUSC", type = "clinical")
  saveRDS(clinical, CLIN_PATH)
  message("✓ Saved clinical cache to: ", CLIN_PATH)
}

# -----------------------------
# Choose gene
# -----------------------------
gene_symbol <- "KRT6A"

# Helper: find gene row index by SYMBOL using rowData
find_gene_row <- function(vsd_obj, gene_sym) {
  rd <- rowData(vsd_obj)
  possible_cols <- c("gene_name", "external_gene_name", "hgnc_symbol", "symbol", "gene")
  hit_cols <- intersect(possible_cols, colnames(rd))
  
  if (length(hit_cols) > 0) {
    for (col in hit_cols) {
      vals <- as.character(rd[[col]])
      idx <- which(toupper(vals) == toupper(gene_sym))
      if (length(idx) > 0) return(idx[1])
    }
  }
  return(NA_integer_)
}

gene_row <- find_gene_row(vsd, gene_symbol)

# If symbol not found, fallback to top DE feature
if (is.na(gene_row)) {
  message(gene_symbol, " not found in rowData(vsd). Falling back to top DE feature...")
  
  de_path <- "results/DE_genes_LUSC.csv"
  if (!file.exists(de_path)) stop("Missing file: ", de_path, " (run 03_differential_expression.R first)")
  
  de <- read.csv(de_path, row.names = 1)
  de <- de[!is.na(de$padj), ]
  de <- de[order(de$padj), ]
  
  top_id <- rownames(de)[1]  # likely Ensembl ID
  message("Top DE feature ID: ", top_id)
  
  gene_row <- which(rownames(expr) == top_id)[1]
  if (is.na(gene_row)) stop("Could not match top DE ID to expression matrix rownames.")
  
  gene_label <- top_id
} else {
  gene_label <- gene_symbol
}

gene_expression <- as.numeric(expr[gene_row, ])

# -----------------------------
# Map sample barcodes -> patient IDs
# -----------------------------
sample_ids <- colnames(expr)
patient_ids <- substr(sample_ids, 1, 12)

surv_df <- data.frame(
  submitter_id = patient_ids,
  expression = gene_expression,
  stringsAsFactors = FALSE
)

# One patient may appear multiple times -> keep first
surv_df <- surv_df[!duplicated(surv_df$submitter_id), ]

# -----------------------------
# Merge expression with clinical
# -----------------------------
if (!"submitter_id" %in% names(clinical)) {
  stop("Clinical table has no 'submitter_id'. Columns are: ", paste(names(clinical), collapse = ", "))
}

merged <- merge(clinical, surv_df, by = "submitter_id")

# -----------------------------
# Define groups: High vs Low expression (median split)
# -----------------------------
median_expr <- median(merged$expression, na.rm = TRUE)
merged$group <- ifelse(merged$expression > median_expr, "High", "Low")

# -----------------------------
# Build survival time + status (robust)
# -----------------------------
death_col  <- intersect(names(merged), c("days_to_death", "days_to_death_demographic"))
follow_col <- intersect(names(merged), c("days_to_last_follow_up", "days_to_last_follow_up_demographic"))

if (length(death_col) == 0 && length(follow_col) == 0) {
  stop("No survival time columns found. Expected days_to_death and/or days_to_last_follow_up.")
}

merged$time <- NA_real_
if (length(death_col) > 0) merged$time <- as.numeric(merged[[death_col[1]]])

if (length(follow_col) > 0) {
  idx <- is.na(merged$time)
  merged$time[idx] <- as.numeric(merged[[follow_col[1]]][idx])
}

if (!"vital_status" %in% names(merged)) {
  stop("No vital_status column found in merged clinical data.")
}

merged$status <- ifelse(tolower(merged$vital_status) %in% c("dead", "deceased"), 1, 0)

# Filter usable rows
merged <- merged[!is.na(merged$time) & merged$time > 0, ]

# -----------------------------
# Sanity checks
# -----------------------------
cat("Gene used:", gene_label, "\n")
cat("Patients after merge:", nrow(merged), "\n")
cat("Events (dead):", sum(merged$status == 1, na.rm = TRUE), "\n")
print(table(merged$group))

if (nrow(merged) < 20) stop("Too few patients with survival info after filtering (n < 20).")

# -----------------------------
# Kaplan–Meier
# -----------------------------
fit <- survfit(Surv(time, status) ~ group, data = merged)

plt <- ggsurvplot(
  fit,
  data = merged,
  pval = TRUE,
  risk.table = TRUE,
  title = paste("TCGA-LUSC survival:", gene_label, "(High vs Low)")
)

# Save BOTH plot + whole figure (plot only is usually enough for README)
out_plot <- paste0("figures/survival_", gsub("[^A-Za-z0-9_\\-]", "_", gene_label), ".png")
ggsave(out_plot, plt$plot, width = 7, height = 5, dpi = 300)

cat("Saved plot:", out_plot, "\n")