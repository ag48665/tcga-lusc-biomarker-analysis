# scripts/05_survival_analysis.R
# Survival (Kaplan–Meier) for TCGA-LUSC using expression (VST) + clinical data

library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)

dir.create("figures", showWarnings = FALSE)

# -----------------------------
# Load data
# -----------------------------
vsd <- readRDS("data/vsd_lusc_tp_nt.rds")
expr <- assay(vsd)

clinical <- GDCquery_clinic(project = "TCGA-LUSC", type = "clinical")

# -----------------------------
# Choose gene (try TP53 first)
# -----------------------------
gene_symbol <- "TP53"   # you can change later, e.g. "MKI67", "CD274", "PDCD1"

# Helper: find gene row index by SYMBOL using rowData, fallback to DE results
find_gene_row <- function(vsd_obj, gene_sym) {
  rd <- rowData(vsd_obj)
  
  # Try common symbol columns in TCGA/Bioc objects
  possible_cols <- c("gene_name", "external_gene_name", "hgnc_symbol", "symbol", "gene")
  hit_col <- intersect(possible_cols, colnames(rd))
  
  if (length(hit_col) > 0) {
    for (col in hit_col) {
      vals <- as.character(rd[[col]])
      idx <- which(toupper(vals) == toupper(gene_sym))
      if (length(idx) > 0) return(idx[1])
    }
  }
  
  # If no symbol mapping exists, return NA
  return(NA_integer_)
}

gene_row <- find_gene_row(vsd, gene_symbol)

# If TP53 not found by symbol, fallback to TOP DE gene from results file
if (is.na(gene_row)) {
  message("TP53 not found as a gene symbol in rowData(vsd). Falling back to top DE gene...")
  
  de_path <- "results/DE_genes_LUSC.csv"
  if (!file.exists(de_path)) stop("Missing file: ", de_path, " (run 03_differential_expression.R first)")
  
  de <- read.csv(de_path, row.names = 1)
  de <- de[!is.na(de$padj), ]
  de <- de[order(de$padj), ]
  
  top_id <- rownames(de)[1]  # likely Ensembl ID
  message("Top DE feature ID: ", top_id)
  
  # Find exact match in expression rownames (works if they are Ensembl IDs)
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

# Merge expression with clinical
if (!"submitter_id" %in% names(clinical)) {
  stop("Clinical table has no 'submitter_id'. Columns are: ", paste(names(clinical), collapse = ", "))
}

merged <- merge(clinical, surv_df, by = "submitter_id")

# -----------------------------
# Define groups: High vs Low expression
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

# Filter usable survival rows
merged <- merged[!is.na(merged$time) & merged$time > 0, ]

# -----------------------------
# Sanity checks
# -----------------------------
cat("Gene used:", gene_label, "\n")
cat("Patients after merge:", nrow(merged), "\n")
cat("Events (dead):", sum(merged$status == 1, na.rm = TRUE), "\n")
print(table(merged$group))

if (nrow(merged) < 20) {
  stop("Too few patients with survival info after filtering (n < 20).")
}

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

out_file <- "figures/survival_gene.png"
ggsave(out_file, plt$plot, width = 7, height = 5, dpi = 300)

cat("Saved plot:", out_file, "\n")

# ---- Clinical cache ----
CLIN_PATH <- "data/clinical_lusc.rds"

if (file.exists(CLIN_PATH)) {
  cat("✓ Using cached clinical data\n")
  clinical <- readRDS(CLIN_PATH)
} else {
  cat("Downloading clinical data from GDC (first run only)...\n")
  clinical <- TCGAbiolinks::GDCquery_clinic(project = "TCGA-LUSC", type = "clinical")
  saveRDS(clinical, CLIN_PATH)
  cat("✓ Saved clinical cache to: ", CLIN_PATH, "\n")
}