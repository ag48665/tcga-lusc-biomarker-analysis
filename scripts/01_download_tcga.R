# =========================================================

# TCGA LUSC DOWNLOAD SCRIPT (STABLE OFFLINE VERSION)

# =========================================================

cat("\n=== TCGA LUSC DATA DOWNLOAD ===\n")

DATA_DIR  <- "data"
DATA_PATH <- file.path(DATA_DIR, "tcga_lusc.rds")

if (!dir.exists(DATA_DIR))
  dir.create(DATA_DIR, recursive = TRUE)

# ---------------------------------------------------------

# 1. LOAD FROM CACHE (OFFLINE MODE)

# ---------------------------------------------------------

if (file.exists(DATA_PATH)) {
  cat("✓ Using cached TCGA dataset (offline mode)\n")
  lusc <- readRDS(DATA_PATH)
  assign("lusc", lusc, envir = .GlobalEnv)
  cat("=== DONE ===\n")
  
}

# ---------------------------------------------------------

# 2. LIBRARIES

# ---------------------------------------------------------

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

# ---------------------------------------------------------

# 3. QUERY

# ---------------------------------------------------------

cat("Creating query...\n")

query <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# ---------------------------------------------------------

# 4. DOWNLOAD

# ---------------------------------------------------------

cat("Downloading data (first run only)...\n")
GDCdownload(query, method = "api")

# ---------------------------------------------------------

# 5. PREPARE (RETRY SAFE)

# ---------------------------------------------------------

cat("Preparing dataset (may retry if GDC is busy)...\n")

success <- FALSE
attempts <- 1

while (!success && attempts <= 5) {
  
  cat("Attempt:", attempts, "\n")
  
  lusc <- tryCatch(
    GDCprepare(query),
    error = function(e) NULL
  )
  
  if (!is.null(lusc)) {
    success <- TRUE
  } else {
    cat("GDC busy... waiting 60 seconds\n")
    Sys.sleep(60)
    attempts <- attempts + 1
  }
}

if (!success)
  stop("GDC API unreachable. Try again later (not your fault).")

# ---------------------------------------------------------

# 6. SAVE OFFLINE CACHE

# ---------------------------------------------------------

saveRDS(lusc, DATA_PATH)
assign("lusc", lusc, envir = .GlobalEnv)

cat("✓ Dataset prepared and cached locally.\n")
cat("From now on the pipeline works offline.\n")
cat("=== DONE ===\n")
