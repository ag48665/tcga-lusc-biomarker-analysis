# ============================
# ENSEMBL -> Gene SYMBOL (ładne nazwy do README)
# ============================
suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

ens_ids <- rownames(mat)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ens_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# fallback: jeśli brak mapowania, zostaw ENSG
gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- ens_ids[is.na(gene_symbols) | gene_symbols == ""]

# usuń duplikaty po mapowaniu (pheatmap lubi unikalne rownames)
dup <- duplicated(gene_symbols)
if (any(dup)) {
  mat <- mat[!dup, , drop = FALSE]
  gene_symbols <- gene_symbols[!dup]
}

rownames(mat) <- gene_symbols

# ============================
# Z-score + clip
# ============================
mat_z <- t(scale(t(mat)))
mat_z[is.na(mat_z)] <- 0
mat_z[mat_z > 2] <- 2
mat_z[mat_z < -2] <- -2