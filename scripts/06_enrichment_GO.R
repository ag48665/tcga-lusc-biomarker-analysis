# scripts/06_enrichment_GO.R

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(DOSE)
  library(ggplot2)
})

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ---- Load DE results ----
res <- read.csv("results/DE_genes_LUSC.csv")

# sanity check
if (!"ensembl_id" %in% colnames(res)) {
  stop("Column 'ensembl_id' not found in DE_genes_LUSC.csv. Check 03_differential_expression.R export.")
}
if (!"padj" %in% colnames(res) || !"log2FoldChange" %in% colnames(res)) {
  stop("Expected columns 'padj' and 'log2FoldChange' not found. Check your DE results format.")
}

# remove NA padj
res <- res[!is.na(res$padj), ]

# significant genes (adjust thresholds if needed)
sig <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]

if (nrow(sig) == 0) {
  stop("No significant genes found with padj < 0.05 & |log2FC| > 1. Try relaxing thresholds.")
}

# ---- Extract clean ENSEMBL IDs from ensembl_id column ----
genes_raw <- sig$ensembl_id

# remove symbol part after "|", if present
genes <- sub("\\|.*$", "", genes_raw)

# remove ENSEMBL version after ".", if present
genes <- sub("\\..*$", "", genes)

# unique + drop empty
genes <- unique(genes)
genes <- genes[genes != "" & !is.na(genes)]

cat("Significant genes (input):", length(genes_raw), "\n")
cat("Unique ENSEMBL IDs (clean):", length(genes), "\n")
cat("Example IDs:", paste(head(genes, 5), collapse = ", "), "\n\n")

# ---- Map ENSEMBL -> ENTREZ ----
gene_map <- bitr(
  genes,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

if (is.null(gene_map) || nrow(gene_map) == 0) {
  stop("bitr() returned 0 mappings. Check that your IDs are ENSEMBL (ENSG...).")
}

cat("Mapped ENSEMBL -> ENTREZ rows:", nrow(gene_map), "\n")
cat("Unique mapped ENSEMBL:", length(unique(gene_map$ENSEMBL)), "\n")
cat("Unique mapped ENTREZ:", length(unique(gene_map$ENTREZID)), "\n\n")

write.csv(gene_map, "results/ensembl_to_entrez_mapping.csv", row.names = FALSE)

entrez_genes <- unique(gene_map$ENTREZID)

# ---- GO enrichment (Biological Process) ----
ego <- enrichGO(
  gene = entrez_genes,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

ego_df <- as.data.frame(ego)
write.csv(ego_df, "results/go_enrichment_BP.csv", row.names = FALSE)

cat("GO terms found:", nrow(ego_df), "\n")

# If no terms found, save a simple message plot instead of a blank image
if (nrow(ego_df) == 0) {
  png("figures/go_enrichment.png", width = 1200, height = 900)
  plot.new()
  text(
    0.5, 0.5,
    "No significant GO BP terms found.\nTry relaxing thresholds (padj/log2FC) or check mappings.",
    cex = 1.2
  )
  dev.off()
  cat("Saved placeholder plot: figures/go_enrichment.png\n")
} else {
  p <- dotplot(ego, showCategory = 15) +
    ggtitle("GO Biological Processes - TCGA LUSC") +
    theme(plot.title = element_text(hjust = 0.5))
  
  png("figures/go_enrichment.png", width = 1200, height = 900)
  print(p)
  dev.off()
  
  cat("Saved plot: figures/go_enrichment.png\n")
}
ego_df <- as.data.frame(ego)
write.csv(ego_df, "results/go_enrichment_BP.csv", row.names = FALSE)
