library(ggplot2)

dir.create("figures", showWarnings = FALSE)

# wczytaj wyniki
res <- read.csv("results/DE_genes_LUSC.csv")

# usuń NA
res <- res[!is.na(res$padj),]

# obliczamy -log10(padj)
res$negLog10Padj <- -log10(res$padj)

# klasyfikacja genów
res$group <- "Not significant"
res$group[res$padj < 0.05 & res$log2FoldChange > 1] <- "Upregulated"
res$group[res$padj < 0.05 & res$log2FoldChange < -1] <- "Downregulated"

# volcano plot
p <- ggplot(res, aes(x = log2FoldChange, y = negLog10Padj, color = group)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_minimal() +
  labs(title = "TCGA-LUSC Differential Expression",
       x = "log2 Fold Change",
       y = "-log10(adjusted p-value)")

print(p)

ggsave("figures/volcano_lusc.png", p, width = 7, height = 5, dpi = 300)
