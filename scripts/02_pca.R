# scripts/02_pca.R
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# wczytaj zapisany obiekt
lusc <- readRDS("data/tcga_lusc.rds")

# sprawdź typy próbek
table(lusc$shortLetterCode)

# wybieramy Tumor (TP) i Normal (NT)
keep <- lusc$shortLetterCode %in% c("TP", "NT")
lusc2 <- lusc[, keep]

# DESeq2 dataset (TYLKO do transformacji — bez DESeq!)
dds <- DESeqDataSet(lusc2, design = ~ shortLetterCode)

# transformacja do PCA
vsd <- vst(dds, blind = FALSE)

# PCA
p <- plotPCA(vsd, intgroup = "shortLetterCode") +
  ggtitle("TCGA-LUSC: Tumor (TP) vs Normal (NT)")

print(p)

ggsave("figures/pca_lusc_tp_vs_nt.png", plot = p, width = 7, height = 5, dpi = 300)

# zapisujemy obiekty na później
saveRDS(dds, "results/dds.rds")
saveRDS(vsd, "results/vsd.rds")