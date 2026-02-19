# scripts/02_pca.R
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)

dir.create("figures", showWarnings = FALSE)

# wczytaj zapisany obiekt
lusc <- readRDS("data/tcga_lusc.rds")

# sprawdź typy próbek
table(lusc$shortLetterCode)

# wybieramy Tumor (TP) i Normal (NT)
keep <- lusc$shortLetterCode %in% c("TP", "NT")
lusc2 <- lusc[, keep]

# DESeq2 dataset
dds <- DESeqDataSet(lusc2, design = ~ shortLetterCode)
dds <- DESeq(dds)

# transformacja do PCA
vsd <- vst(dds, blind = FALSE)

# PCA
p <- plotPCA(vsd, intgroup = "shortLetterCode") + ggtitle("TCGA-LUSC: Tumor (TP) vs Normal (NT)")
print(p)

ggsave("figures/pca_lusc_tp_vs_nt.png", p, width = 7, height = 5, dpi = 300)

# zapisujemy obiekty na później (żeby nie liczyć od nowa)
saveRDS(dds, "data/dds_lusc_tp_nt.rds")
saveRDS(vsd, "data/vsd_lusc_tp_nt.rds")

print(p)

ggsave("figures/pca_lusc_tp_vs_nt.png", plot = p, width = 7, height = 5, dpi = 300)
