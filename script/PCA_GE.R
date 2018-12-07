################################################################################
##                 PCA analysis
################################################################################
library(DESeq2)
library(ggplot2)
library(viridis)

# Working directory
#setwd("")

# Load normalized data
rld.genes <- readRDS("rna/genes.deseq.RDS")

# PCA
pca.genes <- plotPCA(rld.genes, intgroup="stage", returnData = T)
percentVar.genes <- round(100 * attr(pca.genes, "percentVar"))

ggplot(pca.genes, aes(PC1, PC2, color= factor(stage, levels = c("oocyte","2PN","2C","4C","8C","morula","ICM","TE")))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar.genes[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.genes[2],"% variance")) + 
  scale_color_viridis(discrete = T) + labs(color = "stage") +
  theme_bw()
