################################################################################
##                 PCA analysis
################################################################################
library(DESeq2)
library(SummarizedExperiment)

rld.atac <- readRDS("atac/atac.deseq.RDS")

rld.atac.rm.oocyte <- rld.atac[,!colnames(rld.atac) %in% c("Embryo_A_oocyte_1","Embryo_A_oocyte_2")]
pca.atac <- plotPCA(rld.atac.rm.oocyte, intgroup="stage", returnData = T)
percentVar.atac <- round(100 * attr(pca.atac, "percentVar"))

ggplot(pca.atac, aes(PC1, PC2, color= factor(stage, levels = c("2PN","2C","4C","8C","morula","ICM","TE")))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar.atac[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.atac[2],"% variance")) + 
  labs(color = "stage") +
  theme_bw()
