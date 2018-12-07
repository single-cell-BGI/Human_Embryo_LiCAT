####################################################################################
##    Heatmap correlation between replicates of human ESC and Liver LiCAT-seq
####################################################################################

library(DESeq2)
library(SummarizedExperiment)
library(pheatmap)
library(viridis)

## Set working directory
#setwd()

rld.atac <- readRDS("atac/atac.hES.hLiver.deseq.RDS")

sampleDists.atac <- dist(t(assay(rld.atac)))
sampleDistMatrix.atac <- as.matrix(sampleDists.atac)
rownames(sampleDistMatrix.atac) <- rld.atac$celltype
colnames(sampleDistMatrix.atac) <- NULL
pheatmap(sampleDistMatrix.atac,
         clustering_distance_rows=sampleDists.atac,
         clustering_distance_cols=sampleDists.atac,
         color=viridis(100))

pdf("hES.liver.ATAC.correlation.pdf")
pheatmap(cor(assay(rld.atac)),
         color=viridis(100))
dev.off()
