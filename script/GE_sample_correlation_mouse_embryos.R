####################################################################################
##    Heatmap correlation between replicates of mouse embryo LiCAT-seq
####################################################################################

library(DESeq2)
library(pheatmap)

## Working directory
#setwd()

rna.SEobj <- readRDS("rna/rna.mouse.embryo.SEobj.RDS")

sample<- c("mouse4C_1", "mouse4C_2", "mouseMorula_1", "mouseMorula_2")
celltype <- gsub("_\\d","", sample)

sampleDists <- dist(t(assays(rna.SEobj)$FPKM.log2))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix.atac) <- rna.SEobj$celltype
colnames(sampleDistMatrix) <- NULL
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         color=viridis(100))

pheatmap(cor(assays(rna.SEobj)$FPKM.log2),
         color=viridis(100))
