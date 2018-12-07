####################################################################################
##    Heatmap correlation between replicates of mouse embryo LiCAT-seq
####################################################################################

library(DESeq2)
library(pheatmap)

## Working directory
#setwd()

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##        Heatmap correlation between replicates
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
rld.atac <- readRDS("atac/atac.mouse4C.Morula.deseq.RDS")

sampleDists.atac <- dist(t(assay(rld.atac)))
sampleDistMatrix.atac <- as.matrix(sampleDists.atac)
rownames(sampleDistMatrix.atac) <- rld.atac$celltype
colnames(sampleDistMatrix.atac) <- NULL
pheatmap(sampleDistMatrix.atac,
         clustering_distance_rows=sampleDists.atac,
         clustering_distance_cols=sampleDists.atac,
         color=viridis(100))
pheatmap(cor(assay(rld.atac)),
         color=viridis(100))