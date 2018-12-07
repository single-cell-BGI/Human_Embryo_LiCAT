
################################################################################
##    gene expression correlation between duplicates of hESC and liver LiCAT
################################################################################
library(GenomicRanges)
library(SummarizedExperiment)
library(DESeq2)
library(pheatmap)
library(viridis)

rna.count.mat <- readRDS("rna/rna.count.hES.liver.mat.RDS")
hsapiens.genes <- readRDS("rna/hsapiens.genes.RDS")

ref.genes <- hsapiens.genes[hsapiens.genes$gene_id %in% rownames(rna.count.mat)]
rna.count.mat <- rna.count.mat[match(ref.genes$gene_id,rownames(rna.count.mat)),]

rna.SEobj <- SummarizedExperiment(assays =  SimpleList(FPKM = as.matrix(rna.count.mat),
                                                       FPKM.log2 = log2(as.matrix(rna.count.mat) + 1)),
                                  rowRanges = ref.genes)

colData(rna.SEobj)$celltype <- c("hES","hES","hLiver","hLiver")
colnames(rna.SEobj) <- c("hES-1","hES-2","hLiver-1","hLiver-2")



##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##        Heatmap correlation between replicates
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

pdf("hES.liver.licat.RNA.correlation.pdf")
pheatmap(cor(assays(rna.SEobj)$FPKM.log2),
         color=viridis(100))
dev.off()
