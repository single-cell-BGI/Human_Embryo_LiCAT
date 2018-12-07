
#################################################################################################
#             Read count matrix and normalization for chromatin accessibility datasets
#################################################################################################

library(GenomicAlignments)
library(SummarizedExperiment)
library(GenomicRanges)
library(DESeq2)
library(dplyr)
library(reshape2)


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##               Count the read number within each peak                       ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

#setwd("")
bam.files <-list.files(pattern = "sort.rmdup2.bam$")
ref.peak <- read.table("atac/refPeak.bed", header = F, stringsAsFactors = F)
colnames(ref.peak) <- c("chr","start","end","strand","var1","var2")
ref.peak <- makeGRangesFromDataFrame(ref.peak, keep.extra.columns = T)
counts <- summarizeOverlaps(ref.peak,bam.files,mode="Union",inter.feature=TRUE,
fragment=FALSE,singleEnd=TRUE)
saveRDS(counts,"atac/counts.refPeaks.RDS")
##==================================================================================================


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##        Normalization of chromatin accessibility profiles by DESeq2
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

counts.atac <- readRDS("atac/counts.refPeaks.RDS")
colnames(counts.atac) <- gsub(".sort.rmdup2.bam","",colnames(counts.atac))
counts.atac$stage <- sapply(strsplit(colnames(counts.atac),"_"),"[",3)
rownames(counts.atac) <- paste(seqnames(counts.atac), start(counts.atac), end(counts.atac),sep = "_")
counts.atac <- counts.atac[rowSums(assays(counts.atac)$counts > 0) > 2,]

counts.atac.deseq <- DESeqDataSet(counts.atac, ~stage)
rld.atac <- rlogTransformation(counts.atac.deseq)  
saveRDS(rld.atac, "atac/atac.deseq.RDS")
