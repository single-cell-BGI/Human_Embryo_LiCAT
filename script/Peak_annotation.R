#################################################################################################
#                      Peak distribution along the genome
#################################################################################################

library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


setwd("~/Documents/Research/embryo/2_revision/Code/Code_submission")

##---------------------------------------------------------------------------##
##                              Path of peak files                                    
##---------------------------------------------------------------------------##
peakfiles <-  file.path("atac","NarrowPeak","merge",list.files("atac/NarrowPeak/merge"))

peakfiles <- as.list(peakfiles)
names(peakfiles) <- list.files("atac/NarrowPeak/merge")

##---------------------------------------------------------------------------##
##                              Annotate peaks                                    
##---------------------------------------------------------------------------##

peakAnnoList <- lapply(peakfiles, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE)
pdf("distribution of peaks.pdf")
plotAnnoBar(peakAnnoList)
dev.off()