#################################################################################
##               Gene expression data matrix            
#################################################################################
library(GenomicRanges)
library(GenomicAlignments)
library(SummarizedExperiment)
library(DESeq2)
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
organism(edb)

## Change the seqlevels style form Ensembl (default) to UCSC:
seqlevelsStyle(edb) <- "UCSC"

hsapiens.genes <- genes(edb)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##               Count the read number within genes/transcripts             
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

##===============================================================================================
hsapiens.genes <- readRDS("hsapiens.genes.RDS")
hsapiens.genes <- hsapiens.genes[grep("chr\\d{1,2}$",seqnames(hsapiens.genes))]
hsapiens.genes <- keepSeqlevels(hsapiens.genes,paste0("chr",1:22)) 
bam.files <-list.files(pattern = "sort.rmdup.bam$")
count.genes <- summarizeOverlaps(hsapiens.genes,bam.files,mode="Union",inter.feature=TRUE,fragment=FALSE,singleEnd=TRUE)
saveRDS(count.genes,"count.genes.RDS")
##==============================================================================================


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                   Generate DESeqDataSet                     
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

# Genes
count.genes <- readRDS("rna/count.genes.RDS")

colnames(count.genes) <- gsub(".genome.sort.rmdup.bam","",colnames(count.genes))
colData(count.genes)$stage <- sapply(strsplit(colnames(count.genes),"_"),"[",3)

# Filer out genes with low number of reads
count.genes <- count.genes[rowSums(assays(count.genes)$counts > 0) > 2,]

count.genes.deseq <- DESeqDataSet(count.genes, design = ~stage)



##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                 Normalized Count matrix
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
rld.genes <- rlogTransformation(count.genes.deseq) 
saveRDS(rld.genes, "rna/genes.deseq.RDS")

