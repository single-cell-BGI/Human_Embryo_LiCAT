#################################################################################################
#                      Chromatin accessibility level at HCP/ICP/LCP
#################################################################################################
library(GenomicRanges)
library(GenomicAlignments)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(DESeq2)
library(ggplot2)
library(Biostrings)
library(dplyr)
library(reshape2)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#setwd()

stage <- c("Oocyte","Zygote","2-cell","4-cell","8-cell","Morula","ICM","TE")
sample <- paste0(rep(stage, each = 2), "_",rep(c(1,2),8))

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##               Count the read number within gene promoters                  
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Get the locations of all promoters in the hg19 genome
prom <- promoters(TxDb,upstream = 500, downstream = 500) 

# Remove the promoters in sex chromosomes
prom <- prom[seqnames(prom) %in% paste0("chr",1:22)]
saveRDS(prom,"prom.RDS")

# Count number of read within promoter regions
bam.files <-list.files(pattern = "sort.rmdup2.bam$")
ref.peak <- readRDS("prom.RDS")
colnames(ref.peak) <- c("chr","start","end","var1","var2","strand")
ref.peak <- makeGRangesFromDataFrame(ref.peak, keep.extra.columns = T)
counts <- summarizeOverlaps(ref.peak,bam.files,mode="Union",inter.feature=TRUE,fragment=FALSE,singleEnd=TRUE)
##saveRDS(counts,"counts.promoter.RDS")

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##        Normalization of promoter accessibility profiles by DESeq2
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

# Add promoter tags to count matrix
counts.prom <- readRDS("atac/counts.promoter.RDS")
rowData(counts.prom)$tx_id <- prom$tx_id
colnames(counts.prom) <- gsub(".sort.rmdup2.bam","",colnames(counts.prom))
counts.prom$stage <- sapply(strsplit(colnames(counts.prom),"_"),"[",3)
rownames(counts.prom) <- paste(seqnames(counts.prom), start(counts.prom), end(counts.prom),sep = "_")


# construct a deseq datasets
counts.prom.deseq <- DESeqDataSet(counts.prom, ~stage)

# Remove the promoters with low number of read
counts.prom.deseq <- counts.prom.deseq[rowSums(assay(counts.prom.deseq)) > 5,]

# rlog normalization of promoter accessibility
rld.prom <- rlogTransformation(counts.prom.deseq)  
#saveRDS(rld.prom,"atac/prom.deseq.RDS")

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##    Chromatin accessibility level on promoters with different CpG ratio
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

# get sequences of promoter region
prom.view <-Views(BSgenome.Hsapiens.UCSC.hg19,prom)
prom.sequence <- as(prom.view, "DNAStringSet")

# Count the number of "CG" pattern in the promoter regions 
CG.count <- vcountPattern("CG",prom.sequence)
mcols(prom)$CpG.number <- vcountPattern("CG",prom.sequence)

# Divide the promoters into LCP, MCP and HCP based on the rank of number of CpG  
ql <- quantile(mcols(prom)$CpG.number,c(.25,.75))
mcols(prom)$CP <- "ICP"
mcols(prom)$CP[mcols(prom)$CpG.number < ql[1]] <- "LCP"
mcols(prom)$CP[mcols(prom)$CpG.number > ql[2]] <- "HCP"

# Dataframe for ggplot2
rld.prom.df <- as.data.frame(assay(rld.prom))
rld.prom.df$CP <- mcols(prom)$CP[match(rowData(rld.prom)$tx_id,mcols(prom)$tx_id)]

rld.prom.df <- rld.prom.df %>%
  melt(id.vars = "CP",variable.name = "sample",value.name = "expValue") %>%
  mutate(stage = sapply(strsplit(as.character(sample),"_"), "[",3)) 

# plot the promoter accessibility level on promoters with different CpG ratio.
ggplot(rld.prom.df, aes(sample,expValue, fill = factor(CP,levels = c("LCP","ICP","HCP")))) +
  geom_boxplot() + 
  facet_wrap(~factor(stage,levels =
                       c("oocyte","2PN","2C","4C","8C","morula","ICM","TE")),scales="free",nrow = 2)+ xlab("") +
  ylab("Chromatin accessibility level") + 
  labs(fill = "") + 
  ggtitle("Chromatin accessibility") + 
  scale_x_discrete(labels=c("rep1", "rep2")) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle=45, hjust=1))