################################################################################
##                   Add DNAse data
################################################################################

library(DESeq2)
library(AnnotationHub)
library(GenomicRanges)

## set working directory


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##              Promoters divided by different GC content
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
prom <- readRDS("atac/prom.GC.content.HCP.MCP.LCP.gr.RDS")

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                   Add DNAse data
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

epiFiles <- query(ah, c( "dnase", "H1"))
dnase.h1 <- epiFiles[["AH22507"]]

overlap <- findOverlaps(prom, dnase.h1)
prom$dnase <- FALSE; prom$dnase[unique(queryHits(overlap))] <- TRUE


# Dataframe for ggplot2
rld.prom <- readRDS("atac/counts.prom.hES.hLiver.deseq.RDS")
rld.prom.df <- as.data.frame(assay(rld.prom))
rld.prom.df$CP <- mcols(prom)$CP[match(rowData(rld.prom)$gene_id,prom$name)]
rld.prom.df$dnase <- mcols(prom)$dnase[match(rowData(rld.prom)$gene_id,prom$name)]
colnames(rld.prom.df)[1:4] <-  c("hESC-1","hESC-2","hLiver-1","hliver-2")

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                    Plot
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

rld.prom.df <- rld.prom.df %>%
  melt(id.vars = c("CP","dnase"),variable.name = "sample",value.name = "expValue") 


# plot the promoter accessibility level on promoters with different CpG ratio.

p <- ggplot(rld.prom.df, aes(dnase,expValue, fill = factor(dnase))) +
  geom_boxplot() + 
  facet_wrap(~CP + sample)+ 
  xlab("") +
  ylab("Chromatin accessibility level") + 
  scale_fill_viridis(discrete=TRUE) + 
  labs(fill = "DNase") + 
  ggtitle("Chromatin accessibility") + 
  scale_x_discrete(labels=c("rep1", "rep2")) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle=45, hjust=1))

pdf("GC.bias.pdf")
p 
dev.off()
