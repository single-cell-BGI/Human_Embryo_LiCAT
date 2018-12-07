################################################################################
##                 Number of genes deteced
################################################################################
library(dplyr)
library(reshape2)
library(ggplot2)
library(viridis)

#setwd("")

rld.genes <- readRDS("rna/genes.deseq.RDS")

number.genes <- data.frame(
  sample = colnames(rld.genes),
  number = colSums(assay(rld.genes) > 1),
  stages = sapply(strsplit(colnames(rld.genes),"_"), "[",3))

number.genes <- as.data.frame(assay(rld.genes))
number.genes$gene_biotype <- rowData(rld.genes)$gene_biotype

number.genes <- number.genes %>% 
  melt(id.vars = "gene_biotype",variable.name = "sample",value.name = "expValue") %>%
  mutate(stage = sapply(strsplit(as.character(sample),"_"), "[",3)) %>%
  group_by(sample,gene_biotype,stage) %>%
  summarise(number = sum(expValue > 1)) %>%
  ungroup() %>%
  group_by(stage,gene_biotype) %>%
  summarise(mean = mean(number), sd = sd(number)) %>%
  ungroup() %>%
  dplyr::filter(gene_biotype %in% c("protein_coding","antisense","pseudogene","lincRNA"))

pdf("number.of.genes.pdf")
ggplot(number.genes, 
       aes(factor(stage,levels = c("oocyte","2PN","2C","4C","8C","morula","ICM","TE")), 
           mean,
           fill = factor(gene_biotype,levels = c("lincRNA","antisense","pseudogene","protein_coding")))) + geom_bar(position=position_dodge(), stat="identity") + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9)) + 
  ggtitle("Number of genes") + 
  scale_fill_viridis(discrete=TRUE) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle=45, hjust=1)) + 
  xlab("") + 
  ylab("Number of genes") + 
  labs(fill = "")
dev.off()