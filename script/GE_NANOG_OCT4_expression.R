################################################################################
##                 PCA analysis
################################################################################
library(DESeq2)
library(SummarizedExperiment)
library(dplyr)
library(reshape2)

## Set working directory
#setwd()


# Example genes
example_genes <- c("POU5F1","NANOG","GAPDH")
example_genes <- factor(example_genes,levels = example_genes)

# Expression matrix
rna.SEobj <- readRDS("rna/rna.hES.hliver.SEobj.RDS")

example_genes.exp.mat <- as.data.frame(assays(rna.SEobj)$FPKM.log2[match(example_genes,rowData(rna.SEobj)$gene_name),])
example_genes.exp.mat$example_genes <- example_genes

example_genes.exp.mat <- example_genes.exp.mat %>% 
  melt(id.vars = c('example_genes'), variable.name = 'sample', value.name = 'expression')
example_genes.exp.mat$celltype = sapply(strsplit(as.character(example_genes.exp.mat$sample),"-"),"[",1)

# Plot example genes

pdf("pluripotency genes and differentiation genes.pdf")
ggplot(example_genes.exp.mat, aes(sample, expression, fill = celltype)) + 
  geom_col() + 
  facet_wrap(~example_genes,scales="free") +
  xlab("") + 
  labs(fill = "Cell type") +
  theme_bw()
dev.off()