#################################################################################
##              Gene expression dynamics for ERVs                                         
#################################################################################

library(pheatmap)
library(ggplot2)

#save()

# Load rpkm value for gene expression profiles of ERVs
erv.rpkm <- readRDS("rna/erv.rna.rpkm.RDS")

# Entropy

entropy_rela <- sapply(erv.rpkm, function(column) {
  rela <- column/rowSums(erv.rpkm)
})

entropy <- function(a){
  a[a == 0] <- 0.000000001
  b <- -sum(a*log2(a))
  return(b)
}

entropy_value <- apply(entropy_rela,1,entropy)


erv.rna.top <- erv.rpkm[order(entropy_value)[1:50],]

# Heatmap showing genes with highest entropy value
pdf("heatmap.rna.erv.pdf",width = 3,height = 4)
pp <- pheatmap(erv.rna.top,scale = "row", 
               cluster_rows = T, 
               cluster_cols = F, 
               color = colorRampPalette(c("SlateBlue", "white", "red"))(50),
               fontsize_row = 3, clustering_method = "centroid")
dev.off()


# barplot for all top ERVs

colnames(erv.rna.top) <- c("Oocyte","2-PN","2-Cell","4-Cell","8-Cell","Morula","ICM","TE")
erv.rna.top1 <- as.data.frame(t(erv.rna.top))
erv.rna.top1$stage <- factor(rownames(erv.rna.top1),levels = c("Oocyte","2-PN","2-Cell","4-Cell","8-Cell","Morula","ICM","TE"))
erv.rna.top1 <- melt(erv.rna.top1,id = "stage")

pdf("GE.for.tops.ERVs.pdf", width = 12,height = 12)
ggplot(erv.rna.top1, aes(x = stage, y = value)) + 
  geom_col() + 
  facet_wrap(~variable, scales = "free_y",nrow = 8) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 6)) + 
  ylab("RPKM")+ 
  ggtitle("Gene expression")
dev.off()