#################################################################################
##              Chromatin accessibility dynamics for ERVs                                         
#################################################################################

library(pheatmap)
library(ggplot2)

#save()

# Load rpkm value for chromatin accessibility profiles of ERVs
erv.rpkm <- readRDS("atac/erv.atac.rpkm.RDS")

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


erv.atac.top <- erv.rpkm[order(entropy_value)[1:50],]
#erv.atac.top <- erv.atac.top[!rownames(erv.atac.top) %in% c("MLT1M-int","HERV1_LTRb","HERV-Fc1_LTR1","HERV-Fc1_LTR3","HERV-Fc2_LTR","HERV1_LTRc","HERV-Fc1_LTR2","L1P4c"),] 

pdf("heatmap.atac.erv.pdf",width = 3,height = 4)
pp <- pheatmap(erv.atac.top,scale = "row", 
               cluster_rows = T, 
               cluster_cols = F, 
               color = colorRampPalette(c("SlateBlue", "white", "red"))(50),
               fontsize_row = 3, clustering_method = "centroid")
dev.off()


# barplot for all top ERVs

colnames(erv.atac.top) <- c("Oocyte","2-PN","2-Cell","4-Cell","8-Cell","Morula","ICM","TE")
erv.atac.top1 <- as.data.frame(t(erv.atac.top))
erv.atac.top1$stage <- factor(rownames(erv.atac.top1),levels = c("Oocyte","2-PN","2-Cell","4-Cell","8-Cell","Morula","ICM","TE"))
erv.atac.top1 <- melt(erv.atac.top1,id = "stage")
pdf("CA.for.tops.ERVs.pdf", width = 12,height = 12)
ggplot(erv.atac.top1, aes(x = stage, y = value)) + 
  geom_col() + 
  facet_wrap(~variable, scales = "free_y",nrow = 8) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 6)) + 
  ylab("RPKM")+ 
  ggtitle("Chromatin accessibility")
dev.off()