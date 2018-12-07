#################################################################################
##          Known Motif enrichment in peaks gained at each stage                                        
#################################################################################
# The p-value of each motif was caculated by Homer. 
library(ggplot2)

#setwd("")
motif.sig <- readRDS("atac/motif.pval.expression.RDS")

# Select highly enriched and representative motifs
motif.sig.factor <- c("DUX4","Otx2","NFY","Fli1","ETS","GABPA","Elk4","ELF1","Elk1","Sp1","Sp5","GSC","CRX","Pitx1","Klf9","KLF14","KLF3","KLF5","KLF6","Klf4","CTCF","Gata4","GATA3","Gata2","Gata6","Oct4","OCT4-SOX2-TCF-NANOG","TEAD","TEAD2","TEAD4","Jun-AP1","JunB","Fosl2","Bach2","Atf3","GRHL2","EBF1")
motif.sig <- motif.sig[motif.sig$motif %in% motif.sig.factor,]

# Plot motif accessibility and the corresponding TF expression
pdf("motif.pdf", width = 10, height = 4)
ggplot(motif.sig,aes(factor(motif,levels = motif.sig.factor), 
                     factor(group,levels = c("TE","ICM","morula","8Cell","4Cell","2Cell")))) + 
  geom_point(aes(size = logP,color = exp.level))+ 
  labs(size = "accessibility")  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 10)) + 
  ylab("") + 
  xlab("") + 
  ggtitle("Motif") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(terrain.colors(2)))
dev.off()
