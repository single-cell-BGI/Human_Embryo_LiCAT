###############################################################################
##        Gene expression correlation between replicates
###############################################################################
library(DESeq2)
library(SummarizedExperiment)

## working directory
#setwd("")

rld.genes <- readRDS("rna/genes.deseq.RDS")

stage.tag <- c("oocyte","2PN","2C","4C","8C","morula","ICM","TE")

countlist <- lapply(setNames(stage.tag,stage.tag), function(X){
  count <- assay(rld.genes)[,grep(X,colnames(rld.genes))]
  return(count)
})

count.rbind.genes <- as.data.frame(do.call(rbind, countlist))
colnames(count.rbind.genes) <- c("rep1", "rep2")
count.rbind.genes$stage <- rep(stage.tag, each = nrow(rld.genes))

par(mfrow = c(2,4))
smoothScatter(count.rbind.genes[count.rbind.genes$stage == "oocyte",1:2])
r = round(cor(count.rbind.genes[count.rbind.genes$stage == "oocyte",1],count.rbind.genes[count.rbind.genes$stage == "oocyte",2]),2)
legend("top",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(count.rbind.genes[count.rbind.genes$stage == "2PN",1:2])
r = round(cor(count.rbind.genes[count.rbind.genes$stage == "2PN",1],count.rbind.genes[count.rbind.genes$stage == "2PN",2]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(count.rbind.genes[count.rbind.genes$stage == "2C",1:2])
r = round(cor(count.rbind.genes[count.rbind.genes$stage == "2C",1],count.rbind.genes[count.rbind.genes$stage == "2C",2]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(count.rbind.genes[count.rbind.genes$stage == "4C",1:2])
r = round(cor(count.rbind.genes[count.rbind.genes$stage == "4C",1],count.rbind.genes[count.rbind.genes$stage == "4C",2]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(count.rbind.genes[count.rbind.genes$stage == "8C",1:2])
r = round(cor(count.rbind.genes[count.rbind.genes$stage == "8C",1],count.rbind.genes[count.rbind.genes$stage == "8C",2]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(count.rbind.genes[count.rbind.genes$stage == "morula",1:2])
r = round(cor(count.rbind.genes[count.rbind.genes$stage == "morula",1],count.rbind.genes[count.rbind.genes$stage == "morula",2]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(count.rbind.genes[count.rbind.genes$stage == "ICM",1:2])
r = round(cor(count.rbind.genes[count.rbind.genes$stage == "ICM",1],count.rbind.genes[count.rbind.genes$stage == "ICM",2]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(count.rbind.genes[count.rbind.genes$stage == "TE",1:2])
r = round(cor(count.rbind.genes[count.rbind.genes$stage == "TE",1],count.rbind.genes[count.rbind.genes$stage == "TE",2]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")
