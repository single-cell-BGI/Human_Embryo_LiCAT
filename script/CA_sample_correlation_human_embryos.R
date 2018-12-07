################################################################################
##        Smoothscatter between replicates
################################################################################

library(DESeq2)
library(SummarizedExperiment)
rld.atac <- readRDS("atac/atac.deseq.RDS")

stage.tag <- c("oocyte","2PN","2C","4C","8C","morula","ICM","TE")

countlist <- lapply(setNames(stage.tag,stage.tag), function(X){
  count <- assay(rld.atac)[,grep(X,colnames(rld.atac))]
  return(count)
})

count.rbind.atac <- as.data.frame(do.call(rbind, countlist))
colnames(count.rbind.atac) <- c("rep1", "rep2")
count.rbind.atac$stage <- rep(stage.tag, each = nrow(rld.atac))

par(mfrow = c(2,4))

smoothScatter(count.rbind.atac[count.rbind.atac$stage == "oocyte",1:2])
r = round(cor(count.rbind.atac[count.rbind.atac$stage == "oocyte",1],count.rbind.atac[count.rbind.atac$stage == "oocyte",2]),2)
legend("top",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(count.rbind.atac[count.rbind.atac$stage == "2PN",1:2])
r = round(cor(count.rbind.atac[count.rbind.atac$stage == "2PN",1],count.rbind.atac[count.rbind.atac$stage == "2PN",2]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(count.rbind.atac[count.rbind.atac$stage == "2C",1:2])
r = round(cor(count.rbind.atac[count.rbind.atac$stage == "2C",1],count.rbind.atac[count.rbind.atac$stage == "2C",2]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(count.rbind.atac[count.rbind.atac$stage == "4C",1:2])
r = round(cor(count.rbind.atac[count.rbind.atac$stage == "4C",1],count.rbind.atac[count.rbind.atac$stage == "4C",2]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(count.rbind.atac[count.rbind.atac$stage == "8C",1:2])
r = round(cor(count.rbind.atac[count.rbind.atac$stage == "8C",1],count.rbind.atac[count.rbind.atac$stage == "8C",2]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(count.rbind.atac[count.rbind.atac$stage == "morula",1:2])
r = round(cor(count.rbind.atac[count.rbind.atac$stage == "morula",1],count.rbind.atac[count.rbind.atac$stage == "morula",2]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(count.rbind.atac[count.rbind.atac$stage == "ICM",1:2])
r = round(cor(count.rbind.atac[count.rbind.atac$stage == "ICM",1],count.rbind.atac[count.rbind.atac$stage == "ICM",2]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(count.rbind.atac[count.rbind.atac$stage == "TE",1:2])
r = round(cor(count.rbind.atac[count.rbind.atac$stage == "TE",1],count.rbind.atac[count.rbind.atac$stage == "TE",2]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")