################################################################################
##        Genes denisty within each genomic windows
################################################################################

count.2M <- readRDS("atac/counts.2M.RDS")

count.2M <- count.2M[grep("chr\\d{1,2}",rownames(count.2M)),]
assays(count.2M)$normalizedCounts <- sapply(as.data.frame(assays(count.2M)$counts),function(column) column/mean(column))

pdf("geneCount VS read count.pdf",width = 18, height = 5)
par(mfrow = c(2,9))

smoothScatter(rowData(count.2M)$geneNumber, assays(count.2M)$normalizedCounts[,grep("sperm_1",colnames(count.2M))],nrpoints = 0,ylim = c(0,4),xlab = "", ylab = "",main = "Sperm rep1")
r = round(cor(rowData(count.2M)$geneNumber,assays(count.2M)$normalizedCounts[,grep("sperm_1",colnames(count.2M))]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(rowData(count.2M)$geneNumber, assays(count.2M)$normalizedCounts[,grep("oocyte_1",colnames(count.2M))],nrpoints = 0,ylim = c(0,4),xlab = "", ylab = "",main = "Oocyte rep1")
r = round(cor(rowData(count.2M)$geneNumber,assays(count.2M)$normalizedCounts[,grep("oocyte_1",colnames(count.2M))]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(rowData(count.2M)$geneNumber, assays(count.2M)$normalizedCounts[,grep("2PN_1",colnames(count.2M))],nrpoints = 0,ylim = c(0,4),xlab = "", ylab = "",main = "2PN rep1")
r = round(cor(rowData(count.2M)$geneNumber,assays(count.2M)$normalizedCounts[,grep("2PN_1",colnames(count.2M))]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(rowData(count.2M)$geneNumber, assays(count.2M)$normalizedCounts[,grep("2C_1",colnames(count.2M))],nrpoints = 0,ylim = c(0,4),xlab = "", ylab = "",main = "2C rep1")
r = round(cor(rowData(count.2M)$geneNumber,assays(count.2M)$normalizedCounts[,grep("2C_1",colnames(count.2M))]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(rowData(count.2M)$geneNumber, assays(count.2M)$normalizedCounts[,grep("4C_1",colnames(count.2M))],nrpoints = 0,ylim = c(0,4),xlab = "", ylab = "",main = "4C rep1")
r = round(cor(rowData(count.2M)$geneNumber,assays(count.2M)$normalizedCounts[,grep("4C_1",colnames(count.2M))]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(rowData(count.2M)$geneNumber, assays(count.2M)$normalizedCounts[,grep("8C_1",colnames(count.2M))],nrpoints = 0,ylim = c(0,4),xlab = "", ylab = "",main = "8C rep1")
r = round(cor(rowData(count.2M)$geneNumber,assays(count.2M)$normalizedCounts[,grep("8C_1",colnames(count.2M))]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(rowData(count.2M)$geneNumber, assays(count.2M)$normalizedCounts[,grep("morula_1",colnames(count.2M))],nrpoints = 0,ylim = c(0,4),xlab = "", ylab = "",main = "morula rep1")
r = round(cor(rowData(count.2M)$geneNumber,assays(count.2M)$normalizedCounts[,grep("morula_1",colnames(count.2M))]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(rowData(count.2M)$geneNumber, assays(count.2M)$normalizedCounts[,grep("ICM_1",colnames(count.2M))],nrpoints = 0,ylim = c(0,4),xlab = "", ylab = "",main = "ICM rep1")
r = round(cor(rowData(count.2M)$geneNumber,assays(count.2M)$normalizedCounts[,grep("ICM_1",colnames(count.2M))]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(rowData(count.2M)$geneNumber, assays(count.2M)$normalizedCounts[,grep("TE_1",colnames(count.2M))],nrpoints = 0,ylim = c(0,4),xlab = "", ylab = "",main = "TE rep1")
r = round(cor(rowData(count.2M)$geneNumber,assays(count.2M)$normalizedCounts[,grep("TE_1",colnames(count.2M))]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(rowData(count.2M)$geneNumber, assays(count.2M)$normalizedCounts[,grep("sperm_2",colnames(count.2M))],nrpoints = 0,ylim = c(0,4),xlab = "", ylab = "",main = "Sperm rep2")
r = round(cor(rowData(count.2M)$geneNumber,assays(count.2M)$normalizedCounts[,grep("sperm_2",colnames(count.2M))]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(rowData(count.2M)$geneNumber, assays(count.2M)$normalizedCounts[,grep("oocyte_2",colnames(count.2M))],nrpoints = 0,ylim = c(0,4),xlab = "", ylab = "",main = "Oocyte rep2")
r = round(cor(rowData(count.2M)$geneNumber,assays(count.2M)$normalizedCounts[,grep("oocyte_2",colnames(count.2M))]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(rowData(count.2M)$geneNumber, assays(count.2M)$normalizedCounts[,grep("2PN_2",colnames(count.2M))],nrpoints = 0,ylim = c(0,4),xlab = "", ylab = "",main = "2PN rep2")
r = round(cor(rowData(count.2M)$geneNumber,assays(count.2M)$normalizedCounts[,grep("2PN_2",colnames(count.2M))]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(rowData(count.2M)$geneNumber, assays(count.2M)$normalizedCounts[,grep("2C_2",colnames(count.2M))],nrpoints = 0,ylim = c(0,4),xlab = "", ylab = "",main = "2C rep2")
r = round(cor(rowData(count.2M)$geneNumber,assays(count.2M)$normalizedCounts[,grep("2C_2",colnames(count.2M))]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(rowData(count.2M)$geneNumber, assays(count.2M)$normalizedCounts[,grep("4C_2",colnames(count.2M))],nrpoints = 0,ylim = c(0,4),xlab = "", ylab = "",main = "4C rep2")
r = round(cor(rowData(count.2M)$geneNumber,assays(count.2M)$normalizedCounts[,grep("4C_2",colnames(count.2M))]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(rowData(count.2M)$geneNumber, assays(count.2M)$normalizedCounts[,grep("8C_2",colnames(count.2M))],nrpoints = 0,ylim = c(0,4),xlab = "", ylab = "",main = "8C rep2")
r = round(cor(rowData(count.2M)$geneNumber,assays(count.2M)$normalizedCounts[,grep("8C_2",colnames(count.2M))]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(rowData(count.2M)$geneNumber, assays(count.2M)$normalizedCounts[,grep("morula_2",colnames(count.2M))],nrpoints = 0,ylim = c(0,4),xlab = "", ylab = "",main = "morula rep2")
r = round(cor(rowData(count.2M)$geneNumber,assays(count.2M)$normalizedCounts[,grep("morula_2",colnames(count.2M))]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(rowData(count.2M)$geneNumber, assays(count.2M)$normalizedCounts[,grep("ICM_2",colnames(count.2M))],nrpoints = 0,ylim = c(0,4),xlab = "", ylab = "",main = "ICM rep2")
r = round(cor(rowData(count.2M)$geneNumber,assays(count.2M)$normalizedCounts[,grep("ICM_2",colnames(count.2M))]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")

smoothScatter(rowData(count.2M)$geneNumber, assays(count.2M)$normalizedCounts[,grep("TE_2",colnames(count.2M))],nrpoints = 0,ylim = c(0,4),xlab = "", ylab = "",main = "TE rep2")
r = round(cor(rowData(count.2M)$geneNumber,assays(count.2M)$normalizedCounts[,grep("TE_2",colnames(count.2M))]),2)
legend("topleft",legend = paste("Cor =", r ,sep = ""),bty = "n")
dev.off()
