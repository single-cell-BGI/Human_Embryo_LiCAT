#################################################################################
##              DUX4 target gene expression                                         
#################################################################################
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(DESeq2)

#setwd("")


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                 Plot rna expression of example genes           
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

rna.z <- readRDS("rna/rna.mat.zscore.RDS")

rna.z.all <- as.data.frame(rna.z) %>%
  dplyr::select(colnames(rna.z)[grep("oocyte|2PN|2C|4C|8C|morula|ICM|TE",colnames(rna.z))]) %>%
  mutate(gene_name = hsapiens.genes$gene_name[match(rownames(rna.z),hsapiens.genes$gene_id)])


plot_rna_all <- function(genelist) {
  rna.list <- rna.z.all[match(genelist,rna.z.all$gene_name),]
  rna.list <- rna.list %>% 
    melt(id.vars = "gene_name",variable.name = "sample",value.name = "expValue") %>%
    mutate(stage = sapply(strsplit(as.character(sample),"_"),"[",1)) %>%
    group_by(stage,gene_name) %>%
    summarise(mean = mean(expValue), sd = sd(expValue)) %>%
    ungroup()
  
  p <- ggplot(rna.list,aes(factor(stage,levels = c("oocyte","2PN","2C","4C","8C","morula","ICM","TE")),mean,group = 1)) +
    geom_point(color = "Orange") +
    geom_line(color = "Orange") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9),color = "Orange") +
    facet_wrap(~gene_name,scales="free",ncol = 1) +
    theme_classic() +
    xlab("") +
    ylab("Gene expression")
  return(p)
}

genelist <- "DUX4"
genelist1 <- c("KDM4E","ZSCAN4","LEUTX")
genelist2 <- c("DUSP22","ZMYND12","PDZD7")
genelist3 <- c("INO80E","CBX6","BZW2")


pdf("DUX4.rna.pdf",width = 2.5,height = 2)
plot_rna_all(genelist)
dev.off()

pdf("DUX4.target.rna.group1.pdf",width = 2.5,height = 6)
plot_rna_all(genelist1)
dev.off()

pdf("DUX4.target.rna.group2.pdf",width = 2.5,height = 6)
plot_rna_all(genelist2)
dev.off()

pdf("DUX4.target.rna.group3.pdf",width = 2.5,height = 6)
plot_rna_all(genelist3)
dev.off()



##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                 Plot chromatin accessibility of example genes           
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

atac.z <- readRDS("atac/atac.zscore.mat.RDS")

#atac.z <- atac.z[rownames(atac.z) %in% peak2gene.prom$peakID,]

atac.z.all <- as.data.frame(atac.z) %>%
  dplyr::select(colnames(atac.z)[grep("oocyte|2PN|2C|4C|8C|morula|ICM|TE|gene_name",colnames(atac.z))])



plot_atac_all <- function(peak2genelist) {
  atac.list <- atac.z.all[rownames(atac.z.all) %in% peak2genelist$peakID,]
  atac.list <- atac.list %>% 
    mutate(gene_name = peak2genelist$gene[match(rownames(atac.list),peak2genelist$peakID)]) %>%
    melt(id.vars = "gene_name",variable.name = "sample",value.name = "AccValue") %>%
    mutate(stage = sapply(strsplit(as.character(sample),"_"),"[",1)) %>%
    group_by(stage,gene_name) %>%
    summarise(mean = mean(AccValue), sd = sd(AccValue)) %>%
    ungroup()
  
  p <- ggplot(atac.list,aes(factor(stage,levels = c("oocyte","2PN","2C","4C","8C","morula","ICM","TE")),mean,group = 1)) +
    geom_point(color = "ForestGreen") +
    geom_line(color = "ForestGreen") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9),color = "ForestGreen") +
    facet_wrap(~gene_name,scales="free",ncol = 1) +
    theme_classic() +
    xlab("") +
    ylab("Chromatin accessibility")
  return(p)
}

peak2genelist1 <- data.frame(gene = c("ZSCAN4","LEUTX","KDM4E"),
                             peakID = c("chr19_58189292_58190208","chr19_40262355_40263969","chr11_94760410_94761417"))

peak2genelist2 <- data.frame(gene = c("DUSP22","ZMYND12","PDZD7"),
                             peakID = c("chr6_324291_325057","chr1_42915032_42915780","chr10_102789324_102790084"))

peak2genelist3 <- data.frame(gene = c("INO80E","CBX6","BZW2"),
                             peakID = c("chr16_30020207_30020829","chr22_39273268_39273771","chr7_16732217_16732827"))





pdf("DUX4.target.atac.group1.pdf",width = 2.5,height = 6)
plot_atac_all(peak2genelist1)
dev.off()

pdf("DUX4.target.atac.group2.pdf",width = 2.5,height = 6)
plot_atac_all(peak2genelist2)
dev.off()

pdf("DUX4.target.atac.group3.pdf",width = 2.5,height = 6)
plot_atac_all(peak2genelist3)
dev.off()

