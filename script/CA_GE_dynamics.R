#################################################################################
##              CA and GE dynamics during EGA                                         
#################################################################################
# This script analyze the overall CA and GE dynamics of EGA genes and the dynamics of example genes 

library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(pheatmap)
library(DESeq2)
library(viridis)

#setwd("")
hsapiens.genes <- readRDS("rna/hsapiens.genes.RDS")
hsapiens.genes <- hsapiens.genes[seqnames(hsapiens.genes) %in% paste0("chr",1:22)]


## EGA genes

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                         Z-score of all genes                               ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
rna <- readRDS("rna/genes.deseq.RDS")
assays(rna)$deseq <- assay(rna)
assays(rna)$z <- (assay(rna) - rowMeans(assay(rna)))/rowSds(assay(rna))

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                         Select protein-coding genes                        
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
hsapiens.genes.pro.coding <- hsapiens.genes$gene_id[hsapiens.genes$gene_biotype == "protein_coding"]

rna <- rna[rownames(rna) %in% hsapiens.genes.pro.coding,]

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                         RNA expression matrix z-socre                        
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

rna.mat.zscore <- assays(rna)$z
colnames(rna.mat.zscore) <- sapply(strsplit(colnames(rna.mat.zscore),"_R_"), "[",2)
#saveRDS(rna.mat.zscore,"rna/rna.mat.zscore.RDS")


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                  Mean expression of  all genes                             ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
rna.mean <- do.call(cbind,lapply(setNames(unique(as.character(rna$stage)),unique(as.character(rna$stage))), function(X){
  stage.mean <- rowMeans(assays(rna[,grep(X, colnames(rna))])$z)
  return(stage.mean)
}))
rna.mean <- as.data.frame(rna.mean)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                         Select EGA genes                                   ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

# Define a function that find EGA genes
is.ega <- function(gene) {
  if(gene[1] < 0 & sum(gene[3:7] > 1) >= 1) 
    a <- TRUE
  else if(gene[1] > 0 & sum(gene[3:7] >= 2*gene[1]) >0)
    a = TRUE
  else a = FALSE
  return(a)
}

# Keep EGA genes
keep <- apply(rna.mean, 1, is.ega)
rna.ega <- rna.mean[keep,]
#saveRDS(rna.ega,"rna.ega.RDS")

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                mFUZZ dynamics pattern of genes                             ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

#library(Mfuzz)
# mFuzz analysis
#rna.ega.mf <- ExpressionSet(assayData=as.matrix(rna.ega[,3:6]))
#rna.ega.mf.r <- filter.NA(rna.ega.mf, thres=0.5)
#rna.ega.mf.f <- fill.NA(rna.ega.mf.r,mode="mean")
#tmp <- filter.std(rna.ega.mf.f,min.std=0)
#rna.ega.mf.s <- standardise(rna.ega.mf.f)
#cl.rna.ega <- mfuzz(rna.ega.mf.s,c=6,m=mestimate(rna.ega.mf.s))
#mfuzz.plot(rna.ega.mf.s,cl=cl.rna.ega,mfrow=c(1,3),time.labels=seq(0,30,10))
#saveRDS(cl.rna.ega,"rna/cl.rna.ega.RDS")


## Gene expression dynamics of EGA genes

# Filter genes by membership
cl.rna.ega <- readRDS("rna/cl.rna.ega.RDS")
rna.ega <- readRDS("rna/rna.ega.RDS")
membership <- 0.3
rna.ega.hm <- as.data.frame(rna.ega[rowSums(cl.rna.ega$membership > membership) > 0,])

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                Plot different clusters                             
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
rna.ega.hm$gene_name <- hsapiens.genes$gene_name[match(rownames(rna.ega.hm),hsapiens.genes$gene_id)]
rna.ega.hm$cluster <- cl.rna.ega$cluster[match(rownames(rna.ega.hm),names(cl.rna.ega$cluster))]
rownames(rna.ega.hm) <- NULL

rna.ega.hm.m <- rna.ega.hm %>% 
  dplyr::select("2C", "4C", "8C", "morula","gene_name", "cluster") %>%
  melt(id.vars = c("gene_name","cluster"),variable.name = "stage", value.name = "expression")

center <- as.data.frame(cl.rna.ega$centers)
center$cluster <- rownames(center)


center <- center %>%
  melt(id.var = "cluster",variable.name = "stage", value.name = "expression") 



rna.plot <- 
  ggplot(rna.ega.hm.m, aes(stage,expression,group = gene_name)) +
  geom_line(size = 0.01,color="grey") + 
  facet_wrap(~cluster,ncol = 1) + 
  geom_line(data = center,aes(stage,expression,group = 1),color = "orange") + 
  geom_point(data = center,aes(stage,expression,group = 1),color = "orange") +
  ggtitle("EGA gene expression") +
  ylab("Mean gene expression") +
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_classic()

# Number of EGA genes
table(rna.ega.hm$cluster)



## Chromatin accessibility dynamics of EGA genes 


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                         Load peak2gene                                     ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
up_down <- 10e3
peak2gene <- readRDS(paste0("rna/peak2gene_",as.character(up_down/1000),"Kb.RDS"))
peak2gene.prom <- peak2gene[peak2gene$Promoter,]
peak2gene.dist <- peak2gene[!peak2gene$Promoter,]



##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                         atac data matrix                                   ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
#Load atac and rna data matrices
atac <- readRDS("atac/atac.deseq.RDS")
assays(atac)$counts <- assay(atac)
# Z-score of atac
assays(atac)$z <- (assay(atac) - rowMeans(assay(atac)))/rowSds(assay(atac))

atac.mean <- do.call(cbind,lapply(setNames(unique(as.character(atac$stage)),unique(as.character(atac$stage))), function(X){
  stage.mean <- rowMeans(assays(atac[,grep(X, colnames(atac))])$z)
  return(stage.mean)
}))

atac.mean <- as.data.frame(atac.mean)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                           EGA genes promoter peaks                               
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
atac.mean$gene <- peak2gene.prom$gene_name[match(rownames(atac.mean),peak2gene.prom$peakID)]
atac.mean$cluster <- rna.ega.hm$cluster[match(atac.mean$gene,rna.ega.hm$gene_name)]
atac.ega <- atac.mean[!is.na(atac.mean$cluster),]

atac.ega.m <- atac.ega %>%
  dplyr::select("2C", "4C", "8C", "morula","gene", "cluster") %>%
  melt(id.vars = c("gene","cluster"),variable.name = "stage", value.name = "accessibility") %>%
  group_by(cluster,stage) %>%
  summarise(mean_accessibility = mean(accessibility),median_accessibility = median(accessibility)) %>%
  ungroup()

atac.plot <-
  ggplot(atac.ega.m, aes(stage,mean_accessibility,group = 1)) +
  geom_point(color = "ForestGreen") +
  geom_line(color = "ForestGreen") + 
  facet_wrap(~cluster,ncol = 1) + 
  ylab("Mean chromatin accessibility") +
  xlab("") +
  ggtitle("Accessibility (promoter peaks)")+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_classic()



##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                           EGA genes distal peaks                               
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
atac.mean.dist <- atac.mean
atac.mean.dist$gene <- peak2gene.dist$gene_name[match(rownames(atac.mean),peak2gene.dist$peakID)]
atac.mean.dist$cluster <- rna.ega.hm$cluster[match(atac.mean.dist$gene,rna.ega.hm$gene_name)]
atac.dist.ega <- atac.mean.dist[!is.na(atac.mean.dist$cluster),]

table(atac.dist.ega$cluster)

atac.dist.ega.m <- atac.dist.ega %>%
  dplyr::select("2C", "4C", "8C", "morula","gene", "cluster") %>%
  melt(id.vars = c("gene","cluster"),variable.name = "stage", value.name = "accessibility") %>%
  group_by(cluster,stage) %>%
  summarise(mean_accessibility = mean(accessibility),median_accessibility = median(accessibility)) %>%
  ungroup()

atac.dist.plot <-
  ggplot(atac.dist.ega.m, aes(stage,mean_accessibility,group = 1)) +
  geom_point(color = "ForestGreen") +
  geom_line(color = "ForestGreen") + 
  facet_wrap(~cluster,ncol = 1) + 
  ylab("Mean chromatin accessibility") +
  xlab("") +
  ggtitle("Accessibility (distal peaks)")+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_classic()

pdf("atac.rna.assoc.pdf")
grid.arrange(rna.plot,atac.plot,atac.dist.plot,nrow = 1)
dev.off()

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                A table integrates atac and rna for EGA genes                                     
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
atac.ega$peak_type <- "promoter"
atac.dist.ega$peak_type <- "distal"
atac.ega.all <- rbind(atac.ega,atac.dist.ega)
atac.ega.all$Peak_chr <- sapply(strsplit(rownames(atac.ega.all),"_"),"[",1)
atac.ega.all$Peak_start <- sapply(strsplit(rownames(atac.ega.all),"_"),"[",2)
atac.ega.all$Peak_end <- sapply(strsplit(rownames(atac.ega.all),"_"),"[",3)
colnames(atac.ega.all) <- c("Accessibility_2-cell", "Accessibility_Zygote", "Accessibility_4-cell", "Accessibility_8-cell", "Accessibility_ICM", "Accessibility_Morula", "Accessibility_Oocyte", "Accessibility_TE", "Gene", 
                            "Cluster", "Peak_type", "Peak_chr", "Peak_start", "Peak_end")
atac.ega.all <-atac.ega.all[c("Peak_chr", "Peak_start", "Peak_end","Peak_type","Gene","Cluster","Accessibility_2-cell", "Accessibility_4-cell", "Accessibility_8-cell", "Accessibility_Morula")]

rna.ega.hm <- rna.ega.hm[c("2C", "4C", "8C", "morula","gene_name")]
colnames(rna.ega.hm) <- c("Expression_2-cell", "Expression_4-cell", "Expression_8-cell", "Expression_Morula","Gene")

atac_rna.ega_hm <- merge(atac.ega.all,rna.ega.hm,by = "Gene", all.x = T)

#write.table(atac_rna.ega_hm,
#            "Supplementary Table 2 Accessibility and gene expression level of EGA genes.xls",
#            col.names = T,
#            row.names = F,
#            sep = "\t",
#            quote = F)



##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                 Plot rna expression of example genes           
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

rna.z <- readRDS("rna/rna.mat.zscore.RDS")

rna.z.ega <- as.data.frame(rna.z) %>%
  dplyr::select(colnames(rna.z)[grep("2C|4C|8C|morula",colnames(rna.z))]) %>%
  mutate(gene_name = hsapiens.genes$gene_name[match(rownames(rna.z),hsapiens.genes$gene_id)])


plot_rna <- function(genelist) {
  rna.list <- rna.z.ega[match(genelist,rna.z.ega$gene_name),]
  rna.list <- rna.list %>% 
    melt(id.vars = "gene_name",variable.name = "sample",value.name = "expValue") %>%
    mutate(stage = sapply(strsplit(as.character(sample),"_"),"[",1)) %>%
    group_by(stage,gene_name) %>%
    summarise(mean = mean(expValue), sd = sd(expValue)) %>%
    ungroup()
  
  p <- ggplot(rna.list,aes(factor(stage,levels = c("2C","4C","8C","morula")),mean,group = 1)) +
    geom_point(color = "Orange") +
    geom_line(color = "Orange") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9),color = "Orange") +
    facet_wrap(~gene_name,scales="free",ncol = 1) +
    theme_classic() +
    xlab("") +
    ylab("Gene expression")
  return(p)
}



genelist <- c("TEAD4","ARID3B","HDAC2","LIN28A","ZSCAN5B", "SLC30A7")
pdf("example.rna.main.pdf",width = 2.5,height = 12)
plot_rna(genelist)
dev.off()

genelist <- c("WNT5B","CHD3","KLF9","TET1","KDM4E", "WNT9A")
pdf("example.rna.supp.pdf",width = 2.5,height = 12)
plot_rna(genelist)
dev.off()

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                 Plot atac expression of example genes           
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

atac.z <- readRDS("atac/atac.zscore.mat.RDS")

#atac.z <- atac.z[rownames(atac.z) %in% peak2gene.prom$peakID,]

atac.z.ega <- as.data.frame(atac.z) %>%
  dplyr::select(colnames(atac.z)[grep("2C|4C|8C|morula|gene_name",colnames(atac.z))])



plot_atac <- function(peak2genelist) {
  atac.list <- atac.z.ega[rownames(atac.z.ega) %in% peak2genelist$peakID,]
  atac.list <- atac.list %>% 
    mutate(gene_name = peak2genelist$gene[match(rownames(atac.list),peak2genelist$peakID)]) %>%
    melt(id.vars = "gene_name",variable.name = "sample",value.name = "AccValue") %>%
    mutate(stage = sapply(strsplit(as.character(sample),"_"),"[",1)) %>%
    group_by(stage,gene_name) %>%
    summarise(mean = mean(AccValue), sd = sd(AccValue)) %>%
    ungroup()
  
  p <- ggplot(atac.list,aes(factor(stage,levels = c("2C","4C","8C","morula")),mean,group = 1)) +
    geom_point(color = "ForestGreen") +
    geom_line(color = "ForestGreen") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9),color = "ForestGreen") +
    facet_wrap(~gene_name,scales="free",ncol = 1) +
    theme_classic() +
    xlab("") +
    ylab("Chromatin accessibility")
  return(p)
}


peak2genelist <- data.frame(gene =  factor(c("TEAD4",
                                             "ARID3B",
                                             "HDAC2",
                                             "LIN28A",
                                             "ZSCAN5B",
                                             "SLC30A7"),
                                           levels = c("TEAD4",
                                                      "ARID3B",
                                                      "HDAC2",
                                                      "LIN28A",
                                                      "ZSCAN5B",
                                                      "SLC30A7")),
                            peakID = c("chr12_3094752_3095941",
                                       "chr15_74832598_74834318",
                                       "chr6_114332595_114333069",
                                       "chr1_26743372_26744951",
                                       "chr19_56707769_56708623",
                                       "chr1_101437701_101438111"))



pdf("example.atac.main.pdf",width = 2.5,height = 12)
plot_atac(peak2genelist)
dev.off()


peak2genelist <- data.frame(gene =  factor(c("WNT5B",
                                             "CHD3",
                                             "KLF9",
                                             "TET1",
                                             "KDM4E",
                                             "WNT9A"),
                                           levels = c("WNT5B",
                                                      "CHD3",
                                                      "KLF9",
                                                      "TET1",
                                                      "KDM4E",
                                                      "WNT9A")),
                            peakID = c("chr12_1730039_1730536",
                                       "chr17_7816185_7816687",
                                       "chr9_73025421_73028150",
                                       "chr10_70391018_70392130",
                                       "chr11_94760410_94761417",
                                       "chr1_101437701_101438111"))

pdf("example.atac.supp.pdf",width = 2.5,height = 12)
plot_atac(peak2genelist)
dev.off()


## Dynamics of DNMT and TET


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
    facet_wrap(~gene_name,scales="free",nrow = 1) +
    theme_classic() +
    xlab("") +
    ylab("Gene expression")
  return(p)
}

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
    facet_wrap(~gene_name,scales="free",nrow = 1) +
    theme_classic() +
    xlab("") +
    ylab("Chromatin accessibility")
  return(p)
}


genelist4 <- c("TET1","TET2","TET3","DNMT3A","DNMT3B","DNMT1","DNMT3L","TDG")

pdf("DNA methytransferase and demethylase.rna.pdf",width = 14,height = 2)
plot_rna_all(genelist4)
dev.off()


peak2genelist4 <- data.frame(gene = c("DNMT1","DNMT3A","DNMT3B","DNMT3L","TET1","TET2","TET3","TDG"),
                             peakID = c("chr19_10337457_10338533","chr2_25537028_25537528","chr20_31366154_31366887","chr21_45674272_45676265","chr10_70391018_70392130","chr4_106124369_106124879","chr2_74236997_74237861","chr12_104350179_104351613"))

pdf("DNA methytransferase and demethylase.atac.pdf",width = 14,height = 2)
plot_atac_all(peak2genelist4)
dev.off()
