#################################################################################################
#                      Number of peaks at each stage
#################################################################################################

library(GenomicRanges)
library(ggplot2)

setwd("~/Documents/Research/embryo/2_revision/Code/Code_submission")

stage <- c("Sperm","Oocyte","Zygote","2-cell","4-cell","8-cell","Morula","ICM","TE")
sample <- paste0(rep(stage, each = 2), "_",rep(c(1,2),9))

##---------------------------------------------------------------------------##
##                              Import peak files                                    
##---------------------------------------------------------------------------##

# get path of peaks
peak.path <- file.path("atac","bed.atac",list.files("atac/NarrowPeak/",pattern = ".narrowPeak$"))

# import peak data into a list
peak.all <- lapply(peak.path, function(bed) {
  peak<- read.table(bed, header = F,stringsAsFactors = F)
  colnames(peak) <- c("chr","start","end","name","score","strand","fold-change","-log10pvalue","-log10qvalue","summit")
  return(peak)
})

names(peak.all) <- gsub("sort.rmdup2.bam_peaks.","",list.files("atac/bed.atac",pattern = ".narrowPeak$"))

##---------------------------------------------------------------------------##
##                  Make GRangesList from peaks                                   
##---------------------------------------------------------------------------##

peak.all.df <- do.call(rbind,peak.all)
peak.all.df$sample <- sapply(strsplit(rownames(peak.all.df),"\\."),"[",1)

peak.grl <- makeGRangesListFromDataFrame(peak.all.df,split.field = "sample",keep.extra.columns = TRUE)

##---------------------------------------------------------------------------##
##                  Remove peaks in sex chromosomes                                  
##---------------------------------------------------------------------------##

peak.grl.rmXY <- lapply(peak.grl, function(x) x[seqnames(x) %in% paste0("chr",1:22)])

# Re-order the samples 
peak.grl.rmXY <- peak.grl.rmXY[c("Embryo_A_sperm_1","Embryo_A_sperm_2","Embryo_A_oocyte_1", 
                                 "Embryo_A_oocyte_2","Embryo_A_2PN_1", "Embryo_A_2PN_2","Embryo_A_2C_1", "Embryo_A_2C_2",  "Embryo_A_4C_1", "Embryo_A_4C_2", "Embryo_A_8C_1", "Embryo_A_8C_2", 
                                 "Embryo_A_morula_1", "Embryo_A_morula_2", 
                                 "Embryo_A_ICM_1", "Embryo_A_ICM_2", "Embryo_A_TE_1", "Embryo_A_TE_2")]

##---------------------------------------------------------------------------##
##            Find peaks that are present in both replicates                        
##---------------------------------------------------------------------------##
rep.name <- unique(sapply(strsplit(names(peak.grl.rmXY),"_"), "[",3))

peak.merge <- lapply(setNames(rep.name,rep.name), function(X){
  stage.grl <- peak.grl.rmXY[grep(X, names(peak.grl.rmXY))]
  stage.peak <- findOverlaps(stage.grl[[1]],stage.grl[[2]])
  peak1 <- stage.grl[[1]][queryHits(stage.peak)]
  peak2 <- stage.grl[[2]][subjectHits(stage.peak)]
  peak <- reduce(c(peak1, peak2))
  return(peak)
})
names(peak.merge) <- stage

peakfiles <-  file.path("atac","bed.atac","merge",list.files("atac/bed.atac/merge"))

peakfiles <- as.list(peakfiles)
names(peakfiles) <- list.files("atac/bed.atac/merge")

# Number of peaks
elementNROWS(peak.merge)

number.peak.both <- data.frame(stage = stage,number = elementNROWS(peak.merge)
)


pdf("number.of.peaks.in.both.replicates.pdf")
ggplot(number.peak.both, aes(factor(stage,levels = stage), number)) + 
  geom_col() + 
  ggtitle("Number of peaks") + 
  xlab("") +
  ylab("Number of peaks") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle=45, hjust=1)) 
dev.off()

##---------------------------------------------------------------------------##
##                  Gained peaks at each stage                                  
##---------------------------------------------------------------------------##

gain.peak <- lapply(setNames(3:9,stage[3:9]), function(X){
  peak1 <- peak.merge[[X]]
  peak2 <- peak.merge[[X - 1]]
  peak.gain <- setdiff(peak1,peak2)
  return(peak.gain)
})

# Number of peaks
elementNROWS(gain.peak)

##---------------------------------------------------------------------------##
##                  Lost peaks at each stage                                  
##---------------------------------------------------------------------------##

lost.peak <- lapply(setNames(3:9,stage[3:9]), function(X){
  peak1 <- peak.merge[[X]]
  peak2 <- peak.merge[[X - 1]]
  peak.lost <- setdiff(peak2,peak1)
  return(peak.lost)
})

# Number of peaks
elementNROWS(lost.peak)

##---------------------------------------------------------------------------##
##                  Plot number of gained and lost peaks at each stage                        
##---------------------------------------------------------------------------##

number.gain <- data.frame(stage = factor(stage[-(1:2)], levels = stage[-(1:2)]), number = elementNROWS(gain.peak))
number.lost <- data.frame(stage = factor(stage[-(1:2)], levels = stage[-(1:2)]), number = elementNROWS(lost.peak) * (-1))

number.gain.lost <- rbind(number.gain, number.lost)
number.gain.lost$group <- rep(c("Gained","Lost"),each = 7)

pdf("number.of.gained.lost.peaks.pdf")
ggplot(number.gain.lost, aes(stage, number,fill = group)) + 
  geom_col(color = "black") + 
  ggtitle("Number of gained/lost peaks") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle=45, hjust=1)) 
dev.off()
