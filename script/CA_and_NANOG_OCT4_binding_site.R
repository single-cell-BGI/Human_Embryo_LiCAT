################################################################################
##        Overlaps between LiCAT-seq and NANOG/OCT4 binding sites
################################################################################

library(GenomicRanges)
library(AnnotationHub)
library(dplyr)
library(reshape2)
library(ggplot2)

ah <- AnnotationHub()

## working directory
#setwd()


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##       Overlapping peaks with ATAC-seq and NANOG peaks
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

# Load the GrangeList containing peak files from each stage
peak.grl.rmXY <- readRDS("atac/peak.grl.rmXY.hES.Liver.RDS")

# Download NANOG peakfiles from ESCs
epiFiles <- query(ah, c( "Nanog", "chip"))
nanog.h1 <- epiFiles[["AH22662"]]

OverlapWithNanog <- sapply(peak.grl.rmXY, function(gr) sum(countOverlaps(nanog.h1,gr) > 0))

number.nanog.peak <- data.frame(sample = names(OverlapWithNanog),
                                number = OverlapWithNanog,
                                celltype = gsub("_\\d","",names(OverlapWithNanog)))

number.nanog.peak <- number.nanog.peak %>%
  group_by(celltype) %>%
  summarise(sd = sd(number), mean = mean(number)) %>% 
  ungroup()

pdf("nanong.binding.site.pdf")
ggplot(number.nanog.peak, aes(factor(celltype,levels = celltype), mean, fill = celltype)) + 
  geom_col() + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9)) + 
  ggtitle("NANOG binding peaks") + 
  scale_fill_viridis(discrete=TRUE) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle=45, hjust=1)) + 
  xlab("") + 
  ylab("Number of peaks")
dev.off()

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##       Overlapping peaks with ATAC-seq and OCT4 peaks
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

# Download OCT4 peakfiles from ESCs
epiFiles <- query(ah, c( "POU5F1","chip"))
oct4.h1 <- epiFiles[["AH22667"]]

OverlapWithoct4 <- sapply(peak.grl.rmXY, function(gr) sum(countOverlaps(oct4.h1,gr) > 0))

number.oct4.peak <- data.frame(sample = names(OverlapWithoct4),
                               number = OverlapWithoct4,
                               celltype = gsub("_\\d","",names(OverlapWithoct4)))

number.oct4.peak <- number.oct4.peak %>%
  group_by(celltype) %>%
  summarise(sd = sd(number), mean = mean(number)) %>% 
  ungroup()

pdf("OCT4.binding.site.pdf")
ggplot(number.oct4.peak, aes(factor(celltype,levels = celltype), mean, fill = celltype)) + 
  geom_col() + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9)) + 
  ggtitle("OCT4 binding peaks") + 
  scale_fill_viridis(discrete=TRUE) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle=45, hjust=1)) + 
  xlab("") + 
  ylab("Number of peaks")
dev.off()