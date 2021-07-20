#!/usr/bin/env Rscript

### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(tidyverse)
  library(patchwork)

### Load file
  dorthe <- fread("dortheB.kin")

  sc <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
  scA <- data.table(ID1=sc$clone, SCA=sc$SC, medrdA=sc$medrd)
  scB <- data.table(ID2=sc$clone, SCB=sc$SC, medrdB=sc$medrd)

  dortheA <- merge(dorthe, scA, by="ID1")
  dortheAB <- merge(dortheA, scB, by="ID2")

  dortheAB$type <- ifelse(dortheAB$SCA=="A" & dortheAB$SCB=="A", "AvsA",
    ifelse(dortheAB$SCA!="A" & dortheAB$SCB!="A", "OthervsOther", "OthervsA"))

  dortheAB$type <- factor(dortheAB$type, levels=c("AvsA", "OthervsOther", "OthervsA"))

### Graph

King <- ggplot(data=dortheAB[medrdA > 14 & medrdB > 14], aes(x=IBS0, y=Kinship, color=type)) + geom_point() +
  labs(color="Type") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

KingB <- ggplot(data=dortheAB[type!="OthervsA" & medrdA > 14 & medrdB > 14], aes(x=IBS0, y=Kinship, color=type)) + geom_point() +
  labs(color="Type") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


ggsave(King, file="DoertheKing.pdf")
ggsave(KingB, file="DoertheKingB.pdf")
