# Becker et al 2020 (UVA) - Suppl FIG 1 (data provided by Karen)

#################
### libraries ###
#################
library(ggplot2)
library(cowplot)
library(ggbeeswarm)
library(gridExtra)
library(data.table)
library(viridis)
library(dplyr)
library(tidyverse)
library(tidyr)
library(foreach)
library(patchwork)


#######################
### load pheno data ###
#######################

# IBS

# load(file = "output/IBS.long.RData")
load(file = "output/IBS.long_new.RData")
# ibs.long

ibs.long_As <- ibs.long[SC.A == 'A'][SC.B == 'A'][!IBS == 1][, 'IBS']  #[! dist.noTri == 1][, 'IBS']
ibs.long_As[, set := "genetically similar"]
ibs.long_Os <- ibs.long[!SC.A == 'A'][!SC.B == 'A'][!IBS == 1][, 'IBS'] #[! dist.noTri == 1][, 'IBS']
ibs.long_Os[, set := "genetically unique"]

ibs.long_use <- rbind(ibs.long_As,ibs.long_Os)


## FREQUENCY IBS
hist_plot <- ggplot(data = ibs.long_use) +   
  
                geom_histogram(aes(x=IBS), bins = 50, fill='#C0C0C0', color='#000000', alpha=0.2, size=0.5) + 
                
                facet_wrap(~set, nrow=1) +
                
                theme(legend.position="none", 
                      rect = element_rect(fill = "transparent"),
                      panel.grid.major = element_line(colour = "grey70", size=0.25),
                      panel.grid.minor = element_line(colour = "grey90", size=0.1),
                      panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA), 
                      #strip.text.x = element_blank(),
                      #axis.text.x = element_blank(), 
                      #axis.title.x = element_blank(), 
                      #axis.title.y = element_blank(),
                      axis.line = element_line(size = 1),
                      axis.title.x = element_text(size=15,family='Arial'), 
                      axis.title.y = element_text(size=15, family='Arial'),
                      axis.text = element_text(size=15, family='Arial'),
                      strip.text.x = element_text(size = 15, color = "black"),
                      strip.text.y = element_text(size = 15, color = "black"),
                      panel.spacing.x = unit(6, "mm"),
                      panel.spacing.y = unit(6, "mm")) 

hist_plot + labs(x = "IBS", y = "frequency")


# kinship

### Load file
dorthe <- fread("dortheB.kin")

sc <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
scA <- data.table(ID1=sc$clone, SCA=sc$SC, medrdA=sc$medrd)
scB <- data.table(ID2=sc$clone, SCB=sc$SC, medrdB=sc$medrd)

dortheA <- merge(dorthe, scA, by="ID1")
dortheAB <- merge(dortheA, scB, by="ID2")

dortheAB$type <- ifelse(dortheAB$SCA=="A" & dortheAB$SCB=="A", "genetically similar",
                        ifelse(dortheAB$SCA!="A" & dortheAB$SCB!="A", "genetically unique", "OthervsA"))

dortheAB$type <- factor(dortheAB$type, levels=c("genetically similar", "genetically unique", "OthervsA"))
setnames(dortheAB, 'Kinship', 'kinship')

King_plot <- ggplot(data=dortheAB[type!="OthervsA" & medrdA > 14 & medrdB > 14], aes(x=IBS0, y=kinship, color=as.factor(type))) + 
                geom_point(size = 2, pch=19) +
              
                theme(legend.position="none", 
                      rect = element_rect(fill = "transparent"),
                      panel.grid.major = element_line(colour = "grey70", size=0.25),
                      panel.grid.minor = element_line(colour = "grey90", size=0.1),
                      panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA), 
                      #strip.text.x = element_blank(),
                      #axis.text.x = element_blank(), 
                      #axis.title.x = element_blank(), 
                      #axis.title.y = element_blank(),
                      axis.line = element_line(size = 1),
                      axis.title.x = element_text(size=15,family='Arial'), 
                      axis.title.y = element_text(size=15, family='Arial'),
                      axis.text = element_text(size=15, family='Arial'),
                      strip.text.x = element_text(size = 15, color = "black"),
                      strip.text.y = element_text(size = 15, color = "black"),
                      panel.spacing.x = unit(6, "mm"),
                      panel.spacing.y = unit(6, "mm")) 


King_plot + scale_color_manual(values=c("#FF0000","#0000CC"))


####

hist_plot + labs(x = "IBS", y = "frequency")

King_plot + scale_color_manual(values=c("#FF0000","#0000CC"))


patchwork_plots_IBS_kinship <- hist_plot + labs(x = "IBS", y = "frequency") +  
  (King_plot + scale_color_manual(values=c("#FF0000","#0000CC")))   # plot_spacer() +

patchwork_plots_IBS_kinship


