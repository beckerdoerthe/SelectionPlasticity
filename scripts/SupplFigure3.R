# Becker et al - Suppl FIG 3

library(data.table)
library(geomorph)
library(RRPP)
library(reshape2)
library(dplyr)
library(tidyverse)


load(file="data/all_data_final.RData")
all_data_final

### modify here re instar
all_data_use <- all_data_final[, c("cloneid_geno","Geno","SC_unique","group","treatment","instar","i","batch","height","APlength")][instar == 2]
setnames(all_data_use, c('APlength','height'), c('x', 'y'))  

all_data_use_wide <- reshape(all_data_use,   ## modify input data set 
                             idvar=c("Geno", "cloneid_geno", "SC_unique", "group", "instar", "treatment", "batch"),
                             timevar="i",
                             direction="wide")

setkey(all_data_use_wide, Geno)


load(file = "output/ta_out.OandA_I1_16June_fitted.RData")
# load(file = "output/ta_out.OandA_I2_16June_fitted.RData")

ta_out_red_use <- ta_out_red_mod[ta_out_red$Geno %in% all_data_use_wide$Geno][,c(1:1201)]  ## use only 600 data points

# procrustes data to use - "ta_out$fit$LM$fitted"
tmp_data_adj <- ta_out_red_use
names_adj <- as.data.table(ta_out_red_use$Geno)
procrustes_adj <- arrayspecs(tmp_data_adj[,2:ncol(tmp_data_adj)], 600, 2)  ## make 3D array object
dimnames(procrustes_adj)[3] <- names_adj  ## 3D array: number of rows (p), number of columns (k) and number of “sheets” (n), data has p=641, k=2, n= total number of rows.


## PROCRUSTES & TRAJECTORY ANALYSES  
# Step 1: procrustes analysis ----
gpa_trans_adj <- gpagen(procrustes_adj)  #, Proj = TRUE, ProcD = TRUE, curves = NULL, surfaces = NULL, print.progress = TRUE) 

gdf_adj <- geomorph.data.frame(
  Shape = gpa_trans_adj$coords,
  Clone = factor(all_data_use_wide[all_data_use_wide$Geno %in% ta_out_red_use$Geno][,c(1:1201)]$cloneid_geno),
  Predation = factor(all_data_use_wide[all_data_use_wide$Geno %in% ta_out_red_use$Geno][,c(1:1201)]$treatment),
  SC = factor(all_data_use_wide[all_data_use_wide$Geno %in% ta_out_red_use$Geno][,c(1:1201)]$SC_unique),
  Group = factor(all_data_use_wide[all_data_use_wide$Geno %in% ta_out_red_use$Geno][,c(1:1201)]$group),
  Batch = factor(all_data_use_wide[all_data_use_wide$Geno %in% ta_out_red_use$Geno][,c(1:1201)]$batch))


ref<-mshape(gpa_trans_adj$coords)
gp<-gdf_adj$Predation 
gc<-gdf_adj$Clone


# PCA_adj <- gm.prcomp(gpa_trans_adj$coords, groups = gp) 
#save(PCA_adj, file = "PCA_adj_I1.RData")
#save(PCA_adj, file = "PCA_adj_I2.RData")

# load(file = "output/PCA_adj_I1.RData")
load(file = "output/PCA_adj_I2.RData")


## tangent space plots
plotRefToTarget(ref,PCA_adj$shapes$shapes.PC1$min, method = 'TPS')  # links=links
plotRefToTarget(ref,PCA_adj$shapes$shapes.PC1$max, method = 'TPS')  # links=links
plotRefToTarget(ref,PCA_adj$shapes$shapes.PC2$min, method = 'TPS')  # links=links
plotRefToTarget(ref,PCA_adj$shapes$shapes.PC2$max, method = 'TPS')  # links=links
plotRefToTarget(ref,PCA_adj$shapes$shapes.PC3$min, method = 'TPS')  # links=links
plotRefToTarget(ref,PCA_adj$shapes$shapes.PC3$max, method = 'TPS')  # links=links


## scatter plots
PC_data <- as.data.table(PCA_adj$x[,1:3])
PC_data[, Geno:= rownames(PCA_adj$x)]

PC_scatter_plot <- ggplot(data = PC_data, 
                          aes(y=PC3, x=PC2, color=as.factor(gp))) + 
                      geom_point(size = 2) +  
                      theme(legend.position="none", 
                            rect = element_rect(fill = "transparent"),
                            panel.grid.major = element_line(colour = "grey70", size=0.25),
                            panel.grid.minor = element_line(colour = "grey90", size=0.1),
                            panel.background = element_rect(fill = "transparent",colour = NA),
                            plot.background = element_rect(fill = "transparent",colour = NA), 
                            #strip.text.x = element_blank(),
                            #axis.text.x = element_blank(), 
                            axis.title.x = element_blank(),
                            axis.title.y = element_blank(),
                            axis.line = element_line(size = 1),
                            # axis.title.x = element_text(size=30,family='Arial'),
                            # axis.title.y = element_text(size=30, family='Arial'),
                            axis.text = element_text(size=30, family='Arial'),
                            strip.text.x = element_text(size = 30, color = "black"),
                            strip.text.y = element_text(size = 30, color = "black"),
                            panel.spacing.x = unit(6, "mm"),
                            panel.spacing.y = unit(6, "mm")) 
                    
PC_scatter_plot + scale_color_manual(values=c("#000000","#FF0000")) + scale_y_continuous(limits=c(-0.025,0.030), breaks=c(-0.02,0,0.02)) + scale_x_continuous(limits=c(-0.03,0.045), breaks=c(-0.02,0,0.02,0.04))
PC_scatter_plot + scale_color_manual(values=c("#000000","#FF0000")) + scale_y_continuous(limits=c(-0.012,0.015), breaks=c(-0.01,0,0.01)) + scale_x_continuous(limits=c(-0.03,0.045), breaks=c(-0.02,0,0.02,0.04))
PC_scatter_plot + scale_color_manual(values=c("#000000","#FF0000")) + scale_y_continuous(limits=c(-0.012,0.015), breaks=c(-0.01,0,0.01)) + scale_x_continuous(limits=c(-0.025,0.030), breaks=c(-0.02,0,0.02))


