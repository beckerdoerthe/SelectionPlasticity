# Becker et al - Suppl FIG 4 


library(data.table)
library(geomorph)
library(RRPP)
library(reshape2)
library(dplyr)
library(tidyverse)


load(file="data/all_data_final.RData")
all_data_final

### modify here re instar and group - NOTE: only use the induced state for modularity check
load(file="data/all_data_final.RData")
# Os
#all_data_use <- all_data_final[, c("cloneid_geno","Geno","SC_unique","group","treatment","instar","i","batch","height","APlength")][instar == 2][group %in% c('O','ctrl_O')][treatment == 0.5]
# As
all_data_use <- all_data_final[, c("cloneid_geno","Geno","SC_unique","group","treatment","instar","i","batch","height","APlength")][instar == 2][group %in% c('A','ctrl_A')][treatment == 0.5]

setnames(all_data_use, c('APlength','height'), c('x', 'y'))  

all_data_use_wide <- reshape(all_data_use,   ## modify input data set 
                             idvar=c("Geno", "cloneid_geno", "SC_unique", "group", "instar", "treatment", "batch"),
                             timevar="i",
                             direction="wide")

setkey(all_data_use_wide, Geno)


load(file = "ta_out.OandA_I1_16June_fitted.RData")
# load(file = "ta_out.OandA_I2_16June_fitted.RData")

ta_out_red_use <- ta_out_red_mod[ta_out_red$Geno %in% all_data_use_wide$Geno][,c(1:1201)]  ## use only 600 data points

# procrustes data to use - "ta_out$fit$LM$fitted"
tmp_data_adj <- ta_out_red_use
names_adj <- as.data.table(ta_out_red_use$Geno)
procrustes_adj <- arrayspecs(tmp_data_adj[,2:ncol(tmp_data_adj)], 600, 2)  ## make 3D array object
dimnames(procrustes_adj)[3] <- names_adj  ## 3D array: number of rows (p), number of columns (k) and number of “sheets” (n), data has p=641, k=2, n= total number of rows.


## PROCRUSTES & TRAJECTORY ANALYSES  (modified from Andrew's code)
# Step 1: procrustes analysis ----
gpa_trans_adj <- gpagen(procrustes_adj)  #, Proj = TRUE, ProcD = TRUE, curves = NULL, surfaces = NULL, print.progress = TRUE) 

gdf_adj <- geomorph.data.frame(
  Shape = gpa_trans_adj$coords,
  Clone = factor(all_data_use_wide[all_data_use_wide$Geno %in% ta_out_red_use$Geno][,c(1:1201)]$cloneid_geno),
  Predation = factor(all_data_use_wide[all_data_use_wide$Geno %in% ta_out_red_use$Geno][,c(1:1201)]$treatment),
  SC = factor(all_data_use_wide[all_data_use_wide$Geno %in% ta_out_red_use$Geno][,c(1:1201)]$SC_unique),
  Group = factor(all_data_use_wide[all_data_use_wide$Geno %in% ta_out_red_use$Geno][,c(1:1201)]$group),
  Batch = factor(all_data_use_wide[all_data_use_wide$Geno %in% ta_out_red_use$Geno][,c(1:1201)]$batch))



## Model 1
land.gps.body100<-c(rep('a',300),rep('b',300))

MT_body100 <- modularity.test(gpa_trans_adj$coords,
                              land.gps.body100, 
                              CI=FALSE, iter=999)

# save(MT_body100, file = 'Mod1_MT_O_I1.RData')
# save(MT_body100, file = 'Mod1_MT_O_I2.RData')

# save(MT_body100, file = 'Mod1_MT_A_I1.RData')
save(MT_body100, file = 'Mod1_MT_A_I2.RData')



## Model 2
land.gps.body100<-c(rep('a',250),rep('b',350))

MT_body100 <- modularity.test(gpa_trans_adj$coords,
                              land.gps.body100, 
                              CI=FALSE, iter=999)

# save(MT_body100, file = 'Mod2_MT_O_I1.RData')
# save(MT_body100, file = 'Mod2_MT_O_I2.RData')

# save(MT_body100, file = 'Mod2_MT_A_I1.RData')
save(MT_body100, file = 'Mod2_MT_A_I2.RData')



## Model 3
land.gps.body100<-c(rep('a',200),rep('b',400))

MT_body100 <- modularity.test(gpa_trans_adj$coords,
                              land.gps.body100, 
                              CI=FALSE, iter=999)

# save(MT_body100, file = 'Mod3_MT_O_I1.RData')
# save(MT_body100, file = 'Mod3_MT_O_I2.RData')

# save(MT_body100, file = 'Mod3_MT_A_I1.RData')
save(MT_body100, file = 'Mod3_MT_A_I2.RData')



## Model 4
land.gps.body100<-c(rep('a',150),rep('b',450))

MT_body100 <- modularity.test(gpa_trans_adj$coords,
                              land.gps.body100, 
                              CI=FALSE, iter=999)

# save(MT_body100, file = 'Mod4_MT_O_I1.RData')
# save(MT_body100, file = 'Mod4_MT_O_I2.RData')

# save(MT_body100, file = 'Mod4_MT_A_I1.RData')
save(MT_body100, file = 'Mod4_MT_A_I2.RData')



## Model 5
land.gps.body100<-c(rep('a',100),rep('b',500))

MT_body100 <- modularity.test(gpa_trans_adj$coords,
                              land.gps.body100, 
                              CI=FALSE, iter=999)

# save(MT_body100, file = 'Mod5_MT_O_I1.RData')
# save(MT_body100, file = 'Mod5_MT_O_I2.RData')

# save(MT_body100, file = 'Mod5_MT_A_I1.RData')
save(MT_body100, file = 'Mod5_MT_A_I2.RData')



## Model 6
land.gps.body100<-c(rep('a',100),rep('b',200),rep('c',300))

MT_body100 <- modularity.test(gpa_trans_adj$coords,
                              land.gps.body100, 
                              CI=FALSE, iter=999)

# save(MT_body100, file = 'Mod6_MT_O_I1.RData')
# save(MT_body100, file = 'Mod6_MT_O_I2.RData')

# save(MT_body100, file = 'Mod6_MT_A_I1.RData')
save(MT_body100, file = 'Mod6_MT_A_I2.RData')



## Model 7
land.gps.body100<-c(rep('a',100),rep('b',100),rep('c',400))

MT_body100 <- modularity.test(gpa_trans_adj$coords,
                              land.gps.body100, 
                              CI=FALSE, iter=999)

# save(MT_body100, file = 'Mod7_MT_O_I1.RData')
# save(MT_body100, file = 'Mod7_MT_O_I2.RData')

# save(MT_body100, file = 'Mod7_MT_A_I1.RData')
save(MT_body100, file = 'Mod7_MT_A_I2.RData')



## Model 8
land.gps.body100<-c(rep('a',100),rep('b',150),rep('c',350))

MT_body100 <- modularity.test(gpa_trans_adj$coords,
                              land.gps.body100, 
                              CI=FALSE, iter=999)

# save(MT_body100, file = 'Mod8_MT_O_I1.RData')
# save(MT_body100, file = 'Mod8_MT_O_I2.RData')

# save(MT_body100, file = 'Mod8_MT_A_I1.RData')
save(MT_body100, file = 'Mod8_MT_A_I2.RData')



#########

# model 1
load(file = 'Mod1_MT_O_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod1_MT_O_I2.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod1_MT_A_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod1_MT_A_I2.RData')  
summary(MT_body100)
# plot(MT_body100)



# model 2
load(file = 'Mod2_MT_O_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod2_MT_O_I2.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod2_MT_A_I1.RData') 
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod2_MT_A_I2.RData')  
summary(MT_body100)
# plot(MT_body100)



# model 3
load(file = 'Mod3_MT_O_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod3_MT_O_I2.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod3_MT_A_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod3_MT_A_I2.RData')  
summary(MT_body100)
# plot(MT_body100)



# model 4
load(file = 'Mod4_MT_O_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod4_MT_O_I2.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod4_MT_A_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod4_MT_A_I2.RData')  
summary(MT_body100)
# plot(MT_body100)



# model 5
load(file = 'Mod5_MT_O_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod5_MT_O_I2.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod5_MT_A_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod5_MT_A_I2.RData')  
summary(MT_body100)
# plot(MT_body100)


# model 6
load(file = 'Mod6_MT_O_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod6_MT_O_I2.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod6_MT_A_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod6_MT_A_I2.RData')  
summary(MT_body100)
# plot(MT_body100)



# model 7
load(file = 'Mod7_MT_O_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod7_MT_O_I2.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod7_MT_A_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod7_MT_A_I2.RData')  
summary(MT_body100)
# plot(MT_body100)



# model 8
load(file = 'Mod8_MT_O_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod8_MT_O_I2.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod8_MT_A_I1.RData') 
summary(MT_body100)
# plot(MT_body100)

load(file = 'Mod8_MT_A_I2.RData')  
summary(MT_body100)
# plot(MT_body100)





load(file = "data/all_data_final.RData")
# all_data_final

# for trait plots, use one i-th position (e.g., 150) due to data dublication re i-th positions...

## AVERAGE SHAPE 
all_data_final.ag <- all_data_final[i <= 600][, list(height = mean(height)), list(i, treatment, instar_new, SC_group_new)]

ag.line_plot <- ggplot(data = all_data_final.ag[SC_group_new == "cluster O"][instar_new == "instar 1"],  aes(x=i, y=height)) +   
  
                    geom_line(data = all_data_final.ag[treatment == 0.5][SC_group_new == "cluster O"][instar_new == "instar 1"],  aes(x=i, y=height), size = 3, colour = "red") + 
                    
                    # geom_vline(xintercept = 300, linetype="dashed", color = "black", size=2) +
                    # geom_vline(xintercept = 250, linetype="dashed", color = "black", size=2) +
                    # geom_vline(xintercept = 200, linetype="dashed", color = "black", size=2) +
                    # geom_vline(xintercept = 150, linetype="dashed", color = "black", size=2) +
                    # geom_vline(xintercept = 100, linetype="dashed", color = "black", size=2) +
                    # geom_vline(xintercept = c(100, 300), linetype="dashed", color = "black", size=2) +
                    # geom_vline(xintercept = c(100, 250), linetype="dashed", color = "black", size=2) +
                    geom_vline(xintercept = c(100, 200), linetype="dashed", color = "black", size=2) +
  
                    theme(legend.position="none", 
                          rect = element_rect(fill = "transparent"),
                          # panel.grid.major = element_line(colour = "grey70", size=0.25),
                          # panel.grid.minor = element_line(colour = "grey90", size=0.1),
                          panel.background = element_rect(fill = "transparent",colour = NA),
                          plot.background = element_rect(fill = "transparent",colour = NA), 
                          #strip.text.x = element_blank(),
                          #axis.text.x = element_blank(), 
                          #axis.title.x = element_blank(), 
                          #axis.title.y = element_blank(),
                          axis.line = element_line(size = 1),
                          axis.title.x = element_text(size=60,family='Arial'),
                          axis.title.y = element_text(size=60, family='Arial'),
                          axis.text = element_text(size=60, family='Arial'),
                          strip.text.x = element_text(size = 60, color = "black"),
                          strip.text.y = element_text(size = 60, color = "black"),
                          panel.spacing.x = unit(6, "mm"),
                          panel.spacing.y = unit(6, "mm")) 


# tiff(file = "avg_shape_plot.tiff", width = 3200, height = 3200, units = "px", res = 800)
ag.line_plot  + labs(x = "dorsal position", y = "dorsal height (mm)") + scale_x_continuous(limits=c(0,621), breaks=c(0,100,200,300,400,500,600))
# dev.off()


