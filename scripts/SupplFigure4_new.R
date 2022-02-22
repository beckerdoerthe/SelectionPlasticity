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
all_data_use <- all_data_final[, c("cloneid_geno","Geno","SC_unique","group","treatment","instar","i","batch","height","APlength")][instar == 2][group %in% c('O','ctrl_O')][treatment == 0.5][i>=10 &i <=600]
# As
#all_data_use <- all_data_final[, c("cloneid_geno","Geno","SC_unique","group","treatment","instar","i","batch","height","APlength")][instar == 2][group %in% c('A','ctrl_A')][treatment == 0.5][i>=10 &i <=600]

setnames(all_data_use, c('APlength','height'), c('x', 'y'))  

all_data_use_wide <- reshape(all_data_use,   ## modify input data set 
                             idvar=c("Geno", "cloneid_geno", "SC_unique", "group", "instar", "treatment", "batch"),
                             timevar="i",
                             direction="wide")

setkey(all_data_use_wide, Geno)


# load(file = "output/ta_out.OandA_I1_16June_fitted.RData")
load(file = "output/ta_out.OandA_I2_16June_fitted.RData")
## data ist from ....
# load(file="ta_out.OandA_I2_16June.RData")
# 
# ta_out
# 
# ta_out_red <- as.data.table(ta_out$fit$LM$fitted)
# ta_out_red[, Geno := rownames(ta_out$fit$LM$fitted)]
# ta_out_red_mod <- ta_out_red %>% select(Geno, everything())
# 
# save(ta_out_red, ta_out_red_mod, file = "ta_out.OandA_I2_16June_fitted.RData")

ta_out_red_use <- ta_out_red_mod[ta_out_red$Geno %in% all_data_use_wide$Geno][,c(1, 20:1201)]  ## use only 591 data points & GenoIDs that are in the full data set that was used for other analyses

# procrustes data to use - "ta_out$fit$LM$fitted"
tmp_data_adj <- ta_out_red_use
names_adj <- as.data.table(ta_out_red_use$Geno)
procrustes_adj <- arrayspecs(tmp_data_adj[,2:ncol(tmp_data_adj)], 591, 2)  ## make 3D array object
dimnames(procrustes_adj)[3] <- names_adj  ## 3D array: number of rows (p), number of columns (k) and number of “sheets” (n), data has p=641, k=2, n= total number of rows.


## PROCRUSTES & TRAJECTORY ANALYSES - MODULARITY TEST
# Step 1: procrustes analysis ----
gpa_trans_adj <- gpagen(procrustes_adj)  #, Proj = TRUE, ProcD = TRUE, curves = NULL, surfaces = NULL, print.progress = TRUE) 

# gdf_adj <- geomorph.data.frame(
#   Shape = gpa_trans_adj$coords,
#   Clone = factor(all_data_use_wide[all_data_use_wide$Geno %in% ta_out_red_use$Geno]$cloneid_geno),
#   Predation = factor(all_data_use_wide[all_data_use_wide$Geno %in% ta_out_red_use$Geno]$treatment),
#   SC = factor(all_data_use_wide[all_data_use_wide$Geno %in% ta_out_red_use$Geno]$SC_unique),
#   Group = factor(all_data_use_wide[all_data_use_wide$Geno %in% ta_out_red_use$Geno]$group),
#   Batch = factor(all_data_use_wide[all_data_use_wide$Geno %in% ta_out_red_use$Geno]$batch))


## Model 1
land.gps.body_1 <- c(rep('a',291),rep('b',300))

MT_body_1 <- modularity.test(gpa_trans_adj$coords,
                              land.gps.body_1, 
                              CI=FALSE, iter=999)

save(MT_body_1, file = 'output/Mod1_MT_O_I2.RData')


## Model 2
land.gps.body_2 <- c(rep('a',241),rep('b',350))

MT_body_2 <- modularity.test(gpa_trans_adj$coords,
                             land.gps.body_2, 
                             CI=FALSE, iter=999)

save(MT_body_2, file = 'output/Mod2_MT_O_I2.RData')


## Model 3
land.gps.body_3 <- c(rep('a',191),rep('b',400))

MT_body_3 <- modularity.test(gpa_trans_adj$coords,
                             land.gps.body_3, 
                             CI=FALSE, iter=999)

save(MT_body_3, file = 'output/Mod3_MT_O_I2.RData')


## Model 4
land.gps.body_4 <- c(rep('a',141),rep('b',450))

MT_body_4 <- modularity.test(gpa_trans_adj$coords,
                             land.gps.body_4, 
                             CI=FALSE, iter=999)

save(MT_body_4, file = 'output/Mod4_MT_O_I2.RData')


## Model 5
land.gps.body_5 <- c(rep('a',91),rep('b',500))

MT_body_5 <- modularity.test(gpa_trans_adj$coords,
                             land.gps.body_5, 
                             CI=FALSE, iter=999)

save(MT_body_5, file = 'output/Mod5_MT_O_I2.RData')


## Model 6
land.gps.body_6 <- c(rep('a',91),rep('b',200),rep('c',300))

MT_body_6 <- modularity.test(gpa_trans_adj$coords,
                             land.gps.body_6, 
                             CI=FALSE, iter=999)

save(MT_body_6, file = 'output/Mod6_MT_O_I2.RData')


## Model 7
land.gps.body_7 <- c(rep('a',91),rep('b',150),rep('c',350))

MT_body_7 <- modularity.test(gpa_trans_adj$coords,
                             land.gps.body_7, 
                             CI=FALSE, iter=999)

save(MT_body_7, file = 'output/Mod7_MT_O_I2.RData')


## Model 8
land.gps.body_8 <- c(rep('a',91),rep('b',100),rep('c',400))

MT_body_8 <- modularity.test(gpa_trans_adj$coords,
                             land.gps.body_8, 
                             CI=FALSE, iter=999)

save(MT_body_8, file = 'output/Mod8_MT_O_I2.RData')


## Model 9
land.gps.body_9 <- c(rep('a',91),rep('b',50),rep('c',450))

MT_body_9 <- modularity.test(gpa_trans_adj$coords,
                             land.gps.body_9, 
                             CI=FALSE, iter=999)

save(MT_body_9, file = 'output/Mod9_MT_O_I2.RData')


## Model 10
land.gps.body_10 <- c(rep('a',91),rep('b',100),rep('c',200),rep('d',200))

MT_body_10 <- modularity.test(gpa_trans_adj$coords,
                             land.gps.body_10, 
                             CI=FALSE, iter=999)

save(MT_body_10, file = 'output/Mod10_MT_O_I2.RData')


## Model 11
land.gps.body_11 <- c(rep('a',91),rep('b',100),rep('c',150),rep('d',250))

MT_body_11 <- modularity.test(gpa_trans_adj$coords,
                             land.gps.body_11, 
                             CI=FALSE, iter=999)

save(MT_body_11, file = 'output/Mod11_MT_O_I2.RData')


## Model 12
land.gps.body_12 <- c(rep('a',91),rep('b',100),rep('c',100),rep('d',300))

MT_body_12 <- modularity.test(gpa_trans_adj$coords,
                             land.gps.body_12, 
                             CI=FALSE, iter=999)

save(MT_body_12, file = 'output/Mod12_MT_O_I2.RData')


## Model 13
land.gps.body_13 <- c(rep('a',91),rep('b',100),rep('c',50),rep('d', 350))

MT_body_13 <- modularity.test(gpa_trans_adj$coords,
                             land.gps.body_13, 
                             CI=FALSE, iter=999)

save(MT_body_13, file = 'output/Mod13_MT_O_I2.RData')


## Model 14
land.gps.body_14 <- c(rep('a',91),rep('b',100),rep('c',100),rep('d',100),rep('e',100),rep('f',100))

MT_body_14 <- modularity.test(gpa_trans_adj$coords,
                             land.gps.body_14, 
                             CI=FALSE, iter=999)

save(MT_body_14, file = 'output/Mod14_MT_O_I2.RData')


load(file = 'output/Mod1_MT_O_I2.RData')
load(file = 'output/Mod2_MT_O_I2.RData')
load(file = 'output/Mod3_MT_O_I2.RData')
load(file = 'output/Mod4_MT_O_I2.RData')
load(file = 'output/Mod5_MT_O_I2.RData')
load(file = 'output/Mod6_MT_O_I2.RData')
load(file = 'output/Mod7_MT_O_I2.RData')
load(file = 'output/Mod8_MT_O_I2.RData')
load(file = 'output/Mod9_MT_O_I2.RData')
load(file = 'output/Mod10_MT_O_I2.RData')
load(file = 'output/Mod11_MT_O_I2.RData')
load(file = 'output/Mod12_MT_O_I2.RData')
load(file = 'output/Mod13_MT_O_I2.RData')
load(file = 'output/Mod14_MT_O_I2.RData')


#compare models 1-8
model_1.14 <- compare.CR(MT_body_1,
                        MT_body_2,
                        MT_body_3,
                        MT_body_4,
                        MT_body_5,
                        MT_body_6,
                        MT_body_7, 
                        MT_body_8, 
                        MT_body_9,
                        MT_body_10,
                        MT_body_11,
                        MT_body_12,
                        MT_body_13,
                        MT_body_14)

model_1.14NULL <- compare.CR(MT_body_1,
                            MT_body_2,
                            MT_body_3,
                            MT_body_4,
                            MT_body_5,
                            MT_body_6,
                            MT_body_7, 
                            MT_body_8,
                            MT_body_9,
                            MT_body_10,
                            MT_body_11,
                            MT_body_12,
                            MT_body_13,
                            MT_body_14,
                            CR.null = TRUE)




#### old version ####

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


