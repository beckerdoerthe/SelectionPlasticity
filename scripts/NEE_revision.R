## revision NEE


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
library(doMC)
registerDoMC(20)
library(effectsize)
library(rstatix)
library(patchwork)


#######################
### load pheno data ###
#######################

## load all data
load(file = "data/all_data_final.RData")
# all_data_final


## load anova/lm output
load(file = "output/aov_out.RData")
# height_ANOVA[term %in% c('geno','gxe','trt')]

height_ANOVA_wide <- spread(height_ANOVA[term %in% c('geno','gxe','trt')][, -c('lCI_eta','uCI_eta')], term, partial_eta)
height_ANOVA_wide_use <- height_ANOVA_wide[i <= 600]

setnames(height_ANOVA_wide_use, c('geno','gxe','trt'), c('effect_geno','effect_gxe','effect_trt'))

setkey(height_ANOVA_wide_use, i, instar, group)


# ## load H2 estimates
# load(file = "output/H2_O_I1_batch.RData")
# H2_O_I1 <- as.data.table(collected2)
# H2_O_I1[, i:= as.numeric(i)]
# H2_O_I1[, instar:= "instar 1"]
# H2_O_I1[, group:= "O"]
# setkey(H2_O_I1, group,instar,i)
# 
# load(file = "output/H2_O_I2_batch.RData")
# H2_O_I2 <- as.data.table(collected2)
# H2_O_I2[, i:= as.numeric(i)]
# H2_O_I2[, instar:= "instar 2"]
# H2_O_I2[, group:= "O"]
# setkey(H2_O_I2, group,instar,i)
# 
# load(file = "output/H2_A_I1_batch.RData")
# H2_A_I1 <- as.data.table(collected2)
# H2_A_I1[, i:= as.numeric(i)]
# H2_A_I1[, instar:= "instar 1"]
# H2_A_I1[, group:= "A"]
# setkey(H2_A_I1, group,instar,i)
# 
# load(file = "output/H2_A_I2_batch.RData")
# H2_A_I2 <- as.data.table(collected2)
# H2_A_I2[, i:= as.numeric(i)]
# H2_A_I2[, instar:= "instar 2"]
# H2_A_I2[, group:= "A"]
# setkey(H2_A_I2, group,instar,i)
# 
# AandO_H2_all <- rbind(H2_O_I1,
#                       H2_O_I2,
#                       H2_A_I1,
#                       H2_A_I2)
# 
# AandO_H2_all[, group_new := ifelse(AandO_H2_all$group == 'A', 'cluster A', 'cluster O')]
# 
# AandO_H2_all_use <- AandO_H2_all[, -'group']
# setnames(AandO_H2_all_use, 'group_new', 'group')
# 
# AandO_H2_all_wide <- spread(AandO_H2_all_use[label %in% c('C','T')][, -c('stuff.lCI','stuff.mode','stuff.uCI')], label, stuff.mean)
# setnames(AandO_H2_all_wide, c('C', 'T'), c('H2_ctrl', 'H2_trt'))
# 
# setkey(AandO_H2_all_wide, i, instar, group)


## load h2 estimates
load(file="output/O_h2_out_i_14June.RData")

# I1
O_I1_h2_0 <- h2_estimate_I1_0_out[Source %in% c("V(G)/Vp")]
setnames(O_I1_h2_0, c('Variance','SE'),c('h2','SE'))
O_I1_h2_0[, i := as.numeric(i)]
O_I1_h2_0[, instar:= "instar 1"]
O_I1_h2_0[, treatment:= "control"]
O_I1_h2_0[, U_h2 := h2+SE]
O_I1_h2_0[, L_h2 := h2-SE]
setkey(O_I1_h2_0,i)

O_I1_h2_05 <- h2_estimate_I1_05_out[Source %in% c("V(G)/Vp")]
setnames(O_I1_h2_05, c('Variance','SE'),c('h2','SE'))
O_I1_h2_05[, i := as.numeric(i)]
O_I1_h2_05[, instar:= "instar 1"]
O_I1_h2_05[, treatment:= "predation"]
O_I1_h2_05[, U_h2 := h2+SE]
O_I1_h2_05[, L_h2 := h2-SE]
setkey(O_I1_h2_05,i)

O_I2_h2_0 <- h2_estimate_I2_0_out[Source %in% c("V(G)/Vp")]
setnames(O_I2_h2_0, c('Variance','SE'),c('h2','SE'))
O_I2_h2_0[, i := as.numeric(i)]
O_I2_h2_0[, instar:= "instar 2"]
O_I2_h2_0[, treatment:= "control"]
O_I2_h2_0[, U_h2 := h2+SE]
O_I2_h2_0[, L_h2 := h2-SE]
setkey(O_I2_h2_0,i)

O_I2_h2_05 <- h2_estimate_I2_05_out[Source %in% c("V(G)/Vp")]
setnames(O_I2_h2_05, c('Variance','SE'),c('h2','SE'))
O_I2_h2_05[, i := as.numeric(i)]
O_I2_h2_05[, instar:= "instar 2"]
O_I2_h2_05[, treatment:= "predation"]
O_I2_h2_05[, U_h2 := h2+SE]
O_I2_h2_05[, L_h2 := h2-SE]
setkey(O_I2_h2_05,i)

O_h2_all <- rbind(O_I1_h2_0,
                  O_I1_h2_05,
                  O_I2_h2_0,
                  O_I2_h2_05)


O_h2_all_use <- O_h2_all[, c("i","instar", "treatment","h2")][i <= 600]
O_h2_all_use[, group := 'cluster O']

O_h2_all_wide <- spread(O_h2_all_use, treatment, h2)
setnames(O_h2_all_wide, c('control', 'predation'), c('h2_ctrl', 'h2_trt'))

setkey(O_h2_all_wide, i, instar, group)



## load H2 ratio estimates
## Vg/Vm ratio, i.e. log(Va) - log(Vm)

load(file = "output/H2_0_I1_ratio_batch.RData")
H2_ctrl_I1_ratio <- as.data.table(collected2)
H2_ctrl_I1_ratio[, i:= as.numeric(i)]
H2_ctrl_I1_ratio[, instar:= "instar 1"]
H2_ctrl_I1_ratio[, set:= "ctrl"]
setkey(H2_ctrl_I1_ratio, set,instar,i)

load(file = "output/H2_0_I2_ratio_batch.RData")
H2_ctrl_I2_ratio <- as.data.table(collected2)
H2_ctrl_I2_ratio[, i:= as.numeric(i)]
H2_ctrl_I2_ratio[, instar:= "instar 2"]
H2_ctrl_I2_ratio[, set:= "ctrl"]
setkey(H2_ctrl_I2_ratio, set,instar,i)

load(file = "output/H2_05_I1_ratio_batch.RData")
H2_trt_I1_ratio <- as.data.table(collected2)
H2_trt_I1_ratio[, i:= as.numeric(i)]
H2_trt_I1_ratio[, instar:= "instar 1"]
H2_trt_I1_ratio[, set:= "trt"]
setkey(H2_trt_I1_ratio, set,instar,i)

load(file = "output/H2_05_I2_ratio_batch.RData")
H2_trt_I2_ratio <- as.data.table(collected2)
H2_trt_I2_ratio[, i:= as.numeric(i)]
H2_trt_I2_ratio[, instar:= "instar 2"]
H2_trt_I2_ratio[, set:= "trt"]
setkey(H2_trt_I2_ratio, set,instar,i)


H2_ctrl_I1_ratio_use <- H2_ctrl_I1_ratio[label == "RATIO"]
H2_ctrl_I1_ratio_use[, treatment := 0]
H2_ctrl_I1_ratio_use[, instar := 1]
H2_ctrl_I1_ratio_use[, instar_new := ifelse(H2_ctrl_I1_ratio_use$instar == 1, 'instar 1', 'instar 2')]

H2_ctrl_I2_ratio_use <- H2_ctrl_I2_ratio[label == "RATIO"]
H2_ctrl_I2_ratio_use[, treatment := 0]
H2_ctrl_I2_ratio_use[, instar := 2]
H2_ctrl_I2_ratio_use[, instar_new := ifelse(H2_ctrl_I2_ratio_use$instar == 1, 'instar 1', 'instar 2')]

H2_trt_I1_ratio_use <- H2_trt_I1_ratio[label == "RATIO"]
H2_trt_I1_ratio_use[, treatment := 0.5]
H2_trt_I1_ratio_use[, instar := 1]
H2_trt_I1_ratio_use[, instar_new := ifelse(H2_trt_I1_ratio_use$instar == 1, 'instar 1', 'instar 2')]

H2_trt_I2_ratio_use <- H2_trt_I2_ratio[label == "RATIO"]
H2_trt_I2_ratio_use[, treatment := 0.5]
H2_trt_I2_ratio_use[, instar := 2]
H2_trt_I2_ratio_use[, instar_new := ifelse(H2_trt_I2_ratio_use$instar == 1, 'instar 1', 'instar 2')]


H2_ratio_all <- rbind(H2_ctrl_I1_ratio_use,
                      H2_ctrl_I2_ratio_use,
                      H2_trt_I1_ratio_use,
                      H2_trt_I2_ratio_use)

H2_ratio_all[, group := "VAoverVM"]

H2_ratio_all_wide <- spread(H2_ratio_all[, -c('stuff.lCI','stuff.mode','stuff.uCI')][, -c('instar', 'treatment','label')], set, stuff.mean)

setnames(H2_ratio_all_wide, c('instar_new','ctrl','trt'), c('instar','H2_VAoverVM_ctrl','H2_VAoverVM_trt'))

setkey(H2_ratio_all_wide, i, instar, group)


# height_ANOVA_wide_use
# AandO_H2_all_wide
# O_h2_all_wide
# H2_ratio_all_wide


## filter for O data only, yet including the VAoverVM ratio
# add module info and make data long to merge data sets

# effect sizes
effect_sizes <- height_ANOVA_wide_use #[group == 'cluster O']
setkey(effect_sizes, i, instar, group)

effect_sizes[, module := ifelse(effect_sizes$instar == 'instar 1' & effect_sizes$i %in% c(1:300), "module 1", 
                                ifelse(effect_sizes$instar == 'instar 1' & effect_sizes$i %in% c(301:600), 'module 2',
                                              
                                       ifelse(effect_sizes$instar == 'instar 2' & effect_sizes$i %in% c(1:100), 'module 1',
                                              ifelse(effect_sizes$instar == 'instar 2' & effect_sizes$i %in% c(101:200), 'module 2',
                                                     ifelse(effect_sizes$instar == 'instar 2' & effect_sizes$i %in% c(201:600), 'module 3','NA')))))]

effect_sizes_long <- melt(effect_sizes, id.vars=c("i","instar","group","module")) 


# H2 estimates
H2_est <- AandO_H2_all_wide #[group == 'cluster O']
setkey(H2_est, i, instar, group)

H2_est[, module := ifelse(H2_est$instar == 'instar 1' & H2_est$i %in% c(1:300), "module 1", 
                                ifelse(H2_est$instar == 'instar 1' & H2_est$i %in% c(301:600), 'module 2',
                                       
                                       ifelse(H2_est$instar == 'instar 2' & H2_est$i %in% c(1:100), 'module 1',
                                              ifelse(H2_est$instar == 'instar 2' & H2_est$i %in% c(101:200), 'module 2',
                                                     ifelse(H2_est$instar == 'instar 2' & H2_est$i %in% c(201:600), 'module 3','NA')))))]

H2_est_long <- melt(H2_est, id.vars=c("i","instar","group","module")) 


# h2 estimates
h2_est <- O_h2_all_wide
setkey(h2_est, i, instar, group)

h2_est[, module := ifelse(h2_est$instar == 'instar 1' & h2_est$i %in% c(1:300), "module 1", 
                          ifelse(h2_est$instar == 'instar 1' & h2_est$i %in% c(301:600), 'module 2',
                                 
                                 ifelse(h2_est$instar == 'instar 2' & h2_est$i %in% c(1:100), 'module 1',
                                        ifelse(h2_est$instar == 'instar 2' & h2_est$i %in% c(101:200), 'module 2',
                                               ifelse(h2_est$instar == 'instar 2' & h2_est$i %in% c(201:600), 'module 3','NA')))))]

h2_est_long <- melt(h2_est, id.vars=c("i","instar","group","module")) 


# VAoverVM estimates
H2_VAoverVm_est <- H2_ratio_all_wide
setkey(H2_VAoverVm_est, i, instar, group)

H2_VAoverVm_est[, module := ifelse(H2_VAoverVm_est$instar == 'instar 1' & H2_VAoverVm_est$i %in% c(1:300), "module 1", 
                                   ifelse(H2_VAoverVm_est$instar == 'instar 1' & H2_VAoverVm_est$i %in% c(301:600), 'module 2',
                                          
                                          ifelse(H2_VAoverVm_est$instar == 'instar 2' & H2_VAoverVm_est$i %in% c(1:100), 'module 1',
                                                 ifelse(H2_VAoverVm_est$instar == 'instar 2' & H2_VAoverVm_est$i %in% c(101:200), 'module 2',
                                                        ifelse(H2_VAoverVm_est$instar == 'instar 2' & H2_VAoverVm_est$i %in% c(201:600), 'module 3','NA')))))]

H2_VAoverVm_est_long <- melt(H2_VAoverVm_est, id.vars=c("i","instar","group","module")) 


# effect_sizes_long
# h2_est_long
# H2_est_long
# H2_VAoverVm_est_long


effectSize_heritability <- rbind(effect_sizes_long,
                                 h2_est_long,
                                 H2_est_long,
                                 H2_VAoverVm_est_long)

effectSize_heritability[, module := as.factor(module)]

# save(effectSize_heritability, file="effectSize_heritability.RData")
load(file="output/effectSize_heritability.RData")


## plot data
# PREDATION VS CONTROL
## only instar 2 & h2 estimates
box1 <- ggplot(effectSize_heritability[variable %in% c('effect_trt','h2_trt','h2_ctrl')][instar == "instar 2"], aes(y=value, fill=(variable))) + 
                    geom_boxplot(notch=TRUE, outlier.colour="#555555", outlier.shape=8,
                                 outlier.size=2) + 
                    facet_grid(~module) + 
                    theme(legend.position="right", 
                          rect = element_rect(fill = "transparent"),
                          panel.grid.major = element_line(colour = "grey70", size=0.25),
                          panel.grid.minor = element_line(colour = "grey90", size=0.1),
                          panel.background = element_rect(fill = "transparent",colour = NA),
                          plot.background = element_rect(fill = "transparent",colour = NA), 
                          #strip.text.x = element_blank(),
                          axis.text.x = element_blank(), 
                          axis.title.x = element_blank(), 
                          #axis.title.y = element_blank(),
                          axis.line = element_line(size = 1),
                          # axis.title.x = element_text(size=15,family='Arial'), 
                          axis.title.y = element_text(size=15, family='Arial'),
                          axis.text.y = element_text(size=15, family='Arial'),
                          strip.text.x = element_text(size = 15, color = "black"),
                          strip.text.y = element_text(size = 15, color = "black"),
                          panel.spacing.x = unit(6, "mm"),
                          panel.spacing.y = unit(6, "mm")) 

box1 + labs(x = "module", y = expression(Estiamtes~"("~V[a]~OR~plasticity~")")) + scale_fill_manual(values=c("#DCDCDC","#000000","#FF0000"))


# PREDATION, including VAoverVM ratio
box2 <- ggplot(effectSize_heritability[variable %in% c('H2_VAoverVM_trt','H2_VAoverVM_ctrl')], aes(y=value, fill=(variable))) + 
                geom_boxplot(notch=TRUE, outlier.colour="#555555", outlier.shape=8,
                             outlier.size=2) + 
                facet_grid(~module) + 
                theme(legend.position="right", 
                      rect = element_rect(fill = "transparent"),
                      panel.grid.major = element_line(colour = "grey70", size=0.25),
                      panel.grid.minor = element_line(colour = "grey90", size=0.1),
                      panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA), 
                      #strip.text.x = element_blank(),
                      axis.text.x = element_blank(), 
                      axis.title.x = element_blank(), 
                      #axis.title.y = element_blank(),
                      axis.line = element_line(size = 1),
                      # axis.title.x = element_text(size=15,family='Arial'), 
                      axis.title.y = element_text(size=15, family='Arial'),
                      axis.text.y = element_text(size=15, family='Arial'),
                      strip.text.x = element_text(size = 15, color = "black"),
                      strip.text.y = element_text(size = 15, color = "black"),
                      panel.spacing.x = unit(6, "mm"),
                      panel.spacing.y = unit(6, "mm")) 

box2 + labs(x = "module", y = expression(log[10]~"("~V[g]~"/"~V[m]~")")) + scale_fill_manual(values=c("#000000","#FF0000"))

(box1 + labs(x = "module", y = expression(Estiamtes~"("~V[a]~OR~plasticity~")")) + scale_fill_manual(values=c("#DCDCDC","#000000","#FF0000")) + scale_y_continuous(limits=c(0,0.6))) / 
  (box2 + labs(x = "module", y = expression(log[10]~"("~V[g]~"/"~V[m]~")")) + scale_fill_manual(values=c("#000000","#FF0000")) + scale_y_continuous(limits=c(-8.5,10)))



# stats - CTRL vs TRT
effectSize_heritability_stats <- effectSize_heritability[variable %in% c('h2_trt','h2_ctrl')][instar == "instar 2"]
effectSize_heritability_stats_wide <- spread(effectSize_heritability_stats, variable, value)

# module 1
shapiro.test(effectSize_heritability_stats_wide[module == "module 1"]$h2_ctrl)  # not normal 
shapiro.test(effectSize_heritability_stats_wide[module == "module 1"]$h2_trt)  # not normal 

wilcox_mod1 <- wilcox.test(value ~ variable, data = effectSize_heritability_stats[module == "module 1"], paired = FALSE) 
wilcox_mod1

Zscore_mod1 = qnorm(wilcox_mod1$p.value/2)
Zscore_mod1

# (Z = -7.561211, p-value = 3.993e-14)


t.test(effectSize_heritability_stats_wide$h2_trt, effectSize_heritability_stats_wide$h2_ctrl)

# module 2
shapiro.test(effectSize_heritability_stats_wide[module == "module 2"]$h2_ctrl)  # normal 
shapiro.test(effectSize_heritability_stats_wide[module == "module 2"]$h2_trt)  # not normal 

wilcox_mod2 <- wilcox.test(value ~ variable, data = effectSize_heritability_stats[module == "module 2"]) 
wilcox_mod2

Zscore_mod2 = qnorm(wilcox_mod2$p.value/2)
Zscore_mod2

# (Z = -9.075971, p-value < 2.2e-16)


# module 3
shapiro.test(effectSize_heritability_stats_wide[module == "module 3"]$h2_ctrl)  # normal 
shapiro.test(effectSize_heritability_stats_wide[module == "module 3"]$h2_trt)  # not normal 

wilcox_mod3 <- wilcox.test(value ~ variable, data = effectSize_heritability_stats[module == "module 3"]) 
wilcox_mod3

Zscore_mod3 = qnorm(wilcox_mod3$p.value/2)
Zscore_mod3

# (Z = -9.075971, p-value < 2.2e-16)




## correlation along i 

#libraries
library(data.table)
library(tidyverse)
library(tidyr)
library(patchwork)
library(gridExtra)


effectSize_heritability_wide <- spread(effectSize_heritability, variable, value)

# raw data plots 

# effect size 'treatment' & h2 estimates - instar 2 only
p1_trt <- ggplot(effectSize_heritability_wide[group == "cluster O"][instar == "instar 2"], aes(x = effect_trt, y = h2_trt, colour = i))+
            geom_point(size = 1)+
            facet_wrap(~module) + 
            ylim(0, 0.5) +
            theme(legend.position="right", 
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

p1_trt

p1_ctrl <- ggplot(effectSize_heritability_wide[group == "cluster O"][instar == "instar 2"], aes(x = effect_trt, y = h2_ctrl, colour = i))+
                  geom_point(size = 1)+
                  facet_wrap(~module) + 
                  ylim(0, 0.5) +
                  theme(legend.position="right", 
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

p1_ctrl

p1_ctrl / p1_trt


# effect size 'treatment' & H2 estimates - instar 2 only
p2_trt <- ggplot(effectSize_heritability_wide[group == "cluster O"][instar == "instar 2"], aes(x = effect_trt, y = H2_trt, colour = i))+
                geom_point(size = 1)+
                facet_wrap(~module) + 
                ylim(0, 0.5) +
                theme(legend.position="right", 
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

p2_trt

p2_ctrl <- ggplot(effectSize_heritability_wide[group == "cluster O"][instar == "instar 2"], aes(x = effect_trt, y = H2_ctrl, colour = i))+
                  geom_point(size = 1)+
                  facet_wrap(~module) + 
                  ylim(0, 0.5) +
                  theme(legend.position="right", 
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

p2_ctrl

p2_ctrl / p2_trt



p3_trt <- ggplot(effectSize_heritability_wide[group == "cluster A"][instar == "instar 2"], aes(x = effect_trt, y = H2_trt, colour = i))+
                geom_point(size = 1)+
                facet_wrap(~module) + 
                ylim(0, 0.5) +
                theme(legend.position="right", 
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
              
p3_trt

p3_ctrl <- ggplot(effectSize_heritability_wide[group == "cluster A"][instar == "instar 2"], aes(x = effect_trt, y = H2_ctrl, colour = i))+
                  geom_point(size = 1)+
                  facet_wrap(~module) + 
                  ylim(0, 0.5) +
                  theme(legend.position="right", 
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

p3_ctrl

p3_ctrl / p3_trt


(p3_ctrl / p3_trt) | (p2_ctrl / p2_trt)


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
library(doMC)
registerDoMC(20)
library(sjstats)
library(rstatix)
library(robcor)
library(doMC)
registerDoMC(20)
library(lme4)


#######################
### load pheno data ###
#######################

load(file = "data/all_data_final.RData")
# all_data_final

# ignore deme 3 data, i.e., extra replicates (here: only 6 data points among A clones)
all_data_final_mod <- all_data_final[!deme == 3]

# replace C and D reps with A and B (justified due to replicate info only used for within clutch variation...)
all_data_final_mod[, replicate_new := ifelse(all_data_final_mod$replicate == "1C", "1A", 
                                             ifelse(all_data_final_mod$replicate == "1D", "1B", 
                                                    ifelse(all_data_final_mod$replicate == "2C", "2A", 
                                                           ifelse(all_data_final_mod$replicate == "2D", "2B",  all_data_final_mod$replicate))))]

all_data_final_mod[, clutch := ifelse(all_data_final_mod$replicate_new %like% "A", "A", "B")]


## LMER() - extract variance components

# run LMM to estimate the contribution of each level to overall variance
# (within-clutch, between-mother, between-clone, and between-batch) 
# i.e., 'clutch', 'deme', 'cloneID', 'batch_new'

# reduce to one i-th position for animal length, max height, and eye area
all_data_use <- all_data_final_mod[i == 150]

# revision on LMM based on Alan's suggestion: 1|clone, 1|clone/deme, 1|clone/deme/clutch

# animal length

model_length_1 <- lmer(length ~ treatment + (1|cloneid_geno) + (1|batch_new), data = all_data_use)
model_length_1
VarianceComponents_1 <- VarCorr(model_length_1)
VarianceComponents_1

model_length_2 <- lmer(length ~ treatment + (1|cloneid_geno/deme) + (1|batch_new), data = all_data_use)
model_length_2
VarianceComponents_2 <- VarCorr(model_length_2)
VarianceComponents_2

model_length_3 <- lmer(length ~ treatment + (1|cloneid_geno/deme/clutch) + (1|batch_new), data = all_data_use)
model_length_3
VarianceComponents_3 <- VarCorr(model_length_3)
VarianceComponents_3

anova(model_length_1,model_length_2)
anova(model_length_1,model_length_3)


# max dorsal height

model_maxHeight_1 <- lmer(max_height_new ~ treatment + (1|cloneid_geno) + (1|batch_new), data = all_data_use)
model_maxHeight_1
VarianceComponents_1 <- VarCorr(model_maxHeight_1)
VarianceComponents_1

model_maxHeight_2 <- lmer(max_height_new ~ treatment + (1|cloneid_geno/deme) + (1|batch_new), data = all_data_use)
model_maxHeight_2
VarianceComponents_2 <- VarCorr(model_maxHeight_2)
VarianceComponents_2

model_maxHeight_3 <- lmer(max_height_new ~ treatment + (1|cloneid_geno/deme/clutch) + (1|batch_new), data = all_data_use)
model_maxHeight_3
VarianceComponents_3 <- VarCorr(model_maxHeight_3)
VarianceComponents_3

anova(model_maxHeight_1,model_maxHeight_2)
anova(model_maxHeight_1,model_maxHeight_3)





## correlations in modules

trt_h2_instar1 <- effectSize_heritability_wide[group == "cluster O"][, c('i','instar','group','module','effect_trt', 'h2_trt')] %>% 
  dplyr::filter(instar == "instar 1") %>% 
  dplyr::mutate(section = case_when(
    i<=200 ~ "10-200",
    i>100&i<=299 ~ "100-300",
    i>200&i<=399 ~ "200-400",
    i>300&i<=499 ~ "300-500",
    i>400&i<=600 ~ "400-600"
  )) %>% dplyr::tibble()


trt_H2_instar1 <- effectSize_heritability_wide[group == "cluster O"][, c('i','instar','group','module','effect_trt', 'H2_trt')] %>% 
  dplyr::filter(instar == "instar 1") %>% 
  dplyr::mutate(section = case_when(
    i<=200 ~ "10-200",
    i>100&i<=299 ~ "100-300",
    i>200&i<=399 ~ "200-400",
    i>300&i<=499 ~ "300-500",
    i>400&i<=600 ~ "400-600"
  )) %>% dplyr::tibble()


trt_h2_instar2 <- effectSize_heritability_wide[group == "cluster O"][, c('i','instar','group','module','effect_trt', 'h2_trt')] %>% 
  dplyr::filter(instar == "instar 2") %>% 
  dplyr::mutate(section = case_when(
    i<=200 ~ "10-200",
    i>100&i<=299 ~ "100-300",
    i>200&i<=399 ~ "200-400",
    i>300&i<=499 ~ "300-500",
    i>400&i<=600 ~ "400-600"
  )) %>% dplyr::tibble()

trt_H2_instar2 <- effectSize_heritability_wide[group == "cluster O"][, c('i','instar','group','module','effect_trt', 'H2_trt')] %>% 
  dplyr::filter(instar == "instar 2") %>% 
  dplyr::mutate(section = case_when(
    i<=200 ~ "10-200",
    i>100&i<=299 ~ "100-300",
    i>200&i<=399 ~ "200-400",
    i>300&i<=499 ~ "300-500",
    i>400&i<=600 ~ "400-600"
  )) %>% dplyr::tibble()


# trt_h2_instar1
cor_trt_h2_instar1 <- as.data.table(trt_h2_instar1 %>%   
                                      dplyr::group_by(section) %>% 
                                      dplyr::summarise(cor(h2_trt, effect_trt)))

cor_trt_h2_instar1[, position := c(100,200,300,400,500)]

setkey(cor_trt_h2_instar1, position)
setnames(cor_trt_h2_instar1, "cor(h2_trt, effect_trt)", "cor_h2_trt")

plot_trt_h2_instar1 <- ggplot(data=cor_trt_h2_instar1, aes(x=position, y=cor_h2_trt)) +
                              geom_bar(stat="identity") + 
                              # ylab("cor")+xlab("dorsal position") + 
                              ylim(-1, 1) +
                              theme(legend.position="right", 
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

plot_trt_h2_instar1


# trt_H2_instar1

cor_trt_H2_instar1 <- as.data.table(trt_H2_instar1 %>%   
                                      dplyr::group_by(section) %>% 
                                      dplyr::summarise(cor(H2_trt, effect_trt)))

cor_trt_H2_instar1[, position := c(100,200,300,400,500)]

setkey(cor_trt_H2_instar1, position)
setnames(cor_trt_H2_instar1, "cor(H2_trt, effect_trt)", "cor_H2_trt")

plot_trt_H2_instar1 <- ggplot(data=cor_trt_H2_instar1, aes(x=position, y=cor_H2_trt)) +
                              geom_bar(stat="identity") + 
                              # ylab("cor")+xlab("dorsal position") + 
                              ylim(-1, 1) +
                              theme(legend.position="right", 
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

plot_trt_H2_instar1


# trt_h2_instar2

cor_trt_h2_instar2 <- as.data.table(trt_h2_instar2 %>%   ## modify here (!)
                                      dplyr::group_by(section) %>% 
                                      dplyr::summarise(cor(h2_trt, effect_trt)))

cor_trt_h2_instar2[, position := c(100,200,300,400,500)]

setkey(cor_trt_h2_instar2, position)
setnames(cor_trt_h2_instar2, "cor(h2_trt, effect_trt)", "cor_h2_trt")


plot_trt_h2_instar2 <- ggplot(data=cor_trt_h2_instar2, aes(x=position, y=cor_h2_trt)) +
                              geom_bar(stat="identity") + 
                              # ylab("cor")+xlab("dorsal position") + 
                              ylim(-1, 1) +
                              theme(legend.position="right", 
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

plot_trt_h2_instar2


# trt_H2_instar2

cor_trt_H2_instar2 <- as.data.table(trt_H2_instar2 %>%   ## modify here (!)
  dplyr::group_by(section) %>% 
  dplyr::summarise(cor(H2_trt, effect_trt)))

cor_trt_H2_instar2[, position := c(100,200,300,400,500)]

setkey(cor_trt_H2_instar2, position)
setnames(cor_trt_H2_instar2, "cor(H2_trt, effect_trt)", "cor_H2_trt")


plot_trt_H2_instar2 <- ggplot(data=cor_trt_H2_instar2, aes(x=position, y=cor_H2_trt)) +
                              geom_bar(stat="identity") + 
                              # ylab("cor")+xlab("dorsal position") + 
                              ylim(-1, 1) +
                              theme(legend.position="right", 
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

plot_trt_H2_instar2


# (plot_trt_h2_instar1|plot_trt_h2_instar2) / (plot_trt_H2_instar1|plot_trt_H2_instar2)
(p1_trt/p2_trt/p3_trt)/(plot_trt_h2_instar1|plot_trt_h2_instar2) / (plot_trt_H2_instar1|plot_trt_H2_instar2)




###### increase # bins ######

trt_h2_instar1 <- effectSize_heritability_wide[group == "cluster O"][, c('i','instar','group','module','effect_trt', 'h2_trt')] %>% 
  dplyr::filter(instar == "instar 1") %>% 
  dplyr::mutate(section = case_when(
    i<=100 ~ "10-100",
    i>50&i<=150 ~ "50-150",
    i>100&i<=200 ~ "100-200",
    i>150&i<=250 ~ "150-250",
    i>200&i<=300 ~ "200-300",
    i>250&i<=350 ~ "250-350",
    i>300&i<=400 ~ "300-400",
    i>350&i<=450 ~ "350-450",
    i>400&i<=500 ~ "400-500",
    i>450&i<=550 ~ "450-550",
    i>500&i<=600 ~ "500-600"
  )) %>% dplyr::tibble()


trt_H2_instar1 <- effectSize_heritability_wide[group == "cluster O"][, c('i','instar','group','module','effect_trt', 'H2_trt')] %>% 
  dplyr::filter(instar == "instar 1") %>% 
  dplyr::mutate(section = case_when(
    i<=100 ~ "10-100",
    i>50&i<=150 ~ "50-150",
    i>100&i<=200 ~ "100-200",
    i>150&i<=250 ~ "150-250",
    i>200&i<=300 ~ "200-300",
    i>250&i<=350 ~ "250-350",
    i>300&i<=400 ~ "300-400",
    i>350&i<=450 ~ "350-450",
    i>400&i<=500 ~ "400-500",
    i>450&i<=550 ~ "450-550",
    i>500&i<=600 ~ "500-600"
  )) %>% dplyr::tibble()


trt_h2_instar2 <- effectSize_heritability_wide[group == "cluster O"][, c('i','instar','group','module','effect_trt', 'h2_trt')] %>% 
  dplyr::filter(instar == "instar 2") %>% 
  dplyr::mutate(section = case_when(
    i<=100 ~ "10-100",
    i>50&i<=150 ~ "50-150",
    i>100&i<=200 ~ "100-200",
    i>150&i<=250 ~ "150-250",
    i>200&i<=300 ~ "200-300",
    i>250&i<=350 ~ "250-350",
    i>300&i<=400 ~ "300-400",
    i>350&i<=450 ~ "350-450",
    i>400&i<=500 ~ "400-500",
    i>450&i<=550 ~ "450-550",
    i>500&i<=600 ~ "500-600"
  )) %>% dplyr::tibble()

trt_H2_instar2 <- effectSize_heritability_wide[group == "cluster O"][, c('i','instar','group','module','effect_trt', 'H2_trt')] %>% 
  dplyr::filter(instar == "instar 2") %>% 
  dplyr::mutate(section = case_when(
    i<=100 ~ "10-100",
    i>50&i<=150 ~ "50-150",
    i>100&i<=200 ~ "100-200",
    i>150&i<=250 ~ "150-250",
    i>200&i<=300 ~ "200-300",
    i>250&i<=350 ~ "250-350",
    i>300&i<=400 ~ "300-400",
    i>350&i<=450 ~ "350-450",
    i>400&i<=500 ~ "400-500",
    i>450&i<=550 ~ "450-550",
    i>500&i<=600 ~ "500-600"
  )) %>% dplyr::tibble()


# trt_h2_instar1
cor_trt_h2_instar1 <- as.data.table(trt_h2_instar1 %>%   
                                      dplyr::group_by(section) %>% 
                                      dplyr::summarise(cor(h2_trt, effect_trt)))

cor_trt_h2_instar1[, position := c(50,150,200,250,300,350,400,450,500,100,550)]

setkey(cor_trt_h2_instar1, position)
setnames(cor_trt_h2_instar1, "cor(h2_trt, effect_trt)", "cor_h2_trt")

plot_trt_h2_instar1 <- ggplot(data=cor_trt_h2_instar1, aes(x=position, y=cor_h2_trt)) +
  geom_bar(stat="identity") + 
  # ylab("cor")+xlab("dorsal position") + 
  ylim(-1, 1) +
  theme(legend.position="right", 
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

plot_trt_h2_instar1


# trt_H2_instar1

cor_trt_H2_instar1 <- as.data.table(trt_H2_instar1 %>%   
                                      dplyr::group_by(section) %>% 
                                      dplyr::summarise(cor(H2_trt, effect_trt)))

cor_trt_H2_instar1[, position := c(50,150,200,250,300,350,400,450,500,100,550)]

setkey(cor_trt_H2_instar1, position)
setnames(cor_trt_H2_instar1, "cor(H2_trt, effect_trt)", "cor_H2_trt")

plot_trt_H2_instar1 <- ggplot(data=cor_trt_H2_instar1, aes(x=position, y=cor_H2_trt)) +
  geom_bar(stat="identity") + 
  # ylab("cor")+xlab("dorsal position") + 
  ylim(-1, 1) +
  theme(legend.position="right", 
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

plot_trt_H2_instar1


# trt_h2_instar2

cor_trt_h2_instar2 <- as.data.table(trt_h2_instar2 %>%   ## modify here (!)
                                      dplyr::group_by(section) %>% 
                                      dplyr::summarise(cor(h2_trt, effect_trt)))

cor_trt_h2_instar2[, position := c(50,150,200,250,300,350,400,450,500,100,550)]

setkey(cor_trt_h2_instar2, position)
setnames(cor_trt_h2_instar2, "cor(h2_trt, effect_trt)", "cor_h2_trt")


plot_trt_h2_instar2 <- ggplot(data=cor_trt_h2_instar2, aes(x=position, y=cor_h2_trt)) +
  geom_bar(stat="identity") + 
  # ylab("cor")+xlab("dorsal position") + 
  ylim(-1, 1) +
  theme(legend.position="right", 
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

plot_trt_h2_instar2


# trt_H2_instar2

cor_trt_H2_instar2 <- as.data.table(trt_H2_instar2 %>%   ## modify here (!)
                                      dplyr::group_by(section) %>% 
                                      dplyr::summarise(cor(H2_trt, effect_trt)))

cor_trt_H2_instar2[, position := c(50,150,200,250,300,350,400,450,500,100,550)]

setkey(cor_trt_H2_instar2, position)
setnames(cor_trt_H2_instar2, "cor(H2_trt, effect_trt)", "cor_H2_trt")


plot_trt_H2_instar2 <- ggplot(data=cor_trt_H2_instar2, aes(x=position, y=cor_H2_trt)) +
  geom_bar(stat="identity") + 
  # ylab("cor")+xlab("dorsal position") + 
  ylim(-1, 1) +
  theme(legend.position="right", 
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

plot_trt_H2_instar2


# (plot_trt_h2_instar1|plot_trt_h2_instar2) / (plot_trt_H2_instar1|plot_trt_H2_instar2)
(p1_trt/p2_trt/p3_trt)/(plot_trt_h2_instar1|plot_trt_h2_instar2) / (plot_trt_H2_instar1|plot_trt_H2_instar2)


## correlations based on module support from morphometrics analysis

## correlations in modules

trt_h2_instar1 <- effectSize_heritability_wide[group == "cluster O"][, c('i','instar','group','module','effect_trt', 'h2_trt')] %>% 
  dplyr::filter(instar == "instar 1") %>% 
  dplyr::mutate(section = case_when(
    i<=300 ~ "10-300",
    i>300&i<=600 ~ "300-600"
  )) %>% dplyr::tibble()


trt_H2_instar1 <- effectSize_heritability_wide[group == "cluster O"][, c('i','instar','group','module','effect_trt', 'H2_trt')] %>% 
  dplyr::filter(instar == "instar 1") %>% 
  dplyr::mutate(section = case_when(
    i<=300 ~ "10-300",
    i>300&i<=600 ~ "300-600"
  )) %>% dplyr::tibble()


trt_h2_instar2 <- effectSize_heritability_wide[group == "cluster O"][, c('i','instar','group','module','effect_trt', 'h2_trt')] %>% 
  dplyr::filter(instar == "instar 2") %>% 
  dplyr::mutate(section = case_when(
    i<=100 ~ "10-100",
    i>100&i<=300 ~ "100-300",
    i>300&i<=600 ~ "300-600"
  )) %>% dplyr::tibble()

trt_H2_instar2 <- effectSize_heritability_wide[group == "cluster O"][, c('i','instar','group','module','effect_trt', 'H2_trt')] %>% 
  dplyr::filter(instar == "instar 2") %>% 
  dplyr::mutate(section = case_when(
    i<=100 ~ "10-100",
    i>100&i<=300 ~ "100-300",
    i>300&i<=600 ~ "300-600"
  )) %>% dplyr::tibble()


# trt_h2_instar1
cor_trt_h2_instar1 <- as.data.table(trt_h2_instar1 %>%   
                                      dplyr::group_by(section) %>% 
                                      dplyr::summarise(cor(h2_trt, effect_trt)))

cor_trt_h2_instar1[, position := c(200,500)]

setkey(cor_trt_h2_instar1, position)
setnames(cor_trt_h2_instar1, "cor(h2_trt, effect_trt)", "cor_h2_trt")

plot_trt_h2_instar1 <- ggplot(data=cor_trt_h2_instar1, aes(x=position, y=cor_h2_trt)) +
  geom_bar(stat="identity") + 
  # ylab("cor")+xlab("dorsal position") + 
  ylim(-1, 1) +
  theme(legend.position="right", 
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

plot_trt_h2_instar1


# trt_H2_instar1

cor_trt_H2_instar1 <- as.data.table(trt_H2_instar1 %>%   
                                      dplyr::group_by(section) %>% 
                                      dplyr::summarise(cor(H2_trt, effect_trt)))

cor_trt_H2_instar1[, position := c(200,500)]

setkey(cor_trt_H2_instar1, position)
setnames(cor_trt_H2_instar1, "cor(H2_trt, effect_trt)", "cor_H2_trt")

plot_trt_H2_instar1 <- ggplot(data=cor_trt_H2_instar1, aes(x=position, y=cor_H2_trt)) +
  geom_bar(stat="identity") + 
  # ylab("cor")+xlab("dorsal position") + 
  ylim(-1, 1) +
  theme(legend.position="right", 
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

plot_trt_H2_instar1


# trt_h2_instar2

cor_trt_h2_instar2 <- as.data.table(trt_h2_instar2 %>%   ## modify here (!)
                                      dplyr::group_by(section) %>% 
                                      dplyr::summarise(cor(h2_trt, effect_trt)))

cor_trt_h2_instar2[, position := c(50,150,400)]

setkey(cor_trt_h2_instar2, position)
setnames(cor_trt_h2_instar2, "cor(h2_trt, effect_trt)", "cor_h2_trt")


plot_trt_h2_instar2 <- ggplot(data=cor_trt_h2_instar2, aes(x=position, y=cor_h2_trt)) +
  geom_bar(stat="identity") + 
  # ylab("cor")+xlab("dorsal position") + 
  ylim(-1, 1) +
  theme(legend.position="right", 
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

plot_trt_h2_instar2


# trt_H2_instar2

cor_trt_H2_instar2 <- as.data.table(trt_H2_instar2 %>%   ## modify here (!)
                                      dplyr::group_by(section) %>% 
                                      dplyr::summarise(cor(H2_trt, effect_trt)))

cor_trt_H2_instar2[, position := c(50,150,400)]

setkey(cor_trt_H2_instar2, position)
setnames(cor_trt_H2_instar2, "cor(H2_trt, effect_trt)", "cor_H2_trt")


plot_trt_H2_instar2 <- ggplot(data=cor_trt_H2_instar2, aes(x=position, y=cor_H2_trt)) +
  geom_bar(stat="identity") + 
  # ylab("cor")+xlab("dorsal position") + 
  ylim(-1, 1) +
  theme(legend.position="right", 
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

plot_trt_H2_instar2


# (plot_trt_h2_instar1|plot_trt_h2_instar2) / (plot_trt_H2_instar1|plot_trt_H2_instar2)
(p1_trt/p2_trt/p3_trt)/(plot_trt_h2_instar1|plot_trt_h2_instar2) / (plot_trt_H2_instar1|plot_trt_H2_instar2)

