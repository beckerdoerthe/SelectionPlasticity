# Becker et al - FIG 2 & SUPPL FIG 5


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

load(file = "data/all_data_final.RData")
# all_data_final

# for trait plots, use one i-th position (e.g., 150) due to data dublication re i-th positions...

## AVERAGE SHAPE 
all_data_final.ag <- all_data_final[i <= 600][, list(height = mean(height)), list(i, treatment, instar_new, SC_group_new)]

ag.line_plot <- ggplot(data = all_data_final.ag[SC_group_new == "cluster A"],  aes(x=i, y=height)) +   
  
                    geom_line(data = all_data_final.ag[treatment == 0][SC_group_new == "cluster A"],  aes(x=i, y=height), size = 1.5, colour = "black") + 
                    geom_line(data = all_data_final.ag[treatment == 0.5][SC_group_new == "cluster A"],  aes(x=i, y=height), size = 1.5, colour = "red") + 
                    
                    ylim(0, 0.28) +
                    
                    facet_wrap(~instar_new, ncol=1) +
                    
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


# tiff(file = "avg_shape_plot.tiff", width = 3200, height = 3200, units = "px", res = 800)
ag.line_plot  + labs(x = "dorsal position", y = "dorsal height (mm)") + scale_x_continuous(limits=c(0,621), breaks=c(0,300,600)) + scale_color_manual(values=c("#000000","#000000"))
# dev.off()

## FREQUENCY MAX HEIGHT
hist_plot <- ggplot(data = all_data_final[i <= 600][SC_group_new == "cluster A"][max_height_new == height],  aes(x=i)) +   
  
                    geom_histogram(data=all_data_final[i <= 600][SC_group_new == "cluster A"][max_height_new == height][treatment == 0.5], aes(x=i), bins = 50, fill='#ff4c4c', color='#FF0000', alpha=0.2, size=0.5) + 
                    geom_histogram(data=all_data_final[i <= 600][SC_group_new == "cluster A"][max_height_new == height][treatment == 0], aes(x=i), bins = 50, fill='#C0C0C0', color='#000000', alpha=0.2, size=0.5) + 
                    
                    facet_wrap(~instar_new, ncol=1) +
                    
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

# tiff(file = "freq_maxHeight_plot.tiff", width = 3200, height = 3200, units = "px", res = 800)
hist_plot + labs(x = "dorsal position", y = "frequency") + scale_x_continuous(breaks=c(0,150,300))
# dev.off()


# chisq 
chisq.test(table(all_data_final[i <= 300][SC_group_new == "cluster O"][max_height_new == height][instar_new == "instar 1"]$i, all_data_final[i <= 300][SC_group_new == "cluster O"][max_height_new == height][instar_new == "instar 1"]$treatment_new))
# X(df=121, N=411) = 210.32, p-value = 8.926e-07
chisq.test(table(all_data_final[i <= 300][SC_group_new == "cluster O"][max_height_new == height][instar_new == "instar 2"]$i, all_data_final[i <= 300][SC_group_new == "cluster O"][max_height_new == height][instar_new == "instar 2"]$treatment_new))
# X(df=112, N=380) = 222.4, p-value = 2.676e-09       


chisq.test(table(all_data_final[i <= 300][SC_group_new == "cluster A"][max_height_new == height][instar_new == "instar 1"]$i, all_data_final[i <= 300][SC_group_new == "cluster A"][max_height_new == height][instar_new == "instar 1"]$treatment_new))
# X(df=119, N=463) = 284.84, p-value = 1.267e-15
chisq.test(table(all_data_final[i <= 300][SC_group_new == "cluster A"][max_height_new == height][instar_new == "instar 2"]$i, all_data_final[i <= 300][SC_group_new == "cluster A"][max_height_new == height][instar_new == "instar 2"]$treatment_new))
# X(df=106, N=433) = 228.35, p-value = 5.758e-11                


## MAX HEIGHT
dorsal_height_plot <- ggplot(data = all_data_final[i==150][SC_group == "A"],   
                             aes(y=max_height_new, x=as.factor(treatment_new), color=as.factor(treatment_new))) + 
                        geom_beeswarm(size = 0.5, cex=1.5, priority='density', dodge.width=1) +  
                        geom_boxplot(fill="white", width=0.3, outlier.colour='grey', outlier.size = 0.1) +
                        facet_wrap(~instar_new, ncol = 1) +
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


# tiff(file = "dorsal_height_plot.tiff", width = 3200, height = 3200, units = "px", res = 800)
dorsal_height_plot + labs(x = "treatment", y = paste("max dorsal height [mm]")) + scale_color_manual(values=c("#000000","#FF0000"))
# dev.off()


shapiro.test(all_data_final[i <= 600][SC_group_new == "cluster O"][max_height_new == height][instar_new == "instar 1"][treatment_new == "C"]$max_height_new)  # normal
shapiro.test(all_data_final[i <= 600][SC_group_new == "cluster O"][max_height_new == height][instar_new == "instar 1"][treatment_new == "P"]$max_height_new)  # not normal
shapiro.test(all_data_final[i <= 600][SC_group_new == "cluster O"][max_height_new == height][instar_new == "instar 2"][treatment_new == "C"]$max_height_new)  # normal
shapiro.test(all_data_final[i <= 600][SC_group_new == "cluster O"][max_height_new == height][instar_new == "instar 2"][treatment_new == "P"]$max_height_new)  # normal

shapiro.test(all_data_final[i <= 600][SC_group_new == "cluster A"][max_height_new == height][instar_new == "instar 1"][treatment_new == "C"]$max_height_new)  # not normal
shapiro.test(all_data_final[i <= 600][SC_group_new == "cluster A"][max_height_new == height][instar_new == "instar 1"][treatment_new == "P"]$max_height_new)  # normal
shapiro.test(all_data_final[i <= 600][SC_group_new == "cluster A"][max_height_new == height][instar_new == "instar 2"][treatment_new == "C"]$max_height_new)  # normal
shapiro.test(all_data_final[i <= 600][SC_group_new == "cluster A"][max_height_new == height][instar_new == "instar 2"][treatment_new == "P"]$max_height_new)  # normal


# O_I1
wilcox_O_I1 <- wilcox.test(max_height_new ~ treatment_new, data = all_data_final[i <= 600][SC_group_new == "cluster O"][max_height_new == height][instar_new == "instar 1"]) 
wilcox_O_I1

Zscore_O_I1 = qnorm(wilcox_O_I1$p.value/2)
Zscore_O_I1

wilcox_effectsize_O_I1 <- wilcox_effsize(max_height_new ~ treatment_new, data = all_data_final[i <= 600][SC_group_new == "cluster O"][max_height_new == height][instar_new == "instar 1"])
wilcox_effectsize_O_I1

# (Z = -8.95855, p-value < 2.2e-16)


# O_I2
wilcox_O_I2 <- wilcox.test(max_height_new ~ treatment_new, data = all_data_final[i <= 600][SC_group_new == "cluster O"][max_height_new == height][instar_new == "instar 2"]) 
wilcox_O_I2

Zscore_O_I2 = qnorm(wilcox_O_I2$p.value/2)
Zscore_O_I2

wilcox_effectsize_O_I2 <- wilcox_effsize(max_height_new ~ treatment_new, data = all_data_final[i <= 600][SC_group_new == "cluster O"][max_height_new == height][instar_new == "instar 2"])
wilcox_effectsize_O_I2

# (Z = -9.83551, p-value < 2.2e-16)


# A_I1
wilcox_A_I1 <- wilcox.test(max_height_new ~ treatment_new, data = all_data_final[i <= 600][SC_group_new == "cluster A"][max_height_new == height][instar_new == "instar 1"]) 
wilcox_A_I1

Zscore_A_I1 = qnorm(wilcox_A_I1$p.value/2)
Zscore_A_I1

wilcox_effectsize_A_I1 <- wilcox_effsize(max_height_new ~ treatment_new, data = all_data_final[i <= 600][SC_group_new == "cluster A"][max_height_new == height][instar_new == "instar 1"])
wilcox_effectsize_A_I1

# (Z = -8.722, p-value < 2.2e-16)


# A_I2
wilcox_A_I2 <- wilcox.test(max_height_new ~ treatment_new, data = all_data_final[i <= 600][SC_group_new == "cluster A"][max_height_new == height][instar_new == "instar 2"]) 

Zscore_A_I2 = qnorm(wilcox_A_I2$p.value/2)
Zscore_A_I2

wilcox_effectsize_A_I2 <- wilcox_effsize(max_height_new ~ treatment_new, data = all_data_final[i <= 600][SC_group_new == "cluster A"][max_height_new == height][instar_new == "instar 2"])
wilcox_effectsize_A_I2

# (Z = -10.455, p-value < 2.2e-16)


##NECKTEETH
# make proportion/percentage table
nteeth_number <- as.data.table(
                      all_data_final[i == 150] %>% 
                      group_by(SC_group_new, instar_new, treatment_new, nteeth = round(nteeth)) %>%
                      summarize(N = n())) 

nteeth_number[, totalN := case_when(SC_group_new=='cluster O' & instar_new=='instar 1' & treatment_new=='C' ~ sum(nteeth_number[SC_group_new=='cluster O'][instar_new=='instar 1'][treatment_new=='C']$N),
                                    SC_group_new=='cluster O' & instar_new=='instar 1' & treatment_new=='P' ~ sum(nteeth_number[SC_group_new=='cluster O'][instar_new=='instar 1'][treatment_new=='P']$N),
                                    SC_group_new=='cluster O' & instar_new=='instar 2' & treatment_new=='C' ~ sum(nteeth_number[SC_group_new=='cluster O'][instar_new=='instar 2'][treatment_new=='C']$N),
                                    SC_group_new=='cluster O' & instar_new=='instar 2' & treatment_new=='P' ~ sum(nteeth_number[SC_group_new=='cluster O'][instar_new=='instar 2'][treatment_new=='P']$N),
                                    
                                    SC_group_new=='cluster A' & instar_new=='instar 1' & treatment_new=='C' ~ sum(nteeth_number[SC_group_new=='cluster A'][instar_new=='instar 1'][treatment_new=='C']$N),
                                    SC_group_new=='cluster A' & instar_new=='instar 1' & treatment_new=='P' ~ sum(nteeth_number[SC_group_new=='cluster A'][instar_new=='instar 1'][treatment_new=='P']$N),
                                    SC_group_new=='cluster A' & instar_new=='instar 2' & treatment_new=='C' ~ sum(nteeth_number[SC_group_new=='cluster A'][instar_new=='instar 2'][treatment_new=='C']$N),
                                    SC_group_new=='cluster A' & instar_new=='instar 2' & treatment_new=='P' ~ sum(nteeth_number[SC_group_new=='cluster A'][instar_new=='instar 2'][treatment_new=='P']$N))]
 
nteeth_number[, proportion_nteeth := N/totalN]
nteeth_number[, percent_nteeth := proportion_nteeth * 100]    

nteeth_plot <- ggplot(data = nteeth_number[SC_group_new == "cluster A"], 
                         aes(y=percent_nteeth, x=nteeth, fill=as.factor(treatment_new))) + 
                  geom_bar(stat="identity", position = position_dodge2(preserve = "single")) + 
                  facet_wrap(~instar_new, ncol=1) +
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

# tiff(file = "nteeth_plot.tiff", width = 3200, height = 3200, units = "px", res = 800)
nteeth_plot + labs(x = "# of neckteeth", y = "[%]") + scale_x_continuous(limits=c(-0.5,4.5)) + scale_fill_manual(values=c("#000000","#FF0000"))
# dev.off()


chisq.test(table(all_data_final[i == 150][SC_group_new == "cluster O"][instar_new == "instar 1"]$nteeth, all_data_final[i == 150][SC_group_new == "cluster O"][instar_new == "instar 1"]$treatment_new))
# X(df=6, N=412) = 7.8407, p-value = 0.25
chisq.test(table(all_data_final[i == 150][SC_group_new == "cluster O"][instar_new == "instar 2"]$nteeth, all_data_final[i == 150][SC_group_new == "cluster O"][instar_new == "instar 2"]$treatment_new))
# X(df=7, N=381) = 179.83, p-value < 2.2e-16

chisq.test(table(all_data_final[i == 150][SC_group_new == "cluster A"][instar_new == "instar 1"]$nteeth, all_data_final[i == 150][SC_group_new == "cluster A"][instar_new == "instar 1"]$treatment_new))
# X(df=4, N=464) = 10.588, p-value = 0.0316
chisq.test(table(all_data_final[i == 150][SC_group_new == "cluster A"][instar_new == "instar 2"]$nteeth, all_data_final[i == 150][SC_group_new == "cluster A"][instar_new == "instar 2"]$treatment_new))
# X(df=7, N=431) = 242.64, p-value < 2.2e-16



########
########

ag.line_plot  + labs(x = "dorsal position", y = "dorsal height (mm)") + scale_x_continuous(limits=c(0,621), breaks=c(0,300,600)) + scale_color_manual(values=c("#000000","#000000"))

hist_plot + labs(x = "dorsal position", y = "frequency") + scale_x_continuous(breaks=c(0,150,300))

dorsal_height_plot + labs(x = "treatment", y = paste("max dorsal height [mm]")) + scale_color_manual(values=c("#000000","#FF0000"))

nteeth_plot + labs(x = "# of neckteeth", y = "[%]") + scale_x_continuous(limits=c(-0.5,4.5)) + scale_fill_manual(values=c("#000000","#FF0000"))


patchwork_plots_induction <- ag.line_plot + labs(x = "dorsal position", y = "dorsal height (mm)") + scale_x_continuous(limits=c(0,621), breaks=c(0,300,600)) + scale_color_manual(values=c("#000000","#000000")) + 
  hist_plot + labs(x = "dorsal position", y = "frequency") + scale_x_continuous(breaks=c(0,150,300)) + 
  dorsal_height_plot + labs(x = "treatment", y = paste("max dorsal height [mm]")) + scale_color_manual(values=c("#000000","#FF0000")) +
  nteeth_plot + labs(x = "# of neckteeth", y = "[%]") + scale_x_continuous(limits=c(-0.5,4.5)) + scale_fill_manual(values=c("#000000","#FF0000")) +
  plot_layout(ncol=4, widths = c(1,1,1,1))

patchwork_plots_induction




## ANOVA ##

# height_ANOVA <- foreach(group.i=unique(all_data_final$SC_group_new), .errorhandling="remove", .combine="rbind") %do% {
#   
#                   foreach(instar.i=unique(all_data_final$instar_new), .errorhandling="remove", .combine="rbind") %do% {
#                     
#                   tmp_dat_use <- all_data_final[SC_group_new == group.i][instar_new==instar.i]
#                     
#                     # i
#                     foreach(i.i=unique(tmp_dat_use$i), .errorhandling="remove", .combine="rbind") %dopar% {
#                       if(i.i%%10==0) print(i.i)
#                       
#                       tmp <- tmp_dat_use[i==i.i]
#                       
#                       # tmp.lm <- lm(height~cloneid_geno*treatment+batch, tmp)
#                       # tmp.II.aov <- car::Anova(tmp.fit, type = 2)
#                       
#                       tmp.fit <- aov(height~cloneid_geno*treatment+batch, tmp)
# 
#                       etasq_out <- as.data.table(eta_squared(tmp.fit, partial=TRUE, ci = .99))
#                       # etasq_out <- eta_sq(tmp.lm, partial=TRUE, ci.lvl = .99)
#                       
#                       #
#                       data.table(i=i.i,
#                                  instar = instar.i,
#                                  group = group.i,
#                                  
#                                  partial_eta = unlist(c(etasq_out[1,2],
#                                                  etasq_out[2,2],
#                                                  etasq_out[3,2],
#                                                  etasq_out[4,2])),
#                                  
#                                  lCI_eta = unlist(c(etasq_out[1,4],
#                                              etasq_out[2,4],
#                                              etasq_out[3,4],
#                                              etasq_out[4,4])),
#                                  
#                                  uCI_eta = unlist(c(etasq_out[1,5],
#                                              etasq_out[2,5],
#                                              etasq_out[3,5],
#                                              etasq_out[4,5])),             	
#                                  
#                                  term = c("geno", "trt", "batch", "gxe"))  
#                       
#                       
#                     }
#                   }
#                 }
# 
# 
# save(height_ANOVA, file = "aov_out.RData")
load(file = "output/aov_out.RData")


aov_out_plot <- ggplot(data = height_ANOVA[term %in% c('geno','gxe','trt')][group == "cluster O"],     
                       aes(x=i, y=partial_eta, group=term, color=as.factor(term))) +  
                geom_line(size = 2) +
                geom_ribbon(data=height_ANOVA[term %in% c('geno','gxe','trt')][group == "cluster O"], aes(ymin=lCI_eta, ymax=uCI_eta, colour = NA, fill=as.factor(height_ANOVA[term %in% c('geno','gxe','trt')][group == "cluster O"]$term)), alpha=0.3) +
                facet_wrap(~instar, ncol=1) +
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
  
aov_out_plot + labs(x = "dorsal position", y = "effect size") + scale_color_manual(values=c("#0000CC","#A0A0A0","#FF0000")) + scale_fill_manual(values = c("#0000CC","#A0A0A0","#FF0000")) + scale_x_continuous(limits=c(0,600), breaks=c(0,300,600))


## narrow sense heritability 

# narrow sense heritability - O cluster (GCTA output)
# control & treatment

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


h2_instar_i_plot <- ggplot(data = O_h2_all,  aes(x=i, y=h2)) +   
  
                    geom_rect(data = O_h2_all[instar == "instar 1"], aes(xmin = 10, xmax = 300, ymin = -Inf, ymax = Inf),fill = "#F5F5F5", colour="#F5F5F5", alpha=0.1) + 
                    geom_rect(data = O_h2_all[instar == "instar 1"], aes(xmin = 301, xmax = 600, ymin = -Inf, ymax = Inf),fill = "#DCDCDC", colour="#DCDCDC", alpha=0.1) + 
                    
                    geom_rect(data = O_h2_all[instar == "instar 2"], aes(xmin = 10, xmax = 100, ymin = -Inf, ymax = Inf),fill = "#F5F5F5", colour="#F5F5F5", alpha=0.1) + 
                    geom_rect(data = O_h2_all[instar == "instar 2"], aes(xmin = 101, xmax = 200, ymin = -Inf, ymax = Inf),fill = "#E8E8E8", colour="#E8E8E8", alpha=0.1) + 
                    geom_rect(data = O_h2_all[instar == "instar 2"], aes(xmin = 201, xmax = 600, ymin = -Inf, ymax = Inf),fill = "#DCDCDC", colour="#DCDCDC", alpha=0.1) + 
                    
                    geom_ribbon(data=O_h2_all[i <= 600][treatment == 'predation'], aes(ymin=h2-SE, ymax=h2+SE, x=i), fill = "#FF0000", alpha = 0.3)+
                    geom_line(data=O_h2_all[i <= 600][treatment == 'predation'], aes(x=i, y=h2), size = 1, colour = "#FF0000") + 
                    
                    geom_ribbon(data=O_h2_all[i <= 600][treatment == 'control'], aes(ymin=h2-SE, ymax=h2+SE, x=i), fill = "#000000", alpha = 0.3)+
                    geom_line(data=O_h2_all[i <= 600][treatment == 'control'], aes(x=i, y=h2), size = 1, colour = "#000000") + 

                    #xlim(10,600) + 
                    #ylim(0,0.6) +  
  
                    geom_hline(yintercept=0, linetype = "dotted", size = 0.5, colour = "black") +
  
                    facet_wrap(~instar, ncol=1) +
                    
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

h2_instar_i_plot + labs(x = "dorsal position", y = expression(heritability~(h^{2}))) + scale_x_continuous(limits=c(0,600), breaks=c(0,300,600))


shapiro.test(O_h2_all[i >= 100 & i <= 250][instar == "instar 1"][treatment == 'control']$h2)
shapiro.test(O_h2_all[i >= 100 & i <= 250][instar == "instar 2"][treatment == 'predation']$h2)
shapiro.test(O_h2_all[i >= 100 & i <= 250][instar == "instar 1"][treatment == 'control']$h2)
shapiro.test(O_h2_all[i >= 100 & i <= 250][instar == "instar 2"][treatment == 'predation']$h2)


# O_I1_plast
wilcox_O_I1_plast <- wilcox.test(h2 ~ treatment, data =O_h2_all[i >= 100 & i <= 250][instar == "instar 1"]) 
wilcox_O_I1_plast

Zscore_O_I1_plast = qnorm(wilcox_O_I1_plast$p.value/2)
Zscore_O_I1_plast

wilcox_effectsize_O_I1_plast <- wilcox_effsize(h2 ~ treatment, data =O_h2_all[i >= 100 & i <= 250][instar == "instar 1"])
wilcox_effectsize_O_I1_plast

# (Z = -5.071395, p-value = 3.949e-07)


# O_I2_plast
wilcox_O_I2_plast <- wilcox.test(h2 ~ treatment, data =O_h2_all[i >= 100 & i <= 250][instar == "instar 2"]) 
wilcox_O_I2_plast

Zscore_O_I2_plast = qnorm(wilcox_O_I2_plast$p.value/2)
Zscore_O_I2_plast

wilcox_effectsize_O_I2_plast <- wilcox_effsize(h2 ~ treatment, data =O_h2_all[i >= 100 & i <= 250][instar == "instar 2"])
wilcox_effectsize_O_I2_plast

# (Z = -10.40506, p-value < 2.2e-16)



# braod sense heritability - A and O & I1 and I2 (MCMC)
# this is the MCMCglmm output w/ batch as random effect
# control & treatment

load(file = "output/H2_O_I1_batch.RData")
H2_O_I1 <- as.data.table(collected2)
H2_O_I1[, i:= as.numeric(i)]
H2_O_I1[, instar:= "instar 1"]
H2_O_I1[, group:= "O"]
setkey(H2_O_I1, group,instar,i)

load(file = "output/H2_O_I2_batch.RData")
H2_O_I2 <- as.data.table(collected2)
H2_O_I2[, i:= as.numeric(i)]
H2_O_I2[, instar:= "instar 2"]
H2_O_I2[, group:= "O"]
setkey(H2_O_I2, group,instar,i)

load(file = "output/H2_A_I1_batch.RData")
H2_A_I1 <- as.data.table(collected2)
H2_A_I1[, i:= as.numeric(i)]
H2_A_I1[, instar:= "instar 1"]
H2_A_I1[, group:= "A"]
setkey(H2_A_I1, group,instar,i)

load(file = "output/H2_A_I2_batch.RData")
H2_A_I2 <- as.data.table(collected2)
H2_A_I2[, i:= as.numeric(i)]
H2_A_I2[, instar:= "instar 2"]
H2_A_I2[, group:= "A"]
setkey(H2_A_I2, group,instar,i)

AandO_H2_all <- rbind(H2_O_I1,
                      H2_O_I2,
                      H2_A_I1,
                      H2_A_I2)


H2_instar_i_plot <- ggplot(data = AandO_H2_all[i <= 600][label %in% c("C","T")][group == "O"]) +   
  
                    geom_rect(data = AandO_H2_all[instar == "instar 1"], aes(xmin = 10, xmax = 300, ymin = -Inf, ymax = Inf),fill = "#F5F5F5", colour="#F5F5F5", alpha=0.1) + 
                    geom_rect(data = AandO_H2_all[instar == "instar 1"], aes(xmin = 301, xmax = 600, ymin = -Inf, ymax = Inf),fill = "#DCDCDC", colour="#DCDCDC", alpha=0.1) + 
                    
                    geom_rect(data = AandO_H2_all[instar == "instar 2"], aes(xmin = 10, xmax = 100, ymin = -Inf, ymax = Inf),fill = "#F5F5F5", colour="#F5F5F5", alpha=0.1) + 
                    geom_rect(data = AandO_H2_all[instar == "instar 2"], aes(xmin = 101, xmax = 200, ymin = -Inf, ymax = Inf),fill = "#E8E8E8", colour="#E8E8E8", alpha=0.1) + 
                    geom_rect(data = AandO_H2_all[instar == "instar 2"], aes(xmin = 201, xmax = 600, ymin = -Inf, ymax = Inf),fill = "#DCDCDC", colour="#DCDCDC", alpha=0.1) + 
                    
                    # geom_rect(aes(xmin = 10, xmax = 100, ymin = -Inf, ymax = Inf),fill = "#F5F5F5", colour="#F5F5F5", alpha=0.1) + 
                    # geom_rect(aes(xmin = 101, xmax = 200, ymin = -Inf, ymax = Inf),fill = "#E8E8E8", colour="#E8E8E8", alpha=0.1) + 
                    # geom_rect(aes(xmin = 201, xmax = 600, ymin = -Inf, ymax = Inf),fill = "#DCDCDC", colour="#DCDCDC", alpha=0.1) + 
                    
                    geom_ribbon(data=AandO_H2_all[i <= 600][label == 'T'][group == "O"], aes(ymin=stuff.lCI, ymax=stuff.uCI, x=i), fill = "#FF0000", alpha = 0.3)+
                    geom_line(data=AandO_H2_all[i <= 600][label == 'T'][group == "O"], aes(x=i, y=stuff.mean), size = 1, colour = "#FF0000") + 
                    
                    geom_ribbon(data=AandO_H2_all[i <= 600][label == 'C'][group == "O"], aes(ymin=stuff.lCI, ymax=stuff.uCI, x=i), fill = "#000000", alpha = 0.3)+
                    geom_line(data=AandO_H2_all[i <= 600][label == 'C'][group == "O"], aes(x=i, y=stuff.mean), size = 1, colour = "#000000") + 
  
                    geom_hline(yintercept=0, linetype = "dotted", size = 0.5, colour = "black") +
                    
                    #ylim(0,0.58) +   
  
                    facet_wrap(~instar, ncol=1) +
  
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
  
H2_instar_i_plot + labs(x = "dorsal position", y = expression(heritability~(H^{2}))) + scale_x_continuous(limits=c(10,600), breaks=c(0,300,600)) 


shapiro.test(AandO_H2_all[i >= 100 & i <= 300][instar == "instar 1"][group == "O"][label == 'C']$stuff.mean)
shapiro.test(AandO_H2_all[i >= 100 & i <= 300][instar == "instar 2"][group == "O"][label == 'T']$stuff.mean)
shapiro.test(AandO_H2_all[i >= 100 & i <= 300][instar == "instar 1"][group == "A"][label == 'C']$stuff.mean)
shapiro.test(AandO_H2_all[i >= 100 & i <= 300][instar == "instar 2"][group == "A"][label == 'T']$stuff.mean)



# O_I1_plast
wilcox_O_I1_plast <- wilcox.test(stuff.mean ~ label, data = AandO_H2_all[i >= 100 & i <= 250][instar == "instar 1"][label %in% c("C","T")][group == "O"]) 
wilcox_O_I1_plast

Zscore_O_I1_plast = qnorm(wilcox_O_I1_plast$p.value/2)
Zscore_O_I1_plast

wilcox_effectsize_O_I1_plast <- wilcox_effsize(stuff.mean ~ label, data =AandO_H2_all[i >= 100 & i <= 250][instar == "instar 1"][label %in% c("C","T")][group == "O"])
wilcox_effectsize_O_I1_plast

# (Z = -10.99549, p-value < 2.2e-16)


# O_I2_plast
wilcox_O_I2_plast <- wilcox.test(stuff.mean ~ label, data =AandO_H2_all[i >= 100 & i <= 250][instar == "instar 2"][label %in% c("C","T")][group == "O"]) 
wilcox_O_I2_plast

Zscore_O_I2_plast = qnorm(wilcox_O_I2_plast$p.value/2)
Zscore_O_I2_plast

wilcox_effectsize_O_I2_plast <- wilcox_effsize(stuff.mean ~ label, data =AandO_H2_all[i >= 100 & i <= 250][instar == "instar 2"][label %in% c("C","T")][group == "O"])
wilcox_effectsize_O_I2_plast

# (Z = -7.128684, p-value = 1.013e-12)


  
########
########


aov_out_plot + labs(x = "dorsal position", y = "effect size") + scale_color_manual(values=c("#0000CC","#A0A0A0","#FF0000")) + scale_fill_manual(values = c("#0000CC","#A0A0A0","#FF0000")) + scale_x_continuous(limits=c(0,600), breaks=c(0,300,600)) + scale_y_continuous(limits=c(0,0.6), breaks=c(0,0.2,0.4,0.6))

H2_instar_i_plot + labs(x = "dorsal position", y = expression(heritability~(H^{2}))) + scale_x_continuous(limits=c(10,600), breaks=c(0,300,600)) + scale_y_continuous(limits=c(-0.05,0.6), breaks=c(0,0.2,0.4,0.6))

h2_instar_i_plot + labs(x = "dorsal position", y = expression(heritability~(h^{2}))) + scale_x_continuous(limits=c(0,600), breaks=c(0,300,600)) + scale_y_continuous(limits=c(-0.05,0.6), breaks=c(0,0.2,0.4,0.6))


patchwork_plots_induction2 <- aov_out_plot + labs(x = "dorsal position", y = "effect size") + scale_color_manual(values=c("#0000CC","#A0A0A0","#FF0000")) + scale_fill_manual(values = c("#0000CC","#A0A0A0","#FF0000")) + scale_x_continuous(limits=c(0,600), breaks=c(0,300,600)) + scale_y_continuous(limits=c(0,0.6), breaks=c(0,0.2,0.4,0.6)) + 
  H2_instar_i_plot + labs(x = "dorsal position", y = expression(heritability~(H^{2}))) + scale_x_continuous(limits=c(10,600), breaks=c(0,300,600)) + scale_y_continuous(limits=c(-0.05,0.6), breaks=c(0,0.2,0.4,0.6)) +
  h2_instar_i_plot + labs(x = "dorsal position", y = expression(heritability~(h^{2}))) + scale_x_continuous(limits=c(0,600), breaks=c(0,300,600)) + scale_y_continuous(limits=c(-0.05,0.6), breaks=c(0,0.2,0.4,0.6)) +
  # plot_spacer() +
  plot_layout(ncol=3, widths = c(1,1,1))

patchwork_plots_induction2


