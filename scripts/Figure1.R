## Becker et al - FIG 1


#################
### libraries ###
#################

library(data.table)
library(ggplot2)
library(cowplot)
library(ggbeeswarm)
library(grid)
library(gridExtra)
library(viridis)
library(plyr)
library(dplyr)
library(imager)
library(irr)
library(rstatix)
library(patchwork)


#############
### Fig 1 ###
#############

# validation data; includes filebase info (needed to link validation data to pheno data)
load("data/shape_use_21April2020.Rdata")
shape_use
shape_use[, GenoPLUS := paste(cloneid_geno, "_", barcode, "_", treatment, "_", instar, "_", replicate, sep="")]

shape_use[, filebase_cor := substring(shape_use$filebase, 13)]
shape_use[, maxHeight_mm_GUI := max(height_mm, na.rm=T), by=GenoPLUS]
# use only i>=10 and i<=651, and instars 1 and 2
shape_use_small <- shape_use[i %in% c(10:650)][instar %in% c(1,2)][, list(animal_length_mm=mean(animal_length_mm),
                                                                           maxHeight_mm_GUI=mean(maxHeight_mm_GUI), 
                                                                           nteeth=median(nteeth)),
                                                                    list(filebase,filebase_cor,GenoPLUS,treatment,instar,pixel_to_mm)]



load(file="data/GUI_data_plusManual.RData")

# animal length
GUI_data_plusManual_use_length <- GUI_data_plusManual_use[, c("GenoPLUS","treatment","instar","pixel_to_mm","animal_length_mm","animal_length_manual_mm_DB","animal_length_manual_mm_RP","animal_length_manual_mm_LW")][!is.na(animal_length_manual_mm_DB) & !is.na(animal_length_manual_mm_RP) & !is.na(animal_length_manual_mm_LW)]
GUI_data_plusManual_use_length_use <- melt(GUI_data_plusManual_use_length, id.vars=c("GenoPLUS","treatment","instar","pixel_to_mm","animal_length_mm"))
GUI_data_plusManual_use_length_use[, person := ifelse(variable == 'animal_length_manual_mm_DB', "observer 1",
                                                  ifelse(variable == 'animal_length_manual_mm_RP', "observer 2",
                                                         ifelse(variable == 'animal_length_manual_mm_LW', "observer 3", 'NA')))]
setnames(GUI_data_plusManual_use_length_use, 'value', 'animal_length_manual')


length_plot <- ggplot(data=(GUI_data_plusManual_use_length_use)) +
                      geom_point(data = GUI_data_plusManual_use_length_use, aes(x=animal_length_mm, y=animal_length_manual, col=as.factor(person)), size = 2, pch=21) +
                      facet_wrap(~person, ncol=1) +
                      labs(x = expression(animal~length~"[mm]"[~DAPCHA]), y = expression(animal~length~"[mm]"[~MANUAL])) +

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

length_plot + geom_abline(intercept = 0, slope = 1, size=0.8) + scale_color_manual(values=c("#FF0000","#0000CC","#A0A0A0"))

cor_test_obs1 <- cor.test(GUI_data_plusManual_use_length_use[person == "observer 1"]$animal_length_mm, GUI_data_plusManual_use_length_use[person == "observer 1"]$animal_length_manual)
cor_test_obs2 <- cor.test(GUI_data_plusManual_use_length_use[person == "observer 2"]$animal_length_mm, GUI_data_plusManual_use_length_use[person == "observer 2"]$animal_length_manual)
cor_test_obs3 <- cor.test(GUI_data_plusManual_use_length_use[person == "observer 3"]$animal_length_mm, GUI_data_plusManual_use_length_use[person == "observer 3"]$animal_length_manual)

# (r > 0.9971, 95% CI [0.9967, 0.9979], t(707) > 349, p < 2.2e-16)  # all together

# (r = 0.9971245, 95% CI [0.9966681, 0.9975184], t(707) = 349.86, p < 2.2e-16)  # observer 1
# (r = 0.9975097, 95% CI [0.9971144, 0.9978509], t(707) = 376.06, p < 2.2e-16)  # observer 2
# (r = 0.9972577, 95% CI [0.9968224, 0.9976334], t(707) = 358.29, p < 2.2e-16)  # observer 3



# microstage meter
GUI_data_plusManual_use_microstage <- GUI_data_plusManual_use[, c("GenoPLUS","treatment","instar","pixel_to_mm","micro_length_manual_1_px","micro_length_manual_2_px")][!is.na(micro_length_manual_1_px) & !is.na(micro_length_manual_1_px)]
GUI_data_plusManual_use_microstage_use <- melt(GUI_data_plusManual_use_microstage, id.vars=c("GenoPLUS","treatment","instar","pixel_to_mm"))
GUI_data_plusManual_use_microstage_use[, run := ifelse(variable == 'micro_length_manual_1_px', "run 1", "run 2")]
setnames(GUI_data_plusManual_use_microstage_use, 'value', 'micro_length_manual')


microstage_plot <- ggplot(data=(GUI_data_plusManual_use_microstage_use)) +
                        geom_point(data = GUI_data_plusManual_use_microstage_use, 
                                   aes(x=pixel_to_mm, y=micro_length_manual), pch=8, color='grey50', size=3) +
                        facet_wrap(~run, ncol=1) +
                        labs(x = expression(pixel~to~mm[~DAPCHA]), y = expression(pixel~to~mm[~MANUAL])) +

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

microstage_plot + geom_abline(intercept = 0, slope = 1, size=0.8) + scale_color_manual(values=c("#7F7F7F"))


cor_test_run1 <- cor.test(GUI_data_plusManual_use_microstage_use[run == "run 1"]$pixel_to_mm, GUI_data_plusManual_use_microstage_use[run == "run 1"]$micro_length_manual)
cor_test_run2 <- cor.test(GUI_data_plusManual_use_microstage_use[run == "run 2"]$pixel_to_mm, GUI_data_plusManual_use_microstage_use[run == "run 2"]$micro_length_manual)

# (r > 0.9968, 95% CI [0.9953, 0.9979], t(96) = 123.39, p < 2.2e-16)  # run 1
# (r > 0.7319, 95% CI [0.6241, 0.8123], t(96) = 10.524, p < 2.2e-16)  # run 2



# pedestal score
GUI_data_plusManual_use_pedestal <- GUI_data_plusManual_use[, c("GenoPLUS","treatment","instar","pixel_to_mm","maxHeight_mm_GUI","pedestal_score_manual_DB","pedestal_score_manual_RP","pedestal_score_manual_LW")][!is.na(pedestal_score_manual_DB) & !is.na(pedestal_score_manual_RP) & !is.na(pedestal_score_manual_LW)][treatment %in% c(0,0.5)]
GUI_data_plusManual_use_pedestal_use <- melt(GUI_data_plusManual_use_pedestal, id.vars=c("GenoPLUS","treatment","instar","pixel_to_mm","maxHeight_mm_GUI"))
GUI_data_plusManual_use_pedestal_use[, person := ifelse(variable == 'pedestal_score_manual_DB', "observer 1",
                                                  ifelse(variable == 'pedestal_score_manual_RP', "observer 2",
                                                         ifelse(variable == 'pedestal_score_manual_LW', "observer 3", 'NA')))]
setnames(GUI_data_plusManual_use_pedestal_use, 'value', 'pedestal_score_manual')
GUI_data_plusManual_use_pedestal_use[,treatment_new := ifelse(GUI_data_plusManual_use_pedestal_use$treatment == 0, 'control', 'predation')]


treatment_fill <- c(rep("#FFFFFF",16))

pedestal_plot <- ggplot(data=GUI_data_plusManual_use_pedestal_use, 
                        aes(y=maxHeight_mm_GUI, x=as.character(pedestal_score_manual), group=interaction(as.character(pedestal_score_manual),as.factor(treatment)), color=as.factor(person))) +
                    geom_beeswarm(size = 2, pch=21, cex=1.5, priority='density', dodge.width=0.8) +  
                    geom_boxplot(color='black', fill=treatment_fill, outlier.colour='black', outlier.size = 1, width=0.2) + 
                    facet_grid(~person~treatment_new) +
                    labs(y = expression(max~dorsal~height~"[mm]"[~DAPCHA]), x = expression(max~dorsal~height~"[mm]"[~MANUAL])) +
  
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

pedestal_plot + coord_flip() + scale_x_discrete(limits=c("0","30","50")) + scale_y_continuous(limits=c(0.15,0.36), breaks=c(0.15,0.25,0.35)) + scale_color_manual(values=c("#FF0000","#0000CC","#A0A0A0","#FF0000","#0000CC","#A0A0A0"))
 

GUI_data_plusManual_use_pedestal_use_adj <- unique(GUI_data_plusManual_use_pedestal_use)

# data normally distributed? - NO(!) > use wilcox test
shapiro.test(GUI_data_plusManual_use_pedestal_use_adj$maxHeight_mm_GUI)
shapiro.test(GUI_data_plusManual_use_pedestal_use_adj[treatment_new == "control"]$maxHeight_mm_GUI)
shapiro.test(GUI_data_plusManual_use_pedestal_use_adj[treatment_new == "predation"]$maxHeight_mm_GUI)

#t.test(maxHeight_mm_GUI ~ treatment_new, data = GUI_data_plusManual_use_pedestal_use_adj)
wilcox_treatment <- wilcox.test(maxHeight_mm_GUI ~ treatment_new, data = GUI_data_plusManual_use_pedestal_use_adj) 

Zscore = qnorm(wilcox_treatment$p.value/2)
Zscore

#rstat <- abs(Zscore)/sqrt(867) # alternative calc for r
wilcox_effectsize_treatment <- wilcox_effsize(maxHeight_mm_GUI ~ treatment_new, data = GUI_data_plusManual_use_pedestal_use_adj)
wilcox_effectsize_treatment

# (Z = -9.009, r =  0.307, p-value < 2.2e-16)


# use only 2nd instar data here
ggplot(data = GUI_data_plusManual_use_pedestal_use_adj[instar==2], aes(x=maxHeight_mm_GUI, y=pedestal_score_manual)) + geom_point() + geom_smooth(method="lm")
cor_test_pedestal <- cor.test(GUI_data_plusManual_use_pedestal_use_adj$maxHeight_mm_GUI, GUI_data_plusManual_use_pedestal_use_adj$pedestal_score_manual, method = "spearman", exact=F)
# (r = 0.2349, p < 2.2e-16)  

cor_test_pedestal <- cor.test(GUI_data_plusManual_use_pedestal_use_adj[instar==2]$maxHeight_mm_GUI, GUI_data_plusManual_use_pedestal_use_adj[instar==2]$pedestal_score_manual, method = "spearman", exact=F)
# (r = 0.6367, p < 2.2e-16)  


GUI_data_plusManual_use_pedestal_use_wide <- spread(GUI_data_plusManual_use_pedestal_use_adj[, -"variable"], person, pedestal_score_manual)

GUI_data_plusManual_use_pedestal_use_wide_small <- GUI_data_plusManual_use_pedestal_use_wide[, c('observer 1','observer 2','observer 3')]
icc(GUI_data_plusManual_use_pedestal_use_wide_small, model="twoway", type="agreement",unit="single")

GUI_data_plusManual_use_pedestal_use_wide_small <- GUI_data_plusManual_use_pedestal_use_wide[treatment_new == "control"][, c('observer 1','observer 2','observer 3')]
icc(GUI_data_plusManual_use_pedestal_use_wide_small, model="twoway", type="agreement",unit="single")

GUI_data_plusManual_use_pedestal_use_wide_small <- GUI_data_plusManual_use_pedestal_use_wide[treatment_new == "predation"][, c('observer 1','observer 2','observer 3')]
icc(GUI_data_plusManual_use_pedestal_use_wide_small, model="twoway", type="agreement",unit="single")


# (icc = 0.87, 95% CI [0.844, 0.893], F(286,506) = 21.6, p = 1.96e-178)  # all 
# (icc = 0.497, 95% CI [0.405, 0.586], F(156,305) = 4.03, p = 1.44e-25)  # control 
# (icc = 0.817, 95% CI [0.764, 0.861], F(129,251) = 14.6, p = 1.17e-69)  # predation 



#instar assignment

GUI_data_plusManual_use_instar <- GUI_data_plusManual_use[treatment %in% c(0,0.5)][, c("GenoPLUS","treatment","instar","animal_length_mm")][!is.na(animal_length_mm)]
GUI_data_plusManual_use_instar[,treatment_new := ifelse(GUI_data_plusManual_use_instar$treatment == 0, 'control', 'predation')]

group_mean_length <- ddply(GUI_data_plusManual_use_instar, c("instar", "treatment"), summarise, grp.mean=mean(animal_length_mm))


instar_plot <- ggplot(data = GUI_data_plusManual_use_instar, 
                      aes(x=animal_length_mm, color=as.factor(instar))) +
                  geom_density(aes(y= ..count..), size=1) + 
                  #geom_vline(data=group_mean_length, 
                  #           aes(xintercept=grp.mean), linetype="dashed", size = 0.5, colour="black") + 
                  facet_wrap(~treatment_new, ncol=1) +
                  labs(x = expression(animal~length~"[mm]"[~DAPCHA]), y = "density") +

                  theme(legend.position="none", 
                        rect = element_rect(fill = "transparent"),
                        panel.grid.major = element_line(colour = "grey70", size=0.25),
                        panel.grid.minor = element_line(colour = "grey90", size=0.1),
                        panel.background = element_rect(fill = "transparent",colour = NA),
                        plot.background = element_rect(fill = "transparent",colour = NA), 
                        #strip.text.x = element_blank(),
                        axis.text.y = element_blank(), 
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


instar_plot + scale_color_manual(values=c("#7F7F7F","#0000CC","#7F7F7F","#0000CC"))



########
########

length_plot + geom_abline(intercept = 0, slope = 1, size=0.8) + scale_color_manual(values=c("#FF0000","#0000CC","#A0A0A0"))

microstage_plot + geom_abline(intercept = 0, slope = 1, size=0.8) + scale_color_manual(values=c("#7F7F7F"))

pedestal_plot + coord_flip() + scale_x_discrete(limits=c("0","30","50")) + scale_y_continuous(limits=c(0.15,0.31), breaks=c(0.15,0.25,0.35)) + scale_color_manual(values=c("#FF0000","#0000CC","#A0A0A0","#FF0000","#0000CC","#A0A0A0"))

instar_plot + scale_color_manual(values=c("#7F7F7F","#0000CC","#7F7F7F","#0000CC"))



patchwork_plots_induction <- length_plot + geom_abline(intercept = 0, slope = 1, size=0.8) + scale_color_manual(values=c("#FF0000","#0000CC","#A0A0A0")) + 
  microstage_plot + geom_abline(intercept = 0, slope = 1, size=0.8) + scale_color_manual(values=c("#7F7F7F")) + 
  pedestal_plot + coord_flip() + scale_x_discrete(limits=c("0","30","50")) + scale_y_continuous(limits=c(0.15,0.31), breaks=c(0.15,0.25,0.35)) + scale_color_manual(values=c("#FF0000","#0000CC","#A0A0A0","#FF0000","#0000CC","#A0A0A0")) + 
  instar_plot + scale_color_manual(values=c("#7F7F7F","#0000CC","#7F7F7F","#0000CC")) + 
  plot_layout(ncol=4, widths = c(1,1,1.5,1))
                                            
patchwork_plots_induction

