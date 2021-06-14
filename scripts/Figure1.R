## VALIDATION - PHENO DATA MORPHO INDUCTION MIDGE UVA 2017


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


##################
### validation ###
##################

# validation data; includes filebase info (needed to link validation data to pheno data)
load("~/Google Drive/PNAS_Becker_data_scripts/shape_use_21April2020.Rdata")
shape_use
shape_use[, GenoPLUS := paste(cloneid_geno, "_", barcode, "_", treatment, "_", instar, "_", replicate, sep="")]

shape_use[, filebase_cor := substring(shape_use$filebase, 13)]
shape_use[, maxHeight_mm_GUI := max(height_mm, na.rm=T), by=GenoPLUS]
# use only i>=10 and i<=651, and instars 1 and 2
shape_use_small <- shape_use[i %in% c(10:650)][instar %in% c(1,2)][, list(animal_length_mm=mean(animal_length_mm),
                                                                           maxHeight_mm_GUI=mean(maxHeight_mm_GUI), 
                                                                           nteeth=median(nteeth)),
                                                                    list(filebase,filebase_cor,GenoPLUS,treatment,instar,pixel_to_mm)]


# Stage micrometer
# We have two replicates for a random set of 150 'full stage micrometer' from D10
micrometer_1 <- fread("~/Google Drive/PNAS_Becker_data_scripts/validation/validation_fullMicro_150random_1.txt")
micrometer_1 <- micrometer_1[, -2]
names(micrometer_1) <- c('filepath', 'micro_length_manual_1_px', 'initials')
micrometer_1[, filebase_cor := substring(micrometer_1$filepath, 72)]

micrometer_2 <- fread("~/Google Drive/PNAS_Becker_data_scripts/validation/validation_fullMicro_150random_2.txt")
micrometer_2 <- micrometer_2[, -2]
names(micrometer_2) <- c('filepath', 'micro_length_manual_2_px', 'initials')
micrometer_2[, filebase_cor := substring(micrometer_1$filepath, 72)]


# all D10 pics - DB, RP, LW (exclude duplicated entries; use 2nd entry for each pic)
validation_DB <- fread("~/Google Drive/PNAS_Becker_data_scripts/validation/validation_20190311_DB.txt")
names(validation_DB) <- c('filepath', 'pedestal_score_manual_DB', 'animal_length_manual_DB', 'initials')
validation_DB[, filebase_cor := substring(validation_DB$filepath, 68)]
validation_DB[, rank := c(1:dim(validation_DB)[1])]
validation_DB_adj <- validation_DB[order(validation_DB$filebase, -abs(validation_DB$rank)), ] 
validation_DB_adj_nodups <- validation_DB_adj[!duplicated(validation_DB_adj$filebase), ]   
validation_DB_adj_nodups_no1 <- validation_DB_adj_nodups[!pedestal_score_manual_DB == -1]

validation_RP <- fread("~/Google Drive/PNAS_Becker_data_scripts/validation//validation_20190311_RP.txt")
names(validation_RP) <- c('filepath', 'pedestal_score_manual_RP', 'animal_length_manual_RP', 'initials')
validation_RP[,filebase_cor := substring(validation_RP$filepath, 67)]
validation_RP[, rank := c(1:dim(validation_RP)[1])]
validation_RP_adj <- validation_RP[order(validation_RP$filebase, -abs(validation_RP$rank)), ] 
validation_RP_adj_nodups <- validation_RP_adj[!duplicated(validation_RP_adj$filebase), ]  
validation_RP_adj_nodups_no1 <- validation_RP_adj_nodups[!pedestal_score_manual_RP == -1]

validation_LW <- fread("~/Google Drive/PNAS_Becker_data_scripts/validation/validation_20190311_LW.txt")
names(validation_LW) <- c('filepath', 'pedestal_score_manual_LW', 'animal_length_manual_LW', 'initials')
validation_LW[,filebase_cor := substring(validation_LW$filepath, 67)]
validation_LW[, rank := c(1:dim(validation_LW)[1])]
validation_LW_adj <- validation_LW[order(validation_LW$filebase, -abs(validation_LW$rank)), ] 
validation_LW_adj_nodups <- validation_LW_adj[!duplicated(validation_LW_adj$filebase), ]   
validation_LW_adj_nodups_no1 <- validation_LW_adj_nodups[!pedestal_score_manual_LW == -1]


# random 150 D10 pics - DB, RP, LW
validation_random150_DB <- fread("~/Google Drive/PNAS_Becker_data_scripts/validation/validation_20190311_DB_150random.txt")
names(validation_random150_DB) <- c('filepath', 'pedestal_score_manual_random_DB', 'animal_length_manual_random_DB', 'initials')
validation_random150_DB[, filebase_cor := substring(validation_random150_DB$filepath, 68)]
validation_random150_DB_no1 <- validation_random150_DB[!pedestal_score_manual_random_DB == -1]

validation_random150_RP <- fread("~/Google Drive/PNAS_Becker_data_scripts/validation/validation_20190311_RP_150random.txt")
names(validation_random150_RP) <- c('filepath', 'pedestal_score_manual_random_RP', 'animal_length_manual_random_RP', 'initials')
validation_random150_RP[, filebase_cor := substring(validation_random150_RP$filepath, 67)]
validation_random150_RP_no1 <- validation_random150_RP[!pedestal_score_manual_random_RP == -1]

validation_random150_LW <- fread("~/Google Drive/PNAS_Becker_data_scripts/validation/validation_20190311_LW_150random.txt")
names(validation_random150_LW) <- c('filepath', 'pedestal_score_manual_random_LW', 'animal_length_manual_random_LW', 'initials')
validation_random150_LW[, filebase_cor := substring(validation_random150_LW$filepath, 67)]
validation_random150_LW_no1 <- validation_random150_LW[!pedestal_score_manual_random_LW == -1]


# pretty clunky, but works 
GUI_data_plusManual <- merge(shape_use_small,
                             micrometer_1[, c("filebase_cor", "micro_length_manual_1_px")], 
                             by = 'filebase_cor', all.x = T)

GUI_data_plusManual <- merge(GUI_data_plusManual,
                             micrometer_2[, c("filebase_cor", "micro_length_manual_2_px")], 
                             by = 'filebase_cor', all.x = T)

GUI_data_plusManual <- merge(GUI_data_plusManual,
                             validation_DB_adj_nodups_no1[, c("filebase_cor", "pedestal_score_manual_DB", "animal_length_manual_DB")], 
                             by = 'filebase_cor', all.x = T)

GUI_data_plusManual <- merge(GUI_data_plusManual,
                             validation_RP_adj_nodups_no1[, c("filebase_cor", "pedestal_score_manual_RP", "animal_length_manual_RP")], 
                             by = 'filebase_cor', all.x = T)

GUI_data_plusManual <- merge(GUI_data_plusManual,
                             validation_LW_adj_nodups_no1[, c("filebase_cor", "pedestal_score_manual_LW", "animal_length_manual_LW")], 
                             by = 'filebase_cor', all.x = T)

GUI_data_plusManual <- merge(GUI_data_plusManual,
                             validation_random150_DB_no1[, c("filebase_cor", "pedestal_score_manual_random_DB", "animal_length_manual_random_DB")], 
                             by = 'filebase_cor', all.x = T)

GUI_data_plusManual <- merge(GUI_data_plusManual,
                             validation_random150_RP_no1[, c("filebase_cor", "pedestal_score_manual_random_RP", "animal_length_manual_random_RP")], 
                             by = 'filebase_cor', all.x = T)

GUI_data_plusManual <- merge(GUI_data_plusManual,
                             validation_random150_LW_no1[, c("filebase_cor", "pedestal_score_manual_random_LW", "animal_length_manual_random_LW")], 
                             by = 'filebase_cor', all.x = T)


GUI_data_plusManual[, micro_length_manual_avg_px := (micro_length_manual_1_px+micro_length_manual_2_px)/2]
GUI_data_plusManual[, animal_length_manual_mm_DB := animal_length_manual_DB / pixel_to_mm]
GUI_data_plusManual[, animal_length_manual_mm_RP := animal_length_manual_RP / pixel_to_mm]
GUI_data_plusManual[, animal_length_manual_mm_LW := animal_length_manual_LW / pixel_to_mm]
GUI_data_plusManual[, animal_length_manual_mm_random_DB := animal_length_manual_random_DB / pixel_to_mm]
GUI_data_plusManual[, animal_length_manual_mm_random_RP := animal_length_manual_random_RP / pixel_to_mm]
GUI_data_plusManual[, animal_length_manual_mm_random_LW := animal_length_manual_random_LW / pixel_to_mm]
GUI_data_plusManual[, animal_length_manual_mm_avg := (animal_length_manual_mm_DB + animal_length_manual_mm_RP + animal_length_manual_mm_LW)/3]

GUI_data_plusManual[, pedestal_score_manual_avg := (pedestal_score_manual_DB + pedestal_score_manual_RP + pedestal_score_manual_LW)/3]

GUI_data_plusManual_use <- GUI_data_plusManual


#############
### Fig 1 ###
#############


#Dapcha <- load.image("~/Google Drive/PNAS_Becker_data_scripts/DAPCHA.png")
#plot(Dapcha)


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

























###
# supplement



length_plot <- ggplot(data=(GUI_data_plusManual)) +
                  #geom_point(data = GUI_data_plusManual, aes(x=animal_length_manual_mm_random_DB, y=animal_length_manual_mm_RP), pch=22, color='red', bg="red", size=2) +
                  #geom_point(data = GUI_data_plusManual, aes(x=animal_length_manual_mm_RP, y=animal_length_manual_mm_DB), pch=22, color='red', bg="navyblue", size=2) +
                  #geom_point(data = GUI_data_plusManual, aes(x=animal_length_manual_mm_LW, y=animal_length_manual_mm_DB), pch=22, color='red', bg="grey50", size=2) +
                  
                  #geom_point(data = GUI_data_plusManual, aes(x=animal_length_manual_mm_random_RP, y=animal_length_manual_mm_RP), pch=21, color='navyblue', bg="navyblue", size=2) +
                  #geom_point(data = GUI_data_plusManual, aes(x=animal_length_manual_mm_LW, y=animal_length_manual_mm_RP), pch=21, color='navyblue', bg="grey50", size=2) +
                  
                  geom_point(data = GUI_data_plusManual, aes(x=animal_length_manual_mm_random_LW, y=animal_length_manual_mm_LW), pch=24, color='grey50', bg="grey50", size=2) +
                  
                  theme(legend.position="none", 
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_rect(fill = "transparent",colour = NA),
                        plot.background = element_rect(fill = "transparent",colour = NA), 
                        axis.text = element_text(size=20, family='Arial'), 
                        axis.title.x = element_text(size=20,family='Arial'), 
                        axis.title.y = element_text(size=20, family='Arial')) +
                  xlim(0.6,1) +
                  ylim(0.6,1) 

length_plot + geom_abline(intercept = 0, slope = 1, size=0.8) + labs(x='Person 1 (animal_length)', y='Person 2 (animal_length)')



pedestal_plot_DB_0 <- ggplot(data=GUI_data_plusManual[!pedestal_score_manual_DB == 'NA'][treatment==0], aes(y=maxHeight_mm_GUI, x=as.character(pedestal_score_manual_DB))) +
                      geom_beeswarm(size = 1, cex=1.5, priority='density', color='red') +  # dodge.width=0.5
                      geom_boxplot(color='black', fill='white', outlier.colour='black', outlier.size = 1, width=0.2) + 
                      theme(legend.position="none", 
                            panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(),
                            panel.background = element_rect(fill = "transparent",colour = NA),
                            plot.background = element_rect(fill = "transparent",colour = NA), 
                            axis.title.x = element_blank(), 
                            axis.title.y = element_blank(),
                            axis.text.x = element_blank(),
                            # axis.title.x = element_text(size=20,family='Arial'), 
                            # axis.title.y = element_text(size=20, family='Arial'), 
                            axis.text = element_text(size=20, family='Arial')) 

pedestal_plot_DB_0.5 <- ggplot(data=GUI_data_plusManual[!pedestal_score_manual_DB == 'NA'][treatment==0.5], aes(y=maxHeight_mm_GUI, x=as.character(pedestal_score_manual_DB))) +
                        geom_beeswarm(size = 1, cex=1.5, priority='density', color='red') +  # dodge.width=0.5
                        geom_boxplot(color='black', fill='white', outlier.colour='black', outlier.size = 1, width=0.2) + 
                        theme(legend.position="none", 
                              panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank(),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA), 
                              axis.title.x = element_blank(), 
                              axis.title.y = element_blank(),
                              axis.text.x = element_blank(),
                              axis.text.y = element_blank(),
                              # axis.title.x = element_text(size=20,family='Arial'), 
                              # axis.title.y = element_text(size=20, family='Arial'), 
                              axis.text = element_text(size=20, family='Arial')) 

pedestal_plot_RP_0 <- ggplot(data=GUI_data_plusManual[!pedestal_score_manual_RP == 'NA'][treatment==0], aes(y=maxHeight_mm_GUI, x=as.character(pedestal_score_manual_RP))) +
                      geom_beeswarm(size = 1, cex=1.5, priority='density', color='navyblue') +  # dodge.width=0.5
                      geom_boxplot(color='black', fill='white', outlier.colour='black', outlier.size = 1, width=0.2) + 
                      theme(legend.position="none", 
                            panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(),
                            panel.background = element_rect(fill = "transparent",colour = NA),
                            plot.background = element_rect(fill = "transparent",colour = NA), 
                            axis.title.x = element_blank(), 
                            axis.title.y = element_blank(),
                            axis.text.x = element_blank(),
                            # axis.title.x = element_text(size=20,family='Arial'), 
                            # axis.title.y = element_text(size=20, family='Arial'), 
                            axis.text = element_text(size=20, family='Arial')) 

pedestal_plot_RP_0.5 <- ggplot(data=GUI_data_plusManual[!pedestal_score_manual_RP == 'NA'][treatment==0.5], aes(y=maxHeight_mm_GUI, x=as.character(pedestal_score_manual_RP))) +
                        geom_beeswarm(size = 1, cex=1.5, priority='density', color='navyblue') +  # dodge.width=0.5
                        geom_boxplot(color='black', fill='white', outlier.colour='black', outlier.size = 1, width=0.2) + 
                        theme(legend.position="none", 
                              panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank(),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA), 
                              axis.title.x = element_blank(), 
                              axis.title.y = element_blank(),
                              axis.text.x = element_blank(),
                              axis.text.y = element_blank(),
                              # axis.title.x = element_text(size=20,family='Arial'), 
                              # axis.title.y = element_text(size=20, family='Arial'), 
                              axis.text = element_text(size=20, family='Arial'))

pedestal_plot_LW_0 <- ggplot(data=GUI_data_plusManual[!pedestal_score_manual_LW == 'NA'][treatment==0], aes(y=maxHeight_mm_GUI, x=as.character(pedestal_score_manual_LW))) +
                      geom_beeswarm(size = 1, cex=1.5, priority='density', color='grey50') +  # dodge.width=0.5
                      geom_boxplot(color='black', fill='white', outlier.colour='black', outlier.size = 1, width=0.2) + 
                      theme(legend.position="none", 
                            panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(),
                            panel.background = element_rect(fill = "transparent",colour = NA),
                            plot.background = element_rect(fill = "transparent",colour = NA), 
                            axis.title.x = element_blank(), 
                            axis.title.y = element_blank(),
                            # axis.title.x = element_text(size=20,family='Arial'), 
                            # axis.title.y = element_text(size=20, family='Arial'), 
                            axis.text = element_text(size=20, family='Arial')) 

pedestal_plot_LW_0.5 <- ggplot(data=GUI_data_plusManual[!pedestal_score_manual_LW == 'NA'][treatment==0.5], aes(y=maxHeight_mm_GUI, x=as.character(pedestal_score_manual_LW))) +
                        geom_beeswarm(size = 1, cex=1.5, priority='density', color='grey50') +  # dodge.width=0.5
                        geom_boxplot(color='black', fill='white', outlier.colour='black', outlier.size = 1, width=0.2) + 
                        theme(legend.position="none", 
                              panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank(),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA), 
                              axis.title.x = element_blank(), 
                              axis.title.y = element_blank(),
                              axis.text.y = element_blank(),
                              # axis.title.x = element_text(size=20,family='Arial'), 
                              # axis.title.y = element_text(size=20, family='Arial'), 
                              axis.text = element_text(size=20, family='Arial')) 


grid.arrange(pedestal_plot_DB_0 + coord_flip() + scale_x_discrete(limits=c("0","30","50")) + scale_y_continuous(limits=c(0.15,0.37), breaks=c(0.15,0.25,0.35)), 
             pedestal_plot_DB_0.5 + coord_flip() + scale_x_discrete(limits=c("0","30","50")) + scale_y_continuous(limits=c(0.15,0.37), breaks=c(0.15,0.25,0.35)), 
             pedestal_plot_RP_0 + coord_flip() + scale_x_discrete(limits=c("0","30","50")) + scale_y_continuous(limits=c(0.15,0.37), breaks=c(0.15,0.25,0.35)), 
             pedestal_plot_RP_0.5 + coord_flip() + scale_x_discrete(limits=c("0","30","50")) + scale_y_continuous(limits=c(0.15,0.37), breaks=c(0.15,0.25,0.35)), 
             pedestal_plot_LW_0 + coord_flip() + scale_x_discrete(limits=c("0","30","50")) + scale_y_continuous(limits=c(0.15,0.37), breaks=c(0.15,0.25,0.35)), 
             pedestal_plot_LW_0.5 + coord_flip() + scale_x_discrete(limits=c("0","30","50")) + scale_y_continuous(limits=c(0.15,0.37), breaks=c(0.15,0.25,0.35)), 
             ncol=2)
                    






pedestal_overlap <- GUI_data_plusManual[!pedestal_score_manual_DB == 'NA'][!pedestal_score_manual_RP == 'NA'][!pedestal_score_manual_LW == 'NA'][, c('pedestal_score_manual_DB',
                                                                                                                                                     'pedestal_score_manual_RP',
                                                                                                                                                     'pedestal_score_manual_LW',
                                                                                                                                                     'pedestal_score_manual_random_DB',
                                                                                                                                                     'pedestal_score_manual_random_RP',
                                                                                                                                                     'pedestal_score_manual_random_LW', 
                                                                                                                                                     'treatment')]

pedestal_overlap[, DB_DB := ifelse(pedestal_score_manual_DB == pedestal_score_manual_random_DB, '0', 
                              ifelse(pedestal_score_manual_DB > pedestal_score_manual_random_DB, '1',
                              ifelse(pedestal_score_manual_DB < pedestal_score_manual_random_DB, '-1', 'NA')))]

pedestal_overlap[, DB_RP := ifelse(pedestal_score_manual_DB == pedestal_score_manual_RP, '0', 
                            ifelse(pedestal_score_manual_DB > pedestal_score_manual_RP, '1',
                            ifelse(pedestal_score_manual_DB < pedestal_score_manual_RP, '-1', 'NA')))]

pedestal_overlap[, DB_LW := ifelse(pedestal_score_manual_DB == pedestal_score_manual_LW, '0', 
                            ifelse(pedestal_score_manual_DB > pedestal_score_manual_LW, '1',
                            ifelse(pedestal_score_manual_DB < pedestal_score_manual_LW, '-1', 'NA')))]

pedestal_overlap[, RP_RP := ifelse(pedestal_score_manual_RP == pedestal_score_manual_random_RP, '0', 
                            ifelse(pedestal_score_manual_RP > pedestal_score_manual_random_RP, '1',
                            ifelse(pedestal_score_manual_RP < pedestal_score_manual_random_RP, '-1', 'NA')))]

pedestal_overlap[, RP_LW := ifelse(pedestal_score_manual_RP == pedestal_score_manual_LW, '0', 
                            ifelse(pedestal_score_manual_RP > pedestal_score_manual_LW, '1',
                            ifelse(pedestal_score_manual_RP < pedestal_score_manual_LW, '-1', 'NA')))]

pedestal_overlap[, LW_LW := ifelse(pedestal_score_manual_LW == pedestal_score_manual_random_LW, '0', 
                            ifelse(pedestal_score_manual_LW > pedestal_score_manual_random_LW, '1',
                            ifelse(pedestal_score_manual_LW < pedestal_score_manual_random_LW, '-1', 'NA')))]


table(pedestal_overlap[treatment==0]$DB_DB)
table(pedestal_overlap[treatment==0.1]$DB_DB)
table(pedestal_overlap[treatment==0.25]$DB_DB)
table(pedestal_overlap[treatment==0.5]$DB_DB)
table(pedestal_overlap[treatment==1]$DB_DB)

table(pedestal_overlap[treatment==0]$DB_RP)
table(pedestal_overlap[treatment==0.1]$DB_RP)
table(pedestal_overlap[treatment==0.25]$DB_RP)
table(pedestal_overlap[treatment==0.5]$DB_RP)
table(pedestal_overlap[treatment==1]$DB_RP)

table(pedestal_overlap[treatment==0]$DB_LW)
table(pedestal_overlap[treatment==0.1]$DB_LW)
table(pedestal_overlap[treatment==0.25]$DB_LW)
table(pedestal_overlap[treatment==0.5]$DB_LW)
table(pedestal_overlap[treatment==1]$DB_LW)

table(pedestal_overlap[treatment==0]$RP_RP)
table(pedestal_overlap[treatment==0.1]$RP_RP)
table(pedestal_overlap[treatment==0.25]$RP_RP)
table(pedestal_overlap[treatment==0.5]$RP_RP)
table(pedestal_overlap[treatment==1]$RP_RP)

table(pedestal_overlap[treatment==0]$RP_LW)
table(pedestal_overlap[treatment==0.1]$RP_LW)
table(pedestal_overlap[treatment==0.25]$RP_LW)
table(pedestal_overlap[treatment==0.5]$RP_LW)
table(pedestal_overlap[treatment==1]$RP_LW)

table(pedestal_overlap[treatment==0]$LW_LW)
table(pedestal_overlap[treatment==0.1]$LW_LW)
table(pedestal_overlap[treatment==0.25]$LW_LW)
table(pedestal_overlap[treatment==0.5]$LW_LW)
table(pedestal_overlap[treatment==1]$LW_LW)





















# regression function
ggplotRegression <- function(dat, xvar, yvar){
  fml <- paste(yvar, "~", xvar)
  fit <- lm(fml, dat)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point(size = 1) +
    stat_smooth(method = "lm", col = "red", size=0.5) +
    labs(title = paste("adj R2 = ", signif(summary(fit)$adj.r.squared, 3)))
  #                       "Intercept =", signif(fit$coef[[1]],3),
  #                        " Slope =", signif(fit$coef[[2]], 3),
  #                        " P =", signif(summary(fit)$coef[2,4], 3)))
}


### STAGE MICROMETER
# stage micrometer LW vs GUI (use all data, including instar == NA)
micrometer_1 <- ggplotRegression(GUI_data_plusManual[micro_length_manual_1_px>250], "pixel_to_mm", "micro_length_manual_1_px")
micrometer_2 <- ggplotRegression(GUI_data_plusManual, "pixel_to_mm", "micro_length_manual_2_px")
micrometer_avg <- ggplotRegression(GUI_data_plusManual, "pixel_to_mm", "micro_length_manual_avg_px")

micrometer_1_plot <- micrometer_1 + labs(x='GUI', y='manual')
micrometer_2_plot <- micrometer_2 + labs(x='GUI', y='manual')
micrometer_avg_plot <- micrometer_avg + labs(x='GUI', y='manual')

grid.arrange(micrometer_1_plot, micrometer_2_plot, micrometer_avg_plot)

# no real issue w/ outlier in 1st rep: 
#/mnt/spicy_1/daphnia/data 
#display fullMicro_100270_D10_62_ctrl_1A_RigA_20170602T145148.bmp
#>> exclude >> [!filebase_cor == 'D10_62_ctrl_1A_RigA_20170602T145148.bmp']

#### CHECK IN  W/ ALAN  >> leave out 2nd rep from LW


### ANIMAL LENGTH
# animal length - full vs random MANUAL: DBvsDBrandom, RPvsRPrandom, LWvsLWrandom 
animal_length_DBvsDBrandom <- ggplotRegression(GUI_data_plusManual, "animal_length_manual_mm_DB", "animal_length_manual_mm_random_DB")
animal_length_RPvsRPrandom <- ggplotRegression(GUI_data_plusManual, "animal_length_manual_mm_RP", "animal_length_manual_mm_random_RP")
animal_length_LWvsLWrandom <- ggplotRegression(GUI_data_plusManual, "animal_length_manual_mm_LW", "animal_length_manual_mm_random_LW")

animal_length_DBvsDBrandom_plot <- animal_length_DBvsDBrandom + labs(x='length_full', y='length_random', subtitle = "DB_manual")
animal_length_RPvsRPrandom_plot <- animal_length_RPvsRPrandom + labs(x='length_full', y='length_random', subtitle = "RP_manual")
animal_length_LWvsLWrandom_plot <- animal_length_LWvsLWrandom + labs(x='length_full', y='length_random', subtitle = "LW_manual")

grid.arrange(animal_length_DBvsDBrandom_plot, animal_length_RPvsRPrandom_plot, animal_length_LWvsLWrandom_plot)
# There is hardly any variation between full and random sets

# animal length - full vs full MANUAL: DBvsRP, DBvsLW, RPvsLW 
animal_length_DBvsRP <- ggplotRegression(GUI_data_plusManual, "animal_length_manual_mm_DB", "animal_length_manual_mm_RP")
animal_length_DBvsLW <- ggplotRegression(GUI_data_plusManual, "animal_length_manual_mm_DB", "animal_length_manual_mm_LW")
animal_length_RPvsLW <- ggplotRegression(GUI_data_plusManual, "animal_length_manual_mm_RP", "animal_length_manual_mm_LW")

animal_length_DBvsRP_plot <- animal_length_DBvsDBrandom + labs(x='length_DB', y='length_RP')
animal_length_DBvsLW_plot <- animal_length_RPvsRPrandom + labs(x='length_DB', y='length_LW')
animal_length_RPvsLW_plot <- animal_length_LWvsLWrandom + labs(x='length_RP', y='length_LW')

grid.arrange(animal_length_DBvsRP_plot, animal_length_DBvsLW_plot, animal_length_RPvsLW_plot)
# There is hardly any variation between individual modifiers


# animal length - MANUAL vs GUI: average(DB,RP,LW) vs GUI 
animal_length_MANvsGUI <- ggplot(GUI_data_plusManual[instar %in% c(1,2,3)], aes(x = animal_length_manual_mm_avg, y = animal_length_mm, color=as.factor(instar))) +
  geom_point(size = 2) +
  stat_smooth(method = "lm", col = "black", size=1) +
  theme(legend.position="top", axis.text = element_text(size=25))

animal_length_MANvsGUI_plot <- animal_length_MANvsGUI + labs(x='length_manual', y='animal length (mm) GUI')
# Good regression fit when contrasting MANUAL and GUI data sets



### PEDESTAL MAX HEIGHT
GUI_data_plusManual_pedestal <- GUI_data_plusManual[instar %in% c(1,2,3)][, c("filebase", "treatment", "instar", "pedestal_max_height_mm", "pedestal_max_height_mm_adj", 
                                                                              "pedestal_score_manual_DB" , "pedestal_score_manual_RP", "pedestal_score_manual_LW", 
                                                                              "pedestal_score_manual_random_DB", "pedestal_score_manual_random_RP", "pedestal_score_manual_random_LW", 
                                                                              "pedestal_score_manual_avg")][!is.na(pedestal_score_manual_avg)]
GUI_data_plusManual_pedestal_long <- melt(GUI_data_plusManual_pedestal, id.vars=c("filebase", "treatment", "instar"))
colnames(GUI_data_plusManual_pedestal_long) <- c("filebase","treatment","instar","set","pedestal_score")


# pedestal - full vs random MANUAL: DBvsDBrandom, RPvsRPrandom, LWvsLWrandom 
pedestal_score_DBvsDBrandom_plot <- ggplot(data = GUI_data_plusManual_pedestal_long[set %in% c('pedestal_score_manual_DB', 'pedestal_score_manual_random_DB', 
                                                                                               'pedestal_score_manual_RP', 'pedestal_score_manual_random_RP', 
                                                                                               'pedestal_score_manual_LW', 'pedestal_score_manual_random_LW')], 
                                           aes(x=as.factor(treatment), y=pedestal_score, fill = as.factor(treatment))) +
  #geom_boxplot(outlier.size = 0.5, notch=FALSE) +
  geom_jitter(width=0.1,alpha=0.5) +
  #geom_beeswarm(size = 0.5, cex=1.5, priority='density', dodge.width=0.5) +  
  #geom_boxplot(fill="white", outlier.colour='black', outlier.size = 0.1) +  # width=0.3 
  #facet_grid(~treatment) +
  theme(legend.position="top") +
  theme(axis.text.x = element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank())

pedestal_score_RPvsRPrandom_plot <- ggplot(data = GUI_data_plusManual_pedestal_long[set %in% c('pedestal_score_manual_RP', 'pedestal_score_manual_random_RP')], 
                                           aes(x=set, y=pedestal_score, fill = as.factor(set))) +
  geom_boxplot(outlier.size = 0.5, notch=FALSE) +
  geom_jitter(width=0.1,alpha=0.2) +
  facet_grid(~treatment) +
  theme(legend.position="top") +
  theme(axis.text.x = element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank())

pedestal_score_LWvsLWrandom_plot <- ggplot(data = GUI_data_plusManual_pedestal_long[set %in% c('pedestal_score_manual_LW', 'pedestal_score_manual_random_LW')], 
                                           aes(x=set, y=pedestal_score, fill = as.factor(set))) +
  geom_boxplot(outlier.size = 0.5, notch=FALSE) +
  geom_jitter(width=0.1,alpha=0.2) +
  facet_grid(~treatment) +
  theme(legend.position="top") +
  theme(axis.text.x = element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank())

grid.arrange(pedestal_score_DBvsDBrandom_plot, pedestal_score_RPvsRPrandom_plot, pedestal_score_LWvsLWrandom_plot)



# pedstal score - full vs full MANUAL: DBvsRPvsLW
pedestal_score_DBvsRP_plot <- ggplot(data = GUI_data_plusManual_pedestal_long[set %in% c('pedestal_score_manual_DB', 'pedestal_score_manual_RP', 'pedestal_score_manual_LW')], 
                                     aes(x=set, y=pedestal_score, fill = as.factor(set))) +
  geom_boxplot(outlier.size = 0.5, notch=FALSE) +
  geom_jitter(width=0.1,alpha=0.2) +
  facet_grid(~treatment) +
  theme(legend.position="top") +
  theme(axis.text.x = element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank())



# pedestal score - MANUAL vs GUI: average(DB,RP,LW) vs GUI 
pedestal_score_MANvsGUI_plot <- ggplot(data = GUI_data_plusManual_pedestal, 
                                       aes(x=pedestal_score_manual_avg, y=pedestal_max_height_mm_adj, col = as.factor(treatment))) +
  
  geom_boxplot() + 
  theme(legend.position="top") 






length_pond <- ggplot(dap.sum_cor[pond %in% c('DBunk','D8','D10')], 
                      aes(y=animal_length_mm, x=as.factor(treatment), color=as.factor(instar))) + 
  geom_beeswarm(size = 0.5, cex=1.5, priority='density', dodge.width=0.5) +  
  geom_boxplot(fill="black", outlier.colour='black', outlier.size = 0.5) +  # width=0.3 
  facet_grid(~pond) + 
  theme(legend.position="top") 
length_pond + labs(x = "treatment", y = "animal length (mm)")

