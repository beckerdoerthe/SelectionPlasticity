# Becker et al - data input


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
library(sjstats)
library(rstatix)
library(patchwork)


#######################
### load pheno data ###
#######################

load(file = "data/shape_use.ag_filtered_final_02June2020.Rdata")
# final_data

setkey(final_data, cloneid_geno, i)

# reduce control data (choose random 4 Geno IDs for each treatment group in all three ctrl clones each)
# set.seed(123)
# random_ctrl_O <- final_data[group == 'ctrl_O'][, c('cloneid_geno','Geno','treatment')] %>% group_by(cloneid_geno, treatment) %>% sample_n(4)
# ctrlO_clones <- final_data[Geno %in% random_ctrl_O$Geno][, -"max_height"]

## use control clones used for GRM prep (from previous set.seed)
new_controls <- c("April_2017_DBunk_131_101952_0_1_1A", "April_2017_DBunk_131_101989_0_1_1D", "April_2017_DBunk_131_110092_0_1_1A", "April_2017_DBunk_131_112503_0_1_2A", 
                  "April_2017_DBunk_132_101990_0_1_1A", "April_2017_DBunk_132_110098_0_1_2A", "April_2017_DBunk_132_110275_0_1_2B", "April_2017_DBunk_132_111147_0_1_1C", 
                  "April_2017_DBunk_90_101634_0_1_2B", "April_2017_DBunk_90_110581_0_1_2A", "April_2017_DBunk_90_110803_0_1_2B", "April_2017_DBunk_131_101248_0.5_1_1A", 
                  "April_2017_DBunk_131_110405_0.5_1_1A", "April_2017_DBunk_131_111570_0.5_1_1B" , "April_2017_DBunk_132_101270_0.5_1_2A", "April_2017_DBunk_132_102083_0.5_1_1B", 
                  "April_2017_DBunk_132_111579_0.5_1_2A", "April_2017_DBunk_132_112514_0.5_1_1B", "April_2017_DBunk_90_101638_0.5_1_1B", "April_2017_DBunk_90_101641_0.5_1_3A", 
                  "April_2017_DBunk_90_101972_0.5_1_1A", "April_2017_DBunk_90_102089_0.5_1_1D", "April_2017_DBunk_131_101952_0_2_1A", "April_2017_DBunk_131_101989_0_2_1D", 
                  "April_2017_DBunk_131_110092_0_2_1A", "April_2017_DBunk_132_110098_0_2_2A", "April_2017_DBunk_132_110275_0_2_2B", "April_2017_DBunk_132_111147_0_2_1C", 
                  "April_2017_DBunk_90_101634_0_2_2B", "April_2017_DBunk_90_110581_0_2_2A", "April_2017_DBunk_90_110803_0_2_2B", "April_2017_DBunk_131_101248_0.5_2_1A", 
                  "April_2017_DBunk_131_101958_0.5_2_2A", "April_2017_DBunk_131_110405_0.5_2_1A", "April_2017_DBunk_131_111570_0.5_2_1B", "April_2017_DBunk_132_101270_0.5_2_2A", 
                  "April_2017_DBunk_132_102083_0.5_2_1B", "April_2017_DBunk_132_111579_0.5_2_2A", "April_2017_DBunk_90_101638_0.5_2_1B", "April_2017_DBunk_90_101641_0.5_2_3A", 
                  "April_2017_DBunk_90_101972_0.5_2_1A", "April_2017_DBunk_90_102089_0.5_2_1D")

ctrlO_clones <- final_data[GenoPLUS %in% new_controls][, -"max_height"]


set.seed(123)
random_ctrl_A <- final_data[group == 'ctrl_A'][, c('cloneid_geno','Geno','treatment')] %>% group_by(cloneid_geno, treatment) %>% sample_n(4)
ctrlA_clones <- final_data[Geno %in% random_ctrl_A$Geno][, -"max_height"]

O_clones <- final_data[group == "O"][, -"max_height"]
A_clones <- final_data[group == "A"][, -"max_height"]

all_data_Os <- rbind(O_clones,ctrlO_clones)
all_data_As <- rbind(A_clones,ctrlA_clones)

all_data <- rbind(O_clones,ctrlO_clones,A_clones,ctrlA_clones)


# find O clones with highest read depth per cluster
all_data_Os[,revorder:=frankv(medrd,order=-1,ties.method = "first")]  
setkey(all_data_Os, revorder, Geno, i) 
tmp_medrd_Os <- unique(all_data_Os[i==150][,c('cloneid_geno', 'SC_unique', 'medrd')] %>% group_by(SC_unique) %>% slice(1))
medrd_Os <- as.data.table(tmp_medrd_Os[,'cloneid_geno'])

all_data_O.medrd <- all_data[cloneid_geno %in% medrd_Os$cloneid_geno]

# exclude O clones with less than 1 sample per treatment&instar group 
table(all_data_O.medrd$SC_unique, all_data_O.medrd$treatment)

# exclude missing data in ctrl and treatment (CLUNKY CODE, but works!)
tbl_all_data_O <- as.data.table(table(all_data_O.medrd$SC_unique, all_data_O.medrd$treatment, all_data_O.medrd$instar))
setnames(tbl_all_data_O, c('V1','V2','V3','N'), c('SC','treatment','instar','count'))
tbl_all_data_O_wide <- as.data.table(dcast(tbl_all_data_O, SC ~ c(paste0("treatment_", treatment, "_instar_", instar)), fun=mean, value.var='count'))

data_O_wide_use <- all_data_O.medrd[SC_unique %in% tbl_all_data_O_wide[treatment_0.5_instar_1 >= 641 & treatment_0.5_instar_2 >= 641 & treatment_0_instar_1 >= 641 & treatment_0_instar_2 >= 641]$SC]  
Os <- unique(levels(as.factor(data_O_wide_use$cloneid_geno)))  ### 51 unique clones


# exclude A clones with less than 1 sample per treatment&instar group
table(all_data_As$cloneid_geno, all_data_As$treatment, all_data_As$instar)

# exclude missing data in ctrl and treatment (CLUNKY CODE, but works!)
tbl_all_data_A <- as.data.table(table(all_data_As$cloneid_geno, all_data_As$treatment, all_data_As$instar))
setnames(tbl_all_data_A, c('V1','V2','V3','N'), c('cloneid_geno','treatment','instar','count'))
tbl_all_data_A_wide <- as.data.table(dcast(tbl_all_data_A, cloneid_geno ~ c(paste0("treatment_", treatment, "_instar_", instar)), fun=mean, value.var='count'))

data_A_wide_use <- all_data_As[cloneid_geno %in% tbl_all_data_A_wide[treatment_0.5_instar_1 >= 641 & treatment_0.5_instar_2 >= 641 & treatment_0_instar_1 >= 641 & treatment_0_instar_2 >= 641]$cloneid_geno]

# find A clones with highest read depth - ignore here and use all As
data_A_wide_use[,revorder:=frankv(medrd,order=-1,ties.method = "first")]  
setkey(data_A_wide_use, revorder, Geno, i) 
tmp_medrd_As <- unique(data_A_wide_use[i==150][,c('cloneid_geno', 'SC_unique', 'medrd')]$cloneid_geno)
medrd_As <- tmp_medrd_As


# exclude pond "D10"
all_data_final <- all_data[cloneid_geno %in% Os | cloneid_geno %in% medrd_As][! pond == "D10"]
all_data_final[, i:= as.numeric(i)]
all_data_final[, instar_new := ifelse(all_data_final$instar == 1, 'instar 1', 'instar 2')]
all_data_final[, SC_group_new := ifelse(all_data_final$SC_group == "O", 'cluster O', 'cluster A')]
all_data_final[, treatment_new := ifelse(all_data_final$treatment == 0, 'C', 'P')]
all_data_final[, max_height_new := max(height, na.rm=T), by=c('Geno','instar')]

setkey(all_data_final, cloneid_geno, i) 

save(all_data_final, file = "data/all_data_final.RData")



load(file = "data/all_data_final.RData")

BatchInfo <- as.data.table(all_data_final[i == 150][instar == 1][, list(height = mean(height)), list(cloneid_geno, batch, SC_group, treatment)])
setkey(BatchInfo, batch, SC_group, treatment, cloneid_geno)
BatchInfo

write.csv(BatchInfo, file = "BatchInfo.csv", row.names=FALSE)



# Supplemental Table 2 - clone IDs w/ pond and season info
SupplTab <- all_data_final[i == 150][, list(height = mean(height)), list(cloneid_geno, pond, season, SC_group)]
SupplTab[, season_new := ifelse(season == "spring_1_2017", "Spring 2017", 
                           ifelse(season == "spring_2_2017", "Spring 2017", 
                            ifelse(season == "fall_2016", "Fall 2016", 
                             ifelse(season == "spring_2016", "Spring 2016", "NA"))))]

SupplTab_final <- SupplTab[, -c("height","season")]
setnames(SupplTab_final, c("cloneid_geno","pond","SC_group","season_new"), c("clone ID", "pond", "cluster", "season"))
write.csv(SupplTab_final, file="SupplementTable1.csv", row.names=FALSE)


# Supplemental Figure 3 - batch proportions re A vs O and Ctrl vs Predation
batch_AvsO_CtrlvsTrt <- as.data.table(read.delim(file= "batch_info_graph.txt", header = TRUE, sep = "\t"))


Batch_OvsA <- ggplot(batch_AvsO_CtrlvsTrt, aes(fill=as.factor(Cluster), y=NoClonesTested, x=Batch)) + 
                      geom_bar(position="fill", stat="identity") + 
                      theme(legend.position="bottom", 
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

Batch_OvsA + labs(x = "batch", y = paste("proportion")) + scale_fill_manual(values=c("#AFE1AF","#40826D"))


Batch_CtrlvsPred <- ggplot(batch_AvsO_CtrlvsTrt, aes(fill=Treatment, y=NoClonesTested, x=Batch)) + 
                            geom_bar(position="fill", stat="identity") + 
                            theme(legend.position="bottom", 
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

Batch_CtrlvsPred + labs(x = "batch", y = paste("proportion")) + scale_fill_manual(values=c("#000000","#FF0000"))



( Batch_OvsA + labs(x = "batch", y = paste("proportion")) + scale_fill_manual(values=c("#AFE1AF","#40826D")) ) +
( Batch_CtrlvsPred + labs(x = "batch", y = paste("proportion")) + scale_fill_manual(values=c("#000000","#FF0000")) ) + 
  plot_spacer() + 
  plot_spacer() + 
  plot_layout(ncol=4, widths = c(1,1,1,1))



