## Becker et al - shape analysis using GEOMORPH

rm(list=ls()) 

# n number of specimens/individuals
# p number of landmarks
# k number of dimensions

library(data.table)
library(geomorph)
library(reshape2)
library(dplyr)
library(tidyverse)


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

# exclude O clones with less than 2 samples per treatment&instar group 
table(all_data_O.medrd$SC_unique, all_data_O.medrd$treatment)

# exclude missing data in ctrl and treatment (CLUNKY CODE, but works!)
tbl_all_data_O <- as.data.table(table(all_data_O.medrd$SC_unique, all_data_O.medrd$treatment, all_data_O.medrd$instar))
setnames(tbl_all_data_O, c('V1','V2','V3','N'), c('SC','treatment','instar','count'))
tbl_all_data_O_wide <- as.data.table(dcast(tbl_all_data_O, SC ~ c(paste0("treatment_", treatment, "_instar_", instar)), fun=mean, value.var='count'))

data_O_wide_use <- all_data_O.medrd[SC_unique %in% tbl_all_data_O_wide[treatment_0.5_instar_1 > 641 & treatment_0.5_instar_2 > 641 & treatment_0_instar_1 > 641 & treatment_0_instar_2 > 641]$SC]  
Os <- unique(levels(as.factor(data_O_wide_use$cloneid_geno)))  ### 51 unique clones


# exclude A clones with less than 2 samples per treatment&instar group
table(all_data_As$cloneid_geno, all_data_As$treatment, all_data_As$instar)

# exclude missing data in ctrl and treatment (CLUNKY CODE, but works!)
tbl_all_data_A <- as.data.table(table(all_data_As$cloneid_geno, all_data_As$treatment, all_data_As$instar))
setnames(tbl_all_data_A, c('V1','V2','V3','N'), c('cloneid_geno','treatment','instar','count'))
tbl_all_data_A_wide <- as.data.table(dcast(tbl_all_data_A, cloneid_geno ~ c(paste0("treatment_", treatment, "_instar_", instar)), fun=mean, value.var='count'))

data_A_wide_use <- all_data_As[cloneid_geno %in% tbl_all_data_A_wide[treatment_0.5_instar_1 > 641 & treatment_0.5_instar_2 > 641 & treatment_0_instar_1 > 641 & treatment_0_instar_2 > 641]$cloneid_geno]

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

save(all_data_final, file = "all_data_final.RData")


### modify here re instar
all_data_use <- all_data_final[, c("cloneid_geno","Geno","SC_unique","group","treatment","instar","i","batch","height","APlength")][instar == 1]
setnames(all_data_use, c('APlength','height'), c('x', 'y'))  

all_data_use_wide <- reshape(all_data_use,   ## modify input data set 
                             idvar=c("Geno", "cloneid_geno", "SC_unique", "group", "instar", "treatment", "batch"),
                             timevar="i",
                             direction="wide")

setkey(all_data_use_wide, Geno)


# procrustes data to use
tmp_data <- all_data_use_wide
names <- all_data_use_wide[,2]
procrustes <- arrayspecs(tmp_data[,8:ncol(tmp_data)], 641, 2)  ## make 3D array object
dimnames(procrustes)[3] <- names   ## 3D array: number of rows (p), number of columns (k) and number of “sheets” (n), data has p=641, k=2, n= total number of rows.

# links
#links <- cbind(c(2:698), c(3:699))

## PROCRUSTES & TRAJECTORY ANALYSES  (modified from Andrew's code)
# Step 1: procrustes analysis ----
gpa_trans <- gpagen(procrustes)  #, Proj = TRUE, ProcD = TRUE, curves = NULL, surfaces = NULL, print.progress = TRUE) 
str(gpa_trans)
gpa_trans$Csize

# Step 2: create geomorph data frame for procD.lm and trajectory analysis ----
# note that PREDATION is naturally a numeric variable (juju concentration)
# and needs to be coerced to a factor.
gdf <- geomorph.data.frame(
  Shape = gpa_trans$coords,
  Clone = factor(all_data_use_wide$cloneid_geno),
  Predation = factor(all_data_use_wide$treatment),
  Csize = gpa_trans$Csize, 
  SC = factor(all_data_use_wide$SC_unique),
  Group = factor(all_data_use_wide$group),
  Batch = factor(all_data_use_wide$batch))


# Step 3: Visualise the procrustes data ----
pdf("gpa_trans_I1_16June.pdf")
par(mfrow = c(1,1))
plotAllSpecimens(gpa_trans$coords) #, links = links)
dev.off()

jpeg('gpa_trans_I1_16June.jpg')
par(mfrow = c(1,1))
plotAllSpecimens(gpa_trans$coords) #, links = links)
dev.off()


# Step 4: Formal test of whether shape ~ predation risk * genotype ----
# This is the basic first analysis
# Does the effect of treatment (cues) on shape depend on clone?

# procD.lm
mod <- procD.lm(f1 = Shape ~ Clone*Predation+Batch, SS.type="II", data = gdf, seed = 123456, iter = 999, RRPP = TRUE, print.progress = TRUE)
reveal.model.designs(mod)
mod

anova_mod <- anova(mod)  # effect.type = "cohenf"
#attributes(mod)

# manova
manova_mod <- manova.update(mod, print.progress=FALSE, tol=0)
#manova_mod_.001 <- manova.update(mod, print.progress=FALSE, tol=0.001)
#manova_mod_PC10 <- manova.update(mod, print.progress=FALSE, PC.no=10)

manova_mod_out <- summary(manova_mod, test="Pillai")
#summary(manova_mod_.001, test="Pillai")
#summary(manova_mod_PC10, test="Pillai")

save(anova_mod, file='anova.OandA_I1_16June.RData')
save(manova_mod_out, file='manova.OandA_I1_16June.RData')


# Step 5a: Basic Trajectory Analysis ----
# Distance, Direction and Shape of Trajectory

ta_out <- trajectory.analysis(mod, groups = gdf$Clone, traj.pts = gdf$Predation, pca = TRUE, print.progress = TRUE)
ta_out

## PC plot
#TP1 <- plot(ta_out, pch = 21, cex = 1.2, cex.axis = 2, bg = as.numeric(gdf$Predation), col = "gray")
##TP1 <- plot(ta_out, pch = 21, cex = 1.2, cex.axis = 2, bg = as.numeric(gdf$Clone), col = "gray")
#add.trajectories(TP1, traj.bg = 1:nlevels(gdf$Clone),   # superimposes the trajectories 
#                 start.bg = 1:nlevels(gdf$Clone),
#                 end.bg = 1:nlevels(gdf$Clone))

# extract PC scores
#PC_scores_I1 <- TP1$pc.points
#save(PC_scores_I1, file = 'PC_scores_I1.RData')


TA.summary_MD <- summary(ta_out, attribute = "MD")  # Magnitude difference (absolute difference between path distances)
MD_matrix <- TA.summary_MD$summary.table

TA.summary_TC <- summary(ta_out, attribute = "TC")  # angular differences between trajectory principal axes
TC_matrix <- TA.summary_TC$summary.table

save(ta_out, file='ta_out.OandA_I1_16June.RData')
save(MD_matrix, file='MD.OandA_I1_16June.RData')
save(TC_matrix, file='TC.OandA_I1_16June.RData')

