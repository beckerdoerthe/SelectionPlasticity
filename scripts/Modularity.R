# final version # 20 OCT 2020

rm(list=ls()) 

library(data.table)
library(geomorph)
library(RRPP)
library(reshape2)
library(dplyr)
library(tidyverse)

library(ggplot2)
library(cowplot)
library(viridis)
library(plotrix)


# prep pheno data (same as for TA; contrary to grm data, only data included with more than one replicate, i.e. i>641)
load(file = "/mnt/pricey_4/doerthe/GRM_May2020/shape_use.ag_filtered_final_02June2020.Rdata")  
# final_data

setkey(final_data, cloneid_geno, i) 

# reduce control data (choose random 4 Geno IDs for each treatment group in all three ctrl clones)
set.seed(123)
random_ctrl_O <- final_data[group == 'ctrl_O'][, c('cloneid_geno','Geno','treatment')] %>% group_by(cloneid_geno, treatment) %>% sample_n(4)
ctrlO_clones <- final_data[Geno %in% random_ctrl_O$Geno]

set.seed(123)
random_ctrl_A <- final_data[group == 'ctrl_A'][, c('cloneid_geno','Geno','treatment')] %>% group_by(cloneid_geno, treatment) %>% sample_n(4)
ctrlA_clones <- final_data[Geno %in% random_ctrl_A$Geno]

O_clones <- final_data[group == "O"]
A_clones <- final_data[group == "A"]

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
table(all_data_O.medrd$SC_unique, all_data_O.medrd$treatment, all_data_O.medrd$instar)

# exclude missing data in ctrl and treatment (CLUNKY CODE, but works!)
tbl_all_data_O <- as.data.table(table(all_data_O.medrd$SC_unique, all_data_O.medrd$treatment, all_data_O.medrd$instar))
setnames(tbl_all_data_O, c('V1','V2','V3','N'), c('SC','treatment','instar','count'))
tbl_all_data_O_wide <- as.data.table(dcast(tbl_all_data_O, SC ~ c(paste0("treatment_", treatment, "_instar_", instar)), fun=mean, value.var='count'))

data_O_wide_use <- all_data_O.medrd[SC_unique %in% tbl_all_data_O_wide[treatment_0.5_instar_1 > 641 & treatment_0.5_instar_2 > 641 & treatment_0_instar_1 > 641 & treatment_0_instar_2 > 641]$SC]  
Os <- unique(levels(as.factor(data_O_wide_use$cloneid_geno)))  ### 40 unique clones


# exclude A clones with less than 1 sample per treatment&instar group
table(all_data_As$cloneid_geno, all_data_As$treatment, all_data_As$instar)

# exclude missing data in ctrl and treatment (CLUNKY CODE, but works!)
tbl_all_data_A <- as.data.table(table(all_data_As$cloneid_geno, all_data_As$treatment, all_data_As$instar))
setnames(tbl_all_data_A, c('V1','V2','V3','N'), c('cloneid_geno','treatment','instar','count'))
tbl_all_data_A_wide <- as.data.table(dcast(tbl_all_data_A, cloneid_geno ~ c(paste0("treatment_", treatment, "_instar_", instar)), fun=mean, value.var='count'))

data_A_wide_use <- all_data_As[cloneid_geno %in% tbl_all_data_A_wide[treatment_0.5_instar_1 > 641 & treatment_0.5_instar_2 > 641 & treatment_0_instar_1 > 641 & treatment_0_instar_2 > 641]$cloneid_geno]


# find A clones with highest read depth 
data_A_wide_use[,revorder:=frankv(medrd,order=-1,ties.method = "first")]  
setkey(data_A_wide_use, revorder, Geno, i) 
tmp_medrd_As <- unique(data_A_wide_use[i==150][,c('cloneid_geno', 'SC_unique', 'medrd')]$cloneid_geno)
medrd_As <- tmp_medrd_As[1:40]

all_data_final <- all_data[cloneid_geno %in% Os | cloneid_geno %in% medrd_As]


#load TA data - modify input file here (!)
#load(file="ta_out.OandA_I1.RData")
load(file="ta_out.OandA_I2.RData")
ta_out


### modify here re instar and group
# Os
all_data_use <- all_data_final[, c("cloneid_geno","Geno","SC_unique","group","treatment","instar","i","batch","height","APlength")][instar == 2][group %in% c('O','ctrl_O')][treatment == 0.5]
# As
#all_data_use <- all_data_final[, c("cloneid_geno","Geno","SC_unique","group","treatment","instar","i","batch","height","APlength")][instar == 2][group %in% c('A','ctrl_A')][treatment == 0.5]

setnames(all_data_use, c('APlength','height'), c('x', 'y'))  

all_data_use_wide <- reshape(all_data_use,   ## modify input data set 
                             idvar=c("Geno", "cloneid_geno", "SC_unique", "group", "instar", "treatment", "batch"),
                             timevar="i",
                             direction="wide")

setkey(all_data_use_wide, Geno)


ta_out_red <- as.data.table(ta_out$fit$LM$fitted)
ta_out_red[, Geno := rownames(ta_out$fit$LM$fitted)]
ta_out_red_mod <- ta_out_red %>% select(Geno, everything())
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
  Clone = factor(all_data_use_wide$cloneid_geno),
  Predation = factor(all_data_use_wide$treatment),
  SC = factor(all_data_use_wide$SC_unique),
  Group = factor(all_data_use_wide$group),
  Batch = factor(all_data_use_wide$batch))



## (A1)
land.gps.body100<-c(rep('a',300),rep('b',300))

MT_body100 <- modularity.test(gpa_trans_adj$coords,
                              land.gps.body100, 
                              CI=FALSE, iter=999)

# IT_body100 <- integration.test(gpa_trans_adj$coords,
#                                  partition.gp = land.gps.body100,
#                                  iter=20)

save(MT_body100, file = 'O_MT_A1_I2.RData')


# summary(MT_body100)
# plot(MT_body100)

# IT_body100
# plot(IT_body100)


## (A2)
land.gps.body100<-c(rep('a',250),rep('b',350))

MT_body100 <- modularity.test(gpa_trans_adj$coords,
                              land.gps.body100, 
                              CI=FALSE, iter=999)

# IT_body100 <- integration.test(gpa_trans_adj$coords,
#                                  partition.gp = land.gps.body100,
#                                  iter=20)

save(MT_body100, file = 'O_MT_A2_I2.RData')

# summary(MT_body100)
# plot(MT_body100)

# IT_body100
# plot(IT_body100)


## (A3)
land.gps.body100<-c(rep('a',200),rep('b',400))

MT_body100 <- modularity.test(gpa_trans_adj$coords,
                              land.gps.body100, 
                              CI=FALSE, iter=999)

# IT_body100 <- integration.test(gpa_trans_adj$coords,
#                                  partition.gp = land.gps.body100,
#                                  iter=20)

save(MT_body100, file = 'O_MT_A3_I2.RData')

# summary(MT_body100)
# plot(MT_body100)

# IT_body100
# plot(IT_body100)


## (A4)
land.gps.body100<-c(rep('a',150),rep('b',450))

MT_body100 <- modularity.test(gpa_trans_adj$coords,
                              land.gps.body100, 
                              CI=FALSE, iter=999)

# IT_body100 <- integration.test(gpa_trans_adj$coords,
#                                  partition.gp = land.gps.body100,
#                                  iter=20)

save(MT_body100, file = 'O_MT_A4_I2.RData')

# summary(MT_body100)
# plot(MT_body100)

# IT_body100
# plot(IT_body100)


## (A5)
land.gps.body100<-c(rep('a',100),rep('b',500))

MT_body100 <- modularity.test(gpa_trans_adj$coords,
                              land.gps.body100, 
                              CI=FALSE, iter=999)

# IT_body100 <- integration.test(gpa_trans_adj$coords,
#                                  partition.gp = land.gps.body100,
#                                  iter=20)

save(MT_body100, file = 'O_MT_A5_I2.RData')

# summary(MT_body100)
# plot(MT_body100)

# IT_body100
# plot(IT_body100)




## (B1)
land.gps.body100<-c(rep('a',100),rep('b',200),rep('c',300))

MT_body100 <- modularity.test(gpa_trans_adj$coords,
                              land.gps.body100, 
                              CI=FALSE, iter=999)

# IT_body100 <- integration.test(gpa_trans_adj$coords,
#                                  partition.gp = land.gps.body100,
#                                  iter=20)

save(MT_body100, file = 'O_MT_B1_I2.RData')

# summary(MT_body100)
# plot(MT_body100)

# IT_body100
# plot(IT_body100)



## (B2)
land.gps.body100<-c(rep('a',100),rep('b',100),rep('c',400))

MT_body100 <- modularity.test(gpa_trans_adj$coords,
                              land.gps.body100, 
                              CI=FALSE, iter=999)

# IT_body100 <- integration.test(gpa_trans_adj$coords,
#                                  partition.gp = land.gps.body100,
#                                  iter=20)

save(MT_body100, file = 'O_MT_B2_I2.RData')

# summary(MT_body100)
# plot(MT_body100)

# IT_body100
# plot(IT_body100)


## (B3)
land.gps.body100<-c(rep('a',100),rep('b',150),rep('c',350))

MT_body100 <- modularity.test(gpa_trans_adj$coords,
                              land.gps.body100, 
                              CI=FALSE, iter=999)

# IT_body100 <- integration.test(gpa_trans_adj$coords,
#                                  partition.gp = land.gps.body100,
#                                  iter=20)

save(MT_body100, file = 'O_MT_B3_I2.RData')

# summary(MT_body100)
# plot(MT_body100)

# IT_body100
# plot(IT_body100)










###
load(file = 'MT_A1_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'MT_A2_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'MT_A3_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'MT_A4_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'MT_A5_I1.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'MT_B1_I1.RData')  
summary(MT_body100)
MT_body100$CR.mat
# plot(MT_body100)

load(file = 'MT_B2_I1.RData')  
summary(MT_body100)
MT_body100$CR.mat
# plot(MT_body100)

load(file = 'MT_B3_I1.RData')  
summary(MT_body100)
MT_body100$CR.mat
# plot(MT_body100)




load(file = 'MT_A1_I2.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'MT_A2_I2.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'MT_A3_I2.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'MT_A4_I2.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'MT_A5_I2.RData')  
summary(MT_body100)
# plot(MT_body100)

load(file = 'MT_B1_I2.RData')  
summary(MT_body100)
MT_body100$CR.mat
# plot(MT_body100)

load(file = 'MT_B2_I2.RData')  
summary(MT_body100)
MT_body100$CR.mat
# plot(MT_body100)

load(file = 'MT_B3_I2.RData')  
summary(MT_body100)
MT_body100$CR.mat
# plot(MT_body100)
