# Becker et al 2020 (UVA) - FIG 4

rm(list=ls()) 

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


#######################
### load pheno data ###
#######################
load(file = "~/Google Drive/PNAS_Becker_data_scripts/shape_use.ag_filtered_final_02June2020.Rdata")
# final_data

setkey(final_data, cloneid_geno, i)

# reduce control data (choose random 4 Geno IDs for each treatment group in all three ctrl clones)
set.seed(123)
random_ctrl_O <- final_data[group == 'ctrl_O'][, c('cloneid_geno','Geno','treatment')] %>% group_by(cloneid_geno, treatment) %>% sample_n(4)
ctrlO_clones <- final_data[Geno %in% random_ctrl_O$Geno][, -"max_height"]

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


all_data_final <- all_data[cloneid_geno %in% Os | cloneid_geno %in% medrd_As]
all_data_final[, i:= as.numeric(i)]
all_data_final[, instar_new := ifelse(all_data_final$instar == 1, 'instar 1', 'instar 2')]
all_data_final[, SC_group_new := ifelse(all_data_final$SC_group == "O", 'cluster O', 'cluster A')]
all_data_final[, treatment_new := ifelse(all_data_final$treatment == 0, 'C', 'P')]
all_data_final[, max_height_new := max(height, na.rm=T), by=c('Geno','instar')]

setkey(all_data_final, cloneid_geno, i) 

all_data_final[, deme_new := ifelse(all_data_final$deme == 3, 2, all_data_final$deme)]

all_data_final[, replicate_new := ifelse(all_data_final$replicate == "1C", "1A", 
                                   ifelse(all_data_final$replicate == "1D", "1B", 
                                    ifelse(all_data_final$replicate == "2C", "2A", 
                                      ifelse(all_data_final$replicate == "2D", "2B", 
                                       ifelse(all_data_final$replicate == "3A", "2A", 
                                        ifelse(all_data_final$replicate == "3B", "2B", all_data_final$replicate))))))]

all_data_final[, clutch := ifelse(all_data_final$replicate_new %like% "A", "A", "B")]
        
 
## within clutch - i.e. twins born to the same mum
# reduce to one i-th position for animal length, max height, and eye area
all_data_use <- all_data_final[i == 150]

all_data_wide <- as.data.table(dcast(all_data_use, 
                                     SC_group + cloneid_geno + instar + treatment + deme ~ paste0('clutch', clutch), 
                                     fun = mean, 
                                     value.var = "length"))  # change here for max_height_new / length

setkey(all_data_wide, SC_group, instar, treatment)

clone.tab <- all_data_wide[,list(clone=unique(cloneid_geno)), list(SC_group)]


o <- foreach(p=0:1000)%dopar%{
          print(p)
          clone.perm.repA <- clone.tab[,list(cloneid_geno_new=sample(clone), cloneid_geno=clone), list(SC_group)]
          clone.perm.repB <- clone.tab[,list(cloneid_geno_new=sample(clone), cloneid_geno=clone), list(SC_group)]
          
          repA.tmp <- all_data_wide[,c("SC_group", "instar", "treatment", "clutchA", "cloneid_geno"), with=F]
          repB.tmp <- all_data_wide[,c("SC_group", "instar", "treatment", "clutchB", "cloneid_geno"), with=F]
          
          setkey(repA.tmp, SC_group, cloneid_geno)
          setkey(repB.tmp, SC_group, cloneid_geno)
          
          setkey(clone.perm.repA, SC_group, cloneid_geno)
          setkey(clone.perm.repB, SC_group, cloneid_geno)
          
          if(p==0) {
            repA.tmp.perm <- repA.tmp
            repB.tmp.perm <- repB.tmp
            
          }else if(p>0) {
            repA.tmp.perm  <- repA.tmp[J(clone.perm.repA)]
            setnames(repA.tmp.perm, c("cloneid_geno", "cloneid_geno_new"), c("cloneid_geno_orig", "cloneid_geno"))
            
            repB.tmp.perm  <- repB.tmp[J(clone.perm.repB)]
            setnames(repB.tmp.perm, c("cloneid_geno", "cloneid_geno_new"), c("cloneid_geno_orig", "cloneid_geno"))
            
          }
          
          setkey(repA.tmp.perm, SC_group, cloneid_geno, instar, treatment)
          setkey(repB.tmp.perm, SC_group, cloneid_geno, instar, treatment)
          
          perm.m <- merge(repA.tmp.perm, repB.tmp.perm)
          
          #perm.m[,list(cor=cor(rep1A, rep1B, "kendall", use="complete"), 
          #			   cor.p=cor.test(rep1A, rep1B, method="kendall")[[3]][[1]],
          #			   perm=p), 
          #			   list(SC_group, i, instar, treatment)]
          
          perm.m[,list(cor=robcor(clutchA, clutchB), perm=p), list(SC_group, instar, treatment)]
          
        }


o_out <- rbindlist(o)

o.ag <- o_out[,list(cor.mu=median(cor, na.rm=T), cor.lci=quantile(cor, .05, na.rm=T), cor.uci=quantile(cor, .95, na.rm=T)), 
              list(SC_group, instar, treatment, perm=(perm>0))]

# save(o_out, o.ag, file = 'AvsB_length.RData')
# save(o_out, o.ag, file = 'AvsB_maxHeight.RData')


## across demes - i.e. offspring born to different mums but same field isolate
# reduce to one i-th position for animal length, max height, and eye area
all_data_use <- all_data_final[i == 150]

all_data_wide <- as.data.table(dcast(all_data_use, 
                                     SC_group + cloneid_geno + instar + treatment + clutch ~ paste0('deme', deme_new), 
                                     fun = mean, 
                                     value.var = "length"))  # change here for max_height_new / length

setkey(all_data_wide, SC_group, instar, treatment)

clone.tab <- all_data_wide[,list(clone=unique(cloneid_geno)), list(SC_group)]


o <- foreach(p=0:1000)%dopar%{
          print(p)
          clone.perm.deme1 <- clone.tab[,list(cloneid_geno_new=sample(clone), cloneid_geno=clone), list(SC_group)]
          clone.perm.deme2 <- clone.tab[,list(cloneid_geno_new=sample(clone), cloneid_geno=clone), list(SC_group)]
          
          deme1.tmp <- all_data_wide[,c("SC_group", "instar", "treatment", "deme1", "cloneid_geno"), with=F]
          deme2.tmp <- all_data_wide[,c("SC_group", "instar", "treatment", "deme2", "cloneid_geno"), with=F]
          
          setkey(deme1.tmp, SC_group, cloneid_geno)
          setkey(deme2.tmp, SC_group, cloneid_geno)
          
          setkey(clone.perm.deme1, SC_group, cloneid_geno)
          setkey(clone.perm.deme2, SC_group, cloneid_geno)
          
          if(p==0) {
            deme1.tmp.perm <- deme1.tmp
            deme2.tmp.perm <- deme2.tmp
            
          }else if(p>0) {
            deme1.tmp.perm  <- deme1.tmp[J(clone.perm.deme1)]
            setnames(deme1.tmp.perm, c("cloneid_geno", "cloneid_geno_new"), c("cloneid_geno_orig", "cloneid_geno"))
            
            deme2.tmp.perm  <- deme2.tmp[J(clone.perm.deme2)]
            setnames(deme2.tmp.perm, c("cloneid_geno", "cloneid_geno_new"), c("cloneid_geno_orig", "cloneid_geno"))
            
          }
          
          setkey(deme1.tmp.perm, SC_group, cloneid_geno, instar, treatment)
          setkey(deme2.tmp.perm, SC_group, cloneid_geno, instar, treatment)
          
          perm.m <- merge(deme1.tmp.perm, deme2.tmp.perm)
          
          #perm.m[,list(cor=cor(rep1A, rep1B, "kendall", use="complete"), 
          #			   cor.p=cor.test(rep1A, rep1B, method="kendall")[[3]][[1]],
          #			   perm=p), 
          #			   list(SC_group, i, instar, treatment)]
          
          perm.m[,list(cor=robcor(deme1, deme2), perm=p), list(SC_group, instar, treatment)]
          
        }


o_out <- rbindlist(o)

o.ag <- o_out[,list(cor.mu=median(cor, na.rm=T), cor.lci=quantile(cor, .05, na.rm=T), cor.uci=quantile(cor, .95, na.rm=T)), 
              list(SC_group, instar, treatment, perm=(perm>0))]

# save(o_out, o.ag, file = 'd1vsd2_length.RData')
# save(o_out, o.ag, file = 'd1vsd2_maxHeight.RData')



## across As - i.e. randomly draw offspring among As
# reduce to one i-th position for animal length, max height, and eye area
all_data_use <- all_data_final[i == 150]

all_data_wide <- as.data.table(dcast(all_data_use[SC_group == "A"], 
                                     SC_group + cloneid_geno + instar + treatment ~ paste0('SC_', SC_group), 
                                     fun = mean, 
                                     value.var = "length"))  # change here for max_height_new / length

setkey(all_data_wide, SC_group, instar, treatment)

clone.tab <- all_data_wide[,list(clone=unique(cloneid_geno)), list(SC_group)]

all_data_small <- all_data_use[SC_group == "A"][, list(SC_group, cloneid_geno, instar, treatment, max_height_new)]


o <- foreach(b=1:50)%do%{  #
          print(b)
          
          foreach(p=0:1000)%dopar%{  ## increase to 1000
            print(p)
            
              foreach(instar.i = unique(all_data_small$instar), .errorhandling="remove", .combine="rbind") %do% {
                
                foreach(treatment.i = unique(all_data_small$treatment), .errorhandling="remove", .combine="rbind") %do% {
                  
                  clone.perm.d1 <- clone.tab[,list(cloneid_geno_new=sample(clone), cloneid_geno=clone), list(SC_group)]
                  clone.perm.d2 <- clone.tab[,list(cloneid_geno_new=sample(clone), cloneid_geno=clone), list(SC_group)]
                  
                  d1.tmp <- all_data_wide[,c("SC_group", "instar", "treatment", "SC_A", "cloneid_geno"), with=F][instar == instar.i][treatment == treatment.i]
                  d2.tmp <- all_data_wide[,c("SC_group", "instar", "treatment", "SC_A", "cloneid_geno"), with=F][instar == instar.i][treatment == treatment.i]
                  
                  setkey(d1.tmp, SC_group, cloneid_geno)
                  setkey(d2.tmp, SC_group, cloneid_geno)
                  
                  setkey(clone.perm.d1, SC_group, cloneid_geno)
                  setkey(clone.perm.d2, SC_group, cloneid_geno)
                  
                  if(p==0) {
                    d1.tmp.perm <- d1.tmp
                    d2.tmp.perm <- d2.tmp
                    
                  }else if(p>0) {
                    d1.tmp.perm  <- d1.tmp[J(clone.perm.d1[SC_group == "A"])]
                    setnames(d1.tmp.perm, c("cloneid_geno", "cloneid_geno_new"), c("cloneid_geno_orig", "cloneid_geno"))
                    
                    d2.tmp.perm  <- d2.tmp[SC_group == "A"][J(clone.perm.d2[SC_group == "A"])]
                    setnames(d2.tmp.perm, c("cloneid_geno", "cloneid_geno_new"), c("cloneid_geno_orig", "cloneid_geno"))
                    
                  }
                  
                  setkey(d1.tmp.perm, instar, treatment)
                  setkey(d2.tmp.perm, instar, treatment)
                  
                  # only sample one of the two sets
                  d2.tmp.perm_draw <- d2.tmp.perm[SC_group == "A"][sample(.N, nrow(d1.tmp.perm[SC_group == "A"]))]
                  d2.tmp.perm_draw[, SC_A_draw := SC_A]
                  d2.tmp.perm_draw[, cloneid_geno_draw := cloneid_geno]
                  
                  perm.m <- cbind(d1.tmp.perm,d2.tmp.perm_draw[,-"SC_A"])
                  # table(perm.m[,cloneid_geno] != perm.m[,cloneid_geno_draw])  ## remove ?!
                  
                  #perm.m[,list(cor=cor(d1, d2, "kendall", use="complete"), 
                  #			   cor.p=cor.test(d1, d2, method="kendall")[[3]][[1]],
                  #			   perm=p), 
                  #			   list(SC_group, i, instar, treatment)]
                  
                  cor_out <- robcor(perm.m$SC_A, perm.m$SC_A_draw)
                  
                  ## data out
                  data.table(instar = instar.i,
                             treatment = treatment.i,
                             perm = p,
                             boot = b, 
                             cor = cor_out)
                  
                  
                }}}}


o_out <- rbindlist(unlist(o, recursive = FALSE))

o.ag <- o_out[,list(cor.mu=median(cor, na.rm=T), cor.lci=quantile(cor, .05, na.rm=T), cor.uci=quantile(cor, .95, na.rm=T)), 
              list(instar, treatment, boot, perm=(perm>0))]

o.ag.ag <- o_out[,list(cor.mu=median(cor, na.rm=T), cor.lci=quantile(cor, .05, na.rm=T), cor.uci=quantile(cor, .95, na.rm=T)), 
                 list(instar, treatment, perm=(perm>0))]


# save(o_out, o.ag, o.ag.ag, file = 'SC_AvsSC_A_length.RData')

# save(o_out, o.ag, o.ag.ag, file = 'SC_AvsSC_A_maxHeight.RData')



## stats

# within clutch

# length
load(file = '~/Google Drive/PNAS_Becker_data_scripts/AvsB_length.RData')

clutch_length <- o_out[SC_group == 'A'][, -'SC_group']

# I1 control
clutch_length_I1_C <- wilcox.test(clutch_length[instar == 1][treatment == 0][perm > 0]$cor, mu=clutch_length[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clutch_length_I1_C
clutch_Zscore_length_I1_C = qnorm(clutch_length_I1_C$p.value/2)
clutch_Zscore_length_I1_C

# Z = -27.39292, p < 2.2e-16

# I1 predation
clutch_length_I1_P <- wilcox.test(clutch_length[instar == 1][treatment == 0.5][perm > 0]$cor, mu=clutch_length[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_length_I1_P
clutch_Zscore_length_I1_P = qnorm(clutch_length_I1_P$p.value/2)
clutch_Zscore_length_I1_P

# Z = -27.39292, p < 2.2e-16

# I2 control
clutch_length_I2_C <- wilcox.test(clutch_length[instar == 2][treatment == 0][perm > 0]$cor, mu=clutch_length[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clutch_length_I2_C
clutch_Zscore_length_I2_C = qnorm(clutch_length_I2_C$p.value/2)
clutch_Zscore_length_I2_C

# Z = -27.37048, p < 2.2e-16

# I2 predation
clutch_length_I2_C <- wilcox.test(clutch_length[instar == 2][treatment == 0.5][perm > 0]$cor, mu=clutch_length[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_length_I2_C
clutch_Zscore_length_I2_C = qnorm(clutch_length_I2_C$p.value/2)
clutch_Zscore_length_I2_C

# Z = -27.39281, p < 2.2e-16


# max height
load(file = '~/Google Drive/PNAS_Becker_data_scripts/AvsB_maxHeight.RData')

clutch_maxHeight <- o_out[SC_group == 'A'][, -'SC_group']

# I1 control
clutch_maxHeight_I1_C <- wilcox.test(clutch_maxHeight[instar == 1][treatment == 0][perm > 0]$cor, mu=clutch_maxHeight[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clutch_maxHeight_I1_C
clutch_Zscore_maxHeight_I1_C = qnorm(clutch_maxHeight_I1_C$p.value/2)
clutch_Zscore_maxHeight_I1_C

# Z = -25.8415, p < 2.2e-16

# I1 predation
clutch_maxHeight_I1_P <- wilcox.test(clutch_maxHeight[instar == 1][treatment == 0.5][perm > 0]$cor, mu=clutch_maxHeight[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_maxHeight_I1_P
clutch_Zscore_maxHeight_I1_P = qnorm(clutch_maxHeight_I1_P$p.value/2)
clutch_Zscore_maxHeight_I1_P

# Z = -27.38974, p < 2.2e-16

# I2 control
clutch_maxHeight_I2_C <- wilcox.test(clutch_maxHeight[instar == 2][treatment == 0][perm > 0]$cor, mu=clutch_maxHeight[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clutch_maxHeight_I2_C
clutch_Zscore_maxHeight_I2_C = qnorm(clutch_maxHeight_I2_C$p.value/2)
clutch_Zscore_maxHeight_I2_C

# Z = -21.2732, p < 2.2e-16

# I2 predation
clutch_maxHeight_I2_C <- wilcox.test(clutch_maxHeight[instar == 2][treatment == 0.5][perm > 0]$cor, mu=clutch_maxHeight[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_maxHeight_I2_C
clutch_Zscore_maxHeight_I2_C = qnorm(clutch_maxHeight_I2_C$p.value/2)
clutch_Zscore_maxHeight_I2_C

# Z = -27.38679, p < 2.2e-16


# within clone

# length

load(file = '~/Google Drive/PNAS_Becker_data_scripts/d1vsd2_length.RData')

clone_length <- o_out[SC_group == 'A'][, -'SC_group']

# I1 control
clone_length_I1_C <- wilcox.test(clone_length[instar == 1][treatment == 0][perm > 0]$cor, mu=clone_length[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_length_I1_C
clone_Zscore_length_I1_C = qnorm(clone_length_I1_C$p.value/2)
clone_Zscore_length_I1_C

# Z = -27.0619, p < 2.2e-16

# I1 predation
clone_length_I1_P <- wilcox.test(clone_length[instar == 1][treatment == 0.5][perm > 0]$cor, mu=clone_length[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_length_I1_P
clone_Zscore_length_I1_P = qnorm(clone_length_I1_P$p.value/2)
clone_Zscore_length_I1_P

# Z = -27.38722, p < 2.2e-16

# I2 control
clone_length_I2_C <- wilcox.test(clone_length[instar == 2][treatment == 0][perm > 0]$cor, mu=clone_length[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_length_I2_C
clone_Zscore_length_I2_C = qnorm(clone_length_I2_C$p.value/2)
clone_Zscore_length_I2_C

# Z = -26.69225, p < 2.2e-16

# I2 predation
clone_length_I2_C <- wilcox.test(clone_length[instar == 2][treatment == 0.5][perm > 0]$cor, mu=clone_length[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_length_I2_C
clone_Zscore_length_I2_C = qnorm(clone_length_I2_C$p.value/2)
clone_Zscore_length_I2_C

# Z = -15.02519, p < 2.2e-16


# max height
load(file = '~/Google Drive/PNAS_Becker_data_scripts/d1vsd2_maxHeight.RData')

clone_maxHeight <- o_out[SC_group == 'A'][, -'SC_group']

# I1 control
clone_maxHeight_I1_C <- wilcox.test(clone_maxHeight[instar == 1][treatment == 0][perm > 0]$cor, mu=clone_maxHeight[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_maxHeight_I1_C
clone_Zscore_maxHeight_I1_C = qnorm(clone_maxHeight_I1_C$p.value/2)
clone_Zscore_maxHeight_I1_C

# Z = -14.08415, p < 2.2e-16

# I1 predation
clone_maxHeight_I1_P <- wilcox.test(clone_maxHeight[instar == 1][treatment == 0.5][perm > 0]$cor, mu=clone_maxHeight[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_maxHeight_I1_P
clone_Zscore_maxHeight_I1_P = qnorm(clone_maxHeight_I1_P$p.value/2)
clone_Zscore_maxHeight_I1_P

# Z = -27.36971, p < 2.2e-16

# I2 control
clone_maxHeight_I2_C <- wilcox.test(clone_maxHeight[instar == 2][treatment == 0][perm > 0]$cor, mu=clone_maxHeight[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_maxHeight_I2_C
clone_Zscore_maxHeight_I2_C = qnorm(clone_maxHeight_I2_C$p.value/2)
clone_Zscore_maxHeight_I2_C

# Z = -25.59095, p < 2.2e-16

# I2 predation
clone_maxHeight_I2_C <- wilcox.test(clone_maxHeight[instar == 2][treatment == 0.5][perm > 0]$cor, mu=clone_maxHeight[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_maxHeight_I2_C
clone_Zscore_maxHeight_I2_C = qnorm(clone_maxHeight_I2_C$p.value/2)
clone_Zscore_maxHeight_I2_C

# Z = -27.39292, p < 2.2e-16



# A vs A (random draw)

# length

load(file = '~/Google Drive/PNAS_Becker_data_scripts/SC_AvsSC_A_length.RData')

AvsA_length <- o_out
AvsA_length[, permuted := ifelse(AvsA_length$perm == 0, "no", "yes")]

# I1 control
AvsA_length_I1_C <- wilcox.test(cor ~ permuted, data = AvsA_length[instar == 1][treatment == 0]) 
AvsA_length_I1_C
AvsA_Zscore_length_I1_C = qnorm(AvsA_length_I1_C$p.value/2)
AvsA_Zscore_length_I1_C

# Z = -0.9734845, p = 0.3303

# I1 predation
AvsA_length_I1_P <- wilcox.test(cor ~ permuted, data = AvsA_length[instar == 1][treatment == 0.5]) 
AvsA_length_I1_P
AvsA_Zscore_length_I1_P = qnorm(AvsA_length_I1_P$p.value/2)
AvsA_Zscore_length_I1_P

# Z = -0.5066048, p = 0.6124

# I2 control
AvsA_length_I2_C <- wilcox.test(cor ~ permuted, data = AvsA_length[instar == 2][treatment == 0]) 
AvsA_length_I2_C
AvsA_Zscore_length_I2_C = qnorm(AvsA_length_I2_C$p.value/2)
AvsA_Zscore_length_I2_C

# Z = -1.908762, p = 0.05629

# I2 predation
AvsA_length_I2_P <- wilcox.test(cor ~ permuted, data = AvsA_length[instar == 2][treatment == 0.5]) 
AvsA_length_I2_P
AvsA_Zscore_length_I2_P = qnorm(AvsA_length_I2_P$p.value/2)
AvsA_Zscore_length_I2_P

# Z = -0.8004232, p = 0.4235


# max height
load(file = '~/Google Drive/PNAS_Becker_data_scripts/SC_AvsSC_A_maxHeight.RData')

AvsA_maxHeight <- o_out
AvsA_maxHeight[, permuted := ifelse(AvsA_maxHeight$perm == 0, "no", "yes")]


# I1 control
AvsA_maxHeight_I1_C <- wilcox.test(cor ~ permuted, data = AvsA_maxHeight[instar == 1][treatment == 0]) 
AvsA_maxHeight_I1_C
AvsA_Zscore_maxHeight_I1_C = qnorm(AvsA_maxHeight_I1_C$p.value/2)
AvsA_Zscore_maxHeight_I1_C

# Z = -0.8263647, p = 0.4086

# I1 predation
AvsA_maxHeight_I1_P <- wilcox.test(cor ~ permuted, data = AvsA_maxHeight[instar == 1][treatment == 0.5]) 
AvsA_maxHeight_I1_P
AvsA_Zscore_maxHeight_I1_P = qnorm(AvsA_maxHeight_I1_P$p.value/2)
AvsA_Zscore_maxHeight_I1_P

# Z = -0.7613982, p = 0.4464

# I2 control
AvsA_maxHeight_I2_C <- wilcox.test(cor ~ permuted, data = AvsA_maxHeight[instar == 2][treatment == 0]) 
AvsA_maxHeight_I2_C
AvsA_Zscore_maxHeight_I2_C = qnorm(AvsA_maxHeight_I2_C$p.value/2)
AvsA_Zscore_maxHeight_I2_C

# Z = -1.478488, p = 0.1393

# I2 predation
AvsA_maxHeight_I2_P <- wilcox.test(cor ~ permuted, data = AvsA_maxHeight[instar == 2][treatment == 0.5]) 
AvsA_maxHeight_I2_P
AvsA_Zscore_maxHeight_I2_P = qnorm(AvsA_maxHeight_I2_P$p.value/2)
AvsA_Zscore_maxHeight_I2_P

# Z = -0.2268296, p = 0.8206



# within clutch vs among clones

# length

# I1 control
clutch_vs_As_length_I1_C <- wilcox.test(AvsA_length[instar == 1][treatment == 0][perm == 0]$cor, mu=clutch_length[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_length_I1_C
clutch_vs_As_Zscore_length_I1_C = qnorm(clutch_vs_As_length_I1_C$p.value/2)
clutch_vs_As_Zscore_length_I1_C

# Z = -6.149139, p = 7.79e-10

# I1 predation
clutch_vs_As_length_I1_P <- wilcox.test(AvsA_length[instar == 1][treatment == 0.5][perm == 0]$cor, mu=clutch_length[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_length_I1_P
clutch_vs_As_Zscore_length_I1_P = qnorm(clutch_vs_As_length_I1_P$p.value/2)
clutch_vs_As_Zscore_length_I1_P

# Z = -6.129832, p = 8.797e-10

# I2 control
clutch_vs_As_length_I2_C <- wilcox.test(AvsA_length[instar == 2][treatment == 0][perm == 0]$cor, mu=clutch_length[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_length_I2_C
clutch_vs_As_Zscore_length_I2_C = qnorm(clutch_vs_As_length_I2_C$p.value/2)
clutch_vs_As_Zscore_length_I2_C

# Z = -6.091219, p = 1.121e-09

# I2 predation
clutch_vs_As_length_I2_P <- wilcox.test(AvsA_length[instar == 2][treatment == 0.5][perm == 0]$cor, mu=clutch_length[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_length_I2_P
clutch_vs_As_Zscore_length_I2_P = qnorm(clutch_vs_As_length_I2_P$p.value/2)
clutch_vs_As_Zscore_length_I2_P

# Z = -6.023646, p = 1.705e-09


# max height

# I1 control
clutch_vs_As_maxHeight_I1_C <- wilcox.test(AvsA_maxHeight[instar == 1][treatment == 0][perm == 0]$cor, mu=clutch_maxHeight[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_maxHeight_I1_C
clutch_vs_As_Zscore_maxHeight_I1_C = qnorm(clutch_vs_As_maxHeight_I1_C$p.value/2)
clutch_vs_As_Zscore_maxHeight_I1_C

# Z = -4.102643, p = 4.085e-05

# I1 predation
clutch_vs_As_maxHeight_I1_P <- wilcox.test(AvsA_maxHeight[instar == 1][treatment == 0.5][perm == 0]$cor, mu=clutch_maxHeight[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_maxHeight_I1_P
clutch_vs_As_Zscore_maxHeight_I1_P = qnorm(clutch_vs_As_maxHeight_I1_P$p.value/2)
clutch_vs_As_Zscore_maxHeight_I1_P

# Z = -6.091219, p = 1.121e-09

# I2 control
clutch_vs_As_maxHeight_I2_C <- wilcox.test(AvsA_maxHeight[instar == 2][treatment == 0][perm == 0]$cor, mu=clutch_maxHeight[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_maxHeight_I2_C
clutch_vs_As_Zscore_maxHeight_I2_C = qnorm(clutch_vs_As_maxHeight_I2_C$p.value/2)
clutch_vs_As_Zscore_maxHeight_I2_C

# Z = -5.019705, p = 5.175e-07

# I2 predation
clutch_vs_As_maxHeight_I2_P <- wilcox.test(AvsA_maxHeight[instar == 2][treatment == 0.5][perm == 0]$cor, mu=clutch_maxHeight[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_maxHeight_I2_P
clutch_vs_As_Zscore_maxHeight_I2_P = qnorm(clutch_vs_As_maxHeight_I2_P$p.value/2)
clutch_vs_As_Zscore_maxHeight_I2_P

# Z = -6.149139, p = 7.79e-10



# within clone vs among clones

# length

# I1 control
clone_vs_As_length_I1_C <- wilcox.test(AvsA_length[instar == 1][treatment == 0][perm == 0]$cor, mu=clone_length[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_length_I1_C
clone_vs_As_Zscore_length_I1_C = qnorm(clone_vs_As_length_I1_C$p.value/2)
clone_vs_As_Zscore_length_I1_C

# Z = -6.062259, p = 1.342e-09

# I1 predation
clone_vs_As_length_I1_P <- wilcox.test(AvsA_length[instar == 1][treatment == 0.5][perm == 0]$cor, mu=clone_length[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_length_I1_P
clone_vs_As_Zscore_length_I1_P = qnorm(clone_vs_As_length_I1_P$p.value/2)
clone_vs_As_Zscore_length_I1_P

# Z = -5.965726, p = 2.435e-09

# I2 control
clone_vs_As_length_I2_C <- wilcox.test(AvsA_length[instar == 2][treatment == 0][perm == 0]$cor, mu=clone_length[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_length_I2_C
clone_vs_As_Zscore_length_I2_C = qnorm(clone_vs_As_length_I2_C$p.value/2)
clone_vs_As_Zscore_length_I2_C

# Z = -5.309303, p = 1.1e-07

# I2 predation
clone_vs_As_length_I2_P <- wilcox.test(AvsA_length[instar == 2][treatment == 0.5][perm == 0]$cor, mu=clone_length[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_length_I2_P
clone_vs_As_Zscore_length_I2_P = qnorm(clone_vs_As_length_I2_P$p.value/2)
clone_vs_As_Zscore_length_I2_P

# Z = -2.644998, p = 0.008169


# max height

# I1 control
clone_vs_As_maxHeight_I1_C <- wilcox.test(AvsA_maxHeight[instar == 1][treatment == 0][perm == 0]$cor, mu=clone_maxHeight[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_maxHeight_I1_C
clone_vs_As_Zscore_maxHeight_I1_C = qnorm(clone_vs_As_maxHeight_I1_C$p.value/2)
clone_vs_As_Zscore_maxHeight_I1_C

# Z = -1.361112, p = 0.1735

# I1 predation
clone_vs_As_maxHeight_I1_P <- wilcox.test(AvsA_maxHeight[instar == 1][treatment == 0.5][perm == 0]$cor, mu=clone_maxHeight[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_maxHeight_I1_P
clone_vs_As_Zscore_maxHeight_I1_P = qnorm(clone_vs_As_maxHeight_I1_P$p.value/2)
clone_vs_As_Zscore_maxHeight_I1_P

# Z = -5.94642, p = 2.741e-09

# I2 control
clone_vs_As_maxHeight_I2_C <- wilcox.test(AvsA_maxHeight[instar == 2][treatment == 0][perm == 0]$cor, mu=clone_maxHeight[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_maxHeight_I2_C
clone_vs_As_Zscore_maxHeight_I2_C = qnorm(clone_vs_As_maxHeight_I2_C$p.value/2)
clone_vs_As_Zscore_maxHeight_I2_C

# Z = -5.772661, p = 7.803e-09

# I2 predation
clone_vs_As_maxHeight_I2_P <- wilcox.test(AvsA_maxHeight[instar == 2][treatment == 0.5][perm == 0]$cor, mu=clone_maxHeight[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_maxHeight_I2_P
clone_vs_As_Zscore_maxHeight_I2_P = qnorm(clone_vs_As_maxHeight_I2_P$p.value/2)
clone_vs_As_Zscore_maxHeight_I2_P

# Z = -6.149139, p = 7.79e-10





## plots

# within clutch

load(file = '~/Google Drive/PNAS_Becker_data_scripts/AvsB_length.RData')

clutch_length <- o.ag[SC_group == 'A'][, -'SC_group']
clutch_length[, group := "within clutch"]
clutch_length[, set := "length"]


load(file = '~/Google Drive/PNAS_Becker_data_scripts/AvsB_maxHeight.RData')

clutch_maxHeight <- o.ag[SC_group == 'A'][, -'SC_group']
clutch_maxHeight[, group := "within clutch"]
clutch_maxHeight[, set := "max dorsal height"]


# within clone

load(file = '~/Google Drive/PNAS_Becker_data_scripts/d1vsd2_length.RData')

clone_length <- o.ag[SC_group == 'A'][, -'SC_group']
clone_length[, group := "within clone"]
clone_length[, set := "length"]


load(file = '~/Google Drive/PNAS_Becker_data_scripts/d1vsd2_maxHeight.RData')

clone_maxHeight <- o.ag[SC_group == 'A'][, -'SC_group']
clone_maxHeight[, group := "within clone"]
clone_maxHeight[, set := "max dorsal height"]


# A vs A (random draw)

load(file = '~/Google Drive/PNAS_Becker_data_scripts/SC_AvsSC_A_length.RData')

AvsA_length <- o.ag
AvsA_length[, group := "among clones"]
AvsA_length[, set := "length"]


load(file = '~/Google Drive/PNAS_Becker_data_scripts/SC_AvsSC_A_length.ag.RData')

AvsA_length.ag <- o.ag.ag
AvsA_length.ag[, group := "among clones"]
AvsA_length.ag[, set := "length"]


load(file = '~/Google Drive/PNAS_Becker_data_scripts/SC_AvsSC_A_maxHeight.RData')

AvsA_maxHeight <- o.ag
AvsA_maxHeight[, group := "among clones"]
AvsA_maxHeight[, set := "max dorsal height"]


load(file = '~/Google Drive/PNAS_Becker_data_scripts/SC_AvsSC_A_maxHeight.ag.RData')

AvsA_maxHeight.ag <- o.ag.ag
AvsA_maxHeight.ag[, group := "among clones"]
AvsA_maxHeight.ag[, set := "max dorsal height"]


all_R2s <- rbind(clutch_length,
                 clutch_maxHeight,
                 clone_length,
                 clone_maxHeight,
                 #AvsA_length, 
                 AvsA_length.ag,
                 #AvsA_maxHeight, 
                 AvsA_maxHeight.ag)

all_R2s[, group := as.factor(group)]
all_R2s[, instar_new := ifelse(all_R2s$instar == 1, "instar 1", "instar 2")]

all_R2s_use <- all_R2s[set %in% c('length', 'max dorsal height')]


R2_all_plot <- all_R2s_use %>% mutate(group = fct_relevel(group, "within clutch", "within clone", "among clones")) %>%
  
                  ggplot(aes(x=as.factor(set), y=cor.mu)) +  
                
                      geom_pointrange(aes(ymin = cor.lci, ymax = cor.uci, colour = as.factor(treatment+perm), group = as.factor(treatment+perm)), position = position_dodge(width = 0.6), size = 0.7) +	
                      
                      facet_grid(~instar_new ~ group) +
                      
                      ylim(-0.33, 0.6) +

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
                            axis.title.x = element_blank(),
                            axis.title.y = element_text(size=15, family='Arial'),
                            axis.text.x = element_text(size=15, family='Arial', angle = 30, hjust=1), 
                            axis.text.y = element_text(size=15, family='Arial'), 
                            strip.text.x = element_text(size = 15, color = "black"),
                            strip.text.y = element_text(size = 15, color = "black"),
                            panel.spacing.x = unit(6, "mm"),
                            panel.spacing.y = unit(6, "mm")) 

# tiff(file = "R2_supplfig4.tiff")
R2_all_plot + scale_color_manual(values=c("#000000","#FF0000","#D3D3D3","#fed4d2")) + labs(x = "phenotypic trait", y = "correlation coefficient") 
# dev.off()



## Vg/Vm ratio
# log(Va) - log(Vm)

load(file = "~/Google Drive/PNAS_Becker_data_scripts/H2_0_I1_ratio_batch.RData")
H2_ctrl_I1_ratio <- as.data.table(collected2)
H2_ctrl_I1_ratio[, i:= as.numeric(i)]
H2_ctrl_I1_ratio[, instar:= "instar 1"]
H2_ctrl_I1_ratio[, group:= "ctrl"]
setkey(H2_ctrl_I1_ratio, group,instar,i)

load(file = "~/Google Drive/PNAS_Becker_data_scripts/H2_0_I2_ratio_batch.RData")
H2_ctrl_I2_ratio <- as.data.table(collected2)
H2_ctrl_I2_ratio[, i:= as.numeric(i)]
H2_ctrl_I2_ratio[, instar:= "instar 2"]
H2_ctrl_I2_ratio[, group:= "ctrl"]
setkey(H2_ctrl_I2_ratio, group,instar,i)

load(file = "~/Google Drive/PNAS_Becker_data_scripts/H2_05_I1_ratio_batch.RData")
H2_trt_I1_ratio <- as.data.table(collected2)
H2_trt_I1_ratio[, i:= as.numeric(i)]
H2_trt_I1_ratio[, instar:= "instar 1"]
H2_trt_I1_ratio[, group:= "trt"]
setkey(H2_trt_I1_ratio, group,instar,i)

load(file = "~/Google Drive/PNAS_Becker_data_scripts/H2_05_I2_ratio_batch.RData")
H2_trt_I2_ratio <- as.data.table(collected2)
H2_trt_I2_ratio[, i:= as.numeric(i)]
H2_trt_I2_ratio[, instar:= "instar 2"]
H2_trt_I2_ratio[, group:= "trt"]
setkey(H2_trt_I2_ratio, group,instar,i)


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


H2_ration_all <- rbind(H2_ctrl_I1_ratio_use,
                       H2_ctrl_I2_ratio_use,
                       H2_trt_I1_ratio_use,
                       H2_trt_I2_ratio_use)


VaVm_instar_i_plot <- ggplot(data = H2_ration_all[i <= 600]) +  
  
                        geom_rect(aes(xmin = 10, xmax = 300, ymin = -Inf, ymax = Inf),fill = "#F5F5F5", colour="#F5F5F5", alpha=0.1) + 
                        geom_rect(aes(xmin = 301, xmax = 600, ymin = -Inf, ymax = Inf),fill = "#DCDCDC", colour="#DCDCDC", alpha=0.1) + 
                        
                        geom_ribbon(data=H2_ration_all[i <= 600][group == "trt"], aes(ymin=stuff.lCI, ymax=stuff.uCI, x=i), fill = "#FF0000", alpha = 0.3)+
                        geom_line(data=H2_ration_all[i <= 600][group == "trt"], aes(x=i, y=stuff.mean), size = 1, colour = "#FF0000") + 
                        
                        geom_ribbon(data=H2_ration_all[i <= 600][group == "ctrl"], aes(ymin=stuff.lCI, ymax=stuff.uCI, x=i), fill = "#000000", alpha = 0.3)+
                        geom_line(data=H2_ration_all[i <= 600][group == "ctrl"], aes(x=i, y=stuff.mean), size = 1, colour = "#000000") + 
                        
                        geom_hline(yintercept=0, linetype = "dotted", size = 0.5, colour = "black") +
                        
                        #ylim(0,0.58) +   
                        
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

VaVm_instar_i_plot + labs(x = "dorsal position", y = expression(log~"("~V[g]~"/"~V[m]~")")) + scale_x_continuous(limits=c(10,600), breaks=c(0,300,600)) 




########
########


R2_all_plot + scale_color_manual(values=c("#000000","#FF0000","#D3D3D3","#fed4d2")) + labs(x = "phenotypic trait", y = "correlation coefficient") 

VaVm_instar_i_plot + labs(x = "dorsal position", y = expression(log~"("~V[g]~"/"~V[m]~")")) + scale_x_continuous(limits=c(10,600), breaks=c(0,300,600)) 


patchwork_plots_induction4 <- R2_all_plot + scale_color_manual(values=c("#000000","#FF0000","#D3D3D3","#fed4d2")) + labs(x = "phenotypic trait", y = "correlation coefficient") +
  VaVm_instar_i_plot + labs(x = "dorsal position", y = expression(log~"("~V[g]~"/"~V[m]~")")) + scale_x_continuous(limits=c(10,600), breaks=c(0,300,600)) +
  plot_spacer() +
  plot_layout(ncol=3, widths = c(1,1,0.5))

patchwork_plots_induction4














