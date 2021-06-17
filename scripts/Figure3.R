# Becker et al 2020 (UVA) - FIG 4

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
        
 
## within clutch - i.e. twins born to the same mum
# reduce to one i-th position for animal length, max height, and eye area
all_data_use <- all_data_final_mod[i == 150]

# modify here re pheno trait (i.e., length, max_height_new)
all_data_wide <- as.data.table(dcast(all_data_use, 
                                     SC_group + cloneid_geno + instar + treatment + deme ~ paste0('clutch', clutch), 
                                     fun = mean, 
                                     value.var = "max_height_new"))  # change here for max_height_new / length

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
all_data_use <- all_data_final_mod[i == 150]

# modify here re pheno trait (i.e., length, max_height_new)
all_data_wide <- as.data.table(dcast(all_data_use, 
                                     SC_group + cloneid_geno + instar + treatment + clutch ~ paste0('deme', deme), 
                                     fun = mean, 
                                     value.var = "max_height_new"))  # change here for max_height_new / length

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
all_data_use <- all_data_final_mod[i == 150]

# modify here re pheno trait (i.e., length, max_height_new)
all_data_wide <- as.data.table(dcast(all_data_use, 
                                     SC_group + cloneid_geno + instar + treatment ~ paste0('SC_', SC_group), 
                                     fun = mean, 
                                     value.var = "max_height_new"))  # change here for max_height_new / length

setkey(all_data_wide, SC_group, instar, treatment)

clone.tab <- all_data_wide[,list(clone=unique(cloneid_geno)), list(SC_group)]

# reduce data set to A clones only
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
save(o_out, o.ag, o.ag.ag, file = 'SC_AvsSC_A_maxHeight.RData')




## stats

# within clutch

# length
load(file = 'output/AvsB_length.RData')

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

# Z = -27.32505, p < 2.2e-16

# I2 predation
clutch_length_I2_C <- wilcox.test(clutch_length[instar == 2][treatment == 0.5][perm > 0]$cor, mu=clutch_length[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_length_I2_C
clutch_Zscore_length_I2_C = qnorm(clutch_length_I2_C$p.value/2)
clutch_Zscore_length_I2_C

# Z = -27.39292, p < 2.2e-16


# max height
load(file = 'output/AvsB_maxHeight.RData')

clutch_maxHeight <- o_out[SC_group == 'A'][, -'SC_group']

# I1 control
clutch_maxHeight_I1_C <- wilcox.test(clutch_maxHeight[instar == 1][treatment == 0][perm > 0]$cor, mu=clutch_maxHeight[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clutch_maxHeight_I1_C
clutch_Zscore_maxHeight_I1_C = qnorm(clutch_maxHeight_I1_C$p.value/2)
clutch_Zscore_maxHeight_I1_C

# Z = -25.78852, p < 2.2e-16

# I1 predation
clutch_maxHeight_I1_P <- wilcox.test(clutch_maxHeight[instar == 1][treatment == 0.5][perm > 0]$cor, mu=clutch_maxHeight[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_maxHeight_I1_P
clutch_Zscore_maxHeight_I1_P = qnorm(clutch_maxHeight_I1_P$p.value/2)
clutch_Zscore_maxHeight_I1_P

# Z = -27.39292, p < 2.2e-16

# I2 control
clutch_maxHeight_I2_C <- wilcox.test(clutch_maxHeight[instar == 2][treatment == 0][perm > 0]$cor, mu=clutch_maxHeight[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clutch_maxHeight_I2_C
clutch_Zscore_maxHeight_I2_C = qnorm(clutch_maxHeight_I2_C$p.value/2)
clutch_Zscore_maxHeight_I2_C

# Z = -21.86441, p < 2.2e-16

# I2 predation
clutch_maxHeight_I2_C <- wilcox.test(clutch_maxHeight[instar == 2][treatment == 0.5][perm > 0]$cor, mu=clutch_maxHeight[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_maxHeight_I2_C
clutch_Zscore_maxHeight_I2_C = qnorm(clutch_maxHeight_I2_C$p.value/2)
clutch_Zscore_maxHeight_I2_C

# Z = -27.35843, p < 2.2e-16


# within clone

# length

load(file = 'output/d1vsd2_length.RData')

clone_length <- o_out[SC_group == 'A'][, -'SC_group']

# I1 control
clone_length_I1_C <- wilcox.test(clone_length[instar == 1][treatment == 0][perm > 0]$cor, mu=clone_length[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_length_I1_C
clone_Zscore_length_I1_C = qnorm(clone_length_I1_C$p.value/2)
clone_Zscore_length_I1_C

# Z = -26.87691, p < 2.2e-16

# I1 predation
clone_length_I1_P <- wilcox.test(clone_length[instar == 1][treatment == 0.5][perm > 0]$cor, mu=clone_length[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_length_I1_P
clone_Zscore_length_I1_P = qnorm(clone_length_I1_P$p.value/2)
clone_Zscore_length_I1_P

# Z = -27.38482, p < 2.2e-16

# I2 control
clone_length_I2_C <- wilcox.test(clone_length[instar == 2][treatment == 0][perm > 0]$cor, mu=clone_length[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_length_I2_C
clone_Zscore_length_I2_C = qnorm(clone_length_I2_C$p.value/2)
clone_Zscore_length_I2_C

# Z = -25.83406, p < 2.2e-16

# I2 predation
clone_length_I2_C <- wilcox.test(clone_length[instar == 2][treatment == 0.5][perm > 0]$cor, mu=clone_length[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_length_I2_C
clone_Zscore_length_I2_C = qnorm(clone_length_I2_C$p.value/2)
clone_Zscore_length_I2_C

# Z = -23.86439, p < 2.2e-16


# max height
load(file = 'output/d1vsd2_maxHeight.RData')

clone_maxHeight <- o_out[SC_group == 'A'][, -'SC_group']

# I1 control
clone_maxHeight_I1_C <- wilcox.test(clone_maxHeight[instar == 1][treatment == 0][perm > 0]$cor, mu=clone_maxHeight[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_maxHeight_I1_C
clone_Zscore_maxHeight_I1_C = qnorm(clone_maxHeight_I1_C$p.value/2)
clone_Zscore_maxHeight_I1_C

# Z = -10.8457, p < 2.2e-16

# I1 predation
clone_maxHeight_I1_P <- wilcox.test(clone_maxHeight[instar == 1][treatment == 0.5][perm > 0]$cor, mu=clone_maxHeight[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_maxHeight_I1_P
clone_Zscore_maxHeight_I1_P = qnorm(clone_maxHeight_I1_P$p.value/2)
clone_Zscore_maxHeight_I1_P

# Z = -27.35712, p < 2.2e-16

# I2 control
clone_maxHeight_I2_C <- wilcox.test(clone_maxHeight[instar == 2][treatment == 0][perm > 0]$cor, mu=clone_maxHeight[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_maxHeight_I2_C
clone_Zscore_maxHeight_I2_C = qnorm(clone_maxHeight_I2_C$p.value/2)
clone_Zscore_maxHeight_I2_C

# Z = -25.66034, p < 2.2e-16

# I2 predation
clone_maxHeight_I2_C <- wilcox.test(clone_maxHeight[instar == 2][treatment == 0.5][perm > 0]$cor, mu=clone_maxHeight[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_maxHeight_I2_C
clone_Zscore_maxHeight_I2_C = qnorm(clone_maxHeight_I2_C$p.value/2)
clone_Zscore_maxHeight_I2_C

# Z = -27.39292, p < 2.2e-16



# A vs A (random draw)

# length

load(file = 'output/SC_AvsSC_A_length.RData')

AvsA_length <- o_out
AvsA_length[, permuted := ifelse(AvsA_length$perm == 0, "no", "yes")]

# I1 control
AvsA_length_I1_C <- wilcox.test(cor ~ permuted, data = AvsA_length[instar == 1][treatment == 0]) 
AvsA_length_I1_C
AvsA_Zscore_length_I1_C = qnorm(AvsA_length_I1_C$p.value/2)
AvsA_Zscore_length_I1_C

# Z = -0.7049809, p = 0.4808

# I1 predation
AvsA_length_I1_P <- wilcox.test(cor ~ permuted, data = AvsA_length[instar == 1][treatment == 0.5]) 
AvsA_length_I1_P
AvsA_Zscore_length_I1_P = qnorm(AvsA_length_I1_P$p.value/2)
AvsA_Zscore_length_I1_P

# Z = -0.7697712, p = 0.4414

# I2 control
AvsA_length_I2_C <- wilcox.test(cor ~ permuted, data = AvsA_length[instar == 2][treatment == 0]) 
AvsA_length_I2_C
AvsA_Zscore_length_I2_C = qnorm(AvsA_length_I2_C$p.value/2)
AvsA_Zscore_length_I2_C

# Z = -0.06084859, p = 0.9515

# I2 predation
AvsA_length_I2_P <- wilcox.test(cor ~ permuted, data = AvsA_length[instar == 2][treatment == 0.5]) 
AvsA_length_I2_P
AvsA_Zscore_length_I2_P = qnorm(AvsA_length_I2_P$p.value/2)
AvsA_Zscore_length_I2_P

# Z = -0.8034002, p = 0.4217


# max height
load(file = 'output/SC_AvsSC_A_maxHeight.RData')

AvsA_maxHeight <- o_out
AvsA_maxHeight[, permuted := ifelse(AvsA_maxHeight$perm == 0, "no", "yes")]


# I1 control
AvsA_maxHeight_I1_C <- wilcox.test(cor ~ permuted, data = AvsA_maxHeight[instar == 1][treatment == 0]) 
AvsA_maxHeight_I1_C
AvsA_Zscore_maxHeight_I1_C = qnorm(AvsA_maxHeight_I1_C$p.value/2)
AvsA_Zscore_maxHeight_I1_C

# Z = -0.2043841, p = 0.8381

# I1 predation
AvsA_maxHeight_I1_P <- wilcox.test(cor ~ permuted, data = AvsA_maxHeight[instar == 1][treatment == 0.5]) 
AvsA_maxHeight_I1_P
AvsA_Zscore_maxHeight_I1_P = qnorm(AvsA_maxHeight_I1_P$p.value/2)
AvsA_Zscore_maxHeight_I1_P

# Z = -0.1316419, p = 0.8953

# I2 control
AvsA_maxHeight_I2_C <- wilcox.test(cor ~ permuted, data = AvsA_maxHeight[instar == 2][treatment == 0]) 
AvsA_maxHeight_I2_C
AvsA_Zscore_maxHeight_I2_C = qnorm(AvsA_maxHeight_I2_C$p.value/2)
AvsA_Zscore_maxHeight_I2_C

# Z = -0.6432461, p = 0.5201

# I2 predation
AvsA_maxHeight_I2_P <- wilcox.test(cor ~ permuted, data = AvsA_maxHeight[instar == 2][treatment == 0.5]) 
AvsA_maxHeight_I2_P
AvsA_Zscore_maxHeight_I2_P = qnorm(AvsA_maxHeight_I2_P$p.value/2)
AvsA_Zscore_maxHeight_I2_P

# Z = -1.901035, p = 0.0573



# within clutch vs among clones

# length

# I1 control
clutch_vs_As_length_I1_C <- wilcox.test(AvsA_length[instar == 1][treatment == 0][perm == 0]$cor, mu=clutch_length[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_length_I1_C
clutch_vs_As_Zscore_length_I1_C = qnorm(clutch_vs_As_length_I1_C$p.value/2)
clutch_vs_As_Zscore_length_I1_C

# Z = -6.100872, p = 1.055e-09

# I1 predation
clutch_vs_As_length_I1_P <- wilcox.test(AvsA_length[instar == 1][treatment == 0.5][perm == 0]$cor, mu=clutch_length[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_length_I1_P
clutch_vs_As_Zscore_length_I1_P = qnorm(clutch_vs_As_length_I1_P$p.value/2)
clutch_vs_As_Zscore_length_I1_P

# Z = -6.149139, p = 7.79e-10

# I2 control
clutch_vs_As_length_I2_C <- wilcox.test(AvsA_length[instar == 2][treatment == 0][perm == 0]$cor, mu=clutch_length[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_length_I2_C
clutch_vs_As_Zscore_length_I2_C = qnorm(clutch_vs_As_length_I2_C$p.value/2)
clutch_vs_As_Zscore_length_I2_C

# Z = -5.985033, p = 2.163e-09

# I2 predation
clutch_vs_As_length_I2_P <- wilcox.test(AvsA_length[instar == 2][treatment == 0.5][perm == 0]$cor, mu=clutch_length[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_length_I2_P
clutch_vs_As_Zscore_length_I2_P = qnorm(clutch_vs_As_length_I2_P$p.value/2)
clutch_vs_As_Zscore_length_I2_P

# Z = -6.120179, p = 9.347e-10


# max height

# I1 control
clutch_vs_As_maxHeight_I1_C <- wilcox.test(AvsA_maxHeight[instar == 1][treatment == 0][perm == 0]$cor, mu=clutch_maxHeight[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_maxHeight_I1_C
clutch_vs_As_Zscore_maxHeight_I1_C = qnorm(clutch_vs_As_maxHeight_I1_C$p.value/2)
clutch_vs_As_Zscore_maxHeight_I1_C

# Z = -3.61998, p = 0.0002946

# I1 predation
clutch_vs_As_maxHeight_I1_P <- wilcox.test(AvsA_maxHeight[instar == 1][treatment == 0.5][perm == 0]$cor, mu=clutch_maxHeight[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_maxHeight_I1_P
clutch_vs_As_Zscore_maxHeight_I1_P = qnorm(clutch_vs_As_maxHeight_I1_P$p.value/2)
clutch_vs_As_Zscore_maxHeight_I1_P

# Z = -5.85954, p = 4.642e-09

# I2 control
clutch_vs_As_maxHeight_I2_C <- wilcox.test(AvsA_maxHeight[instar == 2][treatment == 0][perm == 0]$cor, mu=clutch_maxHeight[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_maxHeight_I2_C
clutch_vs_As_Zscore_maxHeight_I2_C = qnorm(clutch_vs_As_maxHeight_I2_C$p.value/2)
clutch_vs_As_Zscore_maxHeight_I2_C

# Z = -4.247443, p = 2.162e-05

# I2 predation
clutch_vs_As_maxHeight_I2_P <- wilcox.test(AvsA_maxHeight[instar == 2][treatment == 0.5][perm == 0]$cor, mu=clutch_maxHeight[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_maxHeight_I2_P
clutch_vs_As_Zscore_maxHeight_I2_P = qnorm(clutch_vs_As_maxHeight_I2_P$p.value/2)
clutch_vs_As_Zscore_maxHeight_I2_P

# Z = -5.975379, p = 2.296e-09



# within clone vs among clones

# length

# I1 control
clone_vs_As_length_I1_C <- wilcox.test(AvsA_length[instar == 1][treatment == 0][perm == 0]$cor, mu=clone_length[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_length_I1_C
clone_vs_As_Zscore_length_I1_C = qnorm(clone_vs_As_length_I1_C$p.value/2)
clone_vs_As_Zscore_length_I1_C

# Z = -5.483062, p = 4.18e-08

# I1 predation
clone_vs_As_length_I1_P <- wilcox.test(AvsA_length[instar == 1][treatment == 0.5][perm == 0]$cor, mu=clone_length[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_length_I1_P
clone_vs_As_Zscore_length_I1_P = qnorm(clone_vs_As_length_I1_P$p.value/2)
clone_vs_As_Zscore_length_I1_P

# Z = -6.091219, p = 1.121e-09

# I2 control
clone_vs_As_length_I2_C <- wilcox.test(AvsA_length[instar == 2][treatment == 0][perm == 0]$cor, mu=clone_length[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_length_I2_C
clone_vs_As_Zscore_length_I2_C = qnorm(clone_vs_As_length_I2_C$p.value/2)
clone_vs_As_Zscore_length_I2_C

# Z = -5.251384, p = 1.51e-07

# I2 predation
clone_vs_As_length_I2_P <- wilcox.test(AvsA_length[instar == 2][treatment == 0.5][perm == 0]$cor, mu=clone_length[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_length_I2_P
clone_vs_As_Zscore_length_I2_P = qnorm(clone_vs_As_length_I2_P$p.value/2)
clone_vs_As_Zscore_length_I2_P

# Z = -5.125891, p = 2.961e-07


# max height

# I1 control
clone_vs_As_maxHeight_I1_C <- wilcox.test(AvsA_maxHeight[instar == 1][treatment == 0][perm == 0]$cor, mu=clone_maxHeight[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_maxHeight_I1_C
clone_vs_As_Zscore_maxHeight_I1_C = qnorm(clone_vs_As_maxHeight_I1_C$p.value/2)
clone_vs_As_Zscore_maxHeight_I1_C

# Z = -1.09082, p = 0.2754

# I1 predation
clone_vs_As_maxHeight_I1_P <- wilcox.test(AvsA_maxHeight[instar == 1][treatment == 0.5][perm == 0]$cor, mu=clone_maxHeight[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_maxHeight_I1_P
clone_vs_As_Zscore_maxHeight_I1_P = qnorm(clone_vs_As_maxHeight_I1_P$p.value/2)
clone_vs_As_Zscore_maxHeight_I1_P

# Z = -5.569942, p = 2.548e-08

# I2 control
clone_vs_As_maxHeight_I2_C <- wilcox.test(AvsA_maxHeight[instar == 2][treatment == 0][perm == 0]$cor, mu=clone_maxHeight[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_maxHeight_I2_C
clone_vs_As_Zscore_maxHeight_I2_C = qnorm(clone_vs_As_maxHeight_I2_C$p.value/2)
clone_vs_As_Zscore_maxHeight_I2_C

# Z = -5.21277, p = 1.86e-07

# I2 predation
clone_vs_As_maxHeight_I2_P <- wilcox.test(AvsA_maxHeight[instar == 2][treatment == 0.5][perm == 0]$cor, mu=clone_maxHeight[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_maxHeight_I2_P
clone_vs_As_Zscore_maxHeight_I2_P = qnorm(clone_vs_As_maxHeight_I2_P$p.value/2)
clone_vs_As_Zscore_maxHeight_I2_P

# Z = -6.139485, p = 8.279e-10




## plots

# within clutch

load(file = 'output/AvsB_length.RData')

clutch_length <- o.ag[SC_group == 'A'][, -'SC_group']
clutch_length[, group := "within clutch"]
clutch_length[, set := "length"]


load(file = 'output/AvsB_maxHeight.RData')

clutch_maxHeight <- o.ag[SC_group == 'A'][, -'SC_group']
clutch_maxHeight[, group := "within clutch"]
clutch_maxHeight[, set := "max dorsal height"]


# within clone

load(file = 'output/d1vsd2_length.RData')

clone_length <- o.ag[SC_group == 'A'][, -'SC_group']
clone_length[, group := "within clone"]
clone_length[, set := "length"]


load(file = 'output/d1vsd2_maxHeight.RData')

clone_maxHeight <- o.ag[SC_group == 'A'][, -'SC_group']
clone_maxHeight[, group := "within clone"]
clone_maxHeight[, set := "max dorsal height"]


# A vs A (random draw)

load(file = 'output/SC_AvsSC_A_length.RData')

AvsA_length <- o.ag
AvsA_length[, group := "among clones"]
AvsA_length[, set := "length"]

AvsA_length.ag <- o.ag.ag
AvsA_length.ag[, group := "among clones"]
AvsA_length.ag[, set := "length"]


load(file = 'output/SC_AvsSC_A_maxHeight.RData')

AvsA_maxHeight <- o.ag
AvsA_maxHeight[, group := "among clones"]
AvsA_maxHeight[, set := "max dorsal height"]

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
                      
                      ylim(-0.35, 0.6) +

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

load(file = "output/H2_0_I1_ratio_batch.RData")
H2_ctrl_I1_ratio <- as.data.table(collected2)
H2_ctrl_I1_ratio[, i:= as.numeric(i)]
H2_ctrl_I1_ratio[, instar:= "instar 1"]
H2_ctrl_I1_ratio[, group:= "ctrl"]
setkey(H2_ctrl_I1_ratio, group,instar,i)

load(file = "output/H2_0_I2_ratio_batch.RData")
H2_ctrl_I2_ratio <- as.data.table(collected2)
H2_ctrl_I2_ratio[, i:= as.numeric(i)]
H2_ctrl_I2_ratio[, instar:= "instar 2"]
H2_ctrl_I2_ratio[, group:= "ctrl"]
setkey(H2_ctrl_I2_ratio, group,instar,i)

load(file = "output/H2_05_I1_ratio_batch.RData")
H2_trt_I1_ratio <- as.data.table(collected2)
H2_trt_I1_ratio[, i:= as.numeric(i)]
H2_trt_I1_ratio[, instar:= "instar 1"]
H2_trt_I1_ratio[, group:= "trt"]
setkey(H2_trt_I1_ratio, group,instar,i)

load(file = "output/H2_05_I2_ratio_batch.RData")
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
  
                        geom_rect(aes(xmin = 10, xmax = 100, ymin = -Inf, ymax = Inf),fill = "#F5F5F5", colour="#F5F5F5", alpha=0.1) + 
                        geom_rect(aes(xmin = 101, xmax = 250, ymin = -Inf, ymax = Inf),fill = "#E8E8E8", colour="#E8E8E8", alpha=0.1) + 
                        geom_rect(aes(xmin = 251, xmax = 600, ymin = -Inf, ymax = Inf),fill = "#DCDCDC", colour="#DCDCDC", alpha=0.1) + 
  
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

VaVm_instar_i_plot + labs(x = "dorsal position", y = expression(log[10]~"("~V[g]~"/"~V[m]~")")) + scale_x_continuous(limits=c(10,600), breaks=c(0,300,600)) 




########
########


R2_all_plot + scale_color_manual(values=c("#000000","#FF0000","#D3D3D3","#fed4d2")) + labs(x = "phenotypic trait", y = "correlation coefficient") 

VaVm_instar_i_plot + labs(x = "dorsal position", y = expression(log~"("~V[g]~"/"~V[m]~")")) + scale_x_continuous(limits=c(10,600), breaks=c(0,300,600)) 


patchwork_plots_induction4 <- R2_all_plot + scale_color_manual(values=c("#000000","#FF0000","#D3D3D3","#fed4d2")) + labs(x = "phenotypic trait", y = "correlation coefficient") +
  VaVm_instar_i_plot + labs(x = "dorsal position", y = expression(log[10]~"("~V[g]~"/"~V[m]~")")) + scale_x_continuous(limits=c(10,600), breaks=c(0,300,600)) +
  plot_spacer() +
  plot_layout(ncol=3, widths = c(1,1,0.5))

patchwork_plots_induction4


