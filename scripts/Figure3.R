# Becker et al - FIG 3


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

model_length <- lmer(length ~ treatment + (1|clutch) + (1|deme) + (1|cloneid_geno) + (1|batch_new), data = all_data_use)
VarianceComponents <- VarCorr(model_length)

model_maxHeight <- lmer(max_height_new ~ treatment + (1|clutch) + (1|deme) + (1|cloneid_geno) + (1|batch_new), data = all_data_use)
VarianceComponents <- VarCorr(model_maxHeight)


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
          
          perm.m[,list(cor=robcor(deme1, deme2), perm=p), list(SC_group, instar, treatment)]
          
        }


o_out <- rbindlist(o)

o.ag <- o_out[,list(cor.mu=median(cor, na.rm=T), cor.lci=quantile(cor, .05, na.rm=T), cor.uci=quantile(cor, .95, na.rm=T)), 
              list(SC_group, instar, treatment, perm=(perm>0))]

# save(o_out, o.ag, file = 'd1vsd2_length.RData')
# save(o_out, o.ag, file = 'd1vsd2_maxHeight.RData')



## across As - i.e. randomly draw offspring among As
# reduce to one i-th position for animal length, max height, and eye area
# reduce to A clones only
all_data_use <- all_data_final_mod[i == 150][SC_group == "A"]

# modify here re pheno trait (i.e., length, max_height_new)
all_data_wide <- as.data.table(dcast(all_data_use, 
                                     SC_group + cloneid_geno + instar + treatment ~ paste0('SC_', SC_group), 
                                     fun = mean, 
                                     value.var = "max_height_new"))  # change here for max_height_new / length

setkey(all_data_wide, SC_group, instar, treatment)

clone.tab <- all_data_wide[,list(clone=unique(cloneid_geno)), list(SC_group)]


o <- foreach(b=1:100)%do%{  ## increase to 100
          print(b)
          
          foreach(p=0:1000)%dopar%{  ## increase to 1000
            print(p)
            
              foreach(instar.i = unique(all_data_use$instar), .errorhandling="remove", .combine="rbind") %do% {
                
                foreach(treatment.i = unique(all_data_use$treatment), .errorhandling="remove", .combine="rbind") %do% {
                  
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
                    d1.tmp.perm  <- d1.tmp[J(clone.perm.d1)]
                    setnames(d1.tmp.perm, c("cloneid_geno", "cloneid_geno_new"), c("cloneid_geno_orig", "cloneid_geno"))
                    
                    d2.tmp.perm  <- d2.tmp[J(clone.perm.d2)]
                    setnames(d2.tmp.perm, c("cloneid_geno", "cloneid_geno_new"), c("cloneid_geno_orig", "cloneid_geno"))
                    
                  }
                  
                  setkey(d1.tmp.perm, instar, treatment)
                  setkey(d2.tmp.perm, instar, treatment)
                  
                  # only sample one of the two sets
                  d2.tmp.perm_draw <- d2.tmp.perm[sample(.N, nrow(d1.tmp.perm))]
                  d2.tmp.perm_draw[, SC_A_draw := SC_A]
                  d2.tmp.perm_draw[, cloneid_geno_draw := cloneid_geno]
                  
                  perm.m <- cbind(d1.tmp.perm,d2.tmp.perm_draw[,-"SC_A"])
                  # table(perm.m[,cloneid_geno] != perm.m[,cloneid_geno_draw]) 
                  
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



## As across batches - i.e. randomly paired A's within batch 
# reduce to one i-th position for animal length, max height, and eye area
# reduce to A clones only

all_data_use <- all_data_final_mod[i == 150][SC_group == "A"]

# modify here re pheno trait (i.e., length, max_height_new)
all_data_wide <- as.data.table(dcast(all_data_use, 
                                     SC_group + batch + cloneid_geno + instar + treatment ~ paste0('SC_', SC_group), 
                                     fun = mean, 
                                     value.var = "length"))  # change here for max_height_new / length

setkey(all_data_wide, SC_group, batch, instar, treatment)

clone.tab <- all_data_wide[,list(clone=unique(cloneid_geno)), list(SC_group, batch)]


o <- foreach(b=1:2)%do%{  ## increase to 100
        print(b)
        
        foreach(p=0:2)%dopar%{  ## increase to 1000
          print(p)
          
          foreach(instar.i = unique(all_data_use$instar), .errorhandling="remove", .combine="rbind") %do% {
            
            foreach(treatment.i = unique(all_data_use$treatment), .errorhandling="remove", .combine="rbind") %do% {

              foreach(batch.i = unique(all_data_use$batch), .errorhandling="remove", .combine="rbind") %do% {
                
              clone.perm.d1 <- clone.tab[,list(cloneid_geno_new=sample(clone), cloneid_geno=clone), list(SC_group, batch)]
              clone.perm.d2 <- clone.tab[,list(cloneid_geno_new=sample(clone), cloneid_geno=clone), list(SC_group, batch)]
              
              d1.tmp <- all_data_wide[,c("SC_group", "batch", "instar", "treatment", "SC_A", "cloneid_geno"), with=F][instar == instar.i][treatment == treatment.i][batch == batch.i]
              d2.tmp <- all_data_wide[,c("SC_group", "batch", "instar", "treatment", "SC_A", "cloneid_geno"), with=F][instar == instar.i][treatment == treatment.i][batch == batch.i]
              
              setkey(d1.tmp, SC_group, cloneid_geno)
              setkey(d2.tmp, SC_group, cloneid_geno)
              
              setkey(clone.perm.d1, SC_group, cloneid_geno)
              setkey(clone.perm.d2, SC_group, cloneid_geno)
              
              if(p==0) {
                d1.tmp.perm <- d1.tmp
                d2.tmp.perm <- d2.tmp
                
              }else if(p>0) {
                d1.tmp.perm  <- d1.tmp[J(clone.perm.d1[batch == batch.i][cloneid_geno %in% d1.tmp$cloneid_geno])]
                setnames(d1.tmp.perm, c("cloneid_geno", "cloneid_geno_new"), c("cloneid_geno_orig", "cloneid_geno"))
                
                d2.tmp.perm  <- d2.tmp[J(clone.perm.d2[batch == batch.i][cloneid_geno %in% d1.tmp$cloneid_geno])]
                setnames(d2.tmp.perm, c("cloneid_geno", "cloneid_geno_new"), c("cloneid_geno_orig", "cloneid_geno"))
                
              }
              
              setkey(d1.tmp.perm, instar, treatment)
              setkey(d2.tmp.perm, instar, treatment)
              
              # only sample one of the two sets
              d2.tmp.perm_draw <- d2.tmp.perm[sample(.N, nrow(d1.tmp.perm))]
              d2.tmp.perm_draw[, SC_A_draw := SC_A]
              d2.tmp.perm_draw[, cloneid_geno_draw := cloneid_geno]
              
              perm.m <- cbind(d1.tmp.perm,d2.tmp.perm_draw[,-"SC_A"])
              # table(perm.m[,cloneid_geno] != perm.m[,cloneid_geno_draw])
              
              cor_out <- cor(perm.m$SC_A, perm.m$SC_A_draw)

              ## data out
              data.table(batch = batch.i, 
                         instar = instar.i,
                         treatment = treatment.i,
                         perm = p,
                         boot = b, 
                         cor = cor_out)
              
              
            }}}}}
      

o_out <- rbindlist(unlist(o, recursive = FALSE))

o.ag <- o_out[,list(cor.mu=median(cor, na.rm=T), cor.lci=quantile(cor, .05, na.rm=T), cor.uci=quantile(cor, .95, na.rm=T)), 
              list(instar, treatment, batch, boot, perm=(perm>0))]

o.ag.ag <- o_out[,list(cor.mu=median(cor, na.rm=T), cor.lci=quantile(cor, .05, na.rm=T), cor.uci=quantile(cor, .95, na.rm=T)), 
                 list(instar, treatment, batch, perm=(perm>0))]

o.ag.ag.ag <- o_out[,list(cor.mu=median(cor, na.rm=T), cor.lci=quantile(cor, .05, na.rm=T), cor.uci=quantile(cor, .95, na.rm=T)), 
                 list(instar, treatment, perm=(perm>0))]

# save(o_out, o.ag, o.ag.ag, o.ag.ag.ag, file = 'SC_AvsSC_A_batch_length_RED.RData')
# save(o_out, o.ag, o.ag.ag, o.ag.ag.ag, file = 'SC_AvsSC_A_batch_maxHeight_RED.RData')



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

# Z = -0.1025878, p = 0.9183

# I1 predation
AvsA_length_I1_P <- wilcox.test(cor ~ permuted, data = AvsA_length[instar == 1][treatment == 0.5]) 
AvsA_length_I1_P
AvsA_Zscore_length_I1_P = qnorm(AvsA_length_I1_P$p.value/2)
AvsA_Zscore_length_I1_P

# Z = -0.7975722, p = 0.4251

# I2 control
AvsA_length_I2_C <- wilcox.test(cor ~ permuted, data = AvsA_length[instar == 2][treatment == 0]) 
AvsA_length_I2_C
AvsA_Zscore_length_I2_C = qnorm(AvsA_length_I2_C$p.value/2)
AvsA_Zscore_length_I2_C

# Z = -1.06435, p = 0.2872

# I2 predation
AvsA_length_I2_P <- wilcox.test(cor ~ permuted, data = AvsA_length[instar == 2][treatment == 0.5]) 
AvsA_length_I2_P
AvsA_Zscore_length_I2_P = qnorm(AvsA_length_I2_P$p.value/2)
AvsA_Zscore_length_I2_P

# Z = -0.8084267, p = 0.4188


# max height
load(file = 'output/SC_AvsSC_A_maxHeight.RData')

AvsA_maxHeight <- o_out
AvsA_maxHeight[, permuted := ifelse(AvsA_maxHeight$perm == 0, "no", "yes")]


# I1 control
AvsA_maxHeight_I1_C <- wilcox.test(cor ~ permuted, data = AvsA_maxHeight[instar == 1][treatment == 0]) 
AvsA_maxHeight_I1_C
AvsA_Zscore_maxHeight_I1_C = qnorm(AvsA_maxHeight_I1_C$p.value/2)
AvsA_Zscore_maxHeight_I1_C

# Z = -0.1169427, p = 0.9069

# I1 predation
AvsA_maxHeight_I1_P <- wilcox.test(cor ~ permuted, data = AvsA_maxHeight[instar == 1][treatment == 0.5]) 
AvsA_maxHeight_I1_P
AvsA_Zscore_maxHeight_I1_P = qnorm(AvsA_maxHeight_I1_P$p.value/2)
AvsA_Zscore_maxHeight_I1_P

# Z = -1.161296, p = 0.2455

# I2 control
AvsA_maxHeight_I2_C <- wilcox.test(cor ~ permuted, data = AvsA_maxHeight[instar == 2][treatment == 0]) 
AvsA_maxHeight_I2_C
AvsA_Zscore_maxHeight_I2_C = qnorm(AvsA_maxHeight_I2_C$p.value/2)
AvsA_Zscore_maxHeight_I2_C

# Z = -0.4066517, p = 0.6843

# I2 predation
AvsA_maxHeight_I2_P <- wilcox.test(cor ~ permuted, data = AvsA_maxHeight[instar == 2][treatment == 0.5]) 
AvsA_maxHeight_I2_P
AvsA_Zscore_maxHeight_I2_P = qnorm(AvsA_maxHeight_I2_P$p.value/2)
AvsA_Zscore_maxHeight_I2_P

# Z = -1.811072, p = 0.07013


# A vs A BATCH (random draw)

# length

load(file = 'output/SC_AvsSC_A_batch_length.RData')

AvsA_batch_length <- o_out
AvsA_batch_length[, permuted := ifelse(AvsA_batch_length$perm == 0, "no", "yes")]

# I1 control
AvsA_batch_length_I1_C <- wilcox.test(cor ~ permuted, data = AvsA_batch_length[instar == 1][treatment == 0]) 
AvsA_batch_length_I1_C
AvsA_batch_Zscore_length_I1_C = qnorm(AvsA_batch_length_I1_C$p.value/2)
AvsA_batch_Zscore_length_I1_C

# Z = -0.422821, p = 0.6724

# I1 predation
AvsA_batch_length_I1_P <- wilcox.test(cor ~ permuted, data = AvsA_batch_length[instar == 1][treatment == 0.5]) 
AvsA_batch_length_I1_P
AvsA_batch_Zscore_length_I1_P = qnorm(AvsA_batch_length_I1_P$p.value/2)
AvsA_batch_Zscore_length_I1_P

# Z = -1.350417, p = 0.1769

# I2 control
AvsA_batch_length_I2_C <- wilcox.test(cor ~ permuted, data = AvsA_batch_length[instar == 2][treatment == 0]) 
AvsA_batch_length_I2_C
AvsA_batch_Zscore_length_I2_C = qnorm(AvsA_batch_length_I2_C$p.value/2)
AvsA_batch_Zscore_length_I2_C

# Z = -1.004316, p = 0.3152

# I2 predation
AvsA_batch_length_I2_P <- wilcox.test(cor ~ permuted, data = AvsA_batch_length[instar == 2][treatment == 0.5]) 
AvsA_batch_length_I2_P
AvsA_batch_Zscore_length_I2_P = qnorm(AvsA_batch_length_I2_P$p.value/2)
AvsA_batch_Zscore_length_I2_P

# Z = -0.6890701, p = 0.4908


# max height
load(file = 'output/SC_AvsSC_A_batch_maxHeight.RData')

AvsA_batch_maxHeight <- o_out
AvsA_batch_maxHeight[, permuted := ifelse(AvsA_batch_maxHeight$perm == 0, "no", "yes")]


# I1 control
AvsA_batch_maxHeight_I1_C <- wilcox.test(cor ~ permuted, data = AvsA_batch_maxHeight[instar == 1][treatment == 0]) 
AvsA_batch_maxHeight_I1_C
AvsA_batch_Zscore_maxHeight_I1_C = qnorm(AvsA_batch_maxHeight_I1_C$p.value/2)
AvsA_batch_Zscore_maxHeight_I1_C

# Z = -0.6746085, p = 0.4999

# I1 predation
AvsA_batch_maxHeight_I1_P <- wilcox.test(cor ~ permuted, data = AvsA_batch_maxHeight[instar == 1][treatment == 0.5]) 
AvsA_batch_maxHeight_I1_P
AvsA_batch_Zscore_maxHeight_I1_P = qnorm(AvsA_batch_maxHeight_I1_P$p.value/2)
AvsA_batch_Zscore_maxHeight_I1_P

# Z = -1.497639, p = 0.1342

# I2 control
AvsA_batch_maxHeight_I2_C <- wilcox.test(cor ~ permuted, data = AvsA_batch_maxHeight[instar == 2][treatment == 0]) 
AvsA_batch_maxHeight_I2_C
AvsA_batch_Zscore_maxHeight_I2_C = qnorm(AvsA_batch_maxHeight_I2_C$p.value/2)
AvsA_batch_Zscore_maxHeight_I2_C

# Z = -0.1886146, p = 0.8504

# I2 predation
AvsA_batch_maxHeight_I2_P <- wilcox.test(cor ~ permuted, data = AvsA_batch_maxHeight[instar == 2][treatment == 0.5]) 
AvsA_batch_maxHeight_I2_P
AvsA_batch_Zscore_maxHeight_I2_P = qnorm(AvsA_batch_maxHeight_I2_P$p.value/2)
AvsA_batch_Zscore_maxHeight_I2_P

# Z = -2.806003, p = 0.005016



# within clutch vs among clones

# length

# I1 control
clutch_vs_As_length_I1_C <- wilcox.test(AvsA_length[instar == 1][treatment == 0][perm == 0]$cor, mu=clutch_length[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_length_I1_C
clutch_vs_As_Zscore_length_I1_C = qnorm(clutch_vs_As_length_I1_C$p.value/2)
clutch_vs_As_Zscore_length_I1_C

# Z = -8.680051, p < 2.2e-16

# I1 predation
clutch_vs_As_length_I1_P <- wilcox.test(AvsA_length[instar == 1][treatment == 0.5][perm == 0]$cor, mu=clutch_length[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_length_I1_P
clutch_vs_As_Zscore_length_I1_P = qnorm(clutch_vs_As_length_I1_P$p.value/2)
clutch_vs_As_Zscore_length_I1_P

# Z = -8.676613, p < 2.2e-16

# I2 control
clutch_vs_As_length_I2_C <- wilcox.test(AvsA_length[instar == 2][treatment == 0][perm == 0]$cor, mu=clutch_length[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_length_I2_C
clutch_vs_As_Zscore_length_I2_C = qnorm(clutch_vs_As_length_I2_C$p.value/2)
clutch_vs_As_Zscore_length_I2_C

# Z = -8.37404, p < 2.2e-16

# I2 predation
clutch_vs_As_length_I2_P <- wilcox.test(AvsA_length[instar == 2][treatment == 0.5][perm == 0]$cor, mu=clutch_length[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_length_I2_P
clutch_vs_As_Zscore_length_I2_P = qnorm(clutch_vs_As_length_I2_P$p.value/2)
clutch_vs_As_Zscore_length_I2_P

# Z = -8.676613, p < 2.2e-16


# max height

# I1 control
clutch_vs_As_maxHeight_I1_C <- wilcox.test(AvsA_maxHeight[instar == 1][treatment == 0][perm == 0]$cor, mu=clutch_maxHeight[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_maxHeight_I1_C
clutch_vs_As_Zscore_maxHeight_I1_C = qnorm(clutch_vs_As_maxHeight_I1_C$p.value/2)
clutch_vs_As_Zscore_maxHeight_I1_C

# Z = -5.984404, p = 2.172e-09

# I1 predation
clutch_vs_As_maxHeight_I1_P <- wilcox.test(AvsA_maxHeight[instar == 1][treatment == 0.5][perm == 0]$cor, mu=clutch_maxHeight[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_maxHeight_I1_P
clutch_vs_As_Zscore_maxHeight_I1_P = qnorm(clutch_vs_As_maxHeight_I1_P$p.value/2)
clutch_vs_As_Zscore_maxHeight_I1_P

# Z = -8.676613, p < 2.2e-16

# I2 control
clutch_vs_As_maxHeight_I2_C <- wilcox.test(AvsA_maxHeight[instar == 2][treatment == 0][perm == 0]$cor, mu=clutch_maxHeight[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_maxHeight_I2_C
clutch_vs_As_Zscore_maxHeight_I2_C = qnorm(clutch_vs_As_maxHeight_I2_C$p.value/2)
clutch_vs_As_Zscore_maxHeight_I2_C

# Z = -5.929391, p = 3.041e-09

# I2 predation
clutch_vs_As_maxHeight_I2_P <- wilcox.test(AvsA_maxHeight[instar == 2][treatment == 0.5][perm == 0]$cor, mu=clutch_maxHeight[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clutch_vs_As_maxHeight_I2_P
clutch_vs_As_Zscore_maxHeight_I2_P = qnorm(clutch_vs_As_maxHeight_I2_P$p.value/2)
clutch_vs_As_Zscore_maxHeight_I2_P

# Z = -8.494382, p < 2.2e-16



# within clone vs among clones

# length

# I1 control
clone_vs_As_length_I1_C <- wilcox.test(AvsA_length[instar == 1][treatment == 0][perm == 0]$cor, mu=clone_length[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_length_I1_C
clone_vs_As_Zscore_length_I1_C = qnorm(clone_vs_As_length_I1_C$p.value/2)
clone_vs_As_Zscore_length_I1_C

# Z = -8.542518, p < 2.2e-16

# I1 predation
clone_vs_As_length_I1_P <- wilcox.test(AvsA_length[instar == 1][treatment == 0.5][perm == 0]$cor, mu=clone_length[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_length_I1_P
clone_vs_As_Zscore_length_I1_P = qnorm(clone_vs_As_length_I1_P$p.value/2)
clone_vs_As_Zscore_length_I1_P

# Z = -8.583778, p < 2.2e-16

# I2 control
clone_vs_As_length_I2_C <- wilcox.test(AvsA_length[instar == 2][treatment == 0][perm == 0]$cor, mu=clone_length[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_length_I2_C
clone_vs_As_Zscore_length_I2_C = qnorm(clone_vs_As_length_I2_C$p.value/2)
clone_vs_As_Zscore_length_I2_C

# Z = -7.033093, p = 2.02e-12

# I2 predation
clone_vs_As_length_I2_P <- wilcox.test(AvsA_length[instar == 2][treatment == 0.5][perm == 0]$cor, mu=clone_length[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_length_I2_P
clone_vs_As_Zscore_length_I2_P = qnorm(clone_vs_As_length_I2_P$p.value/2)
clone_vs_As_Zscore_length_I2_P

# Z = -7.531651, p = 5.01e-14


# max height

# I1 control
clone_vs_As_maxHeight_I1_C <- wilcox.test(AvsA_maxHeight[instar == 1][treatment == 0][perm == 0]$cor, mu=clone_maxHeight[instar == 1][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_maxHeight_I1_C
clone_vs_As_Zscore_maxHeight_I1_C = qnorm(clone_vs_As_maxHeight_I1_C$p.value/2)
clone_vs_As_Zscore_maxHeight_I1_C

# Z = -1.720882, p = 0.08527

# I1 predation
clone_vs_As_maxHeight_I1_P <- wilcox.test(AvsA_maxHeight[instar == 1][treatment == 0.5][perm == 0]$cor, mu=clone_maxHeight[instar == 1][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_maxHeight_I1_P
clone_vs_As_Zscore_maxHeight_I1_P = qnorm(clone_vs_As_maxHeight_I1_P$p.value/2)
clone_vs_As_Zscore_maxHeight_I1_P

# Z = -8.631915, p < 2.2e-16

# I2 control
clone_vs_As_maxHeight_I2_C <- wilcox.test(AvsA_maxHeight[instar == 2][treatment == 0][perm == 0]$cor, mu=clone_maxHeight[instar == 2][treatment == 0][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_maxHeight_I2_C
clone_vs_As_Zscore_maxHeight_I2_C = qnorm(clone_vs_As_maxHeight_I2_C$p.value/2)
clone_vs_As_Zscore_maxHeight_I2_C

# Z = -7.294406, p = 3e-13

# I2 predation
clone_vs_As_maxHeight_I2_P <- wilcox.test(AvsA_maxHeight[instar == 2][treatment == 0.5][perm == 0]$cor, mu=clone_maxHeight[instar == 2][treatment == 0.5][perm == 0]$cor, alternative = "two.sided")
clone_vs_As_maxHeight_I2_P
clone_vs_As_Zscore_maxHeight_I2_P = qnorm(clone_vs_As_maxHeight_I2_P$p.value/2)
clone_vs_As_Zscore_maxHeight_I2_P

# Z = -8.676613, p < 2.2e-16




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


# A vs A batches (random draw)

load(file = 'output/SC_AvsSC_A_batch_length.RData')

AvsA_batch_length <- o.ag.ag
AvsA_batch_length[, group := "among batches"]
AvsA_batch_length[, set := "length"]

AvsA_batch_length.ag <- o.ag.ag.ag
AvsA_batch_length.ag[, group := "among batches"]
AvsA_batch_length.ag[, set := "length"]


load(file = 'output/SC_AvsSC_A_batch_maxHeight.RData')

AvsA_batch_maxHeight <- o.ag.ag
AvsA_batch_maxHeight[, group := "among batches"]
AvsA_batch_maxHeight[, set := "max dorsal height"]

AvsA_batch_maxHeight.ag <- o.ag.ag.ag
AvsA_batch_maxHeight.ag[, group := "among batches"]
AvsA_batch_maxHeight.ag[, set := "max dorsal height"]


all_R2s <- rbind(clutch_length,
                 clutch_maxHeight,
                 clone_length,
                 clone_maxHeight,
                 #AvsA_length, 
                 AvsA_length.ag,
                 #AvsA_maxHeight, 
                 AvsA_maxHeight.ag, 
                 #AvsA_batch_length, 
                 AvsA_batch_length.ag,
                 #AvsA_batch_maxHeight, 
                 AvsA_batch_maxHeight.ag)

all_R2s[, group := as.factor(group)]
all_R2s[, instar_new := ifelse(all_R2s$instar == 1, "instar 1", "instar 2")]

all_R2s_use <- all_R2s[set %in% c('length', 'max dorsal height')]


R2_all_plot <- all_R2s_use %>% mutate(group = fct_relevel(group, "within clutch", "within clone", "among clones", "among batches")) %>%
  
                  ggplot(aes(x=as.factor(set), y=cor.mu)) +  
                
                      geom_pointrange(aes(ymin = cor.lci, ymax = cor.uci, colour = as.factor(treatment+perm), group = as.factor(treatment+perm)), position = position_dodge(width = 0.6), size = 0.7) +	
                      
                      facet_grid(~instar_new ~ group) +
                      
                      # ylim(-0.35, 0.6) +

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


H2_ratio_all <- rbind(H2_ctrl_I1_ratio_use,
                       H2_ctrl_I2_ratio_use,
                       H2_trt_I1_ratio_use,
                       H2_trt_I2_ratio_use)


VaVm_instar_i_plot <- ggplot(data = H2_ratio_all[i <= 600]) +  
  
                        # geom_rect(data = H2_ration_all[i <= 600][instar == 1], aes(xmin = 10, xmax = 300, ymin = -Inf, ymax = Inf),fill = "#F5F5F5", colour="#F5F5F5", alpha=0.1) + 
                        # geom_rect(data = H2_ration_all[i <= 600][instar == 1], aes(xmin = 301, xmax = 600, ymin = -Inf, ymax = Inf),fill = "#DCDCDC", colour="#DCDCDC", alpha=0.1) + 
                        # 
                        # geom_rect(data = H2_ration_all[i <= 600][instar == 2], aes(xmin = 10, xmax = 100, ymin = -Inf, ymax = Inf),fill = "#F5F5F5", colour="#F5F5F5", alpha=0.1) + 
                        # geom_rect(data = H2_ration_all[i <= 600][instar == 2], aes(xmin = 101, xmax = 200, ymin = -Inf, ymax = Inf),fill = "#E8E8E8", colour="#E8E8E8", alpha=0.1) + 
                        # geom_rect(data = H2_ration_all[i <= 600][instar == 2], aes(xmin = 201, xmax = 600, ymin = -Inf, ymax = Inf),fill = "#DCDCDC", colour="#DCDCDC", alpha=0.1) + 
                        
                        geom_vline(data = H2_ratio_all[i <= 600][instar_new == "instar 1"], aes(xintercept = 10), linetype="dotted", color = "black", size=1) +
                        geom_vline(data = H2_ratio_all[i <= 600][instar_new == "instar 1"], aes(xintercept = 300), linetype="dotted", color = "black", size=1) +
                        geom_vline(data = H2_ratio_all[i <= 600][instar_new == "instar 1"], aes(xintercept = 600), linetype="dotted", color = "black", size=1) +
                        
                        geom_vline(data = H2_ratio_all[i <= 600][instar_new == "instar 2"], aes(xintercept = 10), linetype="dotted", color = "black", size=1) +
                        geom_vline(data = H2_ratio_all[i <= 600][instar_new == "instar 2"], aes(xintercept = 100), linetype="dotted", color = "black", size=1) +
                        geom_vline(data = H2_ratio_all[i <= 600][instar_new == "instar 2"], aes(xintercept = 200), linetype="dotted", color = "black", size=1) +
                        geom_vline(data = H2_ratio_all[i <= 600][instar_new == "instar 2"], aes(xintercept = 600), linetype="dotted", color = "black", size=1) +
                      
                        geom_ribbon(data=H2_ratio_all[i <= 600][group == "trt"], aes(ymin=stuff.lCI, ymax=stuff.uCI, x=i), fill = "#FF0000", alpha = 0.3)+
                        geom_line(data=H2_ratio_all[i <= 600][group == "trt"], aes(x=i, y=stuff.mean), size = 1, colour = "#FF0000") + 
                        
                        geom_ribbon(data=H2_ratio_all[i <= 600][group == "ctrl"], aes(ymin=stuff.lCI, ymax=stuff.uCI, x=i), fill = "#000000", alpha = 0.3)+
                        geom_line(data=H2_ratio_all[i <= 600][group == "ctrl"], aes(x=i, y=stuff.mean), size = 1, colour = "#000000") + 
                        
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



load(file="output/effectSize_heritability.RData")

effectSize_heritability[, module_new := ifelse(effectSize_heritability$module == 'module 1', 'mod 1',
                                               ifelse(effectSize_heritability$module == 'module 2', 'mod 2',
                                                      ifelse(effectSize_heritability$module == 'module 3', 'mod 3','NA')))]

## correlation
cor(effectSize_heritability[instar == "instar 2"][i <= 100][variable == "effect_trt"]$value,
    effectSize_heritability[instar == "instar 2"][i <= 100][variable == "H2_VAoverVM_trt"]$value)

cor(effectSize_heritability[instar == "instar 2"][i > 100 & i <= 200][variable == "effect_trt"]$value,
    effectSize_heritability[instar == "instar 2"][i > 100 & i <= 200][variable == "H2_VAoverVM_trt"]$value)

cor(effectSize_heritability[instar == "instar 2"][i > 200 & i <= 600][variable == "effect_trt"]$value,
    effectSize_heritability[instar == "instar 2"][i > 200 & i <= 600][variable == "H2_VAoverVM_trt"]$value)



# h2Ratio.aov <- aov(value ~ module*variable, data = effectSize_heritability[instar == "instar 2"][variable %in% c("H2_VAoverVM_trt","H2_VAoverVM_ctrl")])
# summary(h2Ratio.aov)
# TukeyHSD(h2Ratio.aov)
# 
h2Ratio_trt.aov <- aov(value ~ module, data = effectSize_heritability[instar == "instar 2"][variable == "H2_VAoverVM_trt"])
summary(h2Ratio_trt.aov)
TukeyHSD(h2Ratio_trt.aov)

plot(h2Ratio_trt.aov, 1)
plot(h2Ratio_trt.aov, 2)

library(car)
leveneTest(value ~ module, data = effectSize_heritability[instar == "instar 2"][variable == "H2_VAoverVM_trt"])
oneway.test(value ~ module, data = effectSize_heritability[instar == "instar 2"][variable == "H2_VAoverVM_trt"])


# ANOVA assumptions not ok, use Kruskal Wallis test
kruskal.test(value ~ module, data = effectSize_heritability[instar == "instar 2"][variable == "H2_VAoverVM_trt"])
pairwise.wilcox.test(effectSize_heritability[instar == "instar 2"][variable == "H2_VAoverVM_trt"]$value,
                     effectSize_heritability[instar == "instar 2"][variable == "H2_VAoverVM_trt"]$module,
                     p.adjust.method = "BH")

# 
# h2Ratio_ctrl.aov <- aov(value ~ module, data = effectSize_heritability[instar == "instar 2"][variable == "H2_VAoverVM_ctrl"])
# summary(h2Ratio_ctrl.aov)
# TukeyHSD(h2Ratio_ctrl.aov)


#####

shapiro.test(H2_ratio_all[i <= 300][instar_new == "instar 1"][group == "ctrl"]$stuff.mean)  # not normal
shapiro.test(H2_ratio_all[i <= 300][instar_new == "instar 1"][group == "trt"]$stuff.mean)   # not normal

shapiro.test(H2_ratio_all[i > 100 & i <= 200][instar_new == "instar 2"][group == "ctrl"]$stuff.mean)  # normal
shapiro.test(H2_ratio_all[i > 100 & i <= 200][instar_new == "instar 2"][group == "trt"]$stuff.mean)   # normal


# O_I1_plast
wilcox_ratio_I1 <- wilcox.test(stuff.mean ~ group, alternative =  "greater", data = H2_ratio_all[i <= 300][instar_new == "instar 1"][group %in% c("ctrl","trt")]) 
wilcox_ratio_I1

Zscore_ratio_I1 = qnorm(wilcox_ratio_I1$p.value/2)
Zscore_ratio_I1

wilcox_effectsize_ratio_I1 <- wilcox_effsize(stuff.mean ~ group, data = H2_ratio_all[i <= 300][instar_new == "instar 1"][group %in% c("ctrl","trt")])
wilcox_effectsize_ratio_I1

# (Z = -20.72637, p-value < 2.2e-16, effect size = 0.858 - large)


# O_I2_plast
wilcox_ratio_I1 <- wilcox.test(stuff.mean ~ group, alternative =  "greater", data = H2_ratio_all[i > 100 & i <= 200][instar_new == "instar 2"][group %in% c("ctrl","trt")]) 
wilcox_ratio_I1

Zscore_ratio_I1 = qnorm(wilcox_ratio_I1$p.value/2)
Zscore_ratio_I1

wilcox_effectsize_ratio_I1 <- wilcox_effsize(stuff.mean ~ group, data = H2_ratio_all[i > 100 & i <= 200][instar_new == "instar 2"][group %in% c("ctrl","trt")])
wilcox_effectsize_ratio_I1

# (Z = -11.34801, p-value < 2.2e-16, effect size = 0.798 - large)


# test where delta is biggest between control and predation conditions

shapiro.test(effectSize_heritability[instar == "instar 2"][variable == "H2_VAoverVM_ctrl"][module_new == "mod 1"]$value)  # normal
shapiro.test(effectSize_heritability[instar == "instar 2"][variable == "H2_VAoverVM_ctrl"][module_new == "mod 2"]$value)  # normal
shapiro.test(effectSize_heritability[instar == "instar 2"][variable == "H2_VAoverVM_ctrl"][module_new == "mod 3"]$value)  # not normal

shapiro.test(effectSize_heritability[instar == "instar 2"][variable == "H2_VAoverVM_trt"][module_new == "mod 1"]$value)  # normal
shapiro.test(effectSize_heritability[instar == "instar 2"][variable == "H2_VAoverVM_trt"][module_new == "mod 2"]$value)  # normal
shapiro.test(effectSize_heritability[instar == "instar 2"][variable == "H2_VAoverVM_trt"][module_new == "mod 3"]$value)  # not normal


# test for difference in modules

effectSize_heritability[, variable_new := as.character(variable)]


t.test(effectSize_heritability[instar == "instar 2"][variable %in% c("H2_VAoverVM_ctrl", "H2_VAoverVM_trt")][module_new == "mod 1"][variable_new == "H2_VAoverVM_ctrl"]$value,
       effectSize_heritability[instar == "instar 2"][variable %in% c("H2_VAoverVM_ctrl", "H2_VAoverVM_trt")][module_new == "mod 1"][variable_new == "H2_VAoverVM_trt"]$value,
       alternative = "greater")

t.test(effectSize_heritability[instar == "instar 2"][variable %in% c("H2_VAoverVM_ctrl", "H2_VAoverVM_trt")][module_new == "mod 2"][variable_new == "H2_VAoverVM_ctrl"]$value,
       effectSize_heritability[instar == "instar 2"][variable %in% c("H2_VAoverVM_ctrl", "H2_VAoverVM_trt")][module_new == "mod 2"][variable_new == "H2_VAoverVM_trt"]$value,
       alternative = "greater")

t.test(effectSize_heritability[instar == "instar 2"][variable %in% c("H2_VAoverVM_ctrl", "H2_VAoverVM_trt")][module_new == "mod 3"][variable_new == "H2_VAoverVM_ctrl"]$value,
       effectSize_heritability[instar == "instar 2"][variable %in% c("H2_VAoverVM_ctrl", "H2_VAoverVM_trt")][module_new == "mod 3"][variable_new == "H2_VAoverVM_trt"]$value,
       alternative = "greater")


diff_mod1 <- wilcox.test(value ~ variable_new, data = effectSize_heritability[instar == "instar 2"][variable %in% c("H2_VAoverVM_ctrl", "H2_VAoverVM_trt")][module_new == "mod 1"], paired = FALSE, alternative = "greater")
diff_mod1

Zscore_mod1 = qnorm(diff_mod1$p.value/2)
Zscore_mod1

effectsize_diff_mod1 <- wilcox_effsize(value ~ variable_new, data = effectSize_heritability[instar == "instar 2"][variable %in% c("H2_VAoverVM_ctrl", "H2_VAoverVM_trt")][module_new == "mod 1"], paired = FALSE, alternative = "greater")
effectsize_diff_mod1

# (Z = -9.887683, p-value < 2.2e-16; large effect size: 0.728)


diff_mod2 <- wilcox.test(value ~ variable, data = effectSize_heritability[instar == "instar 2"][variable %in% c("H2_VAoverVM_ctrl", "H2_VAoverVM_trt")][module_new == "mod 2"], paired = FALSE, alternative = "greater")
diff_mod2

p_diff_mod2 = qnorm(diff_mod2$p.value/2)
p_diff_mod2

effectsize_diff_mod2 <- wilcox_effsize(value ~ variable_new, data = effectSize_heritability[instar == "instar 2"][variable %in% c("H2_VAoverVM_ctrl", "H2_VAoverVM_trt")][module_new == "mod 2"], paired = FALSE, alternative = "greater")
effectsize_diff_mod2

# (Z = -11.34801, p-value < 2.2e-16; large effect size: 0.798)


diff_mod3 <- wilcox.test(value ~ variable, data = effectSize_heritability[instar == "instar 2"][variable %in% c("H2_VAoverVM_ctrl", "H2_VAoverVM_trt")][module_new == "mod 3"], paired = FALSE, alternative = "greater")
diff_mod3

p_diff_mod3 = qnorm(diff_mod3$p.value/2)
p_diff_mod3

effectsize_diff_mod3 <- wilcox_effsize(value ~ variable_new, data = effectSize_heritability[instar == "instar 2"][variable %in% c("H2_VAoverVM_ctrl", "H2_VAoverVM_trt")][module_new == "mod 3"], paired = FALSE, alternative = "greater")
effectsize_diff_mod3

# (Z = -1.081356e-06, p-value = 1; small effect size: 0.169)



# VAoverVM ratio
box2 <- ggplot(effectSize_heritability[variable %in% c('H2_VAoverVM_trt','H2_VAoverVM_ctrl')][instar == "instar 2"], aes(y=value, fill=(variable))) + 
                geom_boxplot(notch=TRUE, outlier.colour="#555555", outlier.shape=8, outlier.size=2) + 
                facet_grid(~module_new) + 
                theme(legend.position="none", 
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


load(file = "data/all_data_final.RData")
# all_data_final

## AVERAGE SHAPE 
all_data_final.ag <- all_data_final[i <= 600][, list(height = mean(height)), list(i, treatment, instar_new, SC_group_new)]


ag.line_plot_2 <- ggplot(data = all_data_final.ag[treatment == 0.5][instar_new == "instar 2"],  aes(x=i, y=height)) +   
  
                          # geom_rect(aes(xmin = 10, xmax = 100, ymin = -Inf, ymax = Inf),fill = "#F5F5F5", colour="#F5F5F5", alpha=0.1) +
                          # geom_rect(aes(xmin = 101, xmax = 200, ymin = -Inf, ymax = Inf),fill = "#E8E8E8", colour="#E8E8E8", alpha=0.1) +
                          # geom_rect(aes(xmin = 201, xmax = 600, ymin = -Inf, ymax = Inf),fill = "#DCDCDC", colour="#DCDCDC", alpha=0.1) +
                          
                          geom_vline(xintercept = c(10, 100, 200, 600), linetype="dotted", color = "black", size=1) + 
                          
                          geom_line(aes(x=i, y=height), size = 1.5, colour = "red") + 
                          
                          facet_grid(~instar_new) +
                          
                          ylim(0, 0.28) +
                          
                          theme(legend.position="none", 
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
                                axis.text = element_text(size=15, family='Arial'),
                                strip.text.x = element_text(size = 15, color = "black"),
                                strip.text.y = element_text(size = 15, color = "black"),
                                panel.spacing.x = unit(6, "mm"),
                                panel.spacing.y = unit(6, "mm")) 

# tiff(file = "avg_shape_plot.tiff", width = 3200, height = 3200, units = "px", res = 800)
ag.line_plot_2  + labs(x = "dorsal position", y = "dorsal height") + scale_x_continuous(limits=c(10,600), breaks=c(150,300,600)) 
# dev.off()


########
########


R2_all_plot + scale_color_manual(values=c("#000000","#FF0000","#D3D3D3","#fed4d2")) + labs(x = "phenotypic trait", y = "correlation coefficient") 

VaVm_instar_i_plot + labs(x = "dorsal position", y = expression(log~"("~V[g]~"/"~V[m]~")")) + scale_x_continuous(limits=c(10,600), breaks=c(0,300,600)) 

ag.line_plot_2  + labs(x = "dorsal position", y = "dorsal height") + scale_x_continuous(limits=c(10,600), breaks=c(150,300,600)) 

box2 + labs(x = "module", y = expression(log[10]~"("~V[g]~"/"~V[m]~")")) + scale_fill_manual(values=c("#000000","#FF0000"))


# patchwork_plots_induction4 <- ( (R2_all_plot + scale_color_manual(values=c("#000000","#FF0000","#D3D3D3","#fed4d2")) + labs(x = "phenotypic trait", y = "correlation coefficient")) + plot_spacer() ) /
#  ( (VaVm_instar_i_plot + labs(x = "dorsal position", y = expression(log[10]~"("~V[g]~"/"~V[m]~")")) + scale_x_continuous(limits=c(10,600), breaks=c(0,300,600))) + plot_spacer() )  # plot_spacer() +
#   # plot_layout(ncol=3, widths = c(1,1,0.5))
# 
# patchwork_plots_induction4


# patchwork_plots_induction4_A <- R2_all_plot + scale_color_manual(values=c("#000000","#FF0000","#D3D3D3","#fed4d2")) + labs(x = "phenotypic trait", y = "correlation coefficient") +
#   plot_spacer() +
#   plot_layout(ncol=2, widths = c(2,1.8))
# 
# patchwork_plots_induction4_A
#   
#   
# patchwork_plots_induction4_B <- VaVm_instar_i_plot + labs(x = "dorsal position", y = expression(log[10]~"("~V[g]~"/"~V[m]~")")) + scale_x_continuous(limits=c(10,600), breaks=c(0,300,600)) +
#   plot_spacer() +
#   plot_layout(ncol=2, widths = c(2,1.8))
# 
# patchwork_plots_induction4_B
# 



patchwork_plots_induction4 <- R2_all_plot + scale_color_manual(values=c("#000000","#FF0000","#D3D3D3","#fed4d2")) + labs(x = "phenotypic trait", y = "correlation coefficient") +
  plot_spacer() +
  plot_layout(ncol=3, widths = c(2,1,1))  
  
patchwork_plots_induction4


patchwork_plots_induction4.2 <- ( VaVm_instar_i_plot + labs(x = "dorsal position", y = expression(log~"("~V[g]~"/"~V[m]~")")) + scale_x_continuous(limits=c(10,600), breaks=c(0,300,600)) ) |
  ( ag.line_plot_2  + labs(x = "dorsal position", y = "dorsal height") + scale_x_continuous(limits=c(10,600), breaks=c(150,300,600)) ) / 
    ( box2 + labs(x = "module", y = expression(log[10]~"("~V[g]~"/"~V[m]~")")) + scale_fill_manual(values=c("#000000","#FF0000")) ) + 
  plot_layout(nrow=2, heights = c(0.4,1)) |
  plot_spacer() |
  plot_spacer() 

patchwork_plots_induction4.2



  ( plot_spacer() + plot_spacer() ) + plot_layout(ncol=2, widths = c(2,1)) 
  plot_layout(ncol=3, widths = c(2,1,1))  

patchwork_plots_induction4.2




( box2 + plot_spacer() ) / ( box2 + plot_spacer() ) + plot_layout(nrow=2, heights = c(0.4,1)) 


  
