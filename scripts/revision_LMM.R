## variance components of lmer() as alternative to correlation (twin) analysis (Figure 3A)


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


anova(model_length_1,model_length_2)
anova(model_length_1,model_length_3)
anova(model_length_2,model_length_3)


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
anova(model_maxHeight_2,model_maxHeight_3)
