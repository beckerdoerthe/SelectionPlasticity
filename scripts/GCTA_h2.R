#################
### libraries ###
#################
library(ggplot2)
library(cowplot)
library(viridis)
library(data.table)
library(geomorph)
library(reshape2)
library(dplyr)
library(tidyverse)
library(foreach)
library(abind)
library(doMC)
registerDoMC(20)

library(SeqArray)
library(SNPRelate)
library(snpReady)



################################################
## generate phenotype and covfiles for each i ##
################################################

load("data/all_data_final.RData")
all_data_final
setkey(all_data_final, GenoPLUS, i)


pheno_data_use_O_I1_0 <- all_data_final[group %in% c('O','ctrl_O')][instar == 1][treatment == 0]
pheno_data_use_O_I1_0[, id := GenoPLUS]

pheno_data_use_O_I1_05 <- all_data_final[group %in% c('O','ctrl_O')][instar == 1][treatment == 0.5]
pheno_data_use_O_I1_05[, id := GenoPLUS]

pheno_data_use_O_I2_0 <- all_data_final[group %in% c('O','ctrl_O')][instar == 2][treatment == 0]
pheno_data_use_O_I2_0[, id := GenoPLUS]

pheno_data_use_O_I2_05 <- all_data_final[group %in% c('O','ctrl_O')][instar == 2][treatment == 0.5]
pheno_data_use_O_I2_05[, id := GenoPLUS]



foreach(i.i = unique(pheno_data_use_O_I1_0$i), .errorhandling="remove", .combine="rbind") %do% {
  
  phenos <- pheno_data_use_O_I1_0[i == i.i]
  
  write.table(phenos[,c("id", "id", "height")], paste0("/PATH_TO_FILES/I1_0_plink_phenos_i", i.i, ".txt"), quote=F, sep="\t", row.names=F, col.names=F)
  write.table(phenos[,c("id", "id", "batch")], paste0("/PATH_TO_FILES/I1_0_cov_gcta_i", i.i, ".txt"), quote=F, sep="\t", row.names=F, col.names=F)
}


foreach(i.i = unique(pheno_data_use_O_I1_05$i), .errorhandling="remove", .combine="rbind") %do% {
  
  phenos <- pheno_data_use_O_I1_05[i == i.i]
  
  write.table(phenos[,c("id", "id", "height")], paste0("/PATH_TO_FILES/I1_05_plink_phenos_i", i.i, ".txt"), quote=F, sep="\t", row.names=F, col.names=F)
  write.table(phenos[,c("id", "id", "batch")], paste0("/PATH_TO_FILES/I1_05_cov_gcta_i", i.i, ".txt"), quote=F, sep="\t", row.names=F, col.names=F)
}


foreach(i.i = unique(pheno_data_use_O_I2_0$i), .errorhandling="remove", .combine="rbind") %do% {
  
  phenos <- pheno_data_use_O_I2_0[i == i.i]
  
  write.table(phenos[,c("id", "id", "height")], paste0("/PATH_TO_FILES/I2_0_plink_phenos_i", i.i, ".txt"), quote=F, sep="\t", row.names=F, col.names=F)
  write.table(phenos[,c("id", "id", "batch")], paste0("/PATH_TO_FILES/I2_0_cov_gcta_i", i.i, ".txt"), quote=F, sep="\t", row.names=F, col.names=F)
}


foreach(i.i = unique(pheno_data_use_O_I2_05$i), .errorhandling="remove", .combine="rbind") %do% {
  
  phenos <- pheno_data_use_O_I2_05[i == i.i]
  
  write.table(phenos[,c("id", "id", "height")], paste0("/PATH_TO_FILES/I2_05_plink_phenos_i", i.i, ".txt"), quote=F, sep="\t", row.names=F, col.names=F)
  write.table(phenos[,c("id", "id", "batch")], paste0("/PATH_TO_FILES/I2_05_cov_gcta_i", i.i, ".txt"), quote=F, sep="\t", row.names=F, col.names=F)
}



######################################
## calculate heritability estimates ##
######################################

## Os
# along i; batch as covariate & MAF 0.05

for i in {10..650}; do
echo $i

/PATH_TO_GCTA/gcta_1.92.1beta6/gcta64 \
--grm /PATH_TO_FILES/O_I1_0_chrnum_MAF.05 \
--pheno /PATH_TO_FILES/I1_0_plink_phenos_i$i.txt \
--covar /PATH_TO_FILES/I1_0_cov_gcta_i$i.txt \
--reml \
--reml-alg 1 \
--out /PATH_TO_FILES/I1_0_COVbatch_MAF_gcta_i$i \
--thread-num 10

done


for i in {10..650}; do
echo $i

/PATH_TO_GCTA/gcta_1.92.1beta6/gcta64 \
--grm /PATH_TO_FILES/O_I1_05_chrnum_MAF.05 \
--pheno /PATH_TO_FILES/I1_05_plink_phenos_i$i.txt \
--covar /PATH_TO_FILES/I1_05_cov_gcta_i$i.txt \
--reml \
--reml-alg 1 \
--out /PATH_TO_FILES/I1_05_COVbatch_MAF_gcta_i$i \
--thread-num 10

done


for i in {10..650}; do
echo $i

/PATH_TO_GCTA/gcta_1.92.1beta6/gcta64 \
--grm /PATH_TO_FILES/O_I2_0_chrnum_MAF.05 \
--pheno /PATH_TO_FILES/I2_0_plink_phenos_i$i.txt \
--covar /PATH_TO_FILES/I2_0_cov_gcta_i$i.txt \
--reml \
--reml-alg 1 \
--out /PATH_TO_FILES/I2_0_COVbatch_MAF_gcta_i$i \
--thread-num 10

done


for i in {10..650}; do
echo $i

/PATH_TO_GCTA/gcta_1.92.1beta6/gcta64 \
--grm /PATH_TO_FILES/O_I2_05_chrnum_MAF.05 \
--pheno /PATH_TO_FILES/I2_05_plink_phenos_i$i.txt \
--covar /PATH_TO_FILES/I2_05_cov_gcta_i$i.txt \
--reml \
--reml-alg 1 \
--out /PATH_TO_FILES/I2_05_COVbatch_MAF_gcta_i$i \
--thread-num 10

done



## exctract heritability 

rm(list=ls()) 

library(ggplot2)
library(cowplot)
library(viridis)
library(plotrix)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyverse)
library(foreach)
library(abind)
library(doMC)
registerDoMC(20)


# I1_0
files_path_sites_h2 <- unlist(foreach(i=c(10:650))%do%{list.files(path="/PATH_TO_FILES/", pattern=paste0("I1_0_COVbatch_MAF_gcta_i",i,".hsq"), all.files=T)})

h2_estimate_I1_0 <- foreach(file = files_path_sites_h2) %do% {
  
  #load data
  tmp_data <- fread(paste0("/PATH_TO_FILES/", file), fill=T)
  
  tmp_data[, filename := file]              
  tmp_data[, i := gsub(".*_gcta_i(.*)\\..*", "\\1", tmp_data$file)]
  tmp_data[, perm := 0]
  
  as.data.table(tmp_data)
  
}

h2_estimate_I1_0_out <- rbindlist(h2_estimate_I1_0)
h2_estimate_I1_0_out[, group := "I1_0"]


# I1_05
files_path_sites_h2 <- unlist(foreach(i=c(10:650))%do%{list.files(path="/PATH_TO_FILES/", pattern=paste0("I1_05_COVbatch_MAF_gcta_i",i,".hsq"), all.files=T)})

h2_estimate_I1_05 <- foreach(file = files_path_sites_h2) %do% {
  
  #load data
  tmp_data <- fread(paste0("/PATH_TO_FILES/", file), fill=T)
  
  tmp_data[, filename := file]              
  tmp_data[, i := gsub(".*_gcta_i(.*)\\..*", "\\1", tmp_data$file)]
  tmp_data[, perm := 0]
  
  as.data.table(tmp_data)
  
}

h2_estimate_I1_05_out <- rbindlist(h2_estimate_I1_05)
h2_estimate_I1_05_out[, group := "I1_05"]


# I2_0
files_path_sites_h2 <- unlist(foreach(i=c(10:650))%do%{list.files(path="/PATH_TO_FILES/", pattern=paste0("I2_0_COVbatch_MAF_gcta_i",i,".hsq"), all.files=T)})

h2_estimate_I2_0 <- foreach(file = files_path_sites_h2) %do% {
  
  #load data
  tmp_data <- fread(paste0("/PATH_TO_FILES/", file), fill=T)
  
  tmp_data[, filename := file]              
  tmp_data[, i := gsub(".*_gcta_i(.*)\\..*", "\\1", tmp_data$file)]
  tmp_data[, perm := 0]
  
  as.data.table(tmp_data)
  
}

h2_estimate_I2_0_out <- rbindlist(h2_estimate_I2_0)
h2_estimate_I2_0_out[, group := "I2_0"]


# I2_05
files_path_sites_h2 <- unlist(foreach(i=c(10:650))%do%{list.files(path="/PATH_TO_FILES/", pattern=paste0("I2_05_COVbatch_MAF_gcta_i",i,".hsq"), all.files=T)})

h2_estimate_I2_05 <- foreach(file = files_path_sites_h2) %do% {
  
  #load data
  tmp_data <- fread(paste0("/PATH_TO_FILES/", file), fill=T)
  
  tmp_data[, filename := file]              
  tmp_data[, i := gsub(".*_gcta_i(.*)\\..*", "\\1", tmp_data$file)]
  tmp_data[, perm := 0]
  
  as.data.table(tmp_data)
  
}

h2_estimate_I2_05_out <- rbindlist(h2_estimate_I2_05)
h2_estimate_I2_05_out[, group := "I2_05"]


save(h2_estimate_I1_0_out,h2_estimate_I1_05_out,h2_estimate_I2_0_out,h2_estimate_I2_05_out, file='output/O_h2_out_i_14June.RData')


