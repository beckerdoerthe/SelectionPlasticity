# h2 mcmcglmm

library(data.table)
library(tidyverse)
library(MCMCglmm)
library(doMC)
registerDoMC(4)

#data
##########################
### prepare pheno data ###		
##########################


load("data/all_data_final.RData")
head(all_data_final)


## TREATMENT - CONTROL
get_deltah2 <- function(position){
  tmp_c <- filter(all_data_final, i == position, treatment == 0.0, SC_group == "O", instar == 2)
  tmp_t <- filter(all_data_final, i == position, treatment == 0.5, SC_group == "O", instar == 2)
  
  mod_c <- MCMCglmm(height ~ 1, ~cloneid_geno + batch,  
                    data = tmp_c,
                    nitt=65000,thin=50,burnin=15000,
                    verbose = FALSE)
  
  mod_t <- MCMCglmm(height ~ 1, ~cloneid_geno + batch,  
                    data = tmp_t,
                    nitt=65000,thin=50,burnin=15000,
                    verbose = FALSE)
  
  # h2 posteriors for controls
  h2_posts_c <- mod_c$VCV[,"cloneid_geno"]/(mod_c$VCV[,"cloneid_geno"]+mod_c$VCV[,"units"])  
  
  # h2 posteriors for treats
  h2_posts_t <- mod_t$VCV[,"cloneid_geno"]/(mod_t$VCV[,"cloneid_geno"]+mod_t$VCV[,"units"])
  
  # delta h2
  delta_h2_posts <- h2_posts_t - h2_posts_c
  
  # all  modes (consider means?)
  h2_c_mean <- mean(h2_posts_c)
  h2_c_mode <- posterior.mode(h2_posts_c)
  
  h2_t_mean <- mean(h2_posts_t)
  h2_t_mode <- posterior.mode(h2_posts_t)
  
  delta_h2_mean <- mean(delta_h2_posts)
  delta_h2_mode <- posterior.mode(delta_h2_posts)
  
  # all intervals
  h2_CI_c <- HPDinterval(h2_posts_c, 0.95)
  h2_CI_t <- HPDinterval(h2_posts_t, 0.95)
  delta_h2_CI <- HPDinterval(delta_h2_posts, 0.95)
  
  # collect the pieces
  c_modeCI <- data.frame(lCI = h2_CI_c[1],
                         mode = h2_c_mode, 
                         mean = h2_c_mean,
                         uCI = h2_CI_c[2])
  
  t_modeCI <- data.frame(lCI = h2_CI_t[1],
                         mode = h2_t_mode, 
                         mean = h2_t_mean, 
                         uCI = h2_CI_t[2])
  
  delta_modeCI <- data.frame(lCI = delta_h2_CI[1],
                             mode = delta_h2_mode, 
                             mean = delta_h2_mean, 
                             uCI = delta_h2_CI[2])
  
  # return the stuff
  return(data.frame(
    i = position,
    label = c("C", "T", "DELTA"),
    stuff = bind_rows(c_modeCI, t_modeCI, delta_modeCI)))
}

# use on i = 151
# get_deltah2(position = 151)

# using foreach
collected2 <- foreach(i = 10:600, .combine = rbind) %dopar%
  get_deltah2(position = i)


# w/ batch effect
# save(collected2, file = "H2_O_I1_batch.RData")
# save(collected2, file = "H2_O_I2_batch.RData")
# save(collected2, file = "H2_A_I1_batch.RData")
# save(collected2, file = "H2_A_I2_batch.RData")



## VA/Vm
get_deltah2 <- function(position){
  tmp_O <- filter(all_data_final, i == position, treatment == 0.5, SC_group == "O", instar == 2)
  tmp_A <- filter(all_data_final, i == position, treatment == 0.5, SC_group == "A", instar == 2)
  
  mod_O <- MCMCglmm(height ~ 1, ~cloneid_geno + batch,  
                    data = tmp_O,
                    nitt=65000,thin=50,burnin=15000,
                    verbose = FALSE)
  
  mod_A <- MCMCglmm(height ~ 1, ~cloneid_geno + batch,
                    data = tmp_A,
                    nitt=65000,thin=50,burnin=15000,
                    verbose = FALSE)
  
  # h2 posteriors for Os
  h2_posts_O <- mod_O$VCV[,"cloneid_geno"]/(mod_O$VCV[,"cloneid_geno"]+mod_O$VCV[,"units"])  
  
  # h2 posteriors for As
  h2_posts_A <- mod_A$VCV[,"cloneid_geno"]/(mod_A$VCV[,"cloneid_geno"]+mod_A$VCV[,"units"])
  
  # delta h2
  delta_h2_posts <- log10(h2_posts_O) - log10(h2_posts_A)
  
  # all modes & means
  h2_O_mean <- mean(h2_posts_O)
  h2_O_mode <- posterior.mode(h2_posts_O)
  
  h2_A_mean <- mean(h2_posts_A)
  h2_A_mode <- posterior.mode(h2_posts_A)
  
  delta_h2_mean <- mean(delta_h2_posts)
  delta_h2_mode <- posterior.mode(delta_h2_posts)
  
  # all intervals
  h2_CI_O <- HPDinterval(h2_posts_O, 0.95)
  h2_CI_A <- HPDinterval(h2_posts_A, 0.95)
  delta_h2_CI <- HPDinterval(delta_h2_posts, 0.95)
  
  # collect the pieces
  O_modeCI <- data.frame(lCI = h2_CI_O[1],
                         mode = h2_O_mode, 
                         mean = h2_O_mean,
                         uCI = h2_CI_O[2])
  
  A_modeCI <- data.frame(lCI = h2_CI_A[1],
                         mode = h2_A_mode, 
                         mean = h2_A_mean, 
                         uCI = h2_CI_A[2])
  
  delta_modeCI <- data.frame(lCI = delta_h2_CI[1],
                             mode = delta_h2_mode, 
                             mean = delta_h2_mean, 
                             uCI = delta_h2_CI[2])
  
  # return the stuff
  return(data.frame(
    i = position,
    label = c("O", "A", "RATIO"),
    stuff = bind_rows(O_modeCI, A_modeCI, delta_modeCI)))
}


# use on i = 151
# get_deltah2(position = 151)

# using foreach
collected2 <- foreach(i = 10:600, .combine = rbind) %dopar%
  get_deltah2(position = i)


# w/ batch effect
# save(collected2, file = "H2_0_I1_ratio_batch.RData")
# save(collected2, file = "H2_0_I2_ratio_batch.RData")
# save(collected2, file = "H2_05_I1_ratio_batch.RData")
# save(collected2, file = "H2_05_I2_ratio_batch.RData")










# ## treatment effect
# get_deltah2 <- function(position){
#   tmp_O_I1 <- filter(all_data_final, i == position, SC_group == "O", instar == 1)
#   tmp_O_I2 <- filter(all_data_final, i == position, SC_group == "O", instar == 2)
#   tmp_A_I1 <- filter(all_data_final, i == position, SC_group == "A", instar == 1)
#   tmp_A_I2 <- filter(all_data_final, i == position, SC_group == "A", instar == 2)
#   
#   mod_O_I1 <- MCMCglmm(height ~ treatment, ~cloneid_geno + batch, 
#                        data = tmp_O_I1,
#                        nitt=65000,thin=50,burnin=15000,
#                        verbose = FALSE)
#   
#   mod_O_I2 <- MCMCglmm(height ~ treatment, ~cloneid_geno + batch, 
#                        data = tmp_O_I2,
#                        nitt=65000,thin=50,burnin=15000,
#                        verbose = FALSE)
#   
#   mod_A_I1 <- MCMCglmm(height ~ treatment, ~cloneid_geno + batch, 
#                        data = tmp_A_I1,
#                        nitt=65000,thin=50,burnin=15000,
#                        verbose = FALSE)
#   
#   mod_A_I2 <- MCMCglmm(height ~ treatment, ~cloneid_geno + batch, 
#                        data = tmp_A_I2,
#                        nitt=65000,thin=50,burnin=15000,
#                        verbose = FALSE)
#   
#   mod_O_I1_out <- as.data.table(summary(mod_O_I1)$solutions)
#   mod_O_I2_out <- as.data.table(summary(mod_O_I2)$solutions)
#   mod_A_I1_out <- as.data.table(summary(mod_A_I1)$solutions)
#   mod_A_I2_out <- as.data.table(summary(mod_A_I2)$solutions)
#   
#   
#   # collect the pieces
#   mod_O_I1_CI <- data.frame(lCI = mod_O_I1_out$'l-95% CI'[1] + mod_O_I1_out$'l-95% CI'[2],
#                            mean = mod_O_I1_out$'post.mean'[1] + mod_O_I1_out$'post.mean'[2],
#                            uCI = mod_O_I1_out$'u-95% CI'[1] + mod_O_I1_out$'u-95% CI'[2])
# 
#   mod_O_I2_CI <- data.frame(lCI = mod_O_I2_out$'l-95% CI'[1] + mod_O_I2_out$'l-95% CI'[2],
#                             mean = mod_O_I2_out$'post.mean'[1] + mod_O_I2_out$'post.mean'[2],
#                             uCI = mod_O_I2_out$'u-95% CI'[1] + mod_O_I2_out$'u-95% CI'[2])
# 
#   mod_A_I1_CI <- data.frame(lCI = mod_A_I1_out$'l-95% CI'[1] + mod_A_I1_out$'l-95% CI'[2],
#                             mean = mod_A_I1_out$'post.mean'[1] + mod_A_I1_out$'post.mean'[2],
#                             uCI = mod_A_I1_out$'u-95% CI'[1] + mod_A_I1_out$'u-95% CI'[2])
# 
#   mod_A_I2_CI <- data.frame(lCI = mod_A_I2_out$'l-95% CI'[1] + mod_A_I2_out$'l-95% CI'[2],
#                             mean = mod_A_I2_out$'post.mean'[1] + mod_A_I2_out$'post.mean'[2],
#                             uCI = mod_A_I2_out$'u-95% CI'[1] + mod_A_I2_out$'u-95% CI'[2])
#   
#   # return the stuff
#   return(data.frame(
#     i = position,
#     label = c("O_I1", "O_I2", "A_I1", "A_I2"),
#     stuff = bind_rows(mod_O_I1_CI, mod_O_I2_CI, mod_A_I1_CI, mod_A_I2_CI)))
# }
# 
# 
# # use on i = 151
# # get_deltah2(position = 10)
# 
# # using foreach
# collected2 <- foreach(i = 10:650, .combine = rbind) %dopar%
#   get_deltah2(position = i)
# 
# 
# # w/o batch effect
# # save(collected2, file = "H2_treatment.RData")
# 
# # w/ batch effect
# save(collected2, file = "H2_treatment_batch.RData")
# 
# 
# 
# ### traits
# 
# #load("subset_final_data_forAPB.RData")
# 
# load("traits_response_all_OA.Rdata")
# head(traits_response_all_OA)
# 
# 
# # make a function to do the work on position i
# 
# ## here 'position is actually test_group'
# 
# get_deltah2 <- function(position){
#   tmp_c <- filter(traits_response_all_OA, test_group == position, treatment == 0.0, SC_group == "O", instar == 2)
#   tmp_t <- filter(traits_response_all_OA, test_group == position, treatment == 0.5, SC_group == "O", instar == 2)
#   
#   mod_c <- MCMCglmm(response ~ 1, ~cloneid_geno + batch,  
#                     data = tmp_c,
#                     nitt=65000,thin=50,burnin=15000,
#                     verbose = FALSE)
#   
#   mod_t <- MCMCglmm(response ~ 1, ~cloneid_geno + batch,  
#                     data = tmp_t,
#                     nitt=65000,thin=50,burnin=15000,
#                     verbose = FALSE)
#   
#   # h2 posteriors for controls
#   h2_posts_c <- mod_c$VCV[,"cloneid_geno"]/(mod_c$VCV[,"cloneid_geno"]+mod_c$VCV[,"units"])  
#   
#   # h2 posteriors for treats
#   h2_posts_t <- mod_t$VCV[,"cloneid_geno"]/(mod_t$VCV[,"cloneid_geno"]+mod_t$VCV[,"units"])
#   
#   # delta h2
#   delta_h2_posts <- h2_posts_t - h2_posts_c
# 
#   # all  modes (consider means?)
#   h2_c_mean <- mean(h2_posts_c)
#   h2_c_mode <- posterior.mode(h2_posts_c)
# 
#   h2_t_mean <- mean(h2_posts_t)
#   h2_t_mode <- posterior.mode(h2_posts_t)
#   
#   delta_h2_mean <- mean(delta_h2_posts)
#   delta_h2_mode <- posterior.mode(delta_h2_posts)
#   
#   # all intervals
#   h2_CI_c <- HPDinterval(h2_posts_c, 0.95)
#   h2_CI_t <- HPDinterval(h2_posts_t, 0.95)
#   delta_h2_CI <- HPDinterval(delta_h2_posts, 0.95)
#   
#   # collect the pieces
#   c_modeCI <- data.frame(lCI = h2_CI_c[1],
#                          mode = h2_c_mode, 
#                          mean = h2_c_mean,
#                          uCI = h2_CI_c[2])
#   
#   t_modeCI <- data.frame(lCI = h2_CI_t[1],
#                          mode = h2_t_mode, 
#                          mean = h2_t_mean, 
#                          uCI = h2_CI_t[2])
#   
#   delta_modeCI <- data.frame(lCI = delta_h2_CI[1],
#                              mode = delta_h2_mode, 
#                              mean = delta_h2_mean, 
#                              uCI = delta_h2_CI[2])
#   
#   # return the stuff
#   return(data.frame(
#     i = position,
#     label = c("C", "T", "DELTA"),
#     stuff = bind_rows(c_modeCI, t_modeCI, delta_modeCI)))
# }
# 
# 
# # using foreach
# collected2 <- foreach(i = unique(traits_response_all_OA$test_group), .combine = rbind) %dopar%
#   get_deltah2(position = i)
# 
# 
# # save(collected2, file = "H2_traits_O_I1.RData")
# # save(collected2, file = "H2_traits_O_I2.RData")
# # save(collected2, file = "H2_traits_A_I1.RData")
# # save(collected2, file = "H2_traits_A_I2.RData")
# 
# # save(collected2, file = "H2_traits_batch_O_I1.RData")
# # save(collected2, file = "H2_traits_batch_O_I2.RData")
# # save(collected2, file = "H2_traits_batch_A_I1.RData")
# # save(collected2, file = "H2_traits_batch_A_I2.RData")
