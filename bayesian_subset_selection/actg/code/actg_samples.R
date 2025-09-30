# load libraries
library(hdbayes)
library(parallel)
library(matrixStats)
library(MCMCpack)

# load auxiliary functions
source("bayesian_subset_selection/actg/code/aux_scripts/functions.R")

# define data
current_data <- actg036
hist_data <- actg019

# normalise data
current_data$age <- (current_data$age - mean(current_data$age)) /
  (2*sd(current_data$age))
current_data$cd4 <- (current_data$cd4 - mean(current_data$cd4)) /
  (2*sd(current_data$cd4))
hist_data$age <- (hist_data$age - mean(hist_data$age)) /
  (2*sd(hist_data$age))
hist_data$cd4 <- (hist_data$cd4 - mean(hist_data$cd4)) /
  (2*sd(hist_data$cd4))

# set parameters
# beta_pars <- elicit_beta_mean_cv(m0 = 0.5, v0 = 0.005)
delta0 <- 1
lambda0 <- 1
a0_seq <- seq(0, 1, length.out = 21)
data <- list(current_data, hist_data)
family <- binomial(link = "logit")
c0 <- 1
d0 <- 1
# delta0 <- beta_pars$a
# lambda0 <- beta_pars$b
iter_warmup <- 1000
iter_sampling <- 2500

# run complete model
post_samples <- samples_models(data = data,
                               outcome = "outcome",
                               family = family,
                               a0_seq = a0_seq,
                               c0 = c0,
                               d0 = d0,
                               delta0 = delta0,
                               lambda0 = lambda0,
                               iter_warmup = iter_warmup,
                               iter_sampling = iter_sampling,
                               num_cores = 14)

# # run the model for each group
# post_samples_ctrl <- samples_models(data = data_ctrl,
#                                     outcome = "outcome",
#                                     family = family,
#                                     a0_seq = a0_seq,
#                                     c0 = c0,
#                                     d0 = d0,
#                                     delta0 = delta0,
#                                     lambda0 = lambda0,
#                                     iter_warmup = iter_warmup,
#                                     iter_sampling = iter_sampling,
#                                     num_cores = 14)
# 
# post_samples_trt <- samples_models(data = data_trt,
#                                    outcome = "outcome",
#                                    family = family,
#                                    a0_seq = a0_seq,
#                                    c0 = c0,
#                                    d0 = d0,
#                                    delta0 = delta0,
#                                    lambda0 = lambda0,
#                                    iter_warmup = iter_warmup,
#                                    iter_sampling = iter_sampling,
#                                    num_cores = 14)

# calculate marginal means by arm
# mean_models_ctrl <- mean_models_arm(current_data, "outcome", "treatment", family, post_samples_ctrl$post_betam)
# mean_models_trt <- mean_models_arm(current_data, "outcome", "treatment", family, post_samples_trt$post_betam)
mean_models_ctrl <- mean_models_arm(current_data, 
                                    "outcome", 
                                    "treatment", 
                                    family, 
                                    post_samples$post_betam,
                                    0)
mean_models_trt <- mean_models_arm(current_data, 
                                   "outcome", 
                                   "treatment", 
                                   family, 
                                   post_samples$post_betam,
                                   1)
# compute BMA
bma_ctrl <- bma(mean_models_ctrl, post_samples$df_post, 10000)
bma_trt <- bma(mean_models_trt, post_samples$df_post, 10000)


# save results
save(post_samples, file = "bayesian_subset_selection/actg/samples/post_samples.RData")
# save(post_samples_ctrl, file = "bayesian_subset_selection/actg/samples/post_samples_ctrl.RData")
# save(post_samples_trt, file = "bayesian_subset_selection/actg/samples/post_samples_trt.RData")
save(mean_models_ctrl, file = "bayesian_subset_selection/actg/samples/mean_models_ctrl.RData")
save(mean_models_trt, file = "bayesian_subset_selection/actg/samples/mean_models_trt.RData")
save(bma_ctrl, file = "bayesian_subset_selection/actg/samples/bma_ctrl.RData")
save(bma_trt, file = "bayesian_subset_selection/actg/samples/bma_trt.RData")

#-------------------------------------------------------------------------------
# Sample after PSM 
#-------------------------------------------------------------------------------

current_data <- actg036
hist_data <- readRDS("bayesian_subset_selection/actg/data/actg019_after_PSM.rds")

# normalise data
current_data$age <- (current_data$age - mean(current_data$age)) /
  (2*sd(current_data$age))
current_data$cd4 <- (current_data$cd4 - mean(current_data$cd4)) /
  (2*sd(current_data$cd4))
hist_data$age <- (hist_data$age - mean(hist_data$age)) /
  (2*sd(hist_data$age))
hist_data$cd4 <- (hist_data$cd4 - mean(hist_data$cd4)) /
  (2*sd(hist_data$cd4))

# set parameters
delta0 <- 1
lambda0 <- 1
a0_seq <- seq(0, 1, length.out = 21)
data <- list(current_data, hist_data)
family <- binomial(link = "logit")
c0 <- 0.5
d0 <- 0.5
iter_warmup <- 1000
iter_sampling <- 2500

# run complete model
post_samples_after_PSM <- samples_models(data = data,
                               outcome = "outcome",
                               family = family,
                               a0_seq = a0_seq,
                               c0 = c0,
                               d0 = d0,
                               delta0 = delta0,
                               lambda0 = lambda0,
                               iter_warmup = iter_warmup,
                               iter_sampling = iter_sampling,
                               num_cores = 14)
mean_models_ctrl_after_PSM <- mean_models_arm(current_data, 
                                    "outcome", 
                                    "treatment", 
                                    family, 
                                    post_samples_after_PSM$post_betam,
                                    0)
mean_models_trt_after_PSM <- mean_models_arm(current_data, 
                                   "outcome", 
                                   "treatment", 
                                   family, 
                                   post_samples_after_PSM$post_betam,
                                   1)

# compute BMA
bma_ctrl_after_PSM <- bma(mean_models_ctrl_after_PSM, 
                          post_samples_after_PSM$df_post, 10000)
bma_trt_after_PSM <- bma(mean_models_trt_after_PSM, 
                         post_samples_after_PSM$df_post, 10000)

save(post_samples_after_PSM, 
     file = "bayesian_subset_selection/actg/samples/post_samples_after_PSM.RData")
save(mean_models_ctrl_after_PSM,
     file = "bayesian_subset_selection/actg/samples/mean_models_ctrl_after_PSM.RData")
save(mean_models_trt_after_PSM,
     file = "bayesian_subset_selection/actg/samples/mean_models_trt_after_PSM.RData")
save(bma_ctrl_after_PSM,
     file = "bayesian_subset_selection/actg/samples/bma_ctrl_after_PSM.RData")
save(bma_trt_after_PSM,
     file = "bayesian_subset_selection/actg/samples/bma_trt_after_PSM.RData")
