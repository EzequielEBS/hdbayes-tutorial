# load libraries
library(hdbayes)
library(parallel)
library(matrixStats)
library(MCMCpack)

# load auxiliary functions
source("bayesian_subset_selection/actg/code/functions.R")

# define data
current_data <- actg036
hist_data <- actg019

# split data
current_data_ctrl <- current_data[current_data$treatment == 0, 
                                  !(names(current_data) %in% c("treatment"))]
current_data_trt <- current_data[current_data$treatment == 1, 
                                 !(names(current_data) %in% c("treatment"))]

# normalize data
age_stats <- with(current_data,
              c('mean' = mean(age), 'sd' = sd(age)))
cd4_stats <- with(current_data,
              c('mean' = mean(cd4), 'sd' = sd(cd4)))
age_stats_ctrl <- with(current_data_ctrl,
                       c('mean' = mean(age), 'sd' = sd(age)))
cd4_stats_ctrl <- with(current_data_ctrl,
                       c('mean' = mean(cd4), 'sd' = sd(cd4)))
age_stats_trt <- with(current_data_trt,
                      c('mean' = mean(age), 'sd' = sd(age)))
cd4_stats_trt <- with(current_data_trt,
                      c('mean' = mean(cd4), 'sd' = sd(cd4)))
age_stats_hist <- with(hist_data,
                       c('mean' = mean(age), 'sd' = sd(age)))
cd4_stats_hist <- with(hist_data,
                       c('mean' = mean(cd4), 'sd' = sd(cd4)))

current_data$age <- (current_data$age - age_stats['mean']) / age_stats['sd']
current_data$cd4 <- (current_data$cd4 - cd4_stats['mean']) / cd4_stats['sd']
# current_data_ctrl$age <- (current_data_ctrl$age - age_stats['mean']) / age_stats['sd']
# current_data_ctrl$cd4 <- (current_data_ctrl$cd4 - cd4_stats['mean']) / cd4_stats['sd']
# current_data_trt$age <- (current_data_trt$age - age_stats['mean']) / age_stats['sd']
# current_data_trt$cd4 <- (current_data_trt$cd4 - cd4_stats['mean']) / cd4_stats['sd']
current_data_ctrl$age <- (current_data_ctrl$age - age_stats_ctrl['mean']) / age_stats_ctrl['sd']
current_data_ctrl$cd4 <- (current_data_ctrl$cd4 - cd4_stats_ctrl['mean']) / cd4_stats_ctrl['sd']
current_data_trt$age <- (current_data_trt$age - age_stats_trt['mean']) / age_stats_trt['sd']
current_data_trt$cd4 <- (current_data_trt$cd4 - cd4_stats_trt['mean']) / cd4_stats_trt['sd']
hist_data$age <- (hist_data$age - age_stats_hist['mean']) / age_stats_hist['sd']
hist_data$cd4 <- (hist_data$cd4 - cd4_stats_hist['mean']) / cd4_stats_hist['sd']

# set parameters
beta_pars <- elicit_beta_mean_cv(m0 = 0.5, v0 = 0.005)
a0_seq <- seq(0, 1, length.out = 21)
data <- list(current_data, hist_data)
data_ctrl <- list(current_data_ctrl, hist_data)
data_trt <- list(current_data_trt, hist_data)
family <- binomial(link = "logit")
c0 <- 3^0.5
d0 <- 0.5^0.5
delta0 <- beta_pars$a
lambda0 <- beta_pars$b
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

# run the model for each group
post_samples_ctrl <- samples_models(data = data_ctrl,
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

post_samples_trt <- samples_models(data = data_trt,
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

risk_diff <- mean_models_trt - mean_models_ctrl
risk_diff_onlytrt <- mean_models_trt[,2] - mean_models_ctrl[,2]
risk_diff_complete <- mean_models_trt[,15] - mean_models_ctrl[,15]

# compute BMA
# bma_ctrl <- bma(mean_models_ctrl, post_samples_ctrl$df_post, 10000)
# bma_trt <- bma(mean_models_trt, post_samples_trt$df_post, 10000)
bma_ctrl <- bma(mean_models_ctrl, post_samples$df_post, 10000)
bma_trt <- bma(mean_models_trt, post_samples$df_post, 10000)

bma_ard <- bma(risk_diff, post_samples$df_post, 10000)

# save results
save(post_samples, file = "bayesian_subset_selection/actg/data/post_samples.RData")
save(post_samples_ctrl, file = "bayesian_subset_selection/actg/data/post_samples_ctrl.RData")
save(post_samples_trt, file = "bayesian_subset_selection/actg/data/post_samples_trt.RData")
save(mean_models_ctrl, file = "bayesian_subset_selection/actg/data/mean_models_ctrl.RData")
save(mean_models_trt, file = "bayesian_subset_selection/actg/data/mean_models_trt.RData")
save(bma_ctrl, file = "bayesian_subset_selection/actg/data/bma_ctrl.RData")
save(bma_trt, file = "bayesian_subset_selection/actg/data/bma_trt.RData")
