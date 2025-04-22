setwd("C:/Users/Ezequiel/OneDrive - Fundacao Getulio Vargas - FGV/MSc_MAp_CD/hdbayes_project/hdbayes-tutorial/bayesian_subset_selection/actg")

# load libraries
library(hdbayes)
library(parallel)
library(matrixStats)
library(MCMCpack)

# load auxiliary functions
source("code/functions.R")

# normalize data
age_stats <- with(actg036,
                  c('mean' = mean(age), 'sd' = sd(age)))
cd4_stats <- with(actg036,
                  c('mean' = mean(cd4), 'sd' = sd(cd4)))
actg036$age <- ( actg036$age - age_stats['mean'] ) / age_stats['sd']
actg019$age <- ( actg019$age - age_stats['mean'] ) / age_stats['sd']
actg036$cd4 <- ( actg036$cd4 - cd4_stats['mean'] ) / cd4_stats['sd']
actg019$cd4 <- ( actg019$cd4 - cd4_stats['mean'] ) / cd4_stats['sd']

# define data
current_data <- actg036
hist_data <- actg019

# split data
current_data_ctrl <- current_data[current_data$treatment == 0, 
                                !(names(current_data) %in% c("treatment"))]
current_data_trt <- current_data[current_data$treatment == 1, 
                                 !(names(current_data) %in% c("treatment"))]

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
                               num_cores = 10)

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
                                    num_cores = 10)

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
                                   num_cores = 10)

# save results
save(post_samples, file = "data/post_samples.RData")
save(post_samples_ctrl, file = "data/post_samples_ctrl.RData")
save(post_samples_trt, file = "data/post_samples_trt.RData")