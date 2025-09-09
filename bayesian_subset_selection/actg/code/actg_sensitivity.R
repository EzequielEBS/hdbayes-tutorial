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
a0_seq <- seq(0, 1, length.out = 21)
data <- list(current_data, hist_data)
data_ctrl <- list(current_data_ctrl, hist_data)
data_trt <- list(current_data_trt, hist_data)
family <- binomial(link = "logit")



iter_warmup <- 1000
iter_sampling <- 2500

#sensitivity analyses for c0
c0 <- c(5, 10, 50, 100)^(0.5)
beta_pars <- elicit_beta_mean_cv(m0 = 0.5, v0 = 0.008)
delta0 <- beta_pars$a
lambda0 <- beta_pars$b
d0 <- 5^0.5
post_samples_c0 <- lapply(c0, function(c0_val) {
  samples_models(data = data,
                 outcome = "outcome",
                 family = family,
                 a0_seq = a0_seq,
                 c0 = c0_val,
                 d0 = d0,
                 delta0 = delta0,
                 lambda0 = lambda0,
                 iter_warmup = iter_warmup,
                 iter_sampling = iter_sampling,
                 num_cores = 14)
})