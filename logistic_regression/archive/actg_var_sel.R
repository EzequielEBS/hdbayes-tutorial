# load libraries
library(hdbayes)
library(parallel)
library(matrixStats)

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


# construct Beta prior with mean 0.5 and variance 0.005
beta_pars <- elicit_beta_mean_cv(m0 = 0.5, v0 = 0.005)
# define values to compute Z(a_0)
a0_seq <- seq(0, 1, length.out = 21)
# define data
data <- list(current_data, hist_data)
# define family
family <- binomial(link = "logit")
# define values for c_0
c0_seq <- sqrt(c(3, 5 , 10))
# define d_0
d0 <- 0.5^0.5
# define Beta prior parameters
delta0 <- beta_pars$a
lambda0 <- beta_pars$b
# define number of iterations
iter_warmup <- 1000
iter_sampling <- 2500

# run the model for each value of c_0
post_samples_c01 <- samples_models(data = data,
                                   outcome = "outcome",
                                   family = family,
                                   a0_seq = a0_seq,
                                   c0 = c0_seq[1],
                                   d0 = d0,
                                   delta0 = delta0,
                                   lambda0 = lambda0,
                                   iter_warmup = iter_warmup,
                                   iter_sampling = iter_sampling,
                                   num_cores = 10)
post_samples_c02 <- samples_models(data = data,
                                   outcome = "outcome",
                                   family = family,
                                   a0_seq = a0_seq,
                                   c0 = c0_seq[2],
                                   d0 = d0,
                                   delta0 = delta0,
                                   lambda0 = lambda0,
                                   iter_warmup = iter_warmup,
                                   iter_sampling = iter_sampling,
                                   num_cores = 10)
post_samples_c03 <- samples_models(data = data,
                                   outcome = "outcome",
                                   family = family,
                                   a0_seq = a0_seq,
                                   c0 = c0_seq[3],
                                   d0 = d0,
                                   delta0 = delta0,
                                   lambda0 = lambda0,
                                   iter_warmup = iter_warmup,
                                   iter_sampling = iter_sampling,
                                   num_cores = 10)

# print the results
post_samples_c01$df_post
post_samples_c02$df_post
post_samples_c03$df_post
