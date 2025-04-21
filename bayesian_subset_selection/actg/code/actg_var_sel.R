setwd("C:/Users/Ezequiel/OneDrive - Fundacao Getulio Vargas - FGV/Grad MAp FGV/proj_hdbayes/hdbayes-tutorials/bayesian_subset_selection/actg")

# load libraries
library(hdbayes)
library(parallel)
library(matrixStats)

# load auxiliary functions
source("functions.R")

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
