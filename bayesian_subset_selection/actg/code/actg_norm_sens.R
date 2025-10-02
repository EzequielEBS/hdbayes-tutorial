# load libraries
library(hdbayes)
library(parallel)
library(matrixStats)
library(MCMCpack)
library(ggplot2)

# load auxiliary functions
source("bayesian_subset_selection/actg/code/aux_scripts/functions.R")

current_data <- actg036
hist_data <- actg019

current_data$age <- (current_data$age - mean(current_data$age)) /
  (2*sd(current_data$age))
current_data$cd4 <- (current_data$cd4 - mean(current_data$cd4)) /
  (2*sd(current_data$cd4))
hist_data$age <- (hist_data$age - mean(hist_data$age)) /
  (2*sd(hist_data$age))
hist_data$cd4 <- (hist_data$cd4 - mean(hist_data$cd4)) /
  (2*sd(hist_data$cd4))

# set parameters
a0_seq <- seq(0, 1, length.out = 21)
data <- list(current_data, hist_data)
family <- binomial(link = "logit")

iter_warmup <- 1000
iter_sampling <- 2500

#sensitivity analyses for c0
c0 <- c(0.25, 0.5, 1, 2)
d0 <- c(0.25, 0.5, 1, 2)
delta0 <- 1
lambda0 <- 1
post_samples_c0d0 <- lapply(seq_along(c0), function(i) {
  c0_val <- c0[i]
  d0_val <- d0[i]
  samples_models(data = data,
                 outcome = "outcome",
                 family = family,
                 a0_seq = a0_seq,
                 c0 = c0_val,
                 d0 = d0_val,
                 delta0 = delta0,
                 lambda0 = lambda0,
                 iter_warmup = iter_warmup,
                 iter_sampling = iter_sampling,
                 num_cores = 14)
})

save(post_samples_c0d0, 
     file = "bayesian_subset_selection/actg/samples/post_samples_c0d0.RData")

load("bayesian_subset_selection/actg/samples/post_samples_c0d0.RData")

best_model_c0d0 <- do.call(rbind,lapply(post_samples_c0d0, function(post) {
  post$df_post_ord[1,]
}
)
)
best_model_c0d0$c0 <- c(0.25, 0.5, 1, 2)
best_model_c0d0$d0 <- c(0.25, 0.5, 1, 2)
best_model_c0d0$ml <- log(best_model_c0d0$ml)

xtable::xtable(best_model_c0d0, digits = 3)


# After PSM

current_data <- actg036
hist_data <- readRDS("bayesian_subset_selection/actg/data/actg019_after_PSM.rds")

current_data$age <- (current_data$age - mean(current_data$age)) /
  (2*sd(current_data$age))
current_data$cd4 <- (current_data$cd4 - mean(current_data$cd4)) /
  (2*sd(current_data$cd4))
hist_data$age <- (hist_data$age - mean(hist_data$age)) /
  (2*sd(hist_data$age))
hist_data$cd4 <- (hist_data$cd4 - mean(hist_data$cd4)) /
  (2*sd(hist_data$cd4))

# set parameters
a0_seq <- seq(0, 1, length.out = 21)
data <- list(current_data, hist_data)
family <- binomial(link = "logit")

iter_warmup <- 1000
iter_sampling <- 2500

#sensitivity analyses for c0
c0 <- c(0.25, 0.5, 1, 2)
d0 <- c(0.25, 0.5, 1, 2)
delta0 <- 1
lambda0 <- 1
post_samples_c0d0_after_PSM <- lapply(seq_along(c0), function(i) {
  c0_val <- c0[i]
  d0_val <- d0[i]
  samples_models(data = data,
                 outcome = "outcome",
                 family = family,
                 a0_seq = a0_seq,
                 c0 = c0_val,
                 d0 = d0_val,
                 delta0 = delta0,
                 lambda0 = lambda0,
                 iter_warmup = iter_warmup,
                 iter_sampling = iter_sampling,
                 num_cores = 14)
})

save(post_samples_c0d0_after_PSM, 
     file = "bayesian_subset_selection/actg/samples/post_samples_c0d0_after_PSM.RData")

load("bayesian_subset_selection/actg/samples/post_samples_c0d0_after_PSM.RData")

best_model_c0d0_after_PSM <- do.call(rbind,lapply(post_samples_c0d0_after_PSM, function(post) {
  post$df_post_ord[1,]
}
)
)
best_model_c0d0_after_PSM$c0 <- c(0.25, 0.5, 1, 2)
best_model_c0d0_after_PSM$d0 <- c(0.25, 0.5, 1, 2)
best_model_c0d0_after_PSM$ml <- log(best_model_c0d0_after_PSM$ml)



