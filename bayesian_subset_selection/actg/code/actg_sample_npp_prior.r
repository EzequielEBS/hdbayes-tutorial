source("bayesian_subset_selection/actg/code/aux_scripts/get_stan_data.r")

library(hdbayes)
library(parallel)
library(matrixStats)
library(MCMCpack)
library(bayestestR)
library(dplyr)
library(ggplot2)
library(patchwork)

current_data <- actg036
hist_data <- actg019

current_data$treatment <- current_data$treatment - mean(current_data$treatment)
current_data$race <- current_data$race - mean(current_data$race)
hist_data$treatment <- hist_data$treatment - mean(hist_data$treatment)
hist_data$race <- hist_data$race - mean(hist_data$race)

current_data$age <- (current_data$age - mean(current_data$age)) /
  (2*sd(current_data$age))
current_data$cd4 <- (current_data$cd4 - mean(current_data$cd4)) /
  (2*sd(current_data$cd4))
hist_data$age <- (hist_data$age - mean(hist_data$age)) /
  (2*sd(hist_data$age))
hist_data$cd4 <- (hist_data$cd4 - mean(hist_data$cd4)) /
  (2*sd(hist_data$cd4))

data <- list(current_data, hist_data)
family <- binomial(link = "logit")
delta0 <- 1
lambda0 <- 1
iter_warmup <- 1000
iter_sampling <- 2500
formula <- outcome ~ age + treatment + cd4
a0      <- seq(0, 1, length.out = 21) # for demonstration, change it to a large number in practice
## wrapper to obtain log normalizing constant in parallel package
logncfun.wip <- function(a0, ...){
  glm.npp.lognc.wip(
      formula = formula, family = family, a0 = a0, histdata = hist_data,
      ...
    ) 
}

ncores <- detectCores() - 1
a0.lognc.wip <- mclapply(
  X = a0, FUN = logncfun_cauchy, iter_warmup = 1000,
  iter_sampling = 2500, chains = 4,
  mc.cores = ncores
)
a0.lognc.wip <- data.frame( do.call(rbind, a0.lognc.wip) )

standat.wip <- get.stan.data.npp.prior.wip(
  formula        = formula,
  family         = family,
  data.list      = list(hist_data),
  a0.lognc       = a0.lognc.wip$a0,
  lognc          =  matrix(a0.lognc.wip$lognc, ncol = 1),
  offset.list    = NULL,
  disp.mean      = NULL,
  disp.sd        = NULL,
  a0.shape1      = delta0,
  a0.shape2      = lambda0,
  a0.lower       = 0,
  a0.upper       = 1,
  wip            = 1
)

## fit model in cmdstanr
fit.wip <- glm_npp_prior_wip$sample(
  data = standat.wip,
  iter_warmup = 1000, iter_sampling = 2500, chains = 4,
  parallel_chains = 4
)

pars <- c("beta[1]", "beta[2]", "beta[3]", "beta[4]")
draws.wip <- fit.wip$draws(format = 'draws_df', variables = pars) %>%
  select(all_of(pars))

save(draws.wip, file = "bayesian_subset_selection/actg/samples/draws_npp_wip.RData")

## sample normal

logncfun <- function(a0, ...){
  glm.npp.lognc(
      formula = formula, family = family, a0 = a0, histdata = hist_data,
      ...
    )
}

ncores <- detectCores() - 1
a0.lognc.norm <- mclapply(
  X = a0, FUN = logncfun, iter_warmup = 1000,
  iter_sampling = 2500, chains = 4,  beta.sd = 0.5,
  mc.cores = ncores
)
a0.lognc.norm <- data.frame( do.call(rbind, a0.lognc.norm) )

standat.norm <- get.stan.data.npp.prior.wip(
  formula        = formula,
  family         = family,
  data.list      = list(hist_data),
  a0.lognc       = a0.lognc.norm$a0,
  lognc          =  matrix(a0.lognc.norm$lognc, ncol = 1),
  offset.list    = NULL,
  beta.sd        = 0.5,
  disp.mean      = NULL,
  disp.sd        = NULL,
  a0.shape1      = delta0,
  a0.shape2      = lambda0,
  a0.lower       = 0,
  a0.upper       = 1,
  wip            = 0
)

## fit model in cmdstanr
fit.norm <- glm_npp_prior_wip$sample(
  data = standat.norm,
  iter_warmup = 1000, iter_sampling = 2500, chains = 4
)

pars <- c("beta[1]", "beta[2]", "beta[3]", "beta[4]")
draws.norm <- fit.norm$draws(format = 'draws_df', variables = pars) %>%
  select(all_of(pars))

save(draws.norm, file = "bayesian_subset_selection/actg/samples/draws_npp_norm.RData")
