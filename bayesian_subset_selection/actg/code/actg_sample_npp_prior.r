source("bayesian_subset_selection/actg/code/aux_scripts/get_stan_data.r")
source("bayesian_subset_selection/actg/code/aux_scripts/glm_npp_lognc_wip.r")
source("bayesian_subset_selection/actg/code/aux_scripts/glm_npp_wip.r")
source("bayesian_subset_selection/actg/code/aux_scripts/glm_logml_npp_wip.r")
source("bayesian_subset_selection/actg/code/aux_scripts/data_checks.r")
source("bayesian_subset_selection/actg/code/aux_scripts/expfam_loglik.r")

library(hdbayes)
library(parallel)
library(matrixStats)
library(MCMCpack)
library(bayestestR)
library(dplyr)

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
  X = a0, FUN = logncfun.wip, iter_warmup = 1000,
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

c0 <- c(0.25, 0.5, 1, 2)

logncfun <- function(a0, ...){
  glm.npp.lognc(
      formula = formula, family = family, a0 = a0, histdata = hist_data,
      ...
    )
}

a0      <- seq(0, 1, length.out = 21)
ncores <- detectCores() - 1
a0.lognc.c0 <- lapply(c0, function(c0_val) {
  a0.lognc <- mclapply(
    X = a0, FUN = logncfun, iter_warmup = iter_warmup,
    iter_sampling = iter_sampling, chains = 4,  beta.sd = c0_val,
    mc.cores = ncores
  )
  a0.lognc <- data.frame( do.call(rbind, a0.lognc) )
  return(a0.lognc)
  }
)

draws.c0 <- lapply(seq_along(c0), function(i) {
  c0_val <- c0[i]
  standat.norm <- get.stan.data.npp.prior.wip(
    formula        = formula,
    family         = family,
    data.list      = list(hist_data),
    a0.lognc       = a0.lognc.c0[[i]]$a0,
    lognc          =  matrix(a0.lognc.c0[[i]]$a0, ncol = 1),
    offset.list    = NULL,
    beta.sd        = c0_val,
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
    iter_warmup = iter_sampling, iter_sampling = iter_warmup, chains = 4
  )
  
  pars <- c("beta[1]", "beta[2]", "beta[3]", "beta[4]")
  draws.norm <- fit.norm$draws(format = 'draws_df', variables = pars) %>%
    select(all_of(pars))
})

save(draws.c0, file = "bayesian_subset_selection/actg/samples/draws_npp_c0.RData")

## sample for different a0 priors

a0_hyper <- list(c(1,1), c(2,2), c(1,10), c(10,1))
c0 <- 1

draws.a0 <- lapply(seq_along(a0_hyper), function(i) {
  delta0_val <- a0_hyper[[i]][1]
  lambda0_val <- a0_hyper[[i]][2]
  standat.norm <- get.stan.data.npp.prior.wip(
    formula        = formula,
    family         = family,
    data.list      = list(hist_data),
    a0.lognc       = a0.lognc.c0[[3]]$a0,
    lognc          =  matrix(a0.lognc.c0[[3]]$a0, ncol = 1),
    offset.list    = NULL,
    beta.sd        = c0,
    disp.mean      = NULL,
    disp.sd        = NULL,
    a0.shape1      = delta0_val,
    a0.shape2      = lambda0_val,
    a0.lower       = 0,
    a0.upper       = 1,
    wip            = 0
  )
  
  ## fit model in cmdstanr
  fit.norm <- glm_npp_prior_wip$sample(
    data = standat.norm,
    iter_warmup = iter_sampling, iter_sampling = iter_warmup, chains = 4
  )
  
  pars <- c("beta[1]", "beta[2]", "beta[3]", "beta[4]")
  draws.norm <- fit.norm$draws(format = 'draws_df', variables = pars) %>%
    select(all_of(pars))
})

save(draws.a0, file = "bayesian_subset_selection/actg/samples/draws_npp_a0.RData")

#-------------------------------------------------------------------------------
# After PSM
#-------------------------------------------------------------------------------

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
a0.lognc.wip_after_PSM <- mclapply(
  X = a0, FUN = logncfun.wip, iter_warmup = 1000,
  iter_sampling = 2500, chains = 4,
  mc.cores = ncores
)
a0.lognc.wip_after_PSM <- data.frame( do.call(rbind, a0.lognc.wip_after_PSM) )

standat.wip_after_PSM <- get.stan.data.npp.prior.wip(
  formula        = formula,
  family         = family,
  data.list      = list(hist_data),
  a0.lognc       = a0.lognc.wip_after_PSM$a0,
  lognc          =  matrix(a0.lognc.wip_after_PSM$lognc, ncol = 1),
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
fit.wip_after_PSM <- glm_npp_prior_wip$sample(
  data = standat.wip_after_PSM,
  iter_warmup = 1000, iter_sampling = 2500, chains = 4,
  parallel_chains = 4
)

pars <- c("beta[1]", "beta[2]", "beta[3]", "beta[4]")
draws.wip_after_PSM <- fit.wip_after_PSM$draws(format = 'draws_df', variables = pars) %>%
  select(all_of(pars))

save(draws.wip_after_PSM, 
     file = "bayesian_subset_selection/actg/samples/draws_npp_wip_after_PSM.RData")

## sample normal

c0 <- c(0.25, 0.5, 1, 2)

logncfun <- function(a0, ...){
  glm.npp.lognc(
    formula = formula, family = family, a0 = a0, histdata = hist_data,
    ...
  )
}

a0      <- seq(0, 1, length.out = 21)
ncores <- detectCores() - 1
a0.lognc.c0_after_PSM <- lapply(c0, function(c0_val) {
  a0.lognc <- mclapply(
    X = a0, FUN = logncfun, iter_warmup = iter_warmup,
    iter_sampling = iter_sampling, chains = 4,  beta.sd = c0_val,
    mc.cores = ncores
  )
  a0.lognc <- data.frame( do.call(rbind, a0.lognc) )
  return(a0.lognc)
}
)

draws.c0_after_PSM <- lapply(seq_along(c0), function(i) {
  c0_val <- c0[i]
  standat.norm <- get.stan.data.npp.prior.wip(
    formula        = formula,
    family         = family,
    data.list      = list(hist_data),
    a0.lognc       = a0.lognc.c0[[i]]$a0,
    lognc          =  matrix(a0.lognc.c0[[i]]$a0, ncol = 1),
    offset.list    = NULL,
    beta.sd        = c0_val,
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
    iter_warmup = iter_sampling, iter_sampling = iter_warmup, chains = 4
  )
  
  pars <- c("beta[1]", "beta[2]", "beta[3]", "beta[4]")
  draws.norm <- fit.norm$draws(format = 'draws_df', variables = pars) %>%
    select(all_of(pars))
})

save(draws.c0_after_PSM, 
     file = "bayesian_subset_selection/actg/samples/draws_npp_c0_after_PSM.RData")

## sample for different a0 priors

a0_hyper <- list(c(1,1), c(2,2), c(1,10), c(10,1))
c0 <- 1

draws.a0_after_PSM <- lapply(seq_along(a0_hyper), function(i) {
  delta0_val <- a0_hyper[[i]][1]
  lambda0_val <- a0_hyper[[i]][2]
  standat.norm <- get.stan.data.npp.prior.wip(
    formula        = formula,
    family         = family,
    data.list      = list(hist_data),
    a0.lognc       = a0.lognc.c0_after_PSM[[3]]$a0,
    lognc          =  matrix(a0.lognc.c0_after_PSM[[3]]$a0, ncol = 1),
    offset.list    = NULL,
    beta.sd        = c0,
    disp.mean      = NULL,
    disp.sd        = NULL,
    a0.shape1      = delta0_val,
    a0.shape2      = lambda0_val,
    a0.lower       = 0,
    a0.upper       = 1,
    wip            = 0
  )
  
  ## fit model in cmdstanr
  fit.norm <- glm_npp_prior_wip$sample(
    data = standat.norm,
    iter_warmup = iter_sampling, iter_sampling = iter_warmup, chains = 4
  )
  
  pars <- c("beta[1]", "beta[2]", "beta[3]", "beta[4]")
  draws.norm <- fit.norm$draws(format = 'draws_df', variables = pars) %>%
    select(all_of(pars))
})

save(draws.a0_after_PSM, 
     file = "bayesian_subset_selection/actg/samples/draws_npp_a0_after_PSM.RData")
