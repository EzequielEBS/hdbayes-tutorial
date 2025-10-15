#-------------------------------------------------------------------------------
# Sample models with Cauchy prior
#-------------------------------------------------------------------------------

# load libraries and functions
source("bayesian_subset_selection/actg/code/aux_scripts/glm_npp_lognc_wip.r")
source("bayesian_subset_selection/actg/code/aux_scripts/glm_npp_wip.r")
source("bayesian_subset_selection/actg/code/aux_scripts/glm_logml_npp_wip.r")
source("bayesian_subset_selection/actg/code/aux_scripts/glm_pp_wip.r")
source("bayesian_subset_selection/actg/code/aux_scripts/data_checks.r")
source("bayesian_subset_selection/actg/code/aux_scripts/expfam_loglik.r")
source("bayesian_subset_selection/actg/code/aux_scripts/functions.R")

library(hdbayes)
library(parallel)
library(matrixStats)
library(MCMCpack)
library(bayestestR)
library(dplyr)
library(ggplot2)
library(patchwork)

# define data
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
delta0 <- 1
lambda0 <- 1
iter_warmup <- 1000
iter_sampling <- 2500

# create list of formulas
covariates <- colnames(current_data)[colnames(current_data) != "outcome"]
pset <- powerset(covariates)
formulas <- lapply(pset, create_formula, outcome = "outcome")

# get covariates from models
covariates_models <- lapply(pset, get_covariates)

# functions to fit models
logp0 <- function(formula) {
  glm.npp.lognc.wip(
    formula,
    a0 = 1,
    family = family,
    histdata = hist_data,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    chains = 4,
    refresh = 0
  )
}

logncfun.wip <- function(a0, formula) {
  glm.npp.lognc.wip(
    a0 = a0,
    formula = formula,
    family = family, 
    histdata = hist_data, 
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling, 
    chains = 4, 
    refresh = 0
  ) 
}

post_beta <- function(formula) {
  a0.lognc <- lapply(a0_seq, 
                     logncfun.wip, 
                     formula = formula
  )
  a0.lognc <- data.frame(do.call(rbind, a0.lognc))
  fit <- glm.npp.wip(formula, 
                          data = list(current_data,
                                      hist_data),
                          family = family,
                          a0.lognc = a0.lognc$a0,
                          lognc = matrix(a0.lognc$lognc, ncol = 1),
                          a0.shape1 = delta0,
                          a0.shape2 = lambda0,
                          iter_warmup = iter_warmup,
                          iter_sampling = iter_sampling, 
                          chains = 4, 
                          refresh = 0
  )
}


# fit models in parallel
logp0_models <- mclapply(
  X = formulas,
  FUN = logp0,
  mc.cores = 14
)

post_betam  <- mclapply(
  X = formulas,
  FUN = post_beta,
  mc.cores = 14
)

logml_models <- mclapply(
  X = post_betam,
  FUN = glm.logml.npp.wip,
  mc.cores = 14
)

prior_models <- exp(sapply(logp0_models, function(x) x[2]) -
                      logSumExp(sapply(logp0_models, function(x) x[2]))
)

post_models <- post_models(logp0_models, logml_models)

df_post <- data.frame(model = unlist(covariates_models),
                      prior_model = prior_models,
                      ml = exp(unlist(data.frame(do.call(rbind, logml_models))$logml)),
                      post_model = post_models
)

df_post_ord <- df_post[order(df_post$post_model, decreasing = TRUE),]
rownames(df_post_ord) <- 1:nrow(df_post_ord)

post_samples_wip <- list(logp0_models = logp0_models,
     post_betam = post_betam,
     logml_models = logml_models,
     df_post = df_post,
     df_post_ord = df_post_ord)

save(post_samples_wip, 
     file = "bayesian_subset_selection/actg/samples/post_samples_wip.RData")

xtable::xtable(df_post_ord, digits = 3)


#-------------------------------------------------------------------------------
# Sample after PSM 
#-------------------------------------------------------------------------------

current_data <- actg036
hist_data <- readRDS("bayesian_subset_selection/actg/data/actg019_after_PSM.rds")

# normalize data
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
iter_warmup <- 1000
iter_sampling <- 2500

# create list of formulas
covariates <- colnames(current_data)[colnames(current_data) != "outcome"]
pset <- powerset(covariates)
formulas <- lapply(pset, create_formula, outcome = "outcome")

# get covariates from models
covariates_models <- lapply(pset, get_covariates)

logp0 <- function(formula) {
  glm.npp.lognc.wip(
    formula,
    a0 = 1,
    family = family,
    histdata = hist_data,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    chains = 4,
    refresh = 0
  )
}

logncfun.wip <- function(a0, formula) {
  glm.npp.lognc.wip(
    a0 = a0,
    formula = formula,
    family = family, 
    histdata = hist_data, 
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling, 
    chains = 4, 
    refresh = 0
  ) 
}

post_beta <- function(formula) {
  a0.lognc <- lapply(a0_seq, 
                     logncfun.wip, 
                     formula = formula
  )
  a0.lognc <- data.frame(do.call(rbind, a0.lognc))
  fit <- glm.npp.wip(formula, 
                     data = list(current_data,
                                 hist_data),
                     family = family,
                     a0.lognc = a0.lognc$a0,
                     lognc = matrix(a0.lognc$lognc, ncol = 1),
                     a0.shape1 = delta0,
                     a0.shape2 = lambda0,
                     iter_warmup = iter_warmup,
                     iter_sampling = iter_sampling, 
                     chains = 4, 
                     refresh = 0
  )
}



logp0_models <- mclapply(
  X = formulas,
  FUN = logp0,
  mc.cores = 14
)

formula <- formulas[[1]]
a0.lognc <- lapply(a0_seq, 
                   logncfun.wip, 
                   formula = formula
)
a0.lognc <- data.frame(do.call(rbind, a0.lognc))
fit <- glm.npp.wip(formula, 
                   data = list(current_data,
                               hist_data),
                   family = family,
                   a0.lognc = a0.lognc$a0,
                   lognc = matrix(a0.lognc$lognc, ncol = 1),
                   a0.shape1 = delta0,
                   a0.shape2 = lambda0,
                   iter_warmup = iter_warmup,
                   iter_sampling = iter_sampling, 
                   chains = 4, 
                   refresh = 0
)

post_betam  <- mclapply(
  X = formulas,
  FUN = post_beta,
  mc.cores = 14
)

logml_models <- mclapply(
  X = post_betam,
  FUN = glm.logml.npp.wip,
  mc.cores = 14
)

prior_models <- exp(sapply(logp0_models, function(x) x[2]) -
                      logSumExp(sapply(logp0_models, function(x) x[2]))
)

post_models <- post_models(logp0_models, logml_models)

df_post <- data.frame(model = unlist(covariates_models),
                      prior_model = prior_models,
                      ml = exp(unlist(data.frame(do.call(rbind, logml_models))$logml)),
                      post_model = post_models
)

df_post_ord <- df_post[order(df_post$post_model, decreasing = TRUE),]
rownames(df_post_ord) <- 1:nrow(df_post_ord)

post_samples_wip_after_PSM <- list(logp0_models = logp0_models,
                         post_betam = post_betam,
                         logml_models = logml_models,
                         df_post = df_post,
                         df_post_ord = df_post_ord)

save(post_samples_wip_after_PSM, 
     file = "bayesian_subset_selection/actg/samples/post_samples_wip_after_PSM.RData")

xtable::xtable(df_post_ord, digits = 3)
