source("bayesian_subset_selection/actg/code/aux_scripts/glm_npp_lognc_wip.r")
source("bayesian_subset_selection/actg/code/aux_scripts/glm_npp_wip.r")
source("bayesian_subset_selection/actg/code/aux_scripts/glm_logml_npp_wip.r")
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

logncfun <- function(a0,
                     formula){
  if (a0 == 0) {
    pm <- ncol(model.matrix(formula, data = hist_data))
    c(
      'a0'           = a0,
      'lognc'        = -pm * log(pi),
      'min_ess_bulk' = 1000,
      'max_Rhat'     = 0
    )
  }
  else {
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
}

post_beta <- function(formula) {
  a0.lognc <- lapply(a0_seq, 
                     logncfun, 
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
df_post_ord$ml <- log(df_post_ord$ml)

post_samples_cauchy <- list(logp0_models = logp0_models,
     post_betam = post_betam,
     logml_models = logml_models,
     df_post = df_post,
     df_post_ord = df_post_ord)

save(post_samples_cauchy, 
     file = "bayesian_subset_selection/actg/samples/post_samples_cauchy.RData")

xtable::xtable(df_post_ord, digits = 3)

load("bayesian_subset_selection/actg/data/post_samples_d0.RData")

# Plot using ggplot2
ate_cauchy <- plot_ate(post_samples_cauchy)

ate_cauchy + ggtitle("Cauchy prior on coefficients")
ggsave("bayesian_subset_selection/actg/results/figures/ate_cauchy.png",
       ate_cauchy, width = 10, height = 7, units = "in", dpi = 300)

ate_d0 <- lapply(post_samples_d0, function(x) {
  plot_ate(x)
})

all_ate_d0 <- (ate_d0[[1]] + ggtitle(expression(psi[0] == 0.25))+ xlim(-0.1, 0.05) + xlab("")) /
(ate_d0[[2]] + ggtitle(expression(psi[0] == 0.5)) + xlim(-0.1, 0.05) + xlab("")) /
(ate_d0[[3]] + ggtitle(expression(psi[0] == 1)) + xlim(-0.1, 0.05) +xlab("")) /
(ate_d0[[4]] + ggtitle(expression(psi[0] == 2)) + xlim(-0.1, 0.05)) +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')
all_ate_d0

ggsave("bayesian_subset_selection/actg/results/figures/ate_psi0.png",
       all_ate_d0, width = 10, height = 25, units = "in", dpi = 300)
