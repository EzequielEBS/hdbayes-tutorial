setwd("C:/Users/Ezequiel/OneDrive - Fundacao Getulio Vargas - FGV/Grad MAp FGV/proj_hdbayes/hdbayes-tutorials/bayesian_subset_selection")

# load libraries
library(hdbayes)
library(lme4)
library(loo)


# load data
formulas <- readRDS("data/formulas.rds")
covariates_models <- readRDS("data/covariates_models.rds")

# load results
fits_bayesian_model <- readRDS("results/fits_bayesian_model.rds")
fits_glm <- readRDS("results/fits_glm.rds")

# compute metrics

# compute marginal likelihood
marginal_likelihoods <- lapply(fits_bayesian_model, function(x) glm.logml.post(x)$logml)
# compute AIC
AICs <- lapply(fits_glm, function(x) x$aic)
# compute BIC
BICs <- lapply(fits_glm, function(x) BIC(x))

# create data frame
metrics <- data.frame(
  variables = unlist(covariates_models),
  logml = unlist(marginal_likelihoods),
  AIC = unlist(AICs),
  BIC = unlist(BICs)
)

# sort by logml
metrics <- metrics[order(metrics$logml, decreasing = TRUE), ]

# save results
saveRDS(metrics, "results/metrics.rds")
write.csv(metrics, "results/metrics.csv", row.names = FALSE)
