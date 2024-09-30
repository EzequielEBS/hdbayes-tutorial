setwd("C:/Users/Ezequiel/OneDrive - Fundacao Getulio Vargas - FGV/Grad MAp FGV/proj_hdbayes/hdbayes-tutorials/bayesian_subset_selection")

# load libraries
library(hdbayes)
library(tidyverse)
library(matrixStats)

# load data
data <- read_csv("data/chapman.csv")

# set parameters
ncores <- max(1, parallel::detectCores() - 2)

# load formulas
formulas <- readRDS("data/formulas.rds")

# run hdbayes
fits <- lapply(formulas, function(x) 
  glm.post(x, family = binomial(link = "logit"), data = list(data), parallel_chains = ncores))

# save results
saveRDS(fits, "results/fits_bayesian_model.rds")

# compute marginal likelihood
marginal_likelihoods <- lapply(fits, function(x) glm.logml.post(x)$logml)

# save results
saveRDS(marginal_likelihoods, "results/marginal_likelihoods_bayesian_model.rds")

# get_max_id <- function(i) {
#   which(marginal_likelihoods == max(unlist(marginal_likelihoods[lengths(pset) == i])))
# }
# 
# # get best models
# best_models <- lapply(c(1:6), 
#                       function(i) list(formula = formulas[[get_max_id(i)]], 
#                                        logml = marginal_likelihoods[[get_max_id(i)]]))
# 
# df_best_models <- data.frame(matrix(unlist(best_models), ncol = 2, byrow = T))
# colnames(df_best_models) <- c("formula", "logml")
# 
# # compute weights
# df_best_models$weight <- exp(unlist(df_best_models[,c("logml")]) - logSumExp(df_best_models$logml))
