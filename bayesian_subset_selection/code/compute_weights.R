setwd("C:/Users/Ezequiel/OneDrive - Fundacao Getulio Vargas - FGV/Grad MAp FGV/proj_hdbayes/hdbayes-tutorials/bayesian_subset_selection")

# load libraries
library(matrixStats)
library(tidyverse)

# load data
metrics <- read_csv("results/metrics.csv")
  
# compute weights
weights_logml <- exp(unlist(metrics$logml) - logSumExp(unlist(metrics$logml)))
weights_AIC <- unlist(metrics$AIC) / sum(metrics$AIC)
weights_BIC <- unlist(metrics$BIC) / sum(unlist(metrics$BIC))

# create a data frame with the weights
weights <- data.frame(
  variables = metrics$variables,
  logml = weights_logml,
  AIC = weights_AIC,
  BIC = weights_BIC
)

# save the weights
write.csv(weights, "results/weights.csv", row.names = FALSE)
