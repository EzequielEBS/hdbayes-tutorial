setwd("C:/Users/Ezequiel/OneDrive - Fundacao Getulio Vargas - FGV/Grad MAp FGV/proj_hdbayes/hdbayes-tutorials/bayesian_subset_selection")

# load libraries
library(tidyverse)

# load data
data <- read_csv("data/chapman.csv")

# function to compute power set
powerset <- function(x) {
  sets <- lapply(1:(length(x)), function(i) combn(x, i, simplify = F))
  unlist(sets, recursive = F)
}

# function to create formula
create_formula <- function(x) {
  as.formula(paste("CNT ~ ", paste(x, collapse = "+")))
}

# function to get covariates from models
get_covariates <- function(x) {
  paste(x, collapse = ", ")
}

# create list of formulas
covariates <- colnames(data)[1:ncol(data)-1]
pset <- powerset(covariables)
formulas <- lapply(pset, create_formula)

# get covariates from models
covariates_models <- lapply(pset, get_covariates)

# save formulas
saveRDS(formulas, "data/formulas.rds")
saveRDS(covariates_models, "data/covariates_models.rds")
