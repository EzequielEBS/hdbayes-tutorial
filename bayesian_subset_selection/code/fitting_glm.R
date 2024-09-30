setwd("C:/Users/Ezequiel/OneDrive - Fundacao Getulio Vargas - FGV/Grad MAp FGV/proj_hdbayes/hdbayes-tutorials/bayesian_subset_selection")

# load libraries
library(tidyverse)

# load data
data <- read_csv("data/chapman.csv")

# load formulas
formulas <- readRDS("data/formulas.rds")

# run glm
fits <- lapply(formulas, function(x) 
  glm(x, family = binomial(link = "logit"), data = data))

# save fits
saveRDS(fits, "results/fits_glm.rds")
