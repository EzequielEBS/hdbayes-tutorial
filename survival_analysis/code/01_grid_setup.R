###################################################################################
# Set up a grid of model and prior combinations for running the real data analysis
# on the computing cluster.
###################################################################################

library(dplyr)

# Directory to save the grid
save.dir <- 'code'

# Specify models and priors to be evaluated in the analysis
# "PWE" corresponds to the PWEPH model and "CurePWE" to the CurePWEPH model in the manuscript
model <- c("PWE", "CurePWE")
# The analysis compares multiple prior specifications, including the vague prior ("ref"), 
# power prior ("pp"), propensity score–integrated power prior ("psipp"), Bayesian hierarchical model ("bhm"),
# commensurate prior ("cp"), latent exchangeability prior ("leap"), and normalized power prior ("npp").
prior <- c('ref', 'pp', 'psipp', 'bhm', 'cp', 'leap', 'npp')
# Set the number of intervals (J) for the PWEPH and CurePWEPH model
# We use values from 2 to 9. The upper limit (J = 9) ensures that each interval contains at least one event
# in both the stratified current and external data sets.
J      <- 2:9

# We will stratify the current and external data sets by treatment arm for analysis
arm    <- c(0, 1) 

grid   <- expand.grid(
  J = J, arm = arm, prior = prior, model = model
  , stringsAsFactors = FALSE
)
# Save the grid for downstream use 
saveRDS(grid, file = file.path(save.dir, "grid.rds"))
