############################################################################################################
# Propensity score matching (PSM) before dynamic borrowing
# Example: actg019 (external) & actg036 (current)
############################################################################################################# 

library(dbarts)   # BART
library(MatchIt)  # matching
library(hdbayes)  # example data
library(dplyr)
library(posterior)

# Load data
hist <- actg019 # external data
curr <- actg036 # current data

# Center and scale continuous covariates (age and cd4)
hist$age <- scale(hist$age, center = T, scale = T)
hist$cd4 <- scale(hist$cd4, center = T, scale = T)
curr$age <- scale(curr$age, center = T, scale = T)
curr$cd4 <- scale(curr$cd4, center = T, scale = T)

# Stack current and external data sets as a single data set
# and add study label (study = 1 for current; 0 for external)
data_all       <- rbind(curr, hist)
data_all$study <- rep(c(1, 0), times = c(nrow(curr), nrow(hist)))

# Fit probit BART to estimate propensity scores Pr(study = 1 | X)
bart_ps <- bart2(
  formula   = study ~ age + race + cd4,
  data      = data_all,
  n.samples = 100000,
  n.burn    = 50000,
  n.chains  = 4,
  n.threads = 4
)

# Obtain posterior mean of the estimated propensity scores (on the probability scale)
ps_hat <- fitted(bart_ps)  
# all.equal(apply(pnorm(bart_ps$yhat.train), 3, mean), ps_hat) # TRUE

# Take logit transform
logit <- function(p){
  log(p / (1 - p))
}
logit_ps <- logit(ps_hat)

data_all$ps_hat   <- ps_hat
data_all$logit_ps <- logit_ps

# Nearest-neighbor matching on the logit of the propensity scores (PS)
# Caliper is applied on the raw distance scale
m.out <- matchit(
  study ~ logit_ps,
  data = data_all,
  method = "nearest",
  distance = data_all$logit_ps,
  ratio = 1,
  caliper = 0.1,
  std.caliper = F
)
# Keep all current data subjects in the matched data
m.out$weights[m.out$weights == 0 & data_all$study == 1] <- 1
summary(m.out)

# Visual checks
plot(m.out, type = "qq")
plot(m.out, type = "density")

# Extract matched data sets
matched_data <- match.data(m.out)
# Extract matched external data set
hist.matched <- matched_data %>%
  filter(study == 0) %>%
  select(colnames(curr)) %>%
  as.data.frame() %>%
  mutate(
    age = as.numeric(age),
    cd4 = as.numeric(cd4)
  ) 

# Save matched external set for downstream dynamic borrowing analysis
saveRDS(hist.matched, file = "logistic_regression/data/actg019_after_PSM.rds")
