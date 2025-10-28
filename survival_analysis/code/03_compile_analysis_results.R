############################################################################################################
# Compile analysis results obtained from 02_fit_PWE_CurePWE.R
# Goal:
#   - Load per-scenario fits (by model/prior/J/arm)
#   - Compute two-year relapse-free survival probabilities for each arm
#   - Aggregate model-comparison metric (elpd))
#   - Save a single compiled results object for downstream plotting/tables
############################################################################################################# 

library(dplyr)
library(posterior)
library(loo)
library(hdbayes)

# Wrapper functions to compute two-year relapse-free survival probabilities per arm
wrapper.dir <- "code/wrapper"
source(file.path(wrapper.dir, "pwe_get_trt_effect.R"))
source(file.path(wrapper.dir, "curepwe_get_trt_effect.R"))

# Load scenario grid (model / prior / J / arm)
# We compile results by pairing each control-arm scenario with its matching
# treatment-arm scenario (same model, prior, and J). Hence we keep only arm == 0 here.
grid        <- readRDS("code/grid.rds") %>%
  filter( arm == 0 ) %>%
  dplyr::select(!arm)

# Locate fitted result files
results.dir <- "results"
file.list   <- list.files(results.dir, pattern = ".rds")

# Load current trial data (E1690)
curr <- E1690
# Replace 0-day failure times with 0.50 days (converted to years)
curr <- curr %>% mutate(failtime = if_else(failtime == 0, 0.50/365.25, failtime))
# Center and scale age
curr$cage <- scale(curr$age, center = T, scale = T)

surv.trt.list  <- list()
surv.ctl.list  <- list()
surv.diff.list <- list()
elpd.trt.vals  <- vector(length = nrow(grid))
elpd.ctl.vals  <- vector(length = nrow(grid))

for (i in 1:nrow(grid)) {
  grid.id <- grid[i, ]
  prior   <- grid.id$prior
  model   <- grid.id$model

  # Treatment arm
  pattern <- paste0('_', model, "_", prior, "_arm_", 1, '_nintervals_', grid.id$J)
  res.i   <- readRDS(file.path(results.dir, file.list[ grep(pattern, file.list) ]))
  d       <- res.i$draws
  
  # Store elpd value
  elpd.trt.vals[i] <- res.i$res.loo$estimates[1, 1]

  # Predict two-year relapse-free survival probability
  if ( model == "PWE"){
    if ( prior == "psipp" ){
      d.surv <- get.surv.prob.pwe.psipp(
        t = 2,
        post.samples = d
      )

    }else{
      d.surv <- get.surv.prob.pwe(
        t = 2,
        post.samples = d,
        data = curr
      )
    }

  }else{
    if ( prior == "psipp" ){
      d.surv <- get.surv.prob.curepwe.psipp(
        t = 2,
        post.samples = d
      )

    }else{
      d.surv <- get.surv.prob.curepwe(
        t = 2,
        post.samples = d,
        data = curr
      )

    }
  }
  surv.trt.list[[i]] <- d.surv

  # Control arm
  pattern <- paste0('_', model, "_", prior, "_arm_", 0, '_nintervals_', grid.id$J)
  res.i   <- readRDS(file.path(results.dir, file.list[ grep(pattern, file.list) ]))
  d       <- res.i$draws
  
  # Store elpd value
  elpd.ctl.vals[i]  <- res.i$res.loo$estimates[1, 1]

  if ( model == "PWE"){
    if ( prior == "psipp" ){
      d.surv <- get.surv.prob.pwe.psipp(
        t = 2,
        post.samples = d
      )

    }else{
      d.surv <- get.surv.prob.pwe(
        t = 2,
        post.samples = d,
        data = curr
      )
    }

  }else{
    if ( prior == "psipp" ){
      d.surv <- get.surv.prob.curepwe.psipp(
        t = 2,
        post.samples = d
      )

    }else{
      d.surv <- get.surv.prob.curepwe(
        t = 2,
        post.samples = d,
        data = curr
      )

    }
  }
  surv.ctl.list[[i]] <- d.surv
  # Compute difference in two-year relapse-free survival probability
  d.surv.diff         <- surv.trt.list[[i]] - surv.ctl.list[[i]]
  surv.diff.list[[i]] <- d.surv.diff
  # Posterior summaries 
  estim.i <- c(mean = mean(d.surv.diff), sd = sd(d.surv.diff),
               quantile2(d.surv.diff, probs = c(0.5, 0.025, 0.975)),
               prob_greater_0 = mean(d.surv.diff > 0))
  
  prior.i <- prior
  model.i <- paste0(model, " (J = ", grid.id$J, ")")

  if( i == 1 ){
    estim  <- estim.i
    priors <- prior.i
    models <- model.i
  }else{
    estim  <- rbind(estim, estim.i)
    priors <- c(priors, prior.i)
    models <- c(models, model.i)
  }
  print( paste0("######################## Completed iteration ", i, " #########################"))
}

res <- list(
  surv.trt.list  = surv.trt.list,
  surv.ctl.list  = surv.ctl.list,
  surv.diff.list = surv.diff.list,
  elpd.trt.vals  = elpd.trt.vals,
  elpd.ctl.vals  = elpd.ctl.vals,
  estim = estim,
  priors = priors,
  models = models
)

model_order  <- c(paste0("CurePWE (J = ", 2:9, ")"), paste0("PWE (J = ", 2:9, ")"))
# Table of elpd values for treatment arm
elpd.tab.trt <- data.frame(
  model = res$models,
  prior = res$priors,
  elpd  = res$elpd.trt.vals
) %>%
  mutate(model = factor(model, levels = model_order)) %>%
  arrange(prior, model)
rownames(elpd.tab.trt) <- NULL

# Table of elpd values for control arm
elpd.tab.ctl <- data.frame(
  model = res$models,
  prior = res$priors,
  elpd  = res$elpd.ctl.vals
) %>%
  mutate(model = factor(model, levels = model_order)) %>%
  arrange(prior, model)
rownames(elpd.tab.ctl) <- NULL

# Save compiled analysis results
res.all <- c(res,
             list(elpd.tab.trt = elpd.tab.trt,
                  elpd.tab.ctl = elpd.tab.ctl)
)
save(res.all, file = "results/compiled_results/compiled_analysis_results.rds")
