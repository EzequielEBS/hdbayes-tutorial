############################################################################################################
# Analysis of E1690 study (current data) while incorporating information from E1684 study (external data)
# Models: PWEPH ("PWE") and CurePWEPH ("CurePWE") models with the number of intervals (J) for the piecewise
# constant baseline hazards varying from 2 to 9.
# Priors: vague prior ("ref"), power prior ("pp"), propensity score–integrated power prior ("psipp"), 
# Bayesian hierarchical model ("bhm"), commensurate prior ("cp"), latent exchangeability prior ("leap"), 
# and normalized power prior ("npp")
############################################################################################################# 

library(hdbayes)
library(dplyr)
library(survival)
library(loo)
library(posterior)

# MCMC/Sampling parameters
iter_warmup   <- 5000
iter_sampling <- 25000
chains        <- 1

# Load scenario grid and pick up the row corresponding to this SLURM task
grid <- readRDS("code/grid.rds")
## Get scenario for current ID
id <- as.numeric( Sys.getenv('SLURM_ARRAY_TASK_ID') )
if ( is.na(id) )
  id <- 1

# Directory to save results
save.dir <- 'results'
# Extract scenario-specific settings
scen     <- grid[id, ]
prior.id <- scen$prior # prior type (ref/pp/psipp/bhm/cp/leap/npp)
model.id <- scen$model # model type (PWE/CurePWE)
J.id     <- scen$J     # number of intervals for the PWE/CurePWE model
arm.id   <- scen$arm   # treatment arm to analyze (0 = control, 1 = treatment)

# Construct a file name
filename <- file.path(save.dir,
                      paste0('id_', id, '_', model.id, "_", prior.id, "_arm_", arm.id,
                             '_nintervals_', J.id, '.rds'))

# Load and pre-process data
hist <- E1684 # external data
curr <- E1690 # current data

# Replace 0-day failure times with 0.50 days (converted to years)
hist <- hist %>% mutate(failtime = if_else(failtime == 0, 0.50/365.25, failtime)) 
curr <- curr %>% mutate(failtime = if_else(failtime == 0, 0.50/365.25, failtime)) 

# Center and scale age
hist$cage <- scale(hist$age, center = T, scale = T)
curr$cage <- scale(curr$age, center = T, scale = T)

# Analyze one treatment arm at a time (stratified analysis)
curr <- curr %>% filter(treatment == arm.id)
hist <- hist %>% filter(treatment == arm.id)
data.list <- list(curr, hist)

# Main survival formula for PWE/CurePWE (covariates aligned with manuscript)
fmla       <- survival::Surv(failtime, failcens) ~ sex + cage + node_bin
# for psipp, the outcome model is intercept-only;
# the covariates enter via the model for estimating the propensity scores (PS)
fmla.psipp <- survival::Surv(failtime, failcens) ~ 1

# Obtain cut points for discretizing time intervals
# This yields approximately equal numbers of events per interval
nbreaks <- J.id
probs   <- 1:nbreaks / nbreaks
breaks  <- curr %>%
  filter(failcens == 1) %>%
  reframe(quant = quantile(failtime, probs = probs)) %>%
  unlist
breaks  <- c(0, breaks)
breaks[length(breaks)] <- max(10000, 1000 * breaks[length(breaks)])

# Propensity score–integrated power prior (PSIPP) setup (only when prior.id == "psipp")
if( prior.id == "psipp" ){
  # Wrapper function for implementing psipp
  wrapper.dir <- "code/wrapper"
  source(file.path(wrapper.dir, "get_strata_data.R"))
  
  # Covariates used to compute PS and derive strata for borrowing
  ps.formula <- ~ sex + cage + node_bin
  # Number of strata and target total external subjects to borrow
  nStrata    <- 4
  nBorrow    <- nrow(hist)
  
  # Create PS-based strata  
  res.strata <- get.strata.data(
    data = curr, histdata = hist,
    ps_fml_covs     = ps.formula,
    v_arm           = NULL,
    ctl_arm_level   = NULL,
    borrow_ctl_only = FALSE,
    nstrata         = nStrata,
    total_borrow    = nBorrow
  )
  # Optional diagnostics:
  # plot(res.strata$res.psrwe.est, "balance")  # covariate balance check
  
  data.list.PSIPP   <- res.strata$data.list
  strata.list       <- res.strata$strata.list
}

################################ Fit model for the selected scenario #############################
if ( model.id == "PWE" ){
  
  #################################### PWE model ####################################
  
  if ( prior.id == "psipp" ){
    # PSIPP 
    d <- pwe.stratified.pp(
      formula = fmla.psipp,
      data.list = data.list.PSIPP,
      strata.list = strata.list,
      breaks = breaks,
      a0.strata = res.strata$a0.strata,
      beta.mean = 0, beta.sd = 10,
      base.hazard.mean = 0, base.hazard.sd = 10,
      get.loglik = T,
      chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
    )

  }else {
    if ( prior.id == "ref" ){
      # Vague prior (using current data only)
      d <- pwe.post(
        formula = fmla,
        data.list = list(curr),
        breaks = breaks,
        beta.mean = 0, beta.sd = 10,
        base.hazard.mean = 0, base.hazard.sd = 10,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )

    }else if( prior.id == "pp" ){
      # Power prior with fixed a0 = 0.5
      d <- pwe.pp(
        formula = fmla,
        data.list = data.list,
        breaks = breaks,
        a0 = 0.5,
        beta.mean = 0, beta.sd = 10,
        base.hazard.mean = 0, base.hazard.sd = 10,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )
      
    }else if( prior.id == "bhm" ){
      # Bayesian hierarchical model
      d <- pwe.bhm(
        formula = fmla,
        data.list = data.list,
        breaks = breaks,
        meta.mean.mean = 0, meta.mean.sd = 10,
        meta.sd.mean = 0, meta.sd.sd = 0.5,
        base.hazard.mean = 0, base.hazard.sd = 10,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )

    }else if( prior.id == "cp" ){
      # Commensurate prior 
      d <- pwe.commensurate(
        formula = fmla,
        data.list = data.list,
        breaks = breaks,
        beta0.mean = 0, beta0.sd = 10,
        p.spike = 0.1,
        spike.mean = 200, spike.sd = 0.1,
        slab.mean = 0, slab.sd = 5,
        base.hazard.mean = 0, base.hazard.sd = 10,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )

    }else if( prior.id == "npp" ){
      # Normalized power prior 
      # Pre-compute log normalizing constants Z(a0)
      a0 <- seq(0, 1, length.out = 101)

      # Wrapper function to compute log Z(a0) using the external data
      logncfun <- function(a0, ...){
        hdbayes::pwe.npp.lognc(
          formula = fmla, histdata = data.list[[2]], breaks = breaks, a0 = a0,
          beta.mean = 0, beta.sd = 10,
          base.hazard.mean = 0, base.hazard.sd = 10,
          ...
        )
      }
      a0.lognc <- lapply(a0, function(a){
        logncfun(a0 = a,
                 chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                 refresh = 0)
      })
      a0.lognc <- data.frame( do.call(rbind, a0.lognc) )

      d <- pwe.npp(
        formula = fmla,
        data.list = data.list,
        a0.lognc = a0.lognc$a0,
        lognc = a0.lognc$lognc,
        breaks = breaks,
        beta.mean = 0, beta.sd = 10,
        base.hazard.mean = 0, base.hazard.sd = 10,
        a0.shape1 = 1, a0.shape2 = 1,
        a0.lower = 0, a0.upper = 1,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )

    }else if( prior.id == "leap" ){
      # Latent exchangeability prior
      d <- pwe.leap(
        formula = fmla,
        data.list = data.list,
        breaks = breaks,
        K = 2, prob.conc = 1,
        gamma.lower = 0, gamma.upper = 1,
        beta.mean = 0, beta.sd = 10,
        base.hazard.mean = 0, base.hazard.sd = 10,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )

    }else{
      stop("prior.id incorrect")
    }
  }

}else{
  
  ###################################### CurePWE model ######################################
  
  if ( prior.id == "psipp" ){
    # PSIPP
    d <- curepwe.stratified.pp(
      formula = fmla.psipp,
      data.list = data.list.PSIPP,
      strata.list = strata.list,
      breaks = breaks,
      a0.strata = res.strata$a0.strata,
      beta.mean = 0, beta.sd = 10,
      base.hazard.mean = 0, base.hazard.sd = 10,
      logit.pcured.mean = 0, logit.pcured.sd = 3,
      get.loglik = T,
      chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
    )
    
  }else {
    if ( prior.id == "ref" ){
      # Vague prior (using current data only)
      d <- curepwe.post(
        formula = fmla,
        data.list = list(curr),
        breaks = breaks,
        beta.mean = 0, beta.sd = 10,
        base.hazard.mean = 0, base.hazard.sd = 10,
        logit.pcured.mean = 0, logit.pcured.sd = 3,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )
      
    }else if( prior.id == "pp" ){
      # Power prior with a0 = 0.5
      d <- curepwe.pp(
        formula = fmla,
        data.list = data.list,
        breaks = breaks,
        a0 = 0.5,
        beta.mean = 0, beta.sd = 10,
        base.hazard.mean = 0, base.hazard.sd = 10,
        logit.pcured.mean = 0, logit.pcured.sd = 3,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )

    }else if( prior.id == "bhm" ){
      # Bayesian hierarchical model
      d <- curepwe.bhm(
        formula = fmla,
        data.list = data.list,
        breaks = breaks,
        meta.mean.mean = 0, meta.mean.sd = 10,
        meta.sd.mean = 0, meta.sd.sd = 0.5,
        base.hazard.mean = 0, base.hazard.sd = 10,
        logit.pcured.mean = 0, logit.pcured.sd = 3,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )

    }else if( prior.id == "cp" ){
      # Commensurate prior 
      d <- curepwe.commensurate(
        formula = fmla,
        data.list = data.list,
        breaks = breaks,
        beta0.mean = 0, beta0.sd = 10,
        p.spike = 0.1,
        spike.mean = 200, spike.sd = 0.1,
        slab.mean = 0, slab.sd = 5,
        base.hazard.mean = 0, base.hazard.sd = 10,
        logit.pcured.mean = 0, logit.pcured.sd = 3,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )

    }else if( prior.id == "npp" ){
      # Normalized power prior 
      # Pre-compute log normalizing constants Z(a0)
      a0 <- seq(0, 1, length.out = 101)

      # Wrapper function to compute log Z(a0) using the external data
      logncfun <- function(a0, ...){
        hdbayes::curepwe.npp.lognc(
          formula = fmla, histdata = data.list[[2]], breaks = breaks, a0 = a0,
          beta.mean = 0, beta.sd = 10,
          base.hazard.mean = 0, base.hazard.sd = 10,
          logit.pcured.mean = 0, logit.pcured.sd = 3,
          ...
        )
      }
      a0.lognc <- lapply(a0, function(a){
        logncfun(a0 = a,
                 chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                 refresh = 0)
      })
      a0.lognc <- data.frame( do.call(rbind, a0.lognc) )

      d <- curepwe.npp(
        formula = fmla,
        data.list = data.list,
        a0.lognc = a0.lognc$a0,
        lognc = a0.lognc$lognc,
        breaks = breaks,
        beta.mean = 0, beta.sd = 10,
        base.hazard.mean = 0, base.hazard.sd = 10,
        logit.pcured.mean = 0, logit.pcured.sd = 3,
        a0.shape1 = 1, a0.shape2 = 1,
        a0.lower = 0, a0.upper = 1,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )

    }else if( prior.id == "leap" ){
      # Latent exchangeability prior
      d <- curepwe.leap(
        formula = fmla,
        data.list = data.list,
        breaks = breaks,
        K = 2, prob.conc = 1,
        gamma.lower = 0, gamma.upper = 1,
        beta.mean = 0, beta.sd = 10,
        base.hazard.mean = 0, base.hazard.sd = 10,
        logit.pcured.mean = 0, logit.pcured.sd = 3,
        get.loglik = T,
        chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling
      )

    }else{
      stop("prior.id incorrect")
    }
  }

}

################################# Model-comparison metric (elpd) #################################
# Extract pointwise log-likelihood matrix
loglik <- suppressWarnings(
  d %>%
    dplyr::select(starts_with("log_lik")) %>%
    as.matrix
)
# Compute expected log predictive density (elpd) via leave-one-out cross-validation (LOO-CV) method
res.loo <- loo::loo(loglik)

# Save outputs
res <- list(
  'scen'      = scen
  , 'draws'     = d
  , 'res.loo'   = res.loo
)

if ( prior.id == "psipp" ){
  res <- append(
    res,
    list("res.strata" = res.strata)
  )
}
saveRDS(res, filename)
