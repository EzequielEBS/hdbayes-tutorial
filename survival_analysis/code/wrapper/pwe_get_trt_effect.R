library(MCMCpack)

# Estimand: Treatment effect = difference in two-year relapse-free survival probability:
# P(T > 2 | A = 1) - P(T > 2 | A = 0), where A is the treatment indicator (0 = control, 1 = treatment).

#' Predict survival probability at t years for a single arm A, S(t | A), under a PWEPH model (with non-PSIPP priors).
#' 
#' Uses posterior draws from a PWEPH model fit to a single arm (treated or untreated). The treatment indicator
#' is not a covariate here since each arm is fit separately. The arm-specific posterior parameter estimates are
#' applied to the covariates from all subjects in \code{data} to predict S(t | A).
#'
#' @param t                 time at which to evaluate survival probability.
#' @param post.samples      posterior draws from a PWEPH model fit under various priors (other than PSIPP), with 
#'                          an attribute called 'data' which includes the list of variables specified in the data 
#'                          block of the Stan program.
#' @param data              data set with the baseline covariates for all subjects across treatment arms. 
#'
get.surv.prob.pwe <- function(
    t,
    post.samples,
    data
){
  # Extract Stan data
  stan.data = attr(post.samples, 'data')
  J         = stan.data$J
  p         = stan.data$p
  breaks    = stan.data$breaks
  
  # Build covariate matrix X for all subjects
  X         = as.matrix( data[, colnames(stan.data$X1)] )
  # Obtain posterior draws of parameters
  beta      = suppressWarnings(
    as.matrix( post.samples[, colnames(X), drop=F] )
  )
  lambda    = suppressWarnings(
    as.matrix( post.samples[, paste0("basehaz[", 1:J, "]"), drop=F] )
  )
  # Compute linear predictor
  eta       = tcrossprod(beta, X)

  # Find interval index k such that breaks[k] < t <= breaks[k+1]
  interval.id = findInterval(t, breaks, left.open = TRUE)
  if( t == 0 ){
    interval.id = 1
  }

  # Compute cumulative hazard 
  if( J > 1 ){
    # Compute cumulative baseline hazard at each interval
    cumblhaz = apply(lambda, 1, function(x){
      as.numeric( cumsum( x[1:(J-1)] * ( breaks[2:J] - breaks[1:(J-1)] ) ) )
    })
    cumblhaz = matrix(cumblhaz, nrow = J-1)
    cumblhaz = cbind(0, t(cumblhaz))
    
    cumhaz = lambda[, interval.id] * (t - breaks[interval.id]) + cumblhaz[, interval.id]
    
  }else{
    cumhaz = lambda[, interval.id] * (t - breaks[interval.id])
  }
  cumhaz = cumhaz * exp(eta)

  # Subject-specific survival probability at time t for each posterior draw of parameters
  S      = exp( -cumhaz ) # number of posterior draws x number of subjects in `data`
  
  # Use Bayesian bootstrap method 
  # For each posterior draw, average the predicted t-year survival probabilities over subjects
  # with Dirichlet(1, .., 1) weights.
  omega = MCMCpack::rdirichlet(n = nrow(beta), alpha = rep(1, nrow(X)))
  surv  = rowSums( omega * S )

  return(surv)
}

#' Predict S(t | A) for a PWEPH model under PSIPP priors (stratified hazards, no covariates in the outcome model).
#' 
#' Uses posterior draws from a PWEPH model fit under PSIPP, which yields stratum-specific baseline hazards. The
#' outcome model has no covariates.
#'
#' @param t                 time at which to evaluate survival probability.
#' @param post.samples      posterior draws from a PWEPH model fit under PSIPP, with an attribute called 'data' 
#'                          which includes the list of variables specified in the data block of the Stan program.
#'
get.surv.prob.pwe.psipp <- function(
    t,
    post.samples
){
  # Extract Stan data
  stan.data = attr(post.samples, 'data')
  J         = stan.data$J
  p         = stan.data$p
  K         = stan.data$K
  breaks    = stan.data$breaks

  # Posterior draws of stratum-specific baseline hazards
  lambdaMat = suppressWarnings(
    as.matrix( post.samples[, paste0("basehaz", "_stratum_", rep(1:K, each = J), "[", 1:J, "]"), drop=F] )
  )

  # Find interval index k such that breaks[k] < t <= breaks[k+1]
  interval.id = findInterval(t, breaks, left.open = TRUE)
  if( t == 0 ){
    interval.id = 1
  }

  # Compute cumulative hazard per stratum
  if( J > 1 ){
    cumhaz = sapply(1:K, function(k){
      lambda = lambdaMat[, paste0("basehaz", "_stratum_", k, "[", 1:J, "]"), drop = F]

      # Compute cumulative baseline hazard at breakpoints
      cumblhaz = apply(lambda, 1, function(x){
        as.numeric( cumsum( x[1:(J-1)] * ( breaks[2:J] - breaks[1:(J-1)] ) ) )
      })
      cumblhaz = matrix(cumblhaz, nrow = J-1)
      cumblhaz = cbind(0, t(cumblhaz))
      
      return(
        lambda[, interval.id] * (t - breaks[interval.id]) + cumblhaz[, interval.id]
      )
    })

  }else{
    cumhaz = lambdaMat * (t - breaks[interval.id])
  }

  # Stratum-specific survival probability at time t for each posterior draw of parameters
  S    = exp( -cumhaz ) # number of posterior draws x number of strata
  # Average across strata (equal weights)
  surv = rowMeans(S)

  return(surv)
}
