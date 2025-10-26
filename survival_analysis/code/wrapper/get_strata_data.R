# install psrwe package from GitHub
# devtools::install_github("olssol/psrwe")
library(psrwe)

#' Estimate propensity scores (PS), assign subjects into different strata based on PS, and obtain stratum-specific power
#' prior parameter values (a0s) using psrwe package. This wraps `psrwe::psrwe_est()` and `psrwe::psrwe_borrow()`.
#'
#' @param ps_fml_covs     right side of the propensity score (PS) formula involving the covariates, e.g., ~ sex + age.
#' @param v_arm           column name corresponding to arm (treatment v.s. control) assignment.
#' @param borrow_ctl_only whether to borrow historical information on control arm only. If FALSE, then historical information on
#'                        both the treatment and control arm are borrowed, and the strata are formed based on the PS of all
#'                        subjects in current study. Otherwise, only information on control arm will be borrowed, and the strata
#'                        will be formed based on the PS of subjects in the current study control arm. Defaults to TRUE.
#'
get.strata.data = function(
    data, histdata,
    ps_fml_covs,
    ps_method       = c("logistic", "randomforest"),
    v_arm           = "treatment",
    ctl_arm_level   = 0,
    borrow_ctl_only = TRUE,
    nstrata         = 5,
    trim_ab         = c("both", "above", "below", "none"),
    total_borrow    = NULL,
    method          = c("distance", "inverse_distance", "n_current", "n_external")
) {
  # If borrowing only from historical control, restrict historical data to that arm
  if ( borrow_ctl_only ){
    histdata = histdata[ histdata[, v_arm] == ctl_arm_level, ]
  }
  
  n              = nrow(data)     # current data sample size
  n0             = nrow(histdata) # external data sample size
  data.all       = as.data.frame( rbind(data, histdata) )
  data.all$study = rep(c(1, 0), times = c(n, n0))

  # Estimate propensity scores (PS) via psrwe::psrwe_est() function
  data_ps        = psrwe::psrwe_est(
    data = data.all,
    ps_fml = as.formula( paste0("study", ps_fml_covs) ),
    ps_method = ps_method,
    v_grp = "study",
    cur_grp_level = 1,
    v_arm = v_arm,
    ctl_arm_level = ctl_arm_level,
    stra_ctl_only = borrow_ctl_only,
    nstrata = nstrata,
    trim_ab = trim_ab
  )

  # Distribute total number to borrow across strata based on PS via psrwe::psrwe_borrow() function
  if ( is.null(total_borrow) )
    total_borrow = n0
  ps_bor     = psrwe::psrwe_borrow(
    dtaps = data_ps,
    total_borrow = total_borrow,
    method = method
  )
  
  # a0 per stratum (set NA to 0 for strata with no eligible external data)
  a0.strata = ps_bor$Borrow$Alpha
  a0.strata[is.na(a0.strata)] = 0

  # Keep only non-trimmed rows (strata assigned) and carry strata labels
  data.all.strata = ps_bor$data
  data.all.strata = data.all.strata[!is.na(data.all.strata$`_strata_`), ]
  strata          = as.integer( data.all.strata$`_strata_` )

  # Return data frames with the original columns for current/external studies,
  # aligned with their strata assignments
  data.list.new   = list(
    curr = data.all.strata[data.all.strata$study == 1, colnames(data)],
    hist = data.all.strata[data.all.strata$study == 0, colnames(histdata)]
  )
  strata.list     = list(
    curr = strata[ data.all.strata$study == 1 ],
    hist = strata[ data.all.strata$study == 0 ]
  )
  
  res             = list(
    "data.list"        = data.list.new,
    "strata.list"      = strata.list,
    "a0.strata"        = a0.strata,
    "res.psrwe.est"    = data_ps,
    "res.psrwe.borrow" = ps_bor
  )
  return(res)
}
