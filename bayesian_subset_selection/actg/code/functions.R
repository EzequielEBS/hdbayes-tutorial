# function to compute power set
powerset <- function(x) {
  sets <- lapply(1:(length(x)), function(i) combn(x, i, simplify = F))
  unlist(sets, recursive = F)
}

# function to create formula
create_formula <- function(x, outcome) {
  as.formula(paste(paste(outcome, "~", sep = " "), paste(x, collapse = "+")))
}

# function to get covariates from models
get_covariates <- function(x) {
  paste(x, collapse = ", ")
}

elicit_beta_mean_cv <- function(m0, v0 = NULL, cv = 1) {
  if (!is.null(v0)) {
    a <- -(m0 * v0 + m0 ^ 3 - m0 ^ 2) / v0
    b <- ((m0 - 1) * v0 + m0 ^ 3 - 2 * m0 ^ 2 + m0) / v0
  } else{
    a <- -(m0 * (cv * m0) ^ 2 + m0 ^ 3 - m0 ^ 2) / (cv * m0) ^ 2
    b <- ((m0 - 1) * (cv * m0) ^ 2 + m0 ^ 3 - 2 * m0 ^ 2 + m0) / (cv * m0) ^
      2
  }
  if (a < 0 || b < 0) {
    warning("Warning: at least one of the obtained parameters is not valid")
  }
  return(list(a = a, b = b))
}

logp0 <- function(formula,
                  ...) {
  hdbayes::glm.npp.lognc(
    formula = formula,
    ...
  )
}

logncfun <- function(a0,
                     ...){
  hdbayes::glm.npp.lognc(
    a0 = a0,
    ...
  )
}

post_beta <- function(formula, 
                      c0, 
                      data,
                      family,
                      delta0,
                      lambda0,
                      iter_warmup,
                      iter_sampling) {
  a0.lognc <- lapply(a0_seq, 
                     logncfun, 
                     formula = formula, 
                     family = family, 
                     histdata = data[[2]], 
                     iter_warmup = iter_warmup,
                     iter_sampling = iter_sampling, 
                     chains = 4, 
                     refresh = 0,
                     beta.sd = c0
  )
  a0.lognc <- data.frame(do.call(rbind, a0.lognc))
  fit <- hdbayes::glm.npp(formula, 
                          data = data,
                          family = family,
                          a0.lognc = a0.lognc$a0,
                          lognc = matrix(a0.lognc$lognc, ncol = 1),
                          a0.shape1 = delta0,
                          a0.shape2 = lambda0,
                          iter_warmup = iter_warmup,
                          iter_sampling = iter_sampling, 
                          chains = 4, 
                          refresh = 0,
                          beta.sd = c0
  )
}

samples_models <- function(data,
                        outcome,
                        family,
                        a0_seq,
                        c0,
                        d0,
                        delta0,
                        lambda0,
                        iter_warmup,
                        iter_sampling,
                        num_cores = 10) {
  current_data <- data[[1]]
  hist_data <- data[[2]]
  
  # create list of formulas
  covariates <- colnames(current_data)[colnames(current_data) != outcome]
  pset <- powerset(covariates)
  formulas <- lapply(pset, create_formula, outcome = outcome)
  
  # get covariates from models
  covariates_models <- lapply(pset, get_covariates)
  
  cl <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(cl))
  parallel::clusterExport(cl, varlist = c('a0_seq',
                                          'family',
                                          'data',
                                          'c0',
                                          'd0',
                                          'delta0',
                                          'lambda0',
                                          'iter_warmup',
                                          'iter_sampling',
                                          'logncfun'),
                          envir = environment()
  )
  
  logp0_models <- parLapply(
    cl = cl,
    X = formulas,
    fun = logp0,
    family = family,
    histdata = hist_data,
    a0 = 1,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    chains = 4,
    refresh = 0,
    beta.sd = d0
  )

  post_betam  <- parLapply(
    cl = cl,
    X = formulas,
    fun = post_beta,
    c0 = c0, 
    data = data,
    family = family,
    delta0 = delta0,
    lambda0 = lambda0,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling
  )

  logml_models <- parLapply(
    cl = cl,
    X = post_betam,
    fun = hdbayes::glm.logml.npp
  )

  prior_models <- exp(sapply(logp0_models, function(x) x[2]) -
                        logSumExp(sapply(logp0_models, function(x) x[2]))
  )

  post_models <- post_models(logp0_models, logml_models)

  df_post <- data.frame(model = unlist(covariates_models),
                        prior_model = prior_models,
                        ml = exp(unlist(data.frame(do.call(rbind, logml_models))$logml)),
                        post_model = post_models
  )

  df_post_ord <- df_post[order(df_post$post_model, decreasing = TRUE),]
  rownames(df_post_ord) <- 1:nrow(df_post_ord)
  return(list(logp0_models = logp0_models,
              post_betam = post_betam,
              logml_models = logml_models,
              df_post = df_post,
              df_post_ord = df_post_ord))
}

post_models <- function(logp0_models, logml_models) {
  logp0_models <- sapply(logp0_models, function(x) x[2])
  logml_models <- unlist(data.frame(do.call(rbind, logml_models))$logml)
  post_model <- exp(logml_models + logp0_models -
                      logSumExp(logml_models + logp0_models))
}

predict_ <- function(data, outcome, trt, family, beta_post){
  cov_data <- colnames(data)[!(colnames(data) %in% c(outcome, trt))]
  prod <- lapply(beta_post,
                 function(post_samp) {
                   cov <- names(post_samp)[names(post_samp) %in% cov_data]
                   post_samp <- as.data.frame(post_samp[c("(Intercept)", cov)])
                   cov_mat <- cbind(1, data.matrix(data[cov]))
                   cov_mat <- t(cov_mat)
                   beta <- data.matrix(post_samp)
                   beta %*% cov_mat
                 }
  )
  if (family$link == "logit") {
    pred <- lapply(prod, function(df) {
      apply(df, 2, plogis)
    })
  } else if (family$link == "probit") {
    pred <- lapply(prod, function(df) {
      apply(df, 2, pnorm)
    })
  }
  return(pred)
}

mean_models_arm <- function(data, outcome, trt, family, beta_post){
  pred <- predict_(data, outcome, trt, family, beta_post)
  mean_models <- lapply(pred,
                        function(x) {
                          apply(x, 1, function(row){
                            sample_weights <- t(rdirichlet(1, rep(1, nrow(data))))
                            row %*% sample_weights
                          })
                        }
  )
  # mean_models <- lapply(pred, rowMeans)
  return(do.call(cbind, mean_models))
}

mean_arm <- function(post_models, mean_arm_models){
  mean_arm_models %*% post_models
}

bma <- function(x, df_post, n_samples = 1000){
  ks <- sample(1:nrow(df_post), n_samples, replace = TRUE, prob = df_post$post_model)
  samples <- lapply(ks, function(k) {
    id_draw <- sample(1:nrow(x), 1) 
    draw_x <- x[id_draw,k]
  })
  return(do.call(rbind, samples))
}
    
