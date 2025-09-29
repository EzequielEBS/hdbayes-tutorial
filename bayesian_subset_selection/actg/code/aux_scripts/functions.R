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

predict_ <- function(data, outcome, trt, family, beta_post, arm){
  # cov_data <- colnames(data)[!(colnames(data) %in% c(outcome, trt))]
  if (arm == 1) {
    data[[trt]] <- 1
  } else if (arm == 0) {
    data[[trt]] <- 0
  }
  cov_data <- colnames(data)[!(colnames(data) %in% c(outcome))]
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

mean_models_arm <- function(data, outcome, trt, family, beta_post, arm){
  pred <- predict_(data, outcome, trt, family, beta_post, arm)
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

# mean_arm <- function(post_models, mean_arm_models){
#   mean_arm_models %*% post_models
# }

bma <- function(x, df_post, n_samples = 1000){
  ks <- sample(1:nrow(df_post), n_samples, replace = TRUE, prob = df_post$post_model)
  samples <- lapply(ks, function(k) {
    id_draw <- sample(1:nrow(x), 1) 
    draw_x <- x[id_draw,k]
  })
  return(do.call(rbind, samples))
}

plot_ate <- function(samples){
  mean_models_ctrl <- mean_models_arm(current_data, 
                                      "outcome", 
                                      "treatment", 
                                      family, 
                                      samples$post_betam,
                                      0)
  mean_models_trt <- mean_models_arm(current_data, 
                                     "outcome", 
                                     "treatment", 
                                     family, 
                                     samples$post_betam,
                                     1)
  bma_ctrl <- bma(mean_models_ctrl, samples$df_post, 10000)
  bma_trt <- bma(mean_models_trt, samples$df_post, 10000)
  
  df_bma <- data.frame(
    value = bma_trt - bma_ctrl
  )
  
  ci_bma_95 <- bayestestR::ci(df_bma$value, ci = 0.95)
  ci_bma_90 <- bayestestR::ci(df_bma$value, ci = 0.90)
  
  ci_ate_95_lower <- ci_bma_95$CI_low
  ci_ate_95_upper <- ci_bma_95$CI_high
  ci_ate_90_lower <- ci_bma_90$CI_low
  ci_ate_90_upper <- ci_bma_90$CI_high
  
  # Compute density for the entire dataset
  ate_density_data <- density(df_bma$value) # Compute density
  ate_density_df <- data.frame(x = ate_density_data$x, y = ate_density_data$y) # Convert to data frame
  
  # Filter density data for the HDI regions
  ate_density_95 <- ate_density_df %>% filter(x >= ci_ate_95_lower & x <= ci_ate_95_upper)
  ate_density_90 <- ate_density_df %>% filter(x >= ci_ate_90_lower & x <= ci_ate_90_upper)
  
  # Compute the mean, median, and quartiles
  mean_value <- mean(df_bma$value)
  median_value <- median(df_bma$value)
  q1_value <- quantile(df_bma$value, 0.25)
  q3_value <- quantile(df_bma$value, 0.75)
  
  blended_rgb <- round(colMeans(rbind(
    c(135, 206, 235),
    c(70, 130, 180)
  )))
  
  blended_color <- rgb(blended_rgb[1], blended_rgb[2], blended_rgb[3], maxColorValue = 255)
  
  # Define colors for mean, median, and quartiles
  stats_colors <- c("90% \nBCI" = blended_color,
                    
                    "Density" = "skyblue",
                    "ATE = 0" = "black",
                    "OR = 1" = "black"
  )
  
  # Plot using ggplot2
  bma_ate <- ggplot() +
    # Full density curve
    # geom_density(data = df_bma, aes(x = value), color = "skyblue", fill = "skyblue", alpha = 0.5) +
    geom_area(data = ate_density_df %>% filter(x <= ci_ate_90_lower + 
                                                 (ci_ate_90_lower - min(ate_density_df$x))*0.007), 
              aes(x = x, y = y, fill = "Density"), color = "black",
              alpha = 0.7) +
    geom_area(data = ate_density_df %>% filter(x >= ci_ate_90_upper - 
                                                 (max(ate_density_df$x) - ci_ate_90_upper)*0.007), 
              aes(x = x, y = y, fill = "Density"), color = "black",
              alpha = 0.7) +
    
    # # Highlight 95% HDI region
    geom_area(data = ate_density_90, 
              aes(x = x, y = y, fill = "90% \nBCI"), color = "black",
              alpha = 0.7) +
    # 
    # # Highlight 90% HDI region
    # geom_area(data = ate_density_90, aes(x = x, y = y, fill = "90% \nBCI", colour = "Density"), alpha = 0.7) +
    
    # # Add vertical dashed lines for 95% HDI
    # geom_vline(aes(xintercept = ci_95_lower, color = "95% BCI"), linetype = "dashed", size = 1) +
    # geom_vline(aes(xintercept = ci_95_upper, color = "95% BCI"), linetype = "dashed", size = 1) +
    
    # Add vertical dashed lines for 90% BCI
    # geom_vline(aes(xintercept = ci_ate_90_lower, color = "90% BCI"), linetype = "dashed", size = 1) +
    # geom_vline(aes(xintercept = ci_ate_90_upper, color = "90% BCI"), linetype = "dashed", size = 1) +
    
    # Add lines for mean, median, and quartiles
    geom_vline(aes(xintercept = 0, color = "ATE = 0"), linetype = "solid", size = 1) + 
    # geom_vline(aes(xintercept = median_value, color = "Median"), linetype = "solid", size = 1) + 
    # geom_vline(aes(xintercept = q1_value, color = "Q1"), linetype = "dotted", size = 1) + 
    # geom_vline(aes(xintercept = q3_value, color = "Q3"), linetype = "dotted", size = 1) +
    
    # Add a legend for both color and fill
    scale_fill_manual(name = NULL, values = stats_colors, breaks = c("90% \nBCI")) +
    scale_color_manual(name = NULL, values = stats_colors, guide = NULL) +
    
    # Add labels
    labs(title = "",
         x = "Average treatment effect (ATE)", y = "", color = "") +
    
    theme_bw() +
    theme(legend.position = c(.95, .95),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6),
          legend.background = element_rect(fill = "white", color = "black"),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
    ) +
    theme(text = element_text(size = 16),        # Base text size
          axis.title = element_text(size = 18),  # Axis titles
          axis.text = element_text(size = 16),   # Axis tick labels
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16))
  bma_ate
}
