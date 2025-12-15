# ------------------------------------------------------------------------------
# Sample from the posterior distribution of a0 under different priors
# ------------------------------------------------------------------------------

library(hdbayes)
library(parallel)
library(matrixStats)
library(MCMCpack)
library(bayestestR)
library(dplyr)

current_data <- actg036
hist_data <- actg019

current_data$age <- (current_data$age - mean(current_data$age)) /
  (2*sd(current_data$age))
current_data$cd4 <- (current_data$cd4 - mean(current_data$cd4)) /
  (2*sd(current_data$cd4))
hist_data$age <- (hist_data$age - mean(hist_data$age)) /
  (2*sd(hist_data$age))
hist_data$cd4 <- (hist_data$cd4 - mean(hist_data$cd4)) /
  (2*sd(hist_data$cd4))

family <- binomial(link = "logit")
c0 <- c(0.25, 0.5, 1, 2)
a0_hyper <- list(c(1,1), c(2,2), c(1,10), c(10,1))
iter_warmup <- 1000
iter_sampling <- 2500
formula <- outcome ~ age + treatment + cd4
a0      <- seq(0, 1, length.out = 21)
logncfun.wip <- function(a0, ...){
  glm.npp.lognc.wip(
    formula = formula, family = family, a0 = a0, histdata = hist_data,
    ...
  ) 
}

ncores <- detectCores() - 1
a0.lognc.wip <- mclapply(
  X = a0, FUN = logncfun.wip, iter_warmup = 1000,
  iter_sampling = 2500, chains = 4,
  mc.cores = 1
)
a0.lognc.wip <- data.frame( do.call(rbind, a0.lognc.wip) )

draws.a0.post.wip <- sapply(seq_along(a0_hyper), function(i) {
  delta0_val <- a0_hyper[[i]][1]
  lambda0_val <- a0_hyper[[i]][2]
  
  fit <- glm.npp.wip(
    formula = formula,
    family = family,
    data.list = list(data = current_data, histdata = hist_data),
    a0.lognc = a0.lognc.wip$a0,
    lognc = matrix(a0.lognc.wip$lognc, ncol = 1),
    chains = 4, iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    refresh = 0,
    a0.shape1 = delta0_val,
    a0.shape2 = lambda0_val
  )
  
  fit$a0_hist_1
})
save(draws.a0.post.wip, 
     file = "logistic_regression/samples/post_samples_a0_wip.RData")

logncfun <- function(a0, ...){
  glm.npp.lognc(
    formula = formula, family = family, a0 = a0, histdata = hist_data,
    ...
  ) 
}


draws.a0.post.norm <- lapply(seq_along(c0), function(j) {
  a0.lognc <- lapply(
    X = a0, FUN = logncfun, iter_warmup = 1000,
    iter_sampling = 2500, chains = 4,
    beta.sd = c0[j]
  )
  a0.lognc <- data.frame( do.call(rbind, a0.lognc) )
  sapply(seq_along(a0_hyper), function(i) {
    delta0_val <- a0_hyper[[i]][1]
    lambda0_val <- a0_hyper[[i]][2]
    
    fit <- glm.npp(
      formula = formula,
      family = family,
      data.list = list(data = current_data, histdata = hist_data),
      a0.lognc = a0.lognc$a0,
      lognc = matrix(a0.lognc$lognc, ncol = 1),
      chains = 4, iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      refresh = 0,
      a0.shape1 = delta0_val,
      a0.shape2 = lambda0_val,
      beta.sd = c0[j]
    )
    
    fit$a0_hist_1
  })
})

save(draws.a0.post.norm, 
     file = "logistic_regression/samples/post_samples_a0_norm.RData")

# ------------------------------------------------------------------------------
# After PSM
# ------------------------------------------------------------------------------

current_data <- actg036
hist_data <- readRDS("logistic_regression/data/actg019_after_PSM.rds")

current_data$age <- (current_data$age - mean(current_data$age)) /
  (2*sd(current_data$age))
current_data$cd4 <- (current_data$cd4 - mean(current_data$cd4)) /
  (2*sd(current_data$cd4))
hist_data$age <- (hist_data$age - mean(hist_data$age)) /
  (2*sd(hist_data$age))
hist_data$cd4 <- (hist_data$cd4 - mean(hist_data$cd4)) /
  (2*sd(hist_data$cd4))

a0.lognc.wip_after_PSM <- mclapply(
  X = a0, FUN = logncfun.wip, iter_warmup = 1000,
  iter_sampling = 2500, chains = 4,
  mc.cores = 1
)
a0.lognc.wip_after_PSM <- data.frame( do.call(rbind, a0.lognc.wip_after_PSM) )

draws.a0.post.wip_after_PSM <- sapply(seq_along(a0_hyper), function(i) {
  delta0_val <- a0_hyper[[i]][1]
  lambda0_val <- a0_hyper[[i]][2]
  
  fit <- glm.npp.wip(
    formula = formula,
    family = family,
    data.list = list(data = current_data, histdata = hist_data),
    a0.lognc = a0.lognc.wip_after_PSM$a0,
    lognc = matrix(a0.lognc.wip_after_PSM$lognc, ncol = 1),
    chains = 4, iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    refresh = 0,
    a0.shape1 = delta0_val,
    a0.shape2 = lambda0_val
  )
  
  fit$a0_hist_1
})

save(draws.a0.post.wip_after_PSM, 
     file = "logistic_regression/samples/post_samples_a0_wip_after_PSM.RData")

draws.a0.post.norm_after_PSM <- lapply(seq_along(c0), function(j) {
  a0.lognc <- lapply(
    X = a0, FUN = logncfun, iter_warmup = 1000,
    iter_sampling = 2500, chains = 4,
    beta.sd = c0[j]
  )
  a0.lognc <- data.frame( do.call(rbind, a0.lognc) )
  sapply(seq_along(a0_hyper), function(i) {
    delta0_val <- a0_hyper[[i]][1]
    lambda0_val <- a0_hyper[[i]][2]
    
    fit <- glm.npp(
      formula = formula,
      family = family,
      data.list = list(data = current_data, histdata = hist_data),
      a0.lognc = a0.lognc$a0,
      lognc = matrix(a0.lognc$lognc, ncol = 1),
      chains = 4, iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      refresh = 0,
      a0.shape1 = delta0_val,
      a0.shape2 = lambda0_val,
      beta.sd = c0[j]
    )
    
    fit$a0_hist_1
  })
})

save(draws.a0.post.norm_after_PSM, 
     file = "logistic_regression/samples/post_samples_a0_norm_after_PSM.RData")
