#-------------------------------------------------------------------------------
# Prior predictive distributions 
#-------------------------------------------------------------------------------

load("logistic_regression/samples/draws_npp_c0.RData")
load("logistic_regression/samples/draws_npp_c0_after_PSM.RData")
load("logistic_regression/samples/draws_npp_wip.RData")
load("logistic_regression/samples/draws_npp_wip_after_PSM.RData")
load("logistic_regression/samples/draws_npp_a0.RData")
load("logistic_regression/samples/draws_npp_a0_after_PSM.RData")

library(tidyverse)
library(ggridges)
library(hdbayes)
library(parallel)
library(ggplot2)
library(pROC)
library(ggplotify)
library(patchwork)

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

data <- list(current_data, hist_data)
family <- binomial(link = "logit")
formula <- outcome ~ age + treatment + cd4

Xnew        <- stats::model.matrix(formula, current_data)
nsim <- 1000
c0 <- c(0.25, 0.5, 1, 2)
a0_hyper <- list(c(1,1), c(2,2), c(1,10), c(10,1))
d.sub.norm <- lapply(seq_along(c0), function(j){
  draws.c0[[j]][sample(x = seq_len(nrow(draws.c0[[j]])), size = nsim, replace = F), ]
})
d.sub.wip  <- draws.wip[sample(x = seq_len(nrow(draws.wip)), size = nsim, replace = F), ]
d.sub.a0 <- lapply(seq_along(c0), function(j){
  draws.a0[[j]][sample(x = seq_len(nrow(draws.a0[[j]])), size = nsim, replace = F), ]
})

pnew.norm <- lapply(seq_along(c0), function(j) {
  sapply(seq_len(nsim), function(i){
    beta.sim <- as.numeric(d.sub.norm[[j]][i, ])
    p <- binomial('logit')$linkinv(Xnew %*% beta.sim)
    return(p)
  })
})

save(pnew.norm, 
     file = "logistic_regression/samples/ppc_pnew_norm.RData")

pnew.wip <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d.sub.wip[i, ])
  p <- binomial('logit')$linkinv(Xnew %*% beta.sim)
  return(p)
})

save(pnew.wip, 
     file = "logistic_regression/samples/ppc_pnew_wip.RData")

pnew.a0 <- lapply(seq_along(a0_hyper), function(j) {
  sapply(seq_len(nsim), function(i){
    beta.sim <- as.numeric(d.sub.a0[[j]][i, ])
    p <- binomial('logit')$linkinv(Xnew %*% beta.sim)
    return(p)
  })
})

save(pnew.a0, 
     file = "logistic_regression/samples/ppc_pnew_a0.RData")

roc.curve.wip <- roc(current_data$outcome, rowMeans(pnew.wip))
save(roc.curve.wip, 
     file = "logistic_regression/samples/ppc_roc_curve_wip.RData")
# Compute the confidence interval for the AUC
ci.auc.wip <- ci.auc(roc.curve.wip)
save(ci.auc.wip, 
     file = "logistic_regression/samples/ppc_ci_auc_wip.RData")
ci.roc.wip <- ci.se(roc.curve.wip, specificities = seq(0, 1, l = 25))
save(ci.roc.wip, 
     file = "logistic_regression/samples/ppc_ci_roc_wip.RData")

roc.curves.wip <- lapply(1:ncol(pnew.wip), function(i) {
  roc(current_data$outcome, pnew.wip[, i])
})

save(roc.curves.wip, 
     file = "logistic_regression/samples/ppc_roc_curves_wip.RData")

roc.curve.norm <- lapply(seq_along(pnew.norm), function(j) {
  roc(current_data$outcome, rowMeans(pnew.norm[[j]]),
                      ci = TRUE)
  })
save(roc.curve.norm, 
     file = "logistic_regression/samples/ppc_roc_curve_norm.RData")

ci.auc.norm <- lapply(seq_along(roc.curve.norm), function(j) {
  ci.auc(roc.curve.norm[[j]])
})
save(ci.auc.norm, 
     file = "logistic_regression/samples/ppc_ci_auc_norm.RData")

ci.roc.norm <- lapply(seq_along(roc.curve.norm), function(j) {
  ci.se(roc.curve.norm[[j]], specificities = seq(0, 1, l = 25))
})
save(ci.roc.norm, 
     file = "logistic_regression/samples/ppc_ci_roc_norm.RData")

roc.curves.norm <- lapply(seq_along(pnew.norm), function(j) {
  lapply(1:ncol(pnew.norm[[j]]), function(i) {
  roc(current_data$outcome, pnew.norm[[j]][, i])
})
})
save(roc.curves.norm, 
     file = "logistic_regression/samples/ppc_roc_curves_norm.RData")

roc.curve.a0 <- lapply(seq_along(pnew.a0), function(j) {
  roc(current_data$outcome, rowMeans(pnew.a0[[j]]),
      ci = TRUE)
})
save(roc.curve.a0, 
     file = "logistic_regression/samples/ppc_roc_curve_a0.RData")

ci.auc.a0 <- lapply(seq_along(roc.curve.a0), function(j) {
  ci.auc(roc.curve.a0[[j]])
})
save(ci.auc.a0, 
     file = "logistic_regression/samples/ppc_ci_auc_a0.RData")

ci.roc.a0 <- lapply(seq_along(roc.curve.a0), function(j) {
  ci.se(roc.curve.a0[[j]], specificities = seq(0, 1, l = 25))
})
save(ci.roc.a0, 
     file = "logistic_regression/samples/ppc_ci_roc_a0.RData")

roc.curves.a0 <- lapply(seq_along(pnew.a0), function(j) {
  lapply(1:ncol(pnew.a0[[j]]), function(i) {
    roc(current_data$outcome, pnew.a0[[j]][, i])
  })
})
save(roc.curves.a0, 
     file = "logistic_regression/samples/ppc_roc_curves_a0.RData")

#-------------------------------------------------------------------------------
# After PSM
#-------------------------------------------------------------------------------

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

data <- list(current_data, hist_data)
family <- binomial(link = "logit")
formula <- outcome ~ age + treatment + cd4

Xnew        <- stats::model.matrix(formula, current_data)
nsim <- 1000
c0 <- c(0.25, 0.5, 1, 2)
a0_hyper <- list(c(1,1), c(2,2), c(1,10), c(10,1))
d.sub.norm_after_PSM <- lapply(seq_along(c0), function(j){
  draws.c0_after_PSM[[j]][sample(x = seq_len(nrow(draws.c0_after_PSM[[j]])), 
                                 size = nsim, replace = F), ]
})
d.sub.wip_after_PSM  <- draws.wip_after_PSM[sample(x = seq_len(nrow(draws.wip_after_PSM)), 
                                                   size = nsim, replace = F), ]
d.sub.a0_after_PSM <- lapply(seq_along(c0), function(j){
  draws.a0_after_PSM[[j]][sample(x = seq_len(nrow(draws.a0_after_PSM[[j]])), 
                                 size = nsim, replace = F), ]
})

pnew.norm_after_PSM <- lapply(seq_along(c0), function(j) {
  sapply(seq_len(nsim), function(i){
    beta.sim <- as.numeric(d.sub.norm_after_PSM[[j]][i, ])
    p <- binomial('logit')$linkinv(Xnew %*% beta.sim)
    return(p)
  })
})
save(pnew.norm_after_PSM, 
     file = "logistic_regression/samples/ppc_pnew_norm_after_PSM.RData")

pnew.wip_after_PSM <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d.sub.wip_after_PSM[i, ])
  p <- binomial('logit')$linkinv(Xnew %*% beta.sim)
  return(p)
})
save(pnew.wip_after_PSM, 
     file = "logistic_regression/samples/ppc_pnew_wip_after_PSM.RData")

pnew.a0_after_PSM <- lapply(seq_along(a0_hyper), function(j) {
  sapply(seq_len(nsim), function(i){
    beta.sim <- as.numeric(d.sub.a0_after_PSM[[j]][i, ])
    p <- binomial('logit')$linkinv(Xnew %*% beta.sim)
    return(p)
  })
})
save(pnew.a0_after_PSM, 
     file = "logistic_regression/samples/ppc_pnew_a0_after_PSM.RData")


roc.curve.wip_after_PSM <- roc(current_data$outcome, rowMeans(pnew.wip_after_PSM))
save(roc.curve.wip_after_PSM, 
     file = "logistic_regression/samples/ppc_roc_curve_wip_after_PSM.RData")

ci.auc.wip_after_PSM <- ci.auc(roc.curve.wip_after_PSM)
save(ci.auc.wip_after_PSM, 
     file = "logistic_regression/samples/ppc_ci_auc_wip_after_PSM.RData")

ci.roc.wip_after_PSM <- ci.se(roc.curve.wip_after_PSM, specificities = seq(0, 1, l = 25))
save(ci.roc.wip_after_PSM, 
     file = "logistic_regression/samples/ppc_ci_roc_wip_after_PSM.RData")

roc.curves.wip_after_PSM <- lapply(1:ncol(pnew.wip_after_PSM), function(i) {
  roc(current_data$outcome, pnew.wip_after_PSM[, i])
})
save(roc.curves.wip_after_PSM, 
     file = "logistic_regression/samples/ppc_roc_curves_wip_after_PSM.RData")

roc.curve.norm_after_PSM <- lapply(seq_along(pnew.norm_after_PSM), function(j) {
  roc(current_data$outcome, rowMeans(pnew.norm_after_PSM[[j]]),
      ci = TRUE)
})
save(roc.curve.norm_after_PSM, 
     file = "logistic_regression/samples/ppc_roc_curve_norm_after_PSM.RData")

ci.auc.norm_after_PSM <- lapply(seq_along(roc.curve.norm_after_PSM), function(j) {
  ci.auc(roc.curve.norm_after_PSM[[j]])
})
save(ci.auc.norm_after_PSM, 
     file = "logistic_regression/samples/ppc_ci_auc_norm_after_PSM.RData")

ci.roc.norm_after_PSM <- lapply(seq_along(roc.curve.norm_after_PSM), function(j) {
  ci.se(roc.curve.norm_after_PSM[[j]], specificities = seq(0, 1, l = 25))
})
save(ci.roc.norm_after_PSM, 
     file = "logistic_regression/samples/ppc_ci_roc_norm_after_PSM.RData")

roc.curves.norm_after_PSM <- lapply(seq_along(pnew.norm_after_PSM), function(j) {
  lapply(1:ncol(pnew.norm_after_PSM[[j]]), function(i) {
    roc(current_data$outcome, pnew.norm_after_PSM[[j]][, i])
  })
})
save(roc.curves.norm_after_PSM, 
     file = "logistic_regression/samples/ppc_roc_curves_norm_after_PSM.RData")

roc.curve.a0_after_PSM <- lapply(seq_along(pnew.a0_after_PSM), function(j) {
  roc(current_data$outcome, rowMeans(pnew.a0_after_PSM[[j]]),
      ci = TRUE)
})
save(roc.curve.a0_after_PSM, 
     file = "logistic_regression/samples/ppc_roc_curve_a0_after_PSM.RData")

ci.auc.a0_after_PSM <- lapply(seq_along(roc.curve.a0_after_PSM), function(j) {
  ci.auc(roc.curve.a0_after_PSM[[j]])
})
save(ci.auc.a0_after_PSM, 
     file = "logistic_regression/samples/ppc_ci_auc_a0_after_PSM.RData")

ci.roc.a0_after_PSM <- lapply(seq_along(roc.curve.a0_after_PSM), function(j) {
  ci.se(roc.curve.a0_after_PSM[[j]], specificities = seq(0, 1, l = 25))
})
save(ci.roc.a0_after_PSM, 
     file = "logistic_regression/samples/ppc_ci_roc_a0_after_PSM.RData")

roc.curves.a0_after_PSM <- lapply(seq_along(pnew.a0_after_PSM), function(j) {
  lapply(1:ncol(pnew.a0_after_PSM[[j]]), function(i) {
    roc(current_data$outcome, pnew.a0_after_PSM[[j]][, i])
  })
})
save(roc.curves.a0_after_PSM, 
     file = "logistic_regression/samples/ppc_roc_curves_a0_after_PSM.RData")
