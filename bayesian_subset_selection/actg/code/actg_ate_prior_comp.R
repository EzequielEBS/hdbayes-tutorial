load("bayesian_subset_selection/actg/samples/post_samples_c0d0.RData")
load("bayesian_subset_selection/actg/samples/post_samples_wip.RData")
load("bayesian_subset_selection/actg/samples/post_samples_c0d0_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/post_samples_wip_after_PSM.RData")

library(ggplot2)
library(MCMCpack)
library(bayestestR)
library(patchwork)
library(hdbayes)
library(dplyr)
library(tidyverse)

source("bayesian_subset_selection/actg/code/aux_scripts/functions.R")

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
ate_wip <- plot_ate(post_samples_wip)

c0 <- c(0.25, 0.5, 1, 2)
# d0 <- c(0.25, 0.5, 1, 2)
# m0_v0 <- list(c(1, 1), c(2, 2), c(1, 10), c(10, 1))

# ates_d0 <- lapply(post_samples_d0, plot_ate)
# ates_c0 <- lapply(post_samples_c0, plot_ate)
# ates_m0v0 <- lapply(post_samples_m0v0, plot_ate)
ates_c0d0 <- lapply(post_samples_c0d0, plot_ate)


ate_wip_norm <- 
  (ate_wip + 
     ggtitle("Cauchy") + 
     xlim(-.1, .1) + xlab("")) /
  (ates_c0d0[[1]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[1])^2 * I[p^"(m)"]*")")) +
     xlim(-.1, .1) + xlab("")) /
  (ates_c0d0[[2]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[2])^2 * I[p^"(m)"]*")")) + 
     xlim(-.1, .1) + xlab("")) /
  (ates_c0d0[[3]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[3])^2 * I[p^"(m)"]*")")) + 
     xlim(-.1, .1) +xlab("")) /
  (ates_c0d0[[4]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[4])^2 * I[p^"(m)"]*")")) + 
     xlim(-.1, .1)) & 
  theme(legend.position = 'none')
ate_wip_norm

ggsave("bayesian_subset_selection/actg/results/figures/ate_wip_norm.png",
       ate_wip_norm, width = 7, height = 20, units = "in", dpi = 300)

#-------------------------------------------------------------------------------
# After PSM
#-------------------------------------------------------------------------------

current_data <- actg036
hist_data <- readRDS("bayesian_subset_selection/actg/data/actg019_after_PSM.rds")

current_data$age <- (current_data$age - mean(current_data$age)) /
  (2*sd(current_data$age))
current_data$cd4 <- (current_data$cd4 - mean(current_data$cd4)) /
  (2*sd(current_data$cd4))
hist_data$age <- (hist_data$age - mean(hist_data$age)) /
  (2*sd(hist_data$age))
hist_data$cd4 <- (hist_data$cd4 - mean(hist_data$cd4)) /
  (2*sd(hist_data$cd4))


family <- binomial(link = "logit")
ate_wip_after_PSM <- plot_ate(post_samples_wip_after_PSM)

c0 <- c(0.25, 0.5, 1, 2)
ates_c0d0_after_PSM <- lapply(post_samples_c0d0_after_PSM, plot_ate)


ate_wip_norm_after_PSM <- 
  (ate_wip_after_PSM + 
     ggtitle("Cauchy") + 
     xlim(-.2, .2) + xlab("")) /
  (ates_c0d0_after_PSM[[1]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[1])^2 * I[p^"(m)"]*")")) +
     xlim(-.2, .2) + xlab("")) /
  (ates_c0d0_after_PSM[[2]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[2])^2 * I[p^"(m)"]*")")) + 
     xlim(-.2, .2) + xlab("")) /
  (ates_c0d0_after_PSM[[3]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[3])^2 * I[p^"(m)"]*")")) + 
     xlim(-.2, .2) +xlab("")) /
  (ates_c0d0_after_PSM[[4]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[4])^2 * I[p^"(m)"]*")")) + 
     xlim(-.2, .2)) &
  theme(legend.position = 'bottom')
ate_wip_norm_after_PSM

ggsave("bayesian_subset_selection/actg/results/figures/ate_wip_norm_after_PSM.png",
       ate_wip_norm_after_PSM, width = 7, height = 20, units = "in", dpi = 300)

#-------------------------------------------------------------------------------
# Compare ATE distributions
#-------------------------------------------------------------------------------

models_before_after_PSM <- 
  (models_wip_norm) |
  (models_wip_norm_after_PSM) +
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom',
        legend.justification = "right",
        legend.box.just = "right")
models_before_after_PSM <- models_before_after_PSM + 
  plot_annotation(
    title = "Comparison of models before and after PSM",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
  )

ates_wip_norm_after_PSM <- 
  (ate_wip_norm) |
  (ate_wip_norm_after_PSM) +
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom',
        legend.justification = "right",
        legend.box.just = "right")

ates_wip_norm_after_PSM <- ates_wip_norm_after_PSM +
  plot_annotation(
    title = "Comparison of ATE distributions before and after PSM",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
  )
ates_wip_norm_after_PSM
ggsave("bayesian_subset_selection/actg/results/figures/ates_wip_norm_before_after_PSM.png",
       ates_wip_norm_after_PSM, width = 12, height = 14, units = "in", dpi = 300)
