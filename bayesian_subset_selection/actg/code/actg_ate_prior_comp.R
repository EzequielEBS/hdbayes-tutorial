load("bayesian_subset_selection/actg/samples/post_samples_c0d0.RData")
load("bayesian_subset_selection/actg/samples/post_samples_wip.RData")
# load("bayesian_subset_selection/actg/samples/post_samples_norm.RData")

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

current_data$treatment <- current_data$treatment - mean(current_data$treatment)
current_data$race <- current_data$race - mean(current_data$race)
hist_data$treatment <- hist_data$treatment - mean(hist_data$treatment)
hist_data$race <- hist_data$race - mean(hist_data$race)

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
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[1])^2 * I[p]*")")) +
     xlim(-.1, .1) + xlab("")) /
  (ates_c0d0[[2]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[2])^2 * I[p]*")")) + 
     xlim(-.1, .1) + xlab("")) /
  (ates_c0d0[[3]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[3])^2 * I[p]*")")) + 
     xlim(-.1, .1) +xlab("")) /
  (ates_c0d0[[4]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[4])^2 * I[p]*")")) + 
     xlim(-.1, .1)) +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ate_wip_norm

ggsave("bayesian_subset_selection/actg/results/figures/ate_wip_norm.png",
       ate_wip_norm, width = 7, height = 20, units = "in", dpi = 300)

# ate_c0 <-
#   (ates_c0[[1]] + ggtitle(bquote(tau[0] == .(c0[1]))) + xlim(-.1, .1) + xlab("")) /
#   (ates_c0[[2]] + ggtitle(bquote(tau[0] == .(c0[2]))) + xlim(-.1, .1) + xlab("")) /
#   (ates_c0[[3]] + ggtitle(bquote(tau[0] == .(c0[3]))) + xlim(-.1, .1) +xlab("")) /
#   (ates_c0[[4]] + ggtitle(bquote(tau[0] == .(c0[4]))) + xlim(-.1, .1)) +
#   plot_layout(guides = "collect") & theme(legend.position = 'bottom')
# ate_c0
# 
# ggsave("bayesian_subset_selection/actg/results/figures/ate_tau0.png",
#        ate_c0, width = 7, height = 20, units = "in", dpi = 300)

# ate_m0v0 <-
#   (ates_m0v0[[1]] + 
#    ggtitle(bquote((delta[0]*', '*lambda[0]) == (.(m0_v0[[1]][1])*', '*.(m0_v0[[1]][2])))) +
#    xlim(-.1, .1) + xlab("")) /
#   (ates_m0v0[[2]] +
#    ggtitle(bquote((delta[0]*', '*lambda[0]) == (.(m0_v0[[2]][1])*', '*.(m0_v0[[2]][2])))) +
#    xlim(-.1, .1) + xlab("")) /
#   (ates_m0v0[[3]] +
#    ggtitle(bquote((delta[0]*', '*lambda[0]) == (.(m0_v0[[3]][1])*', '*.(m0_v0[[3]][2])))) +
#    xlim(-.1, .1) +xlab("")) /
#   (ates_m0v0[[4]] +
#    ggtitle(bquote((delta[0]*', '*lambda[0]) == (.(m0_v0[[4]][1])*', '*.(m0_v0[[4]][2])))) +
#    xlim(-.1, .1)) +
#   plot_layout(guides = "collect") & theme(legend.position = 'bottom')
# ate_m0v0
  

# ggsave("bayesian_subset_selection/actg/results/figures/ate_m0v0.png",
#        ate_m0v0, width = 7, height = 20, units = "in", dpi = 300)

