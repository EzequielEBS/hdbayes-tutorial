load("bayesian_subset_selection/actg/samples/post_samples_d0.RData") # psi0
load("bayesian_subset_selection/actg/samples/post_samples_c0.RData") # tau0
load("bayesian_subset_selection/actg/samples/post_samples_m0v0.RData")
load("bayesian_subset_selection/actg/samples/post_samples_cauchy.RData")
load("bayesian_subset_selection/actg/samples/post_samples_norm.RData")

library(ggplot2)
library(MCMCpack)
library(bayestestR)

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

ate_wip <- plot_ate(post_samples_cauchy)
ate_norm <- plot_ate(post_samples)

ate_wip_norm <- (ate_wip + ggtitle("Cauchy prior on coefficients")) /
  (ate_norm + ggtitle("Normal prior on coefficients"))

ate_wip + ggtitle("Cauchy prior on coefficients")
ggsave("bayesian_subset_selection/actg/results/figures/ate_cauchy.png",
       ate_cauchy, width = 10, height = 7, units = "in", dpi = 300)

ate_d0 <- lapply(post_samples_d0, function(x) {
  plot_ate(x)
})

all_ate_d0 <- (ate_d0[[1]] + ggtitle(expression(psi[0] == 0.25))+ xlim(-0.1, 0.05) + xlab("")) /
  (ate_d0[[2]] + ggtitle(expression(psi[0] == 0.5)) + xlim(-0.1, 0.05) + xlab("")) /
  (ate_d0[[3]] + ggtitle(expression(psi[0] == 1)) + xlim(-0.1, 0.05) +xlab("")) /
  (ate_d0[[4]] + ggtitle(expression(psi[0] == 2)) + xlim(-0.1, 0.05)) +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')
all_ate_d0

ggsave("bayesian_subset_selection/actg/results/figures/ate_psi0.png",
       all_ate_d0, width = 10, height = 25, units = "in", dpi = 300)
