#---------------------------------------------------------------------------------------
# Make ATE plots under different priors
#---------------------------------------------------------------------------------------

# Load data
load("bayesian_subset_selection/actg/samples/post_samples_c0d0.RData")
load("bayesian_subset_selection/actg/samples/post_samples_wip.RData")
load("bayesian_subset_selection/actg/samples/post_samples_c0d0_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/post_samples_wip_after_PSM.RData")

# load libraries
library(ggplot2)
library(MCMCpack)
library(bayestestR)
library(patchwork)
library(hdbayes)
library(dplyr)
library(tidyverse)

# load functions
source("bayesian_subset_selection/actg/code/aux_scripts/functions.R")

# define data
current_data <- actg036
hist_data <- actg019

# normalize age and cd4
current_data$age <- (current_data$age - mean(current_data$age)) /
  (2*sd(current_data$age))
current_data$cd4 <- (current_data$cd4 - mean(current_data$cd4)) /
  (2*sd(current_data$cd4))
hist_data$age <- (hist_data$age - mean(hist_data$age)) /
  (2*sd(hist_data$age))
hist_data$cd4 <- (hist_data$cd4 - mean(hist_data$cd4)) /
  (2*sd(hist_data$cd4))

# define hyperparameters
c0 <- c(0.25, 0.5, 1, 2)
family <- binomial(link = "logit")

# create ATE plots
ate_wip <- plot_ate(post_samples_wip)
ates_c0d0 <- lapply(post_samples_c0d0, plot_ate)


ate_wip_norm <- 
  (ate_wip + 
     ggtitle("Cauchy") + 
     xlim(-.1, .1) + xlab("") + 
     theme(legend.position = 'none')) +
  (ates_c0d0[[1]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[1])^2 * I[p^"(m)"]*")")) +
     xlim(-.1, .1) + xlab("") + 
     theme(legend.position = 'none')) +
  (ates_c0d0[[2]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[2])^2 * I[p^"(m)"]*")")) + 
     xlim(-.1, .1) + xlab("") + 
     theme(legend.position = 'none')) +
  (ates_c0d0[[3]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[3])^2 * I[p^"(m)"]*")")) + 
     xlim(-.1, .1) +
     theme(legend.position = 'none')
   ) +
  (ates_c0d0[[4]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[4])^2 * I[p^"(m)"]*")")) + 
     xlim(-.1, .1) + 
     theme(legend.position = c(0.95, 0.95),
           legend.background = element_rect(fill = "white", color = "black"))) +
  plot_layout(ncol = 2) 
ate_wip_norm

ggsave("bayesian_subset_selection/actg/results/figures/ate_wip_norm.png",
       ate_wip_norm, 
       width = 14, height = 15, units = "in", dpi = 300)

#-------------------------------------------------------------------------------
# After PSM
#-------------------------------------------------------------------------------

current_data <- actg036
hist_data <- readRDS("bayesian_subset_selection/actg/data/actg019_after_PSM2.rds")

current_data$age <- (current_data$age - mean(current_data$age)) /
  (2*sd(current_data$age))
current_data$cd4 <- (current_data$cd4 - mean(current_data$cd4)) /
  (2*sd(current_data$cd4))
hist_data$age <- (hist_data$age - mean(hist_data$age)) /
  (2*sd(hist_data$age))
hist_data$cd4 <- (hist_data$cd4 - mean(hist_data$cd4)) /
  (2*sd(hist_data$cd4))

ate_wip_after_PSM <- plot_ate(post_samples_wip_after_PSM)
ates_c0d0_after_PSM <- lapply(post_samples_c0d0_after_PSM, plot_ate)


ate_wip_norm_after_PSM <- 
  (ate_wip_after_PSM + 
     ggtitle("Cauchy") + 
     xlim(-.2, .2) + xlab("") + 
     theme(legend.position = 'none')) +
  (ates_c0d0_after_PSM[[1]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[1])^2 * I[p^"(m)"]*")")) +
     xlim(-.2, .2) + xlab("") + 
     theme(legend.position = 'none')) +
  (ates_c0d0_after_PSM[[2]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[2])^2 * I[p^"(m)"]*")")) + 
     xlim(-.2, .2) + xlab("") + 
     theme(legend.position = 'none')) +
  (ates_c0d0_after_PSM[[3]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[3])^2 * I[p^"(m)"]*")")) + 
     xlim(-.2, .2) + 
     theme(legend.position = 'none')
     ) +
  (ates_c0d0_after_PSM[[4]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[4])^2 * I[p^"(m)"]*")")) + 
     xlim(-.2, .2) +
     theme(legend.position = c(0.95, 0.95),
           legend.background = element_rect(fill = "white", color = "black"))) +
  plot_layout(ncol = 2)
ate_wip_norm_after_PSM

ggsave("bayesian_subset_selection/actg/results/figures/ate_wip_norm_after_PSM.png",
       ate_wip_norm_after_PSM, width = 14, height = 15, units = "in", dpi = 300)

#-------------------------------------------------------------------------------
# Compare ATE distributions
#-------------------------------------------------------------------------------

ates_wip_norm_after_PSM <- 
  ((ate_wip + 
      ggtitle("Before PSM", subtitle = "Cauchy") + 
      xlim(-.1, .1) + xlab("")) +
     (ates_c0d0[[1]] + 
        ggtitle(label = NULL,
                subtitle = bquote("Normal" * "(" * 0 * "," * .(c0[1])^2 * I[p^"(m)"]*")")) +
        xlim(-.1, .1) + xlab("")) +
     (ates_c0d0[[2]] + 
        ggtitle(label = NULL,
                subtitle = bquote("Normal" * "(" * 0 * "," * .(c0[2])^2 * I[p^"(m)"]*")")) + 
        xlim(-.1, .1) + xlab("")) +
     (ates_c0d0[[3]] + 
        ggtitle(label = NULL,
                subtitle = bquote("Normal" * "(" * 0 * "," * .(c0[3])^2 * I[p^"(m)"]*")")) + 
        xlim(-.1, .1) +
        xlab("") 
     ) +
     (ates_c0d0[[4]] + 
        ggtitle(label = NULL,
                subtitle = bquote("Normal" * "(" * 0 * "," * .(c0[4])^2 * I[p^"(m)"]*")")) + 
        xlim(-.1, .1)) +
     plot_layout(ncol = 1) & 
     theme(legend.position = 'none') 
   ) |
  ((((ate_wip_after_PSM + 
      ggtitle("After PSM", subtitle = "Cauchy") + 
      xlim(-.2, .2) + xlab("")) +
     (ates_c0d0_after_PSM[[1]] + 
        ggtitle(label = NULL,
                subtitle = bquote("Normal" * "(" * 0 * "," * .(c0[1])^2 * I[p^"(m)"]*")")) +
        xlim(-.2, .2) + xlab("")) +
     (ates_c0d0_after_PSM[[2]] + 
        ggtitle(label = NULL,
                subtitle = bquote("Normal" * "(" * 0 * "," * .(c0[2])^2 * I[p^"(m)"]*")")) + 
        xlim(-.2, .2) + xlab("")) +
     (ates_c0d0_after_PSM[[3]] + 
        ggtitle(label = NULL,
                subtitle = bquote("Normal" * "(" * 0 * "," * .(c0[3])^2 * I[p^"(m)"]*")")) + 
        xlim(-.2, .2) +
        xlab("")
     )) & 
      theme(legend.position = 'none') ) +
     (ates_c0d0_after_PSM[[4]] + 
        ggtitle(label = NULL,
                subtitle = bquote("Normal" * "(" * 0 * "," * .(c0[4])^2 * I[p^"(m)"]*")")) + 
        xlim(-.2, .2) +
        theme(legend.position = c(0.95, 0.95),
              legend.background = element_rect(fill = "white", color = "black"))) +
     plot_layout(ncol = 1)
   )
ates_wip_norm_after_PSM
ggsave("bayesian_subset_selection/actg/results/figures/ates_wip_norm_before_after_PSM.png",
       ates_wip_norm_after_PSM, width = 14, height = 20, units = "in", dpi = 300)
