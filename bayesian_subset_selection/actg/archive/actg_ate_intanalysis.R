library(ggplot2)
library(patchwork)
library(reshape2)
library(MCMCpack)
library(broom)
library(tidyverse)
library(marginaleffects)
library(hdbayes)

# load auxiliary functions
source("code/functions.R")

# load samples
load("data/post_samples.RData")
load("data/post_samples_ctrl.RData")
load("data/post_samples_trt.RData")

# define data
current_data <- actg036
hist_data <- actg019

# split data
current_data_ctrl <- current_data[current_data$treatment == 0, 
                                  !(names(current_data) %in% c("treatment"))]
current_data_trt <- current_data[current_data$treatment == 1, 
                                 !(names(current_data) %in% c("treatment"))]

# normalize data
age_stats_ctrl <- with(current_data_ctrl,
                       c('mean' = mean(age), 'sd' = sd(age)))
cd4_stats_ctrl <- with(current_data_ctrl,
                       c('mean' = mean(cd4), 'sd' = sd(cd4)))
age_stats_trt <- with(current_data_trt,
                       c('mean' = mean(age), 'sd' = sd(age)))
cd4_stats_trt <- with(current_data_trt,
                       c('mean' = mean(cd4), 'sd' = sd(cd4)))
age_stats_hist <- with(hist_data,
                       c('mean' = mean(age), 'sd' = sd(age)))
cd4_stats_hist <- with(hist_data,
                       c('mean' = mean(cd4), 'sd' = sd(cd4)))

current_data_ctrl$age <- (current_data_ctrl$age - age_stats_ctrl['mean']) / age_stats_ctrl['sd']
current_data_ctrl$cd4 <- (current_data_ctrl$cd4 - cd4_stats_ctrl['mean']) / cd4_stats_ctrl['sd']
current_data_trt$age <- (current_data_trt$age - age_stats_trt['mean']) / age_stats_trt['sd']
current_data_trt$cd4 <- (current_data_trt$cd4 - cd4_stats_trt['mean']) / cd4_stats_trt['sd']
hist_data$age <- (hist_data$age - age_stats_hist['mean']) / age_stats_hist['sd']
hist_data$cd4 <- (hist_data$cd4 - cd4_stats_hist['mean']) / cd4_stats_hist['sd']

# empirical ate
mean_ctrl <- mean(current_data_ctrl$outcome)
mean_trt <- mean(current_data_trt$outcome)
raw_ate <- mean_trt - mean_ctrl

# set parameters
data <- list(current_data, hist_data)
data_ctrl <- list(current_data_ctrl, hist_data)
data_trt <- list(current_data_trt, hist_data)
family <- binomial(link = "logit")

# REWRITE
# calculate ATE
mean_models_ctrl <- mean_models_arm(current_data_ctrl, "outcome", "treatment", family, post_samples_ctrl$post_betam)
mean_models_trt <- mean_models_arm(current_data_trt, "outcome", "treatment", family, post_samples_trt$post_betam)
mean_ctrl <- mean_arm(post_samples_ctrl$df_post$post_model, mean_models_ctrl)
mean_trt <- mean_arm(post_samples_trt$df_post$post_model, mean_models_trt)
ate <- mean_trt - mean_ctrl

########################################################################################
# Plot the ate distribution
########################################################################################

# Create a data frame with the ATE values
df_ate <- data.frame(
  ate = ate
)

# Plot the histogram
plot_ate <- ggplot(df_ate, aes(x = ate)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Posterior distribution of the average treatment effect (ATE)", x = "ATE", y = "") +
  theme_gray()

plot_ate

df_or <- data.frame(
  or = mean_ctrl / (1-mean_ctrl) / 
    (mean_trt / (1-mean_trt))
)

# Plot the histogram
plot_or <- ggplot(df_or, aes(x = or)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Posterior distribution of the odds ratio", x = "OR", y = "") +
  theme_gray()

plot_or

# Save plot
ggsave("results/figures/posterior_distribution_or.png",
       plot_or, width = 10, height = 7, units = "in", dpi = 300)

########################################################################################
# Plot posterior distribution of model covariates
########################################################################################


# Create a filter to separate models with and without treatment effect
filter_models_trt <- unlist(lapply(post_samples$post_betam, function(df) {
  if ("treatment" %in% colnames(df)) TRUE else FALSE
}))

# Extract covariates and posterior samples for models with treatment effect
cov_models_trt <- post_samples$df_post[filter_models_trt, "model"]
post_beta_models_trt <- post_samples$post_betam[filter_models_trt]

# Create individual plots for models with treatment effect
plots_post_trt <- lapply(seq_along(post_beta_models_trt), function(i) {
  ggplot(exp(post_beta_models_trt[[i]]), aes(x = treatment)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
    labs(title = paste("Model:", cov_models_trt[i]), x = "", y = "") +
    xlim(c(0, 1.5)) +
    theme_gray() +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black")
})

# Combine all plots using patchwork for models with treatment effect
combined_plot_post_trt <- Reduce(`+`, plots_post_trt) + 
  plot_annotation(title = "Posterior distribution of the odds ratio of treatment effect") +
  theme_gray()
combined_plot_post_trt

# Save plot
ggsave("results/figures/posterior_distribution_or_treatment_effect.png",
       combined_plot_post_trt, width = 15, height = 10, units = "in", dpi = 300)

# Extract covariates and posterior samples for models without treatment effect
cov_models_wtrt <- post_samples$df_post[!filter_models_trt, "model"]
post_beta_models_wtrt <- post_samples$post_betam[!filter_models_trt]

# create individual plots for models without treatment effect
plots_post_wtrt <- lapply(seq_along(post_beta_models_wtrt), function(i) {
  covariates <- colnames(post_beta_models_wtrt[[i]])[!(colnames(post_beta_models_wtrt[[i]]) %in% 
                                                         c("lp__", 
                                                           "a0_hist_1", 
                                                           "logit_a0s[1]",
                                                           ".chain", ".iteration", ".draw"))]
  plots_i <- lapply(seq_along(covariates), function(j) {
    ggplot(exp(post_beta_models_wtrt[[i]]), aes(x = !!sym(covariates[j]))) +
      geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
      labs(title = covariates[j], x = "", y = "") +
      theme_gray() +
      geom_vline(xintercept = 1, linetype = "dashed", color = "black")
  })
  wrap_elements(Reduce(`+`, plots_i) + plot_layout(ncol = length(plots_i)) + 
                  plot_annotation(title = paste("Model:", cov_models_wtrt[i])) +
                  theme_gray())
})

# Combine all plots using patchwork for models without treatment effect
combined_plot_post_wtrt <- Reduce(`+`, plots_post_wtrt) + plot_layout(ncol = 2) +
  plot_annotation(title = "Posterior distribution of the odds ratio (models without treatment effect)") +
  theme_gray()
combined_plot_post_wtrt

# Save plot
ggsave("results/figures/posterior_distribution_or_without_treatment_effect.png",
       combined_plot_post_wtrt, width = 20, height = 10, units = "in", dpi = 300)

########################################################################################
# Plot means by arm
########################################################################################


# Convert matrices to dataframes and reshape them
df_mean_models_ctrl <- as.data.frame(mean_models_ctrl)
df_mean_models_trt <- as.data.frame(mean_models_trt)

models <- paste0("Covariates: ", post_samples_ctrl$df_post$model)
titles_post_ctrl <- paste("Posterior (ctrl): ", format(round(post_samples_ctrl$df_post$post_model, digits = 3), 
                                                       nsmall = 3))
titles_post_trt <- paste("Posterior (trt): ", format(round(post_samples_trt$df_post$post_model, digits = 3),
                                                     nsmall = 3))
titles_post <- paste(models, titles_post_ctrl, titles_post_trt, sep = "\n")

df_odds_ratio <- df_mean_models_ctrl / (1-df_mean_models_ctrl) / 
  (df_mean_models_trt / (1-df_mean_models_trt)) 
colnames(df_odds_ratio) <- titles_post
df_odds_ratio <- melt(df_odds_ratio)


# Add identifiers for matrices
df_mean_models_ctrl$arm <- "crtl"
df_mean_models_trt$arm <- "trt"

# Combine the data
df_mean_models <- rbind(
  melt(df_mean_models_ctrl, id.vars = "arm"),  # Convert to long format
  melt(df_mean_models_trt, id.vars = "arm")
)

df_mean_models$variable <- titles_post[match(df_mean_models$variable, 
                                             paste0("V", 1:(ncol(df_mean_models_ctrl)-1)), 
                                             nomatch = 0)]

# Plot histograms
plot_means_arm <- ggplot(df_mean_models, aes(x = value, fill = arm)) +
  geom_density(alpha = 0.6, position = "identity", color = "black") +
  facet_wrap(~ factor(variable, levels = titles_post), 
             scales = "free",
             ncol = 2) +  # Separate plots for each column
  labs(
    title = "",
    x = "Average effect",
    y = "",
    fill = "Arm"
  ) +
  scale_fill_manual(values = c("blue", "red")) +
  xlim(range(df_mean_models$value)) +
  theme_bw() +
  theme(
    legend.position = c(0.95, 0.05),   # (x, y) coordinates — bottom right
    legend.justification = c(1, 0),
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "black"),
  )

plot_means_arm

# Save plot
ggsave("results/figures/posterior_distribution_means_by_arm.png",
       plot_means_arm, width = 15, height = 20, units = "in", dpi = 300, scale = 0.5)

# plot odds ratio
plot_or_by_model <- ggplot(df_odds_ratio, aes(x = value)) +
  geom_density(alpha = 0.6, position = "identity", color = "skyblue", fill = "skyblue") +
  geom_vline(aes(xintercept = 1), linetype = "dotted", size = 1) +
  facet_wrap(~ factor(variable, levels = titles_post),, 
             scales = "free",
             ncol = 2) +  # Separate plots for each column
  labs(
    title = "",
    x = "Odds ratio",
    y = ""
  ) +
  xlim(range(df_odds_ratio$value)) +
  theme_bw() 

plot_or_by_model

# Save plot
ggsave("results/figures/posterior_distribution_or_by_model.png",
       plot_or_by_model, width = 15, height = 20, units = "in", dpi = 300, scale = 0.5)

########################################################################################
# Plots of the posterior distribution of beta by arm
########################################################################################


# Extract posterior samples of beta for each arm
post_beta_ctrl <- post_samples_ctrl$post_betam
post_beta_trt <- post_samples_trt$post_betam

# Combine posterior samples of beta for each arm
post_beta_comb <- lapply(seq_along(post_beta_ctrl), function(i) {
  post_beta_ctrl_i <- as.data.frame(post_beta_ctrl[[i]])
  post_beta_trt_i <- as.data.frame(post_beta_trt[[i]])
  post_beta_ctrl_i$arm <- "crtl"
  post_beta_trt_i$arm <- "trt"
  rbind(post_beta_ctrl_i, post_beta_trt_i)
})

plots_post_beta <- lapply(seq_along(post_beta_comb), function(i) {
  covariates <- colnames(post_beta_comb[[i]])[!(colnames(post_beta_comb[[i]]) %in% 
                                                  c("lp__", 
                                                    "a0_hist_1", 
                                                    "logit_a0s[1]",
                                                    ".chain",
                                                    ".iteration", 
                                                    ".draw",
                                                    "arm",
                                                    "(Intercept)"))]
  plots_i <- lapply(seq_along(covariates), function(j) {
    ggplot(post_beta_comb[[i]], aes(x = !!sym(covariates[j]), fill = arm)) +
      geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
      labs(title = covariates[j], x = "", y = "") +
      scale_fill_manual(values = c("blue", "red")) +
      theme_gray()
  })
  wrap_elements(Reduce(`/`, plots_i) + 
                  #plot_layout(ncol = length(plots_i),
                  #            widths = rep(1, length(plots_i))) + 
                  plot_annotation(title = paste("Model:", post_samples_ctrl$df_post$model[i])) +
                  theme_gray())
})

len_cov <- lapply(seq_along(post_beta_comb), function(i) {
  covariates <- colnames(post_beta_comb[[i]])[!(colnames(post_beta_comb[[i]]) %in% 
                                                  c("lp__", 
                                                    "a0_hist_1", 
                                                    "logit_a0s[1]",
                                                    ".chain",
                                                    ".iteration", 
                                                    ".draw",
                                                    "arm",
                                                    "(Intercept)"))]
  length(covariates)
})

# Separate plots
plots_post_beta1 <- plots_post_beta[len_cov == 1]
plots_post_beta2 <- plots_post_beta[len_cov == 2]
plots_post_beta3 <- plots_post_beta[len_cov == 3]

combined_plot_post_beta_1 <- Reduce(`+`, plots_post_beta1) +
  plot_annotation(title = "Posterior distribution of covariate effects by arm - Part 1") +
  theme_gray()
combined_plot_post_beta_2 <- Reduce(`+`, plots_post_beta2) +
  plot_annotation(title = "Posterior distribution of covariate effects by arm - Part 2") +
  theme_gray()
combined_plot_post_beta_3 <- Reduce(`+`, plots_post_beta3) +
  plot_annotation(title = "Posterior distribution of covariate effects by arm - Part 3") +
  theme_gray()

combined_plot_post_beta_1
combined_plot_post_beta_2
combined_plot_post_beta_3

# Save plot
ggsave("results/figures/posterior_distribution_cov_eff_by_arm1.png",
       combined_plot_post_beta_1, width = 20, height = 7, units = "in", dpi = 300)
ggsave("results/figures/posterior_distribution_cov_eff_by_arm2.png",
       combined_plot_post_beta_2, width = 20, height = 14, units = "in", dpi = 300)
ggsave("results/figures/posterior_distribution_cov_eff_by_arm3.png",
       combined_plot_post_beta_3, width = 7, height = 21, units = "in", dpi = 300)
