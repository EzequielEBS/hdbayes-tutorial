setwd("C:/Users/Ezequiel/OneDrive - Fundacao Getulio Vargas - FGV/Grad MAp FGV/proj_hdbayes/hdbayes-tutorials/bayesian_subset_selection/actg")

library(ggplot2)
library(patchwork)
library(reshape2)
library(MCMCpack)

# load auxiliary functions
source("functions.R")

# load samples
load("data/post_samples.RData")
load("data/post_samples_ctrl.RData")
load("data/post_samples_trt.RData")

# normalize data
age_stats <- with(actg036,
                  c('mean' = mean(age), 'sd' = sd(age)))
cd4_stats <- with(actg036,
                  c('mean' = mean(cd4), 'sd' = sd(cd4)))
actg036$age <- ( actg036$age - age_stats['mean'] ) / age_stats['sd']
actg019$age <- ( actg019$age - age_stats['mean'] ) / age_stats['sd']
actg036$cd4 <- ( actg036$cd4 - cd4_stats['mean'] ) / cd4_stats['sd']
actg019$cd4 <- ( actg019$cd4 - cd4_stats['mean'] ) / cd4_stats['sd']

# define data
current_data <- actg036
hist_data <- actg019

# split data
current_data_ctrl <- current_data[current_data$treatment == 0, 
                                  !(names(current_data) %in% c("treatment"))]
current_data_trt <- current_data[current_data$treatment == 1, 
                                 !(names(current_data) %in% c("treatment"))]

# set parameters
data <- list(current_data, hist_data)
data_ctrl <- list(current_data_ctrl, hist_data)
data_trt <- list(current_data_trt, hist_data)
family <- binomial(link = "logit")

# calculate ATE
mean_models_ctrl <- mean_models_arm(current_data_ctrl, "outcome", family, post_samples_ctrl$post_betam)
mean_models_trt <- mean_models_arm(current_data_trt, "outcome", family, post_samples_trt$post_betam)
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

# Save plot
ggsave("figures/posterior_distribution_ate.png",
       plot_ate, width = 10, height = 7, units = "in", dpi = 300)

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
  ggplot(post_beta_models_trt[[i]], aes(x = treatment)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
    labs(title = paste("Model:", cov_models_trt[i]), x = "", y = "") +
    xlim(c(-2.5, 1.5)) +
    theme_gray()
})

# Combine all plots using patchwork for models with treatment effect
combined_plot_post_trt <- Reduce(`+`, plots_post_trt) + 
  plot_annotation(title = "Posterior distribution of treatment effect in models") +
  theme_gray()
combined_plot_post_trt

# Save plot
ggsave("figures/posterior_distribution_treatment_effect.png",
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
    ggplot(post_beta_models_wtrt[[i]], aes(x = !!sym(covariates[j]))) +
      geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
      labs(title = covariates[j], x = "", y = "") +
      theme_gray()
  })
  wrap_elements(Reduce(`+`, plots_i) + plot_layout(ncol = length(plots_i)) + 
                  plot_annotation(title = paste("Model:", cov_models_wtrt[i])) +
                  theme_gray())
})

# Combine all plots using patchwork for models without treatment effect
combined_plot_post_wtrt <- Reduce(`+`, plots_post_wtrt) + plot_layout(ncol = 2) +
  plot_annotation(title = "Posterior distribution of model covariates without treatment effect") +
  theme_gray()
combined_plot_post_wtrt

# Save plot
ggsave("figures/posterior_distribution_covariates_without_treatment_effect.png",
       combined_plot_post_wtrt, width = 20, height = 10, units = "in", dpi = 300)

########################################################################################
# Plot means by arm
########################################################################################


# Convert matrices to dataframes and reshape them
df_mean_models_ctrl <- as.data.frame(mean_models_ctrl)
df_mean_models_trt <- as.data.frame(mean_models_trt)

# Add identifiers for matrices
df_mean_models_ctrl$arm <- "crtl"
df_mean_models_trt$arm <- "trt"

# Combine the data
df_mean_models <- rbind(
  melt(df_mean_models_ctrl, id.vars = "arm"),  # Convert to long format
  melt(df_mean_models_trt, id.vars = "arm")
)

models <- paste0("Model: ", post_samples_ctrl$df_post$model)
titles_post_ctrl <- paste("Posterior (ctrl): ", format(round(post_samples_ctrl$df_post$post_model, digits = 3), 
                                                       nsmall = 3))
titles_post_trt <- paste("Posterior (trt): ", format(round(post_samples_trt$df_post$post_model, digits = 3),
                                                     nsmall = 3))
titles_post <- paste(models, titles_post_ctrl, titles_post_trt, sep = "\n")


df_mean_models$variable <- titles_post[match(df_mean_models$variable, 
                                             paste0("V", 1:(ncol(df_mean_models_ctrl)-1)), 
                                             nomatch = 0)]

# Plot histograms
plot_means_arm <- ggplot(df_mean_models, aes(x = value, fill = arm)) +
  geom_histogram(alpha = 0.6, bins = 30, position = "identity") +
  facet_wrap(~ variable, 
             scales = "free",
             ncol = 2) +  # Separate plots for each column
  labs(
    title = "Posterior distribution of marginal means by arm",
    x = "",
    y = ""
  ) +
  scale_fill_manual(values = c("blue", "red")) +
  xlim(range(df_mean_models$value)) +
  theme_gray()

# Save plot
ggsave("figures/posterior_distribution_means_by_arm.png",
       plot_means_arm, width = 15, height = 20, units = "in", dpi = 300, scale = 0.5)

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
  plot_annotation(title = "Posterior distribution of model covariates by arm - Part 1") +
  theme_gray()
combined_plot_post_beta_2 <- Reduce(`+`, plots_post_beta2) +
  plot_annotation(title = "Posterior distribution of model covariates by arm - Part 2") +
  theme_gray()
combined_plot_post_beta_3 <- Reduce(`+`, plots_post_beta3) +
  plot_annotation(title = "Posterior distribution of model covariates by arm - Part 3") +
  theme_gray()

combined_plot_post_beta_1
combined_plot_post_beta_2
combined_plot_post_beta_3

# Save plot
ggsave("figures/posterior_distribution_covariates_by_arm1.png",
       combined_plot_post_beta_1, width = 20, height = 7, units = "in", dpi = 300)
ggsave("figures/posterior_distribution_covariates_by_arm2.png",
       combined_plot_post_beta_2, width = 20, height = 14, units = "in", dpi = 300)
ggsave("figures/posterior_distribution_covariates_by_arm3.png",
       combined_plot_post_beta_3, width = 7, height = 21, units = "in", dpi = 300)
