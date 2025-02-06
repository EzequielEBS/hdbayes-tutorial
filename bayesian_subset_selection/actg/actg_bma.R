setwd("C:/Users/Ezequiel/OneDrive - Fundacao Getulio Vargas - FGV/Grad MAp FGV/proj_hdbayes/hdbayes-tutorials/bayesian_subset_selection/actg")

# load libraries
library(bayestestR)
library(ggplot2)
library(dplyr)
library(MCMCpack)
library(bayestestR)
library(dplyr)

# load auxiliary functions
source("functions.R")

# load samples
load("post_samples.RData")
load("post_samples_ctrl.RData")
load("post_samples_trt.RData")

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

# calculate marginal means by arm
mean_models_ctrl <- mean_models_arm(current_data_ctrl, "outcome", family, post_samples_ctrl$post_betam)
mean_models_trt <- mean_models_arm(current_data_trt, "outcome", family, post_samples_trt$post_betam)

bma_ctrl <- bma(mean_models_ctrl, post_samples_ctrl$df_post, 10000)
bma_trt <- bma(mean_models_trt, post_samples_trt$df_post, 10000)

df_bma_arm <- data.frame(
  value = c(bma_ctrl, bma_trt),
  group = c(rep("ctrl", length(bma_ctrl)), rep("trt", length(bma_trt)))
)

# Plot the histograms
bma_marg_means <- ggplot(df_bma_arm, aes(x = value, fill = group)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(title = "", x = "Marginal Mean", y = "") +
  theme_gray() +
  scale_fill_manual(values = c("blue", "red"))

# Save the plot
ggsave("figures/bma_marg_means.png", bma_marg_means, width = 8, height = 6, units = "in", dpi = 300)

df_bma <- data.frame(
  value = bma_trt - bma_ctrl
)

ci_hdi_bma_95 <- ci(df_bma$value, method = "HDI", ci = 0.95)
ci_hdi_bma_90 <- ci(df_bma$value, method = "HDI", ci = 0.90)

hdi_95_lower <- ci_hdi_bma_95$CI_low
hdi_95_upper <- ci_hdi_bma_95$CI_high
hdi_90_lower <- ci_hdi_bma_90$CI_low
hdi_90_upper <- ci_hdi_bma_90$CI_high

# Compute density for the entire dataset
density_data <- density(df_bma$value) # Compute density
density_df <- data.frame(x = density_data$x, y = density_data$y) # Convert to data frame

# Filter density data for the HDI regions
density_95 <- density_df %>% filter(x >= hdi_95_lower & x <= hdi_95_upper)
density_90 <- density_df %>% filter(x >= hdi_90_lower & x <= hdi_90_upper)

# Compute the mean, median, and quartiles
mean_value <- mean(df_bma$value)
median_value <- median(df_bma$value)
q1_value <- quantile(df_bma$value, 0.25)
q3_value <- quantile(df_bma$value, 0.75)

# Define colors for mean, median, and quartiles
stats_colors <- c("90%" = "beige", 
                  "95%" = "bisque4", 
                  "Mean" = "blueviolet",
                  "Median" = "#83E22B",
                  "Q1" = "#FFA500",
                  "Q3" = "#FFA500")
                  
                  

# Plot using ggplot2
bma_ate <- ggplot() +
  # Full density curve
  geom_line(data = density_df, aes(x = x, y = y), color = "skyblue", size = 1) +
  
  # Highlight 95% HDI region
  geom_area(data = density_95, aes(x = x, y = y, fill = "95%"), alpha = 0.3) +
  
  # Highlight 90% HDI region
  geom_area(data = density_90, aes(x = x, y = y, fill = "90%"), alpha = 0.3) +
  
  # Add vertical dashed lines for 95% HDI
  # geom_vline(aes(xintercept = hdi_95_lower, color = "95% HDI"), linetype = "dashed", size = 1) +
  # geom_vline(aes(xintercept = hdi_95_upper, color = "95% HDI"), linetype = "dashed", size = 1) +
  # 
  # # Add vertical dashed lines for 90% HDI
  # geom_vline(aes(xintercept = hdi_90_lower, color = "90% HDI"), linetype = "dashed", size = 1) +
  # geom_vline(aes(xintercept = hdi_90_upper, color = "90% HDI"), linetype = "dashed", size = 1) +
  
  # Add lines for mean, median, and quartiles
  geom_vline(aes(xintercept = mean_value, color = "Mean"), linetype = "solid", size = 1) + 
  geom_vline(aes(xintercept = median_value, color = "Median"), linetype = "solid", size = 1) + 
  geom_vline(aes(xintercept = q1_value, color = "Q1"), linetype = "dotted", size = 1) + 
  geom_vline(aes(xintercept = q3_value, color = "Q3"), linetype = "dotted", size = 1) +
  
  # Add a legend for both color and fill
  scale_fill_manual(values = stats_colors) +
  scale_color_manual(values = stats_colors) +
  
  # Add labels
  labs(title = "",
       x = "Average treatment effect (ATE)", y = "", fill = "BCI", color = "Statistics") +
  
  theme_gray()

# save the plot
ggsave("figures/bma_ate.png", bma_ate, width = 8, height = 6, units = "in", dpi = 300)

# summary of the BMA

df_bma_summary <- data.frame(
  statistic = c("Mean", "Median", "Q1", "Q3", "90% BCI", "95% BCI"),
  value = c(format(round(mean_value, 3), nsmall = 3),
            format(round(median_value, 3), nsmall = 3),
            format(round(q1_value, 3), nsmall = 3),
            format(round(q3_value, 3), nsmall = 3),
            paste0("[", format(round(hdi_90_lower, 3), nsmall = 3), ", ", format(round(hdi_90_upper, 3), nsmall = 3), "]"),
            paste0("[", format(round(hdi_95_lower, 3), nsmall = 3), ", ", format(round(hdi_95_upper, 3), nsmall = 3), "]")
  )
)

print(xtable::xtable(df_bma_summary, caption = "BMA summary (ATE)"), 
      type = "latex", include.rownames = FALSE)
