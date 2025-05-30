setwd("C:/Users/Ezequiel/OneDrive - Fundacao Getulio Vargas - FGV/MSc_MAp_CD/hdbayes_project/hdbayes-tutorial/bayesian_subset_selection/actg")

# load libraries
library(bayestestR)
library(ggplot2)
library(dplyr)
library(MCMCpack)
library(patchwork)

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

# set parameters
data <- list(current_data, hist_data)
data_ctrl <- list(current_data_ctrl, hist_data)
data_trt <- list(current_data_trt, hist_data)
family <- binomial(link = "logit")

# calculate marginal means by arm
mean_models_ctrl <- mean_models_arm(current_data_ctrl, "outcome", "treatment", family, post_samples_ctrl$post_betam)
mean_models_trt <- mean_models_arm(current_data_trt, "outcome", "treatment", family, post_samples_trt$post_betam)

bma_ctrl <- bma(mean_models_ctrl, post_samples_ctrl$df_post, 10000)
bma_trt <- bma(mean_models_trt, post_samples_trt$df_post, 10000)

df_bma_arm <- data.frame(
  value = c(bma_ctrl, bma_trt),
  group = c(rep("ctrl", length(bma_ctrl)), rep("trt", length(bma_trt)))
)

df_or <- data.frame(or = bma_trt / (1-bma_trt) / bma_ctrl / (1-bma_ctrl))

ci_or_90 <- ci(df_or$or, ci = 0.90)

ci_or_90_lower <- ci_or_90$CI_low
ci_or_90_upper <- ci_or_90$CI_high

# Compute density for the entire dataset
or_density_data <- density(df_or$or) # Compute density
or_density_df <- data.frame(x = or_density_data$x, y = or_density_data$y) # Convert to data frame

# Filter density data for the HDI regions
or_density_90 <- or_density_df %>% filter(x >= ci_or_90_lower & x <= ci_or_90_upper)

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

# Plot the density of the odds ratio
plot_or <- ggplot() +
  # geom_density(data = df_or %>% filter(or <= ci_or_90_lower), aes(x = or), 
  #              color = "skyblue", fill = "skyblue", alpha = 0.5) +
  geom_area(data = or_density_df %>% filter(x <= ci_or_90_lower + 0.005), 
            aes(x = x, y = y, fill = "Density", colour = "Density"), 
            alpha = 0.7) +
  geom_area(data = or_density_df %>% filter(x >= ci_or_90_upper - 0.005), 
            aes(x = x, y = y, fill = "Density", colour = "Density"), 
            alpha = 0.7) +
  geom_area(data = or_density_90, 
            aes(x = x, y = y, fill = "90% \nBCI", colour = "Density"), 
            alpha = 0.7) +
  geom_vline(aes(xintercept = 1, colour = "OR = 1"), linetype = "dotted", size = 1)  +
  labs(title = "",
       x = "Odds ratio", y = "") +
  # Add a legend for both color and fill
  scale_fill_manual(name = NULL, values = stats_colors, breaks = c("90% \nBCI")) +
  scale_color_manual(name = NULL, values = stats_colors, guide = NULL) +
  # scale_color_manual(name = NULL, values = stats_colors) +
  theme_minimal() +
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
  )

plot_or

# Save plot
ggsave("results/figures/posterior_distribution_or.png",
       plot_or, width = 10, height = 7, units = "in", dpi = 300)

# Plot the histograms
bma_marg_means <- ggplot(df_bma_arm, aes(x = value, fill = group)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(title = "", x = "Marginal Mean", y = "") +
  theme_gray() +
  scale_fill_manual(values = c("blue", "red"))

bma_marg_means

# Save the plot
ggsave("results/figures/bma_marg_means.png", bma_marg_means, width = 8, height = 6, units = "in", dpi = 300)

df_bma <- data.frame(
  value = bma_trt - bma_ctrl
)

ci_bma_95 <- ci(df_bma$value, ci = 0.95)
ci_bma_90 <- ci(df_bma$value, ci = 0.90)

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

# Plot using ggplot2
bma_ate <- ggplot() +
  # Full density curve
  # geom_density(data = df_bma, aes(x = value), color = "skyblue", fill = "skyblue", alpha = 0.5) +
  geom_area(data = ate_density_df %>% filter(x <= ci_ate_90_lower + 0.0001), 
            aes(x = x, y = y, fill = "Density", colour = "Density"), 
            alpha = 0.7) +
  geom_area(data = ate_density_df %>% filter(x >= ci_ate_90_upper - 0.0003), 
            aes(x = x, y = y, fill = "Density", colour = "Density"), 
            alpha = 0.7) +
  
  # # Highlight 95% HDI region
  # geom_area(data = density_95, aes(x = x, y = y, fill = "95%"), alpha = 0.3) +
  # 
  # # Highlight 90% HDI region
  geom_area(data = ate_density_90, aes(x = x, y = y, fill = "90% \nBCI", colour = "Density"), alpha = 0.7) +
  
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
  
  theme_minimal() +
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
  )

bma_ate

# save the plot
ggsave("results/figures/bma_ate.png", bma_ate, width = 8, height = 6, units = "in", dpi = 300)

# put or and bma_ate together
bma_ate_or <- bma_ate + plot_or + plot_layout(ncol = 2, guides = 'collect')
bma_ate_or

# save the plot
ggsave("results/figures/bma_ate_or.png", bma_ate_or, width = 16, height = 6, units = "in", dpi = 300)

# summary of the BMA
df_bma_summary <- data.frame(
  statistic = c("Mean", "Median", "Q1", "Q3", "90% BCI", "95% BCI"),
  value = c(format(round(mean_value, 3), nsmall = 3),
            format(round(median_value, 3), nsmall = 3),
            format(round(q1_value, 3), nsmall = 3),
            format(round(q3_value, 3), nsmall = 3),
            paste0("[", format(round(ci_90_lower, 3), nsmall = 3), ", ", format(round(ci_90_upper, 3), nsmall = 3), "]"),
            paste0("[", format(round(ci_95_lower, 3), nsmall = 3), ", ", format(round(ci_95_upper, 3), nsmall = 3), "]")
  )
)

print(xtable::xtable(df_bma_summary, caption = "BMA summary (ATE)"), 
      type = "latex", include.rownames = FALSE)
