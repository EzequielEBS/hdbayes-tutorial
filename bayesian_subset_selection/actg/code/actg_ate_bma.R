#---------------------------------------------------------------------------------------
# This script generates plots for the average treatment effect (ATE) and odds ratio (OR)
# using Bayesian Model Averaging (BMA) results from the ACTG dataset.
#---------------------------------------------------------------------------------------

# load libraries
library(bayestestR)
library(ggplot2)
library(dplyr)
library(MCMCpack)
library(patchwork)
library(reshape2)
library(hdbayes)

# load auxiliary functions
source("bayesian_subset_selection/actg/code/aux_scripts/functions.R")

# load samples
load("bayesian_subset_selection/actg/samples/post_samples.RData")
load("bayesian_subset_selection/actg/samples/mean_models_ctrl.RData")
load("bayesian_subset_selection/actg/samples/mean_models_trt.RData")
load("bayesian_subset_selection/actg/samples/bma_ctrl.RData")
load("bayesian_subset_selection/actg/samples/bma_trt.RData")

# Define colors
blended_rgb <- round(colMeans(rbind(
  c(135, 206, 235),
  c(70, 130, 180)
)))
blended_color <- rgb(blended_rgb[1], blended_rgb[2], blended_rgb[3], maxColorValue = 255)
stats_colors <- c("90% \nBCI" = blended_color,
                  "Density" = "skyblue",
                  "ATE = 0" = "black",
                  "OR = 1" = "black"
)

# Create a data frame for the ATE
df_bma <- data.frame(
  value = bma_trt - bma_ctrl
)
# Compute credibility intervals
ci_bma_95 <- bayestestR::ci(df_bma$value, ci = 0.95)
ci_bma_90 <- bayestestR::ci(df_bma$value, ci = 0.90)
ci_ate_95_lower <- ci_bma_95$CI_low
ci_ate_95_upper <- ci_bma_95$CI_high
ci_ate_90_lower <- ci_bma_90$CI_low
ci_ate_90_upper <- ci_bma_90$CI_high

# Compute density for the entire dataset
ate_density_data <- density(df_bma$value)
ate_density_df <- data.frame(x = ate_density_data$x, y = ate_density_data$y) 

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
  geom_area(data = ate_density_df %>% filter(x <= ci_ate_90_lower + 0.0003), 
            aes(x = x, y = y, fill = "Density"), color = "black",
            alpha = 0.7) +
  geom_area(data = ate_density_df %>% filter(x >= ci_ate_90_upper - 0.0003), 
            aes(x = x, y = y, fill = "Density"), color = "black",
            alpha = 0.7) +
  geom_area(data = ate_density_90, 
            aes(x = x, y = y, fill = "90% \nBCI"), color = "black",
            alpha = 0.7) +
  geom_vline(aes(xintercept = 0, color = "ATE = 0"), linetype = "solid", size = 1) + 
  scale_fill_manual(name = NULL, values = stats_colors, breaks = c("90% \nBCI")) +
  scale_color_manual(name = NULL, values = stats_colors, guide = NULL) +
  labs(title = "",
       x = "Average treatment effect (ATE)", y = "", color = "") +
  theme_bw() +
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
  ) +
  theme(text = element_text(size = 16),        
        axis.title = element_text(size = 18),  
        axis.text = element_text(size = 16),   
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16))
# save the plot
ggsave("bayesian_subset_selection/actg/results/figures/bma_ate.png",
       bma_ate, width = 8, height = 6, units = "in", dpi = 300)

# create BMA data frames
df_bma_arm <- data.frame(
  value = c(bma_ctrl, bma_trt),
  group = c(rep("ctrl", length(bma_ctrl)), rep("trt", length(bma_trt)))
)
# create OR data frame
df_or <- data.frame(or = bma_trt / (1-bma_trt) / bma_ctrl / (1-bma_ctrl))

# compute credibility intervals for the odds ratio
ci_or_90 <- bayestestR::ci(df_or$or, ci = 0.90)
ci_or_90_lower <- ci_or_90$CI_low
ci_or_90_upper <- ci_or_90$CI_high

# Compute density for the entire dataset
or_density_data <- density(df_or$or)
or_density_df <- data.frame(x = or_density_data$x, y = or_density_data$y)

# Filter density data for the HDI regions
or_density_90 <- or_density_df %>% filter(x >= ci_or_90_lower & x <= ci_or_90_upper)

# Define colors for mean, median, and quartiles
stats_colors <- c("90% \nBCI" = blended_color,
                  "Density" = "skyblue",
                  "ATE = 0" = "black",
                  "OR = 1" = "black"
                  )

# Plot the density of the odds ratio
plot_or <- ggplot() +
  geom_area(data = or_density_df %>% filter(x <= ci_or_90_lower + 0.01), 
            aes(x = x, y = y, fill = "Density"), color = "black",
            alpha = 0.7) +
  geom_area(data = or_density_df %>% filter(x >= ci_or_90_upper - 0.01), 
            aes(x = x, y = y, fill = "Density"), color = "black",
            alpha = 0.7) +
  geom_area(data = or_density_90, 
            aes(x = x, y = y, fill = "90% \nBCI"), color = "black",
            alpha = 0.7) +
  geom_vline(aes(xintercept = 1, colour = "OR = 1"), linetype = "dotted", size = 1)  +
  labs(title = "",
       x = "Odds ratio", y = "") +
  scale_fill_manual(name = NULL, values = stats_colors, breaks = c("90% \nBCI")) +
  scale_color_manual(name = NULL, values = stats_colors, guide = NULL) +
  theme_bw() +
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
  ) +
  theme(text = element_text(size = 16),        
        axis.title = element_text(size = 18),  
        axis.text = element_text(size = 16),   
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16))
# Save plot
ggsave("bayesian_subset_selection/actg/results/figures/posterior_distribution_or.png",
       plot_or, width = 10, height = 7, units = "in", dpi = 300)

# put or and bma_ate together
bma_ate_or <- (bma_ate + theme(legend.position = "none")) + 
  plot_or + plot_layout(ncol = 2) &
  theme(
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 14)
  )
# save the plot
ggsave("bayesian_subset_selection/actg/results/figures/bma_ate_or.png",
       bma_ate_or, width = 14, height = 8, units = "in", dpi = 300)

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


#---------------------------------------------------------------------------------------
# Plot means by arm
#---------------------------------------------------------------------------------------


# Convert matrices to dataframes and reshape them
df_mean_models_ctrl <- as.data.frame(mean_models_ctrl)
df_mean_models_trt <- as.data.frame(mean_models_trt)
# Get model titles
models <- post_samples$df_post$model
titles_post <- models
# Add identifiers for matrices
df_mean_models_ctrl$arm <- "Ctrl"
df_mean_models_trt$arm <- "Trt"
# Combine the data
df_mean_models <- rbind(
  melt(df_mean_models_ctrl, id.vars = "arm"), 
  melt(df_mean_models_trt, id.vars = "arm")
)
# Map variable names to titles
df_mean_models$variable <- titles_post[match(df_mean_models$variable, 
                                             paste0("V", 1:(ncol(df_mean_models_ctrl)-1)), 
                                             nomatch = 0)]
# Create annotation data frame for probabilities
ann <- data.frame(
  variable = titles_post, 
  x = rep(0.12, length(titles_post)), 
  y = rep(35, length(titles_post)), 
  label = paste("Prob:", 
                format(round(post_samples$df_post$post_model, digits = 3), nsmall = 3))
)

# Plot
plot_means_arm <- ggplot(df_mean_models, aes(x = value, fill = arm)) +
  geom_density(alpha = 0.6, position = "identity", color = "black") +
  facet_wrap(~ factor(variable, levels = titles_post), 
             scales = "fixed",
             ncol = 3) +
  geom_label(
    data = ann,
    aes(x, y, label = label),
    inherit.aes = FALSE,
    fill = "white", 
    color = "black",
    hjust = 0
  ) +
  labs(
    title = "",
    x = "Average effect",
    y = "",
    fill = "Arm"
  ) +
  scale_fill_manual(values = c("#66A8D0", "#D08E66")) +
  xlim(floor(range(df_mean_models$value)*1e3)/1e3) +
  theme_bw() +
  theme(
    legend.position = c(0.97, 0.05),  
    legend.justification = c(1, 0),
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "black"),
  ) +
  theme(text = element_text(size = 12),        
        axis.title = element_text(size = 14),  
        axis.text = element_text(size = 12),   
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 11))

# Plot top 6 models
plot_means_arm_top6 <- ggplot(df_mean_models %>%
                                filter(variable %in% c("age, treatment, cd4", 
                                                       "age, treatment, race, cd4",
                                                       "treatment, cd4",
                                                       "treatment, race, cd4",
                                                       "age, cd4",
                                                       "age, race, cd4"
                                )), 
                              aes(x = value, fill = arm)) +
  geom_density(alpha = 0.7, position = "identity", color = "black") +
  facet_wrap(~ factor(variable,
                      levels = c("age, treatment, cd4", 
                                 "age, treatment, race, cd4",
                                 "treatment, cd4",
                                 "treatment, race, cd4",
                                 "age, cd4",
                                 "age, race, cd4"
                                 )
                      ), 
             scales = "fixed",
             ncol = 3) +
  geom_label(
    data = ann %>%
      filter(variable %in% c("age, treatment, cd4", 
                             "age, treatment, race, cd4",
                             "treatment, cd4",
                             "treatment, race, cd4",
                             "age, cd4",
                             "age, race, cd4"
      )),
    aes(x, y, label = label),
    inherit.aes = FALSE,
    fill = "white",   
    color = "black",  
    hjust = 0,
    size = 6
  ) +
  labs(
    title = "",
    x = "Average effect",
    y = "",
    fill = "Arm"
  ) +
  scale_fill_manual(values = c("#66A8D0", "#D08E66")) +
  xlim(floor(range(df_mean_models$value)*1e3)/1e3) +
  theme_bw() +
  theme(
    legend.position = c(0.97, 0.05),   
    legend.justification = c(1, 0),
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "black"),
  ) +
  theme(text = element_text(size = 11),        
        axis.title = element_text(size = 14),  
        axis.text = element_text(size = 12),   
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 16))

# Save plots
ggsave("bayesian_subset_selection/actg/results/figures/posterior_distribution_means_by_arm.png",
       plot_means_arm, width = 11, height = 14, units = "in", dpi = 300)
ggsave("bayesian_subset_selection/actg/results/figures/posterior_distribution_means_by_arm_top6.png",
       plot_means_arm_top6, width = 14, height = 8, units = "in", dpi = 300)

#---------------------------------------------------------------------------------------
# Plot or by arm
#---------------------------------------------------------------------------------------

# Compute odds ratios for each model
df_mean_models_or <- as.data.frame(
  do.call(cbind,
  lapply(1:(ncol(df_mean_models_trt)-1), function(j) {
  or <- df_mean_models_trt[,j]/(1-df_mean_models_trt[,j]) / 
    df_mean_models_ctrl[,j]/(1-df_mean_models_ctrl[,j])
})
)
)
# Reshape the data
df_mean_models_or <- melt(df_mean_models_or)
# Map variable names to titles
df_mean_models_or$variable <- titles_post[match(df_mean_models_or$variable, 
                               paste0("V", 1:(ncol(df_mean_models_trt)-1)), 
                               nomatch = 0)]
# Create annotation data frame for probabilities
ann_or <- df_mean_models_or %>%
  group_by(variable) %>%
  summarise(
    x = 1.8,          
    y = max(density(value)$y) - max(density(value)$y)/5
  )
ann_or <- merge(ann_or, 
      data.frame(label = paste("Prob:", 
                               format(round(post_samples$df_post$post_model, digits = 3), nsmall = 3)),
                 variable = post_samples$df_post$model),
      by = "variable")
# Plot
plot_or_models <- ggplot(df_mean_models_or, aes(x = value)) +
  geom_density(alpha = 0.6, position = "identity", fill = "#66A8D0") +
  geom_vline(aes(xintercept = 1), linetype = "dotted", size = 1) +
  facet_wrap(~ factor(variable, levels = titles_post), 
             scales = "free_y",
             ncol = 3) +  
  geom_label(
    data = ann_or,
    aes(x, y, label = label),
    inherit.aes = FALSE,
    fill = "white",   
    color = "black",  
    hjust = 0
  ) +
  labs(
    title = "",
    x = "Odds ratio",
    y = ""
  ) +
  xlim(floor(range(df_mean_models_or$value)*1e3)/1e3) +
  theme_bw() +
  theme(
    legend.position = c(0.95, 0.05),   
    legend.justification = c(1, 0),
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "black"),
  ) +
  theme(text = element_text(size = 12),        
        axis.title = element_text(size = 14),  
        axis.text = element_text(size = 12),   
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 11))
# Save plot
ggsave("bayesian_subset_selection/actg/results/figures/posterior_distribution_or_by_model.png",
       plot_or_models, width = 11, height = 14, units = "in", dpi = 300)
