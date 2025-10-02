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
load("bayesian_subset_selection/actg/samples/post_samples_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/mean_models_ctrl_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/mean_models_trt_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/bma_ctrl_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/bma_trt_after_PSM.RData")

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

df_bma_after_PSM <- data.frame(
  value = bma_trt_after_PSM - bma_ctrl_after_PSM
)

ci_bma_90_after_PSM <- bayestestR::ci(df_bma_after_PSM$value, ci = 0.90)
ci_ate_90_lower_after_PSM <- ci_bma_90_after_PSM$CI_low
ci_ate_90_upper_after_PSM <- ci_bma_90_after_PSM$CI_high

# Compute density for the entire dataset
ate_density_data_after_PSM <- density(df_bma_after_PSM$value) # Compute density
ate_density_df_after_PSM <- data.frame(x = ate_density_data_after_PSM$x, 
                                       y = ate_density_data_after_PSM$y) # Convert to data frame

# Filter density data for the HDI regions
ate_density_90_after_PSM <- ate_density_df_after_PSM %>% 
  filter(x >= ci_ate_90_lower_after_PSM & x <= ci_ate_90_upper_after_PSM)

# Plot using ggplot2
bma_ate_after_PSM <- ggplot() +
  geom_area(data = ate_density_df_after_PSM %>% 
              filter(x <= ci_ate_90_lower_after_PSM + 
                       (ci_ate_90_lower_after_PSM - 
                          min(ate_density_df_after_PSM$x))*0.007), 
            aes(x = x, y = y, fill = "Density"), color = "black",
            alpha = 0.7) +
  geom_area(data = ate_density_df_after_PSM %>% 
              filter(x >= ci_ate_90_upper_after_PSM - 
                       (max(ate_density_df_after_PSM$x) - 
                          ci_ate_90_upper_after_PSM)*0.007), 
            aes(x = x, y = y, fill = "Density"), color = "black",
            alpha = 0.7) +
  
  # # Highlight 95% HDI region
  geom_area(data = ate_density_90_after_PSM, 
            aes(x = x, y = y, fill = "90% \nBCI"), color = "black",
            alpha = 0.7) +
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
  
  theme_bw() +
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
  ) +
  theme(text = element_text(size = 16),        # Base text size
        axis.title = element_text(size = 18),  # Axis titles
        axis.text = element_text(size = 16),   # Axis tick labels
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16))

bma_ate_after_PSM

# save the plot
ggsave("bayesian_subset_selection/actg/results/figures/bma_ate_after_PSM.png",
       bma_ate_after_PSM, width = 8, height = 6, units = "in", dpi = 300)

# create BMA data frames
df_bma_arm_after_PSM <- data.frame(
  value = c(bma_ctrl_after_PSM, bma_trt_after_PSM),
  group = c(rep("ctrl", length(bma_ctrl_after_PSM)), 
            rep("trt", length(bma_trt_after_PSM)))
)

df_or_after_PSM <- data.frame(or = bma_trt_after_PSM / (1-bma_trt_after_PSM) / 
                      bma_ctrl_after_PSM / (1-bma_ctrl_after_PSM))

# compute credible intervals for the odds ratio
ci_or_after_PSM_90 <- bayestestR::ci(df_or_after_PSM$or, ci = 0.90)
ci_or_after_PSM_90_lower <- ci_or_after_PSM_90$CI_low
ci_or_after_PSM_90_upper <- ci_or_after_PSM_90$CI_high

# Compute density for the entire dataset
or_density_data_after_PSM <- density(df_or_after_PSM$or) # Compute density
or_density_df_after_PSM <- data.frame(x = or_density_data_after_PSM$x, y = or_density_data_after_PSM$y) # Convert to data frame

# Filter density data for the HDI regions
or_density_90_after_PSM <- or_density_df_after_PSM %>% filter(x >= ci_or_after_PSM_90_lower & x <= ci_or_after_PSM_90_upper)

# Plot the density of the odds ratio
plot_or_after_PSM <- ggplot() +
  # geom_density(data = or_density_df,
  #              aes(x = x, y = y),
  #              color = "black") +
  geom_area(data = or_density_df_after_PSM %>% 
              filter(x <= ci_or_after_PSM_90_lower + 
                       (ci_or_after_PSM_90_lower - 
                        min(or_density_df_after_PSM$x))*0.01), 
            aes(x = x, y = y, fill = "Density"), color = "black",
            alpha = 0.7) +
  geom_area(data = or_density_df_after_PSM %>% 
              filter(x >= ci_or_after_PSM_90_upper - 
                       (max(or_density_df_after_PSM$x) - 
                          ci_or_after_PSM_90_upper)*0.005), 
            aes(x = x, y = y, fill = "Density"), color = "black",
            alpha = 0.7) +
  geom_area(data = or_density_90_after_PSM, 
            aes(x = x, y = y, fill = "90% \nBCI"), color = "black",
            alpha = 0.7) +
  geom_vline(aes(xintercept = 1, colour = "OR = 1"), linetype = "dotted", size = 1)  +
  labs(title = "",
       x = "Odds ratio", y = "") +
  # Add a legend for both color and fill
  scale_fill_manual(name = NULL, values = stats_colors, breaks = c("90% \nBCI")) +
  scale_color_manual(name = NULL, values = stats_colors, guide = NULL) +
  # scale_color_manual(name = NULL, values = stats_colors) +
  theme_bw() +
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
  ) +
  theme(text = element_text(size = 16),        # Base text size
        axis.title = element_text(size = 18),  # Axis titles
        axis.text = element_text(size = 16),   # Axis tick labels
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16))

plot_or_after_PSM

# Save plot
ggsave("bayesian_subset_selection/actg/results/figures/posterior_distribution_or_after_PSM.png",
       plot_or_after_PSM, width = 10, height = 7, units = "in", dpi = 300)

# put or and bma_ate together
bma_ate_or_after_PSM <- (bma_ate_after_PSM + theme(legend.position = "none")) + 
  plot_or_after_PSM + plot_layout(ncol = 2) &
  theme(
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 14)
  )
bma_ate_or_after_PSM

# save the plot
ggsave("bayesian_subset_selection/actg/results/figures/bma_ate_or_after_PSM.png",
       bma_ate_or_after_PSM, width = 14, height = 8, units = "in", dpi = 300)

########################################################################################
# Plot means by arm
########################################################################################


# Convert matrices to dataframes and reshape them
df_mean_models_ctrl_after_PSM <- as.data.frame(mean_models_ctrl_after_PSM)
df_mean_models_trt_after_PSM <- as.data.frame(mean_models_trt_after_PSM)

# models <- paste0("Covariates: ", post_samples_ctrl$df_post$model)
# titles_post_ctrl <- paste("Ctrl:", format(round(post_samples_ctrl$df_post$post_model, digits = 3), 
#                                                        nsmall = 3))
# titles_post_trt <- paste("Trt:", format(round(post_samples_trt$df_post$post_model, digits = 3),
#                                                      nsmall = 3))
# models <- paste0("Covariates: ", post_samples$df_post$model)
models_after_PSM <- post_samples_after_PSM$df_post$model
# titles_post_ctrl <- paste("Ctrl:", format(round(post_samples$df_post$post_model, digits = 3), 
#                                           nsmall = 3))
# titles_post_trt <- paste("Trt:", format(round(post_samples$df_post$post_model, digits = 3),
#                                         nsmall = 3))
# titles_post <- paste(models, titles_post_ctrl, titles_post_trt, sep = " | ")
# titles_post <- paste(models, 
#                      paste("Prob:", 
#                            format(round(post_samples$df_post$post_model, digits = 3), nsmall = 3)),
#                      sep = " | ")
titles_post_after_PSM <- models_after_PSM

# Add identifiers for matrices
df_mean_models_ctrl_after_PSM$arm <- "Ctrl"
df_mean_models_trt_after_PSM$arm <- "Trt"

# Combine the data
df_mean_models_after_PSM <- rbind(
  melt(df_mean_models_ctrl_after_PSM, id.vars = "arm"),  # Convert to long format
  melt(df_mean_models_trt_after_PSM, id.vars = "arm")
)

df_mean_models_after_PSM$variable <- titles_post_after_PSM[match(df_mean_models_after_PSM$variable, 
                                             paste0("V", 1:(ncol(df_mean_models_ctrl_after_PSM)-1)), 
                                             nomatch = 0)]

ann_after_PSM <- data.frame(
  variable = titles_post_after_PSM, 
  x = rep(0.12, length(titles_post_after_PSM)), # position inside plot
  y = rep(35, length(titles_post_after_PSM)), # position inside plot
  label = paste("Prob:", 
                format(round(post_samples_after_PSM$df_post$post_model, digits = 3), nsmall = 3))
)

# Plot histograms
plot_means_arm_after_PSM <- ggplot(df_mean_models_after_PSM, aes(x = value, fill = arm)) +
  geom_density(alpha = 0.6, position = "identity", color = "black") +
  facet_wrap(~ factor(variable, levels = titles_post_after_PSM), 
             scales = "fixed",
             ncol = 3) +  # Separate plots for each column
  geom_label(
    data = ann_after_PSM,
    aes(x, y, label = label),
    inherit.aes = FALSE,
    fill = "white",   # box color
    color = "black",    # text color
    hjust = 0
  ) +
  labs(
    title = "",
    x = "Average effect",
    y = "",
    fill = "Arm"
  ) +
  scale_fill_manual(values = c("#66A8D0", "#D08E66")) +
  xlim(floor(range(df_mean_models_after_PSM$value)*1e3)/1e3) +
  theme_bw() +
  theme(
    legend.position = c(0.97, 0.05),   # (x, y) coordinates — bottom right
    legend.justification = c(1, 0),
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "black"),
  ) +
  theme(text = element_text(size = 12),        # Base text size
        axis.title = element_text(size = 14),  # Axis titles
        axis.text = element_text(size = 12),   # Axis tick labels
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 11))

plot_means_arm_after_PSM

# plot_means_arm_top6 <- ggplot(df_mean_models %>%
#                                 filter(variable %in% c("age, treatment, cd4", 
#                                                        "age, treatment, race, cd4",
#                                                        "treatment, cd4",
#                                                        "treatment, race, cd4",
#                                                        "age, cd4",
#                                                        "age, race, cd4"
#                                 )), 
#                               aes(x = value, fill = arm)) +
#   geom_density(alpha = 0.7, position = "identity", color = "black") +
#   facet_wrap(~ factor(variable,
#                       levels = c("age, treatment, cd4", 
#                                  "age, treatment, race, cd4",
#                                  "treatment, cd4",
#                                  "treatment, race, cd4",
#                                  "age, cd4",
#                                  "age, race, cd4"
#                       )
#   ), 
#   scales = "fixed",
#   ncol = 3) +  # Separate plots for each column
#   geom_label(
#     data = ann %>%
#       filter(variable %in% c("age, treatment, cd4", 
#                              "age, treatment, race, cd4",
#                              "treatment, cd4",
#                              "treatment, race, cd4",
#                              "age, cd4",
#                              "age, race, cd4"
#       )),
#     aes(x, y, label = label),
#     inherit.aes = FALSE,
#     fill = "white",   # box color
#     color = "black",    # text color
#     hjust = 0,
#     size = 6
#   ) +
#   labs(
#     title = "",
#     x = "Average effect",
#     y = "",
#     fill = "Arm"
#   ) +
#   scale_fill_manual(values = c("#66A8D0", "#D08E66")) +
#   xlim(floor(range(df_mean_models$value)*1e3)/1e3) +
#   theme_bw() +
#   theme(
#     legend.position = c(0.97, 0.05),   # (x, y) coordinates — bottom right
#     legend.justification = c(1, 0),
#     legend.title = element_text(size = 10),
#     legend.background = element_rect(fill = "white", color = "black"),
#   ) +
#   theme(text = element_text(size = 11),        # Base text size
#         axis.title = element_text(size = 14),  # Axis titles
#         axis.text = element_text(size = 12),   # Axis tick labels
#         legend.title = element_text(size = 14),
#         legend.text = element_text(size = 12),
#         strip.text = element_text(size = 16))
# 
# plot_means_arm_top6

# Save plot
ggsave("bayesian_subset_selection/actg/results/figures/posterior_distribution_means_by_arm_after_PSM.png",
       plot_means_arm_after_PSM, width = 11, height = 14, units = "in", dpi = 300)
# Save plot
# ggsave("bayesian_subset_selection/actg/results/figures/posterior_distribution_means_by_arm_top6.png",
#        plot_means_arm_top6, width = 14, height = 8, units = "in", dpi = 300)

########################################################################################
# Plot or by arm
########################################################################################

df_mean_models_or_after_PSM <- as.data.frame(
  do.call(cbind,
          lapply(1:(ncol(df_mean_models_trt_after_PSM)-1), function(j) {
            or <- df_mean_models_trt_after_PSM[,j]/(1-df_mean_models_trt_after_PSM[,j]) / 
              df_mean_models_ctrl_after_PSM[,j]/(1-df_mean_models_ctrl_after_PSM[,j])
          })
  )
)

df_mean_models_or_after_PSM <- melt(df_mean_models_or_after_PSM)
df_mean_models_or_after_PSM$variable <- titles_post_after_PSM[match(df_mean_models_or_after_PSM$variable, 
                                                paste0("V", 1:(ncol(df_mean_models_trt_after_PSM)-1)), 
                                                nomatch = 0)]

ann_or_after_PSM <- df_mean_models_or_after_PSM %>%
  group_by(variable) %>%
  summarise(
    x = 1.8,          # fixed x
    y = max(density(value)$y) - max(density(value)$y)/5# y depends on max value of the facet
  )
ann_or_after_PSM <- merge(ann_or_after_PSM, 
                data.frame(label = paste("Prob:", 
                                         format(round(post_samples_after_PSM$df_post$post_model, digits = 3), nsmall = 3)),
                           variable = post_samples_after_PSM$df_post$model),
                by = "variable")

plot_or_models_after_PSM <- ggplot(df_mean_models_or_after_PSM, aes(x = value)) +
  geom_density(alpha = 0.6, position = "identity", fill = "#66A8D0") +
  geom_vline(aes(xintercept = 1), linetype = "dotted", size = 1) +
  facet_wrap(~ factor(variable, levels = titles_post_after_PSM), 
             scales = "free_y",
             ncol = 3) +  # Separate plots for each column
  geom_label(
    data = ann_or_after_PSM,
    aes(x, y, label = label),
    inherit.aes = FALSE,
    fill = "white",   # box color
    color = "black",    # text color
    hjust = 0
  ) +
  labs(
    title = "",
    x = "Odds ratio",
    y = ""
  ) +
  xlim(floor(range(df_mean_models_or_after_PSM$value)*1e3)/1e3) +
  theme_bw() +
  theme(
    legend.position = c(0.95, 0.05),   # (x, y) coordinates — bottom right
    legend.justification = c(1, 0),
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "black"),
  ) +
  theme(text = element_text(size = 12),        # Base text size
        axis.title = element_text(size = 14),  # Axis titles
        axis.text = element_text(size = 12),   # Axis tick labels
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 11))

plot_or_models_after_PSM

ggsave("bayesian_subset_selection/actg/results/figures/posterior_distribution_or_by_model_after_PSM.png",
       plot_or_models_after_PSM, width = 11, height = 14, units = "in", dpi = 300)
