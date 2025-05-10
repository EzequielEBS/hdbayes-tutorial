setwd("C:/Users/Ezequiel/OneDrive - Fundacao Getulio Vargas - FGV/Grad MAp FGV/proj_hdbayes/hdbayes-tutorials/bayesian_subset_selection")

# load libraries
library(tidyverse)
library(ggplot2)
library(viridis)
library(reshape2)

# load data
metrics <- read_csv("results/metrics.csv")
weights <- read_csv("results/weights.csv")

# barplot of weights
weights_reshape <- melt(weights)

barplot_weights <- ggplot(weights_reshape, aes(x = variables, y = value)) +
  geom_bar(stat = "identity", aes(fill = variable), position = "dodge") +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Weights") +
  xlab("Variables") +
  ylab("Weights") +
  labs(fill = "")

# plot weights
plot_weight_aic_vs_logml <- ggplot(weights, aes(x=AIC, y=logml)) + 
  geom_point() +
  ggtitle("Weights: AIC vs LogML") +
  xlab("Weight: AIC") +
  ylab("Weight: LogML") 
  # xlim(0.01485, 0.01675)
plot_weight_bic_vs_logml <- ggplot(weights, aes(x=BIC, y=logml)) + 
  geom_point() +
  ggtitle("Weights: BIC vs LogML") +
  xlab("Weight: BIC") +
  ylab("Weight: LogML") 
  # xlim(0.01485, 0.01675)

# plot metrics
plot_aic_vs_logml <- ggplot(metrics, aes(x=AIC, y=logml)) + 
  geom_point() +
  ggtitle("AIC vs LogML") +
  xlab("AIC") +
  ylab("LogML") +
  xlim(140, 173)

plot_bic_vs_logml <- ggplot(metrics, aes(x=BIC, y=logml)) +
  geom_point() +
  ggtitle("BIC vs LogML") +
  xlab("BIC") +
  ylab("LogML") +
  xlim(140, 173)

# save plots
ggsave("results/figures/weights_aic_vs_logml.png", plot_weight_aic_vs_logml)
ggsave("results/figures/weights_bic_vs_logml.png", plot_weight_bic_vs_logml)
ggsave("results/figures/aic_vs_logml.png", plot_aic_vs_logml)
ggsave("results/figures/bic_vs_logml.png", plot_bic_vs_logml)
