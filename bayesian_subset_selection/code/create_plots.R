setwd("C:/Users/Ezequiel/OneDrive - Fundacao Getulio Vargas - FGV/Grad MAp FGV/proj_hdbayes/hdbayes-tutorials/bayesian_subset_selection")

# load libraries
library(tidyverse)
library(ggplot2)
library(viridis)

# load data
metrics <- read_csv("results/metrics.csv")
weights <- read_csv("results/weights.csv")

head(metrics)
head(weights)

# plot weights

plot_aic_vs_logml <- ggplot(weights, aes(x=AIC, y=logml)) + 
  geom_point() +
  ggtitle("Weights: AIC vs LogML") +
  xlab("Weight: AIC") +
  ylab("Weight: LogML") +
  xlim(0.01485, 0.01675)
plot_bic_vs_logml <- ggplot(weights, aes(x=BIC, y=logml)) + 
  geom_point() +
  ggtitle("Weights: BIC vs LogML") +
  xlab("Weight: BIC") +
  ylab("Weight: LogML") +
  xlim(0.01485, 0.01675)


# save plots
ggsave("results/figures/weights_aic_vs_logml.png", plot_aic_vs_logml)
ggsave("results/figures/weights_bic_vs_logml.png", plot_bic_vs_logml)
