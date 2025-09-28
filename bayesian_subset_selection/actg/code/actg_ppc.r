load("bayesian_subset_selection/actg/samples/draws_npp_norm.RData")
load("bayesian_subset_selection/actg/samples/draws_npp_wip.RData")

library(tidyverse)
library(ggridges)
library(hdbayes)
library(parallel)
library(ggplot2)

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

data <- list(current_data, hist_data)
family <- binomial(link = "logit")
formula <- outcome ~ age + treatment + cd4

Xnew        <- stats::model.matrix(formula, current_data)

## for demonstration, only generate 200 data sets
nsim <- 1000
d.sub.norm <- draws.norm[sample(x = seq_len(nrow(draws.norm)), size = nsim, replace = F), ]
d.sub.wip  <- draws.wip[sample(x = seq_len(nrow(draws.wip)), size = nsim, replace = F), ]

## generate new outcome y
pnew.norm <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d.sub.norm[i, ])
  p <- binomial('logit')$linkinv(Xnew %*% beta.sim)
  return(p)
})
pnew.wip <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d.sub.wip[i, ])
  p <- binomial('logit')$linkinv(Xnew %*% beta.sim)
  return(p)
})
ynew.norm <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d.sub.norm[i, ])
  p <- binomial('logit')$linkinv(Xnew %*% beta.sim)
  return(rbinom(length(p), 1, p))
})
ynew.wip <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d.sub.wip[i, ])
  p <- binomial('logit')$linkinv(Xnew %*% beta.sim)
  return(rbinom(length(p), 1, p))
})

counts_sim_norm <- colSums(ynew.norm)
counts_sim_wip  <- colSums(ynew.wip)
counts_obs <- sum(hist_data$outcome)

i <- 1
annotations <- data.frame(
  x = c(3, 5, 7.5),
  y = c("A", "B", "C"),
  label = c("Peak A", "Peak B", "Peak C")
)

plots_pnew <- mclapply(1:nrow(pnew.wip), function(i){
  df_pnew <- data.frame(
    norm = pnew.norm[i, ],
    wip = pnew.wip[i, ]
  ) %>%
    pivot_longer(cols = everything(), names_to = "prior", values_to = "prob") %>%
    mutate(prior_label = case_when(
      prior == "wip" ~ "beta[0] %~% Cauchy(0,10) *','* beta[2:p] %~% Cauchy(0,2.5)",
      prior == "norm" ~ "beta %~% Normal(0,0.5^2*I[p])"
    ))
  
  df_pnew$prior_label <- ifelse(
    df_pnew$prior == "wip",
    "atop(beta[0] %~% Cauchy(0,10), beta[1:(p-1)] %~% Cauchy(0,2.5))",
    "beta %~% Normal(0, 0.5^2 * I[p])"
  )
  
  blended_rgb <- round(colMeans(rbind(
    c(135, 206, 235),
    c(70, 130, 180)
  )))
  
  blended_color <- rgb(blended_rgb[1], blended_rgb[2], blended_rgb[3], maxColorValue = 255)
  
  ggplot(df_pnew, aes(x = prob, y = prior_label)) +
    geom_density_ridges(alpha=0.6, fill = blended_color)+
    theme_ridges() +
    annotate(
      "label",
      x = Inf,            # X position of the box
      y = Inf, 
      hjust = 1.5,  # Slightly inside the plot
      vjust = 1.5,  # Slightly below the top
      label = paste0("y = ", current_data$outcome[i]),  # Dynamic label
      color = "black",      # Text color
      fill = "white",       # Box background color
      label.size = 0.4,     # Border thickness of the box
      size = 3
    ) +
    scale_y_discrete(labels = function(x) parse(text = x)) +
    theme_minimal() +
    theme(legend.position = "none",
          # legend.position = c(0.8, 0.85),
          legend.background = element_rect(fill = "white", color = "black")) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
    ) +
    xlim(c(0,1)) +
    theme(text = element_text(size = 10),        # Base text size
          axis.title = element_text(size = 10),  # Axis titles
          axis.text = element_text(size = 10),   # Axis tick labels
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 10)) +
    labs(
      x = "Probability of event",
      y = "",
    )
  filename <- paste0("bayesian_subset_selection/actg/results/figures/ppc/ppc_pnew_", i, ".png")
  ggsave(filename, width = 6, height = 4, units = "in", dpi = 300)
},
mc.cores = 15)

df_ynew <- data.frame(
  norm = counts_sim_norm,
  wip = counts_sim_wip
) %>%
  pivot_longer(cols = everything(), names_to = "prior", values_to = "count") %>%
  mutate(prior_label = case_when(
    prior == "wip" ~ "beta[0] %~% Cauchy(0,10) *','* beta[2:p] %~% Cauchy(0,2.5)",
    prior == "norm" ~ "beta %~% Normal(0,0.5^2*I[p])"
  ))


df_ynew$prior_label <- ifelse(
  df_ynew$prior == "wip",
  "atop(beta[0] %~% Cauchy(0,10), beta[1:(p-1)] %~% Cauchy(0,2.5))",
  "beta %~% Normal(0, 0.5^2 * I[p])"
)

plot_ynew <- ggplot(df_ynew, aes(x = count, y = prior_label)) +
  geom_density_ridges(alpha=0.6, stat="binline", bins=30)+
  theme_ridges() +
  
  # vertical lines (optional: same as dataset)
  geom_vline(aes(xintercept = counts_obs),
             color = "red", linewidth = 1) +
  scale_y_discrete(labels = function(x) parse(text = x)) +
  theme_minimal() +
  theme(legend.position = "none",
        # legend.position = c(0.8, 0.85),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
  ) +
  theme(text = element_text(size = 11),        # Base text size
        axis.title = element_text(size = 16),  # Axis titles
        axis.text = element_text(size = 16),   # Axis tick labels
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 16)) +
  labs(
    x = "Number of events",
    y = "",
  )

plot_ynew
