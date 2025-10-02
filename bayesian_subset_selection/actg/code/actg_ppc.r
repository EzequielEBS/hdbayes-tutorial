load("bayesian_subset_selection/actg/samples/draws_npp_norm.RData")
load("bayesian_subset_selection/actg/samples/draws_npp_wip.RData")

library(tidyverse)
library(ggridges)
library(hdbayes)
library(parallel)
library(ggplot2)
library(pROC)
library(ggplotify)

current_data <- actg036
hist_data <- actg019

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


plots_pnew <- mclapply(1:nrow(pnew.wip), function(i){
  df_pnew <- data.frame(
    norm = pnew.norm[i, ],
    wip = pnew.wip[i, ]
  ) %>%
    pivot_longer(cols = everything(), names_to = "prior", values_to = "prob") %>%
    mutate(prior_label = case_when(
      prior == "wip" ~ "Cauchy",
      prior == "norm" ~ "Normal"
    ))
  
  df_pnew$prior_label <- factor(df_pnew$prior_label, levels = c("Normal", "Cauchy"))
  
  blended_rgb <- round(colMeans(rbind(
    c(135, 206, 235),
    c(70, 130, 180)
  )))
  
  blended_color <- rgb(blended_rgb[1], blended_rgb[2], blended_rgb[3], maxColorValue = 255)
  
  plot <- ggplot(df_pnew, aes(x = prob, y = prior_label)) +
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
    theme_bw() +
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
  ggsave(filename, plot = plot, width = 6, height = 4, units = "in", dpi = 300)
  return(plot)
},
mc.cores = 15)

pnew_13_21 <- (plots_pnew[[13]] + 
                 scale_y_discrete(expand = c(0, 0)) +
                 scale_x_continuous(expand = c(0, 0)) +
                 xlim(c(0,1)) ) + 
              (plots_pnew[[21]] + 
                 scale_y_discrete(expand = c(0, 0)) +
                 scale_x_continuous(expand = c(0, 0)) +
                 xlim(c(0,1)) &
                 theme(axis.text.y = element_blank())) +
  plot_layout(guides = "collect") 
pnew_13_21

ggsave("bayesian_subset_selection/actg/results/figures/ppc/ppc_pnew_13_21.png",
       pnew_13_21, width = 10, height = 4, units = "in", dpi = 300)

pnew_2_10 <- (plots_pnew[[2]] + 
                 scale_y_discrete(expand = c(0, 0)) +
                 scale_x_continuous(expand = c(0, 0)) +
                 xlim(c(0,1)) ) + 
  (plots_pnew[[10]] + 
     scale_y_discrete(expand = c(0, 0)) +
     scale_x_continuous(expand = c(0, 0)) +
     xlim(c(0,1)) &
     theme(axis.text.y = element_blank())) +
  plot_layout(guides = "collect")
pnew_2_10

ggsave("bayesian_subset_selection/actg/results/figures/ppc/ppc_pnew_2_10.png",
       pnew_2_10, width = 10, height = 4, units = "in", dpi = 300)

roc.curve.wip <- roc(current_data$outcome, rowMeans(pnew.wip))
# Compute the confidence interval for the AUC
ci.auc.wip <- ci.auc(roc.curve.wip)
ci.roc.wip <- ci.se(roc.curve.wip, specificities = seq(0, 1, l = 25))
# Optionally, you can display the AUC with its CI directly on the plot
plot_roc.wip <- as.ggplot(function() {
  plot(roc.curve.wip,
       col = blended_color,
       main = "Cauchy",
       print.auc = FALSE)
  plot(ci.roc.wip, type = "shape", col = rgb(0.2, 0.4, 0.6, 0.2))
  text(0.6, 0.2, 
       labels = sprintf("AUC = %.3f (%.3f-%.3f)", 
                        roc.curve.wip$auc, ci.auc.wip[1], ci.auc.wip[3]),
       cex = 1.2)
})

roc.curves.wip <- lapply(1:ncol(pnew.wip), function(i) {
  roc(current_data$outcome, pnew.wip[, i])
})

plot_roc.curves.wip <- as.ggplot(function() {
  plot(roc.curves.wip[[1]], col = blended_color, main = "Cauchy")
  for(i in 2:length(roc.curves.wip)) {
    plot(roc.curves.wip[[i]], col = blended_color, add = TRUE)
  }
})

roc.curve.norm <- roc(current_data$outcome, rowMeans(pnew.norm),
                      ci = TRUE)
ci.auc.norm <- ci.auc(roc.curve.norm)
ci.roc.norm <- ci.se(roc.curve.norm, specificities = seq(0, 1, l = 25))

plot_roc.norm <- as.ggplot(function() {
  plot(roc.curve.norm, col = blended_color, main = "Normal", print.auc = FALSE)
  plot(ci.roc.norm, type = "shape", col = rgb(0.2, 0.4, 0.6, 0.2))
  text(0.6, 0.2, 
       labels = sprintf("AUC = %.3f (%.3f-%.3f)", 
                        roc.curve.norm$auc, ci.auc.norm[1], ci.auc.norm[3]),
       cex = 1.2)
})

roc.curves.norm <- lapply(1:ncol(pnew.norm), function(i) {
  roc(current_data$outcome, pnew.norm[, i])
})

plot_roc.curves.norm <- as.ggplot(function() {
  plot(roc.curves.norm[[1]], col = blended_color, main = "Normal")
  for(i in 2:length(roc.curves.norm)) {
    plot(roc.curves.norm[[i]], col = blended_color, add = TRUE)
  }
})

roc_wip_norm <- (plot_roc.wip) + (plot_roc.norm)
roc_wip_norm
ggsave("bayesian_subset_selection/actg/results/figures/ppc/ppc_roc_wip_norm.png",
       roc_wip_norm, width = 15, height = 8, units = "in", dpi = 300)

roc_all_wip_norm <- (plot_roc.curves.wip) + (plot_roc.curves.norm)
roc_all_wip_norm
ggsave("bayesian_subset_selection/actg/results/figures/ppc/ppc_roc_all_wip_norm.png",
       roc_all_wip_norm, width = 15, height = 8, units = "in", dpi = 300)


#-------------------------------------------------------------------------------
# After PSM
#-------------------------------------------------------------------------------

load("bayesian_subset_selection/actg/samples/draws_npp_norm_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/draws_npp_wip_after_PSM.RData")

current_data <- actg036
hist_data <- readRDS("bayesian_subset_selection/actg/data/actg019_after_PSM.rds")

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
d.sub.norm_after_PSM <- 
  draws.norm_after_PSM[sample(x = seq_len(nrow(draws.norm_after_PSM)),
                                          size = nsim, replace = F), ]
d.sub.wip_after_PSM  <- 
  draws.wip_after_PSM[sample(x = seq_len(nrow(draws.wip_after_PSM)), 
                             size = nsim, replace = F), ]

## generate new outcome y
pnew.norm_after_PSM <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d.sub.norm_after_PSM[i, ])
  p <- binomial('logit')$linkinv(Xnew %*% beta.sim)
  return(p)
})
pnew.wip_after_PSM <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d.sub.wip_after_PSM[i, ])
  p <- binomial('logit')$linkinv(Xnew %*% beta.sim)
  return(p)
})


plots_pnew_after_PSM <- mclapply(1:nrow(pnew.wip_after_PSM), function(i){
  df_pnew <- data.frame(
    norm = pnew.norm_after_PSM[i, ],
    wip = pnew.wip_after_PSM[i, ]
  ) %>%
    pivot_longer(cols = everything(), names_to = "prior", values_to = "prob") %>%
    mutate(prior_label = case_when(
      prior == "wip" ~ "Cauchy",
      prior == "norm" ~ "Normal"
    ))
  
  df_pnew$prior_label <- factor(df_pnew$prior_label, levels = c("Normal", "Cauchy"))
  
  blended_rgb <- round(colMeans(rbind(
    c(135, 206, 235),
    c(70, 130, 180)
  )))
  
  blended_color <- rgb(blended_rgb[1], blended_rgb[2], blended_rgb[3], maxColorValue = 255)
  
  plot <- ggplot(df_pnew, aes(x = prob, y = prior_label)) +
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
    theme_bw() +
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
  filename <- 
    paste0("bayesian_subset_selection/actg/results/figures/ppc/ppc_pnew_", i, 
           "_after_PSM.png")
  ggsave(filename, plot = plot, width = 6, height = 4, units = "in", dpi = 300)
  return(plot)
},
mc.cores = 15)

pnew_13_21_after_PSM <- (plots_pnew_after_PSM[[13]] + 
                 scale_y_discrete(expand = c(0, 0)) +
                 scale_x_continuous(expand = c(0, 0)) +
                 xlim(c(0,1)) ) + 
  (plots_pnew_after_PSM[[21]] + 
     scale_y_discrete(expand = c(0, 0)) +
     scale_x_continuous(expand = c(0, 0)) +
     xlim(c(0,1)) &
     theme(axis.text.y = element_blank())) +
  plot_layout(guides = "collect") 
pnew_13_21_after_PSM

ggsave("bayesian_subset_selection/actg/results/figures/ppc/ppc_pnew_13_21_after_PSM.png",
       pnew_13_21_after_PSM, width = 10, height = 4, units = "in", dpi = 300)

pnew_2_10_after_PSM <- (plots_pnew_after_PSM[[2]] + 
                scale_y_discrete(expand = c(0, 0)) +
                scale_x_continuous(expand = c(0, 0)) +
                xlim(c(0,1)) ) + 
  (plots_pnew_after_PSM[[10]] + 
     scale_y_discrete(expand = c(0, 0)) +
     scale_x_continuous(expand = c(0, 0)) +
     xlim(c(0,1)) &
     theme(axis.text.y = element_blank())) +
  plot_layout(guides = "collect")
pnew_2_10_after_PSM

ggsave("bayesian_subset_selection/actg/results/figures/ppc/ppc_pnew_2_10_after_PSM.png",
       pnew_2_10_after_PSM, width = 10, height = 4, units = "in", dpi = 300)

roc.curve.wip_after_PSM<- roc(current_data$outcome, rowMeans(pnew.wip_after_PSM))
# Compute the confidence interval for the AUC
ci.auc.wip_after_PSM <- ci.auc(roc.curve.wip_after_PSM)
ci.roc.wip_after_PSM <- ci.se(roc.curve.wip_after_PSM, specificities = seq(0, 1, l = 25))
# Optionally, you can display the AUC with its CI directly on the plot
plot_roc.wip_after_PSM <- as.ggplot(function() {
  plot(roc.curve.wip_after_PSM,
       col = blended_color,
       main = "Cauchy",
       print.auc = FALSE)
  plot(ci.roc.wip, type = "shape", col = rgb(0.2, 0.4, 0.6, 0.2))
  text(0.6, 0.2, 
       labels = sprintf("AUC = %.3f (%.3f-%.3f)", 
                        roc.curve.wip_after_PSM$auc, 
                        ci.auc.wip_after_PSM[1], 
                        ci.auc.wip_after_PSM[3]),
       cex = 1.2)
})

roc.curves.wip_after_PSM <- lapply(1:ncol(pnew.wip_after_PSM), function(i) {
  roc(current_data$outcome, pnew.wip_after_PSM[, i])
})

plot_roc.curves.wip_after_PSM <- as.ggplot(function() {
  plot(roc.curves.wip_after_PSM[[1]], col = blended_color, main = "Cauchy")
  for(i in 2:length(roc.curves.wip_after_PSM)) {
    plot(roc.curves.wip_after_PSM[[i]], col = blended_color, add = TRUE)
  }
})

roc.curve.norm_after_PSM <- roc(current_data$outcome, rowMeans(pnew.norm_after_PSM),
                      ci = TRUE)
ci.auc.norm_after_PSM <- ci.auc(roc.curve.norm_after_PSM)
ci.roc.norm_after_PSM <- ci.se(roc.curve.norm_after_PSM, specificities = seq(0, 1, l = 25))

plot_roc.norm_after_PSM <- as.ggplot(function() {
  plot(roc.curve.norm_after_PSM, col = blended_color, main = "Normal", print.auc = FALSE)
  plot(ci.roc.norm_after_PSM, type = "shape", col = rgb(0.2, 0.4, 0.6, 0.2))
  text(0.6, 0.2, 
       labels = sprintf("AUC = %.3f (%.3f-%.3f)", 
                        roc.curve.norm_after_PSM$auc, 
                        ci.auc.norm_after_PSM[1], 
                        ci.auc.norm_after_PSM[3]),
       cex = 1.2)
})

roc.curves.norm_after_PSM <- lapply(1:ncol(pnew.norm_after_PSM), function(i) {
  roc(current_data$outcome, pnew.norm_after_PSM[, i])
})

plot_roc.curves.norm_after_PSM <- as.ggplot(function() {
  plot(roc.curves.norm_after_PSM[[1]], col = blended_color, main = "Normal")
  for(i in 2:length(roc.curves.norm_after_PSM)) {
    plot(roc.curves.norm_after_PSM[[i]], col = blended_color, add = TRUE)
  }
})

roc_wip_norm_after_PSM <- (plot_roc.wip_after_PSM) + (plot_roc.norm_after_PSM)
roc_wip_norm_after_PSM
ggsave("bayesian_subset_selection/actg/results/figures/ppc/ppc_roc_wip_norm_after_PSM.png",
       roc_wip_norm_after_PSM, width = 15, height = 8, units = "in", dpi = 300)

roc_all_wip_norm_after_PSM <- (plot_roc.curves.wip_after_PSM) + 
  (plot_roc.curves.norm_after_PSM)
roc_all_wip_norm_after_PSM
ggsave("bayesian_subset_selection/actg/results/figures/ppc/ppc_roc_all_wip_norm_after_PSM.png",
       roc_all_wip_norm_after_PSM, width = 15, height = 8, units = "in", dpi = 300)
