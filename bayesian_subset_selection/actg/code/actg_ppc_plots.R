# load files
load("bayesian_subset_selection/actg/samples/ppc_pnew_wip.RData")
load("bayesian_subset_selection/actg/samples/ppc_pnew_wip_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/ppc_pnew_norm.RData")
load("bayesian_subset_selection/actg/samples/ppc_pnew_norm_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/ppc_pnew_a0.RData")
load("bayesian_subset_selection/actg/samples/ppc_pnew_a0_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/ppc_roc_curve_wip.RData")
load("bayesian_subset_selection/actg/samples/ppc_roc_curve_wip_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/ppc_roc_curve_norm.RData")
load("bayesian_subset_selection/actg/samples/ppc_roc_curve_norm_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/ppc_roc_curve_a0.RData")
load("bayesian_subset_selection/actg/samples/ppc_roc_curve_a0_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/ppc_roc_curves_wip.RData")
load("bayesian_subset_selection/actg/samples/ppc_roc_curves_wip_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/ppc_roc_curves_norm.RData")
load("bayesian_subset_selection/actg/samples/ppc_roc_curves_norm_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/ppc_roc_curves_a0.RData")
load("bayesian_subset_selection/actg/samples/ppc_roc_curves_a0_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/ppc_ci_auc_wip.RData")
load("bayesian_subset_selection/actg/samples/ppc_ci_auc_wip_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/ppc_ci_auc_norm.RData")
load("bayesian_subset_selection/actg/samples/ppc_ci_auc_norm_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/ppc_ci_auc_a0.RData")
load("bayesian_subset_selection/actg/samples/ppc_ci_auc_a0_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/ppc_ci_roc_wip.RData")
load("bayesian_subset_selection/actg/samples/ppc_ci_roc_wip_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/ppc_ci_roc_norm.RData")
load("bayesian_subset_selection/actg/samples/ppc_ci_roc_norm_after_PSM.RData")
load("bayesian_subset_selection/actg/samples/ppc_ci_roc_a0.RData")
load("bayesian_subset_selection/actg/samples/ppc_ci_roc_a0_after_PSM.RData")
load("bayesian_subset_selection/actg/results/figures/ppc/plots_pnew.RData")
load("bayesian_subset_selection/actg/results/figures/ppc/plots_pnew_after_PSM.RData")
load("bayesian_subset_selection/actg/results/figures/ppc_a0/plots_pnew_a0.RData")
load("bayesian_subset_selection/actg/results/figures/ppc_a0/plots_pnew_a0_after_PSM.RData")

# load libraries
library(ggplot2)
library(RColorBrewer)
library(grid)

# load data
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

c0 <- c(0.25, 0.5, 1, 2)
a0_hyper <- list(c(1,1), c(2,2), c(1,10), c(10,1))

#-------------------------------------------------------------------------------
# Posterior predictive distributions
#-------------------------------------------------------------------------------

plots_pnew <- mclapply(1:nrow(pnew.wip), function(i){
  plot <- ggplot() +
    geom_density(data = data.frame(y = pnew.wip[i, ]), 
                 aes(x = y, color = "wip"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = pnew.norm[[1]][i, ]),
                 aes(x = y, color = "c01"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = pnew.norm[[2]][i, ]),
                 aes(x = y, color = "c02"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = pnew.norm[[3]][i, ]),
                 aes(x = y, color = "c03"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = pnew.norm[[4]][i, ]), 
                 aes(x = y, color = "c04"),
                 linewidth = 1) + 
    scale_color_manual(
      name = NULL,
      values = c(
        "wip" = "black",
        "c01"    = "#66A8D0",
        "c02"    = "#D06673",
        "c03" = "#7f7f7f",
        "c04" = "#D0C366"
      ),
      labels = c(
        "wip" = expression(Cauchy),
        "c01" = bquote(Normal(0, .(c0[1])^2 * I[p^(m)])),
        "c02" = bquote(Normal(0, .(c0[2])^2 * I[p^(m)])),
        "c03" = bquote(Normal(0, .(c0[3])^2 * I[p^(m)])),
        "c04" = bquote(Normal(0, .(c0[4])^2 * I[p^(m)]))
      ),
      breaks = c("wip",
                 "c01",
                 "c02",
                 "c03",
                 "c04"
                 )
    ) +
    annotation_custom(
      grid::textGrob(
        label = paste0("Obs outcome: ", current_data$outcome[i]),
        x = unit(0.81, "npc"),   # 81% from left
        y = unit(0.96, "npc"),   # 73% from bottom
        gp = gpar(
          col = "black",         # text color
          fontsize = 10,
          fill = "white"         # background fill
        )
      )
    ) +
    theme_bw() +
    theme(legend.position = c(0.81, 0.73),
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
mc.cores = 4)
save(plots_pnew, 
     file = "bayesian_subset_selection/actg/results/figures/ppc/plots_pnew.RData")


plots_pnew_after_PSM <- mclapply(1:nrow(pnew.wip_after_PSM), function(i){
  plot <- ggplot() +
    geom_density(data = data.frame(y = pnew.wip_after_PSM[i, ]), 
                 aes(x = y, color = "wip"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = pnew.norm_after_PSM[[1]][i, ]),
                 aes(x = y, color = "c01"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = pnew.norm_after_PSM[[2]][i, ]),
                 aes(x = y, color = "c02"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = pnew.norm_after_PSM[[3]][i, ]),
                 aes(x = y, color = "c03"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = pnew.norm_after_PSM[[4]][i, ]), 
                 aes(x = y, color = "c04"),
                 linewidth = 1) + 
    scale_color_manual(
      name = NULL,
      values = c(
        "wip" = "black",
        "c01"    = "#66A8D0",
        "c02"    = "#D06673",
        "c03" = "#7f7f7f",
        "c04" = "#D0C366"
      ),
      labels = c(
        "wip" = expression(Cauchy),
        "c01" = bquote(Normal(0, .(c0[1])^2 * I[p^(m)])),
        "c02" = bquote(Normal(0, .(c0[2])^2 * I[p^(m)])),
        "c03" = bquote(Normal(0, .(c0[3])^2 * I[p^(m)])),
        "c04" = bquote(Normal(0, .(c0[4])^2 * I[p^(m)]))
      ),
      breaks = c("wip",
                 "c01",
                 "c02",
                 "c03",
                 "c04"
      )
    ) +
  annotation_custom(
    grid::textGrob(
      label = paste0("Obs outcome: ", current_data$outcome[i]),
      x = unit(0.81, "npc"),   # 81% from left
      y = unit(0.96, "npc"),   # 73% from bottom
      gp = gpar(
        col = "black",         # text color
        fontsize = 10,
        fill = "white"         # background fill
      )
    )
  ) +
    theme_bw() +
    theme(legend.position = c(0.81, 0.73),
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
  filename <- paste0("bayesian_subset_selection/actg/results/figures/ppc/ppc_pnew_", 
                     i, "_after_PSM.png")
  ggsave(filename, plot = plot, width = 6, height = 4, units = "in", dpi = 300)
  return(plot)
},
mc.cores = 8)
save(plots_pnew_after_PSM,
     file = "bayesian_subset_selection/actg/results/figures/ppc/plots_pnew_after_PSM.RData")


plots_pnew_a0 <- mclapply(1:nrow(pnew.wip), function(i){
  plot <- ggplot() +
    geom_density(data = data.frame(y = pnew.a0[[1]][i, ]),
                 aes(x = y, color = "a01"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = pnew.a0[[2]][i, ]),
                 aes(x = y, color = "a02"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = pnew.a0[[3]][i, ]),
                 aes(x = y, color = "a03"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = pnew.a0[[4]][i, ]), 
                 aes(x = y, color = "a04"),
                 linewidth = 1) + 
    scale_color_manual(
      name = NULL,
      values = c(
        "a01" = brewer.pal(4, "Dark2")[1],
        "a02" = brewer.pal(4, "Dark2")[2],
        "a03" = brewer.pal(4, "Dark2")[3],
        "a04" = brewer.pal(4, "Dark2")[4]
      ),
      labels = c(
        "a01" = bquote("Beta"*"("*.(a0_hyper[[1]][1])*", "*.(a0_hyper[[1]][2])*")"),
        "a02" = bquote("Beta"*"("*.(a0_hyper[[2]][1])*", "*.(a0_hyper[[2]][2])*")"),
        "a03" = bquote("Beta"*"("*.(a0_hyper[[3]][1])*", "*.(a0_hyper[[3]][2])*")"),
        "a04" = bquote("Beta"*"("*.(a0_hyper[[4]][1])*", "*.(a0_hyper[[4]][2])*")")
      ),
      breaks = c("a01",
                 "a02",
                 "a03",
                 "a04"
      )
    ) +
    annotation_custom(
      grid::textGrob(
        label = paste0("Obs outcome: ", current_data$outcome[i]),
        x = unit(0.81, "npc"),   # 81% from left
        y = unit(0.96, "npc"),   # 73% from bottom
        gp = gpar(
          col = "black",         # text color
          fontsize = 10,
          fill = "white"         # background fill
        )
      )
    ) +
    theme_bw() +
    theme(legend.position = c(0.81, 0.73),
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
  filename <- paste0("bayesian_subset_selection/actg/results/figures/ppc_a0/ppc_pnew_", i, ".png")
  ggsave(filename, plot = plot, width = 6, height = 4, units = "in", dpi = 300)
  return(plot)
},
mc.cores = 1)
save(plots_pnew_a0, 
     file = "bayesian_subset_selection/actg/results/figures/ppc_a0/plots_pnew_a0.RData")


plots_pnew_a0_after_PSM <- mclapply(1:nrow(pnew.wip_after_PSM), function(i){
  plot <- ggplot() +
    geom_density(data = data.frame(y = pnew.a0_after_PSM[[1]][i, ]),
                 aes(x = y, color = "a01"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = pnew.a0_after_PSM[[2]][i, ]),
                 aes(x = y, color = "a02"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = pnew.a0_after_PSM[[3]][i, ]),
                 aes(x = y, color = "a03"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = pnew.a0_after_PSM[[4]][i, ]), 
                 aes(x = y, color = "a04"),
                 linewidth = 1) + 
    scale_color_manual(
      name = NULL,
      values = c(
        "a01" = brewer.pal(4, "Dark2")[1],
        "a02" = brewer.pal(4, "Dark2")[2],
        "a03" = brewer.pal(4, "Dark2")[3],
        "a04" = brewer.pal(4, "Dark2")[4]
      ),
      labels = c(
        "a01" = bquote("Beta"*"("*.(a0_hyper[[1]][1])*", "*.(a0_hyper[[1]][2])*")"),
        "a02" = bquote("Beta"*"("*.(a0_hyper[[2]][1])*", "*.(a0_hyper[[2]][2])*")"),
        "a03" = bquote("Beta"*"("*.(a0_hyper[[3]][1])*", "*.(a0_hyper[[3]][2])*")"),
        "a04" = bquote("Beta"*"("*.(a0_hyper[[4]][1])*", "*.(a0_hyper[[4]][2])*")")
      ),
      breaks = c("a01",
                 "a02",
                 "a03",
                 "a04"
      )
    ) +
    annotation_custom(
      grid::textGrob(
        label = paste0("Obs outcome: ", current_data$outcome[i]),
        x = unit(0.81, "npc"),   # 81% from left
        y = unit(0.96, "npc"),   # 73% from bottom
        gp = gpar(
          col = "black",         # text color
          fontsize = 10,
          fill = "white"         # background fill
        )
      )
    ) +
    theme_bw() +
    theme(legend.position = c(0.81, 0.73),
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
  filename <- paste0("bayesian_subset_selection/actg/results/figures/ppc_a0/ppc_pnew_", 
                     i, "_after_PSM.png")
  ggsave(filename, plot = plot, width = 6, height = 4, units = "in", dpi = 300)
  return(plot)
},
mc.cores = 1)

save(plots_pnew_a0_after_PSM,
     file = "bayesian_subset_selection/actg/results/figures/ppc_a0/plots_pnew_a0_after_PSM.RData")


#-------------------------------------------------------------------------------
# ROC curves
#-------------------------------------------------------------------------------

library(plotROC)
library(pROC)


plot_roc.wip <- ggplot(
  data.frame(
    out = current_data$outcome,
    prob = rowMeans(pnew.wip)
), aes(d = out, m = prob)) + 
  geom_ribbon(
    data = data.frame(
      fpr = 1 - as.numeric(rownames(ci.roc.wip)),
      tpr_lower = ci.roc.wip[,1],
      tpr_upper = ci.roc.wip[,3]
    ),
    aes(x = fpr, ymin = tpr_lower, ymax = tpr_upper),
    fill = rgb(0.2, 0.4, 0.6, 0.2),
    inherit.aes = FALSE
  ) +
  geom_roc(labels = FALSE, color = "#66A8D0") + 
  annotate("label",
           x = 0.9, y = 0,
           label = sprintf("AUC = %.3f (%.3f–%.3f)",
                           auc(roc.curve.wip),
                           ci.auc(roc.curve.wip)[1],
                           ci.auc(roc.curve.wip)[3]),
           size = 3, fill = "white", color = "black") +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  theme_bw() +
  theme(legend.position = c(0.81, 0.73),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
  ) +
  theme(text = element_text(size = 10),        # Base text size
        axis.title = element_text(size = 10),  # Axis titles
        axis.text = element_text(size = 10),   # Axis tick labels
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 10)) +
  style_roc()

ggsave("bayesian_subset_selection/actg/results/figures/ppc/ppc_roc_wip.png",
       plot = plot_roc.wip, width = 7.5, height = 4, units = "in", dpi = 300)


plot_roc.wip_after_PSM <- ggplot(
  data.frame(
    out = current_data$outcome,
    prob = rowMeans(pnew.wip_after_PSM)
  ), aes(d = out, m = prob)) + 
  geom_ribbon(
    data = data.frame(
      fpr = 1 - as.numeric(rownames(ci.roc.wip_after_PSM)),
      tpr_lower = ci.roc.wip_after_PSM[,1],
      tpr_upper = ci.roc.wip_after_PSM[,3]
    ),
    aes(x = fpr, ymin = tpr_lower, ymax = tpr_upper),
    fill = rgb(0.2, 0.4, 0.6, 0.2),
    inherit.aes = FALSE
  ) +
  geom_roc(labels = FALSE, color = "#66A8D0") + 
  annotate("label",
           x = 0.9, y = 0,
           label = sprintf("AUC = %.3f (%.3f–%.3f)",
                           auc(roc.curve.wip_after_PSM),
                           ci.auc(roc.curve.wip_after_PSM)[1],
                           ci.auc(roc.curve.wip_after_PSM)[3]),
           size = 3, fill = "white", color = "black") +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  theme_bw() +
  theme(legend.position = c(0.81, 0.73),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
  ) +
  theme(text = element_text(size = 10),        # Base text size
        axis.title = element_text(size = 10),  # Axis titles
        axis.text = element_text(size = 10),   # Axis tick labels
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 10)) +
  style_roc()
ggsave("bayesian_subset_selection/actg/results/figures/ppc/ppc_roc_wip_after_PSM.png",
       plot = plot_roc.wip_after_PSM, width = 7.5, height = 4, units = "in", dpi = 300)

plot_roc.norm <- lapply(seq_along(roc.curve.norm), function(j) {
  ggplot(
    data.frame(
      out = current_data$outcome,
      prob = rowMeans(pnew.norm[[j]])
    ), aes(d = out, m = prob)) + 
    geom_ribbon(
      data = data.frame(
        fpr = 1 - as.numeric(rownames(ci.roc.norm[[j]])),
        tpr_lower = ci.roc.norm[[j]][,1],
        tpr_upper = ci.roc.norm[[j]][,3]
      ),
      aes(x = fpr, ymin = tpr_lower, ymax = tpr_upper),
      fill = rgb(0.2, 0.4, 0.6, 0.2),
      inherit.aes = FALSE
    ) +
    geom_roc(labels = FALSE, color = "#66A8D0") + 
    annotate("label",
             x = 0.9, y = 0,
             label = sprintf("AUC = %.3f (%.3f–%.3f)",
                             auc(roc.curve.norm[[j]]),
                             ci.auc(roc.curve.norm[[j]])[1],
                             ci.auc(roc.curve.norm[[j]])[3]),
             size = 3, fill = "white", color = "black") +
    theme_bw() +
    theme(legend.position = c(0.81, 0.73),
          legend.background = element_rect(fill = "white", color = "black")) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
    ) +
    theme(text = element_text(size = 10),        # Base text size
          axis.title = element_text(size = 10),  # Axis titles
          axis.text = element_text(size = 10),   # Axis tick labels
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 10)) +
    style_roc()
})

plot_roc.norm_after_PSM <- lapply(seq_along(roc.curve.norm_after_PSM), function(j) {
  ggplot(
    data.frame(
      out = current_data$outcome,
      prob = rowMeans(pnew.norm_after_PSM[[j]])
    ), aes(d = out, m = prob)) + 
    geom_ribbon(
      data = data.frame(
        fpr = 1 - as.numeric(rownames(ci.roc.norm_after_PSM[[j]])),
        tpr_lower = ci.roc.norm_after_PSM[[j]][,1],
        tpr_upper = ci.roc.norm_after_PSM[[j]][,3]
      ),
      aes(x = fpr, ymin = tpr_lower, ymax = tpr_upper),
      fill = rgb(0.2, 0.4, 0.6, 0.2),
      inherit.aes = FALSE
    ) +
    geom_roc(labels = FALSE, color = "#66A8D0") + 
    annotate("label",
             x = 0.9, y = 0,
             label = sprintf("AUC = %.3f (%.3f–%.3f)",
                             auc(roc.curve.norm_after_PSM[[j]]),
                             ci.auc(roc.curve.norm_after_PSM[[j]])[1],
                             ci.auc(roc.curve.norm_after_PSM[[j]])[3]),
             size = 3, fill = "white", color = "black") +
    theme_bw() +
    theme(legend.position = c(0.81, 0.73),
          legend.background = element_rect(fill = "white", color = "black")) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
    ) +
    theme(text = element_text(size = 10),        # Base text size
          axis.title = element_text(size = 10),  # Axis titles
          axis.text = element_text(size = 10),   # Axis tick labels
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 10)) +
    style_roc()
})

plot_roc.a0 <- lapply(seq_along(roc.curve.a0), function(j) {
  ggplot(
    data.frame(
      out = current_data$outcome,
      prob = rowMeans(pnew.a0[[j]])
    ), aes(d = out, m = prob)) + 
    geom_ribbon(
      data = data.frame(
        fpr = 1 - as.numeric(rownames(ci.roc.a0[[j]])),
        tpr_lower = ci.roc.a0[[j]][,1],
        tpr_upper = ci.roc.a0[[j]][,3]
      ),
      aes(x = fpr, ymin = tpr_lower, ymax = tpr_upper),
      fill = rgb(0.2, 0.4, 0.6, 0.2),
      inherit.aes = FALSE
    ) +
    geom_roc(labels = FALSE, color = "#66A8D0") + 
    annotate("label",
             x = 0.9, y = 0,
             label = sprintf("AUC = %.3f (%.3f–%.3f)",
                             auc(roc.curve.a0[[j]]),
                             ci.auc(roc.curve.a0[[j]])[1],
                             ci.auc(roc.curve.a0[[j]])[3]),
             size = 3, fill = "white", color = "black") +
    theme_bw() +
    theme(legend.position = c(0.81, 0.73),
          legend.background = element_rect(fill = "white", color = "black")) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
    ) +
    theme(text = element_text(size = 10),        # Base text size
          axis.title = element_text(size = 10),  # Axis titles
          axis.text = element_text(size = 10),   # Axis tick labels
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 10)) +
    style_roc()
})

plot_roc.a0_after_PSM <- lapply(seq_along(roc.curve.a0_after_PSM), function(j) {
  ggplot(
    data.frame(
      out = current_data$outcome,
      prob = rowMeans(pnew.a0_after_PSM[[j]])
    ), aes(d = out, m = prob)) + 
    geom_ribbon(
      data = data.frame(
        fpr = 1 - as.numeric(rownames(ci.roc.a0_after_PSM[[j]])),
        tpr_lower = ci.roc.a0_after_PSM[[j]][,1],
        tpr_upper = ci.roc.a0_after_PSM[[j]][,3]
      ),
      aes(x = fpr, ymin = tpr_lower, ymax = tpr_upper),
      fill = rgb(0.2, 0.4, 0.6, 0.2),
      inherit.aes = FALSE
    ) +
    geom_roc(labels = FALSE, color = "#66A8D0") + 
    annotate("label",
             x = 0.9, y = 0,
             label = sprintf("AUC = %.3f (%.3f–%.3f)",
                             pROC::auc(roc.curve.a0_after_PSM[[j]]),
                             ci.auc(roc.curve.a0_after_PSM[[j]])[1],
                             ci.auc(roc.curve.a0_after_PSM[[j]])[3]),
             size = 3, fill = "white", color = "black") +
    theme_bw() +
    theme(legend.position = c(0.81, 0.73),
          legend.background = element_rect(fill = "white", color = "black")) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
    ) +
    theme(text = element_text(size = 10),        # Base text size
          axis.title = element_text(size = 10),  # Axis titles
          axis.text = element_text(size = 10),   # Axis tick labels
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 10)) +
    style_roc()
})


plot_roc.curves.wip <- ggplot(
  melt_roc(
    cbind(
      data.frame(
        out = current_data$outcome
      ),
      as.data.frame(pnew.wip)
    ),
    d = "out",
    m = paste0("V", 1:(ncol(pnew.wip)))
  ), 
  aes(d = D, m = M, fill = name)) + 
  geom_roc(labels = F, color = "#66A8D0") + 
  theme_bw() +
  theme(legend.position = c(0.81, 0.73),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
  ) +
  theme(text = element_text(size = 10),        # Base text size
        axis.title = element_text(size = 10),  # Axis titles
        axis.text = element_text(size = 10),   # Axis tick labels
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 10)) +
  style_roc()

plot_roc.curves.wip_after_PSM <- ggplot(
  melt_roc(
    cbind(
      data.frame(
        out = current_data$outcome
      ),
      as.data.frame(pnew.wip_after_PSM)
    ),
    d = "out",
    m = paste0("V", 1:(ncol(pnew.wip_after_PSM)))
  ), 
  aes(d = D, m = M, fill = name)) + 
  geom_roc(labels = F, color = "#66A8D0") + 
  theme_bw() +
  theme(legend.position = c(0.81, 0.73),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
  ) +
  theme(text = element_text(size = 10),        # Base text size
        axis.title = element_text(size = 10),  # Axis titles
        axis.text = element_text(size = 10),   # Axis tick labels
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 10)) +
  style_roc()

plot_roc.curves.norm <- lapply(seq_along(roc.curves.norm), function(j) {
  ggplot(
    melt_roc(
      cbind(
        data.frame(
          out = current_data$outcome
        ),
        as.data.frame(pnew.norm[[j]])
      ),
      d = "out",
      m = paste0("V", 1:(ncol(pnew.norm[[j]])))
    ), 
    aes(d = D, m = M, fill = name)) + 
    geom_roc(labels = F, color = "#66A8D0") + 
    theme_bw() +
    theme(legend.position = c(0.81, 0.73),
          legend.background = element_rect(fill = "white", color = "black")) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
    ) +
    theme(text = element_text(size = 10),        # Base text size
          axis.title = element_text(size = 10),  # Axis titles
          axis.text = element_text(size = 10),   # Axis tick labels
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 10)) +
    style_roc()
})

plot_roc.curves.norm_after_PSM <- lapply(seq_along(roc.curves.norm_after_PSM), function(j) {
  ggplot(
    melt_roc(
      cbind(
        data.frame(
          out = current_data$outcome
        ),
        as.data.frame(pnew.norm_after_PSM[[j]])
      ),
      d = "out",
      m = paste0("V", 1:(ncol(pnew.norm_after_PSM[[j]])))
    ), 
    aes(d = D, m = M, fill = name)) + 
    geom_roc(labels = F, color = "#66A8D0") + 
    theme_bw() +
    theme(legend.position = c(0.81, 0.73),
          legend.background = element_rect(fill = "white", color = "black")) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
    ) +
    theme(text = element_text(size = 10),        # Base text size
          axis.title = element_text(size = 10),  # Axis titles
          axis.text = element_text(size = 10),   # Axis tick labels
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 10)) +
    style_roc()
})

plot_roc.curves.a0 <- lapply(seq_along(roc.curves.a0), function(j) {
  ggplot(
    melt_roc(
      cbind(
        data.frame(
          out = current_data$outcome
        ),
        as.data.frame(pnew.a0[[j]])
      ),
      d = "out",
      m = paste0("V", 1:(ncol(pnew.a0[[j]])))
    ), 
    aes(d = D, m = M, fill = name)) + 
    geom_roc(labels = F, color = "#66A8D0") + 
    theme_bw() +
    theme(legend.position = c(0.81, 0.73),
          legend.background = element_rect(fill = "white", color = "black")) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
    ) +
    theme(text = element_text(size = 10),        # Base text size
          axis.title = element_text(size = 10),  # Axis titles
          axis.text = element_text(size = 10),   # Axis tick labels
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 10)) +
    style_roc()
})

plot_roc.curves.a0_after_PSM <- lapply(seq_along(roc.curves.a0_after_PSM), function(j) {
  ggplot(
    melt_roc(
      cbind(
        data.frame(
          out = current_data$outcome
        ),
        as.data.frame(pnew.a0_after_PSM[[j]])
      ),
      d = "out",
      m = paste0("V", 1:(ncol(pnew.a0_after_PSM[[j]])))
    ), 
    aes(d = D, m = M, fill = name)) + 
    geom_roc(labels = F, color = "#66A8D0") + 
    theme_bw() +
    theme(legend.position = c(0.81, 0.73),
          legend.background = element_rect(fill = "white", color = "black")) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
    ) +
    theme(text = element_text(size = 10),        # Base text size
          axis.title = element_text(size = 10),  # Axis titles
          axis.text = element_text(size = 10),   # Axis tick labels
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 10)) +
    style_roc()
})


## Combine plots


roc_wip_norm <- 
  (plot_roc.wip +
     ggtitle("Cauchy") +
     style_roc(xlab = "")) +
  (plot_roc.norm[[1]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[1])^2 * I[p^"(m)"]*")")) +
     style_roc(xlab = "")) +
  (plot_roc.norm[[2]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[2])^2 * I[p^"(m)"]*")")) + 
     style_roc(xlab = "")) +
  (plot_roc.norm[[3]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[3])^2 * I[p^"(m)"]*")")) + 
     style_roc(xlab = "")
   ) +
  (plot_roc.norm[[4]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[4])^2 * I[p^"(m)"]*")")) 
  ) +
  plot_layout(ncol = 2)
roc_wip_norm
ggsave("bayesian_subset_selection/actg/results/figures/ppc/ppc_roc_wip_norm.png",
       roc_wip_norm, width = 14, height = 10.8, units = "in", dpi = 300)

roc_a0 <- 
  (plot_roc.a0[[1]] + 
     ggtitle(bquote("Beta"*"("*.(a0_hyper[[1]][1])*", "*.(a0_hyper[[1]][2])*")")) +
     style_roc(xlab = "")) +
  (plot_roc.a0[[2]] + 
     ggtitle(bquote("Beta"*"("*.(a0_hyper[[2]][1])*", "*.(a0_hyper[[2]][2])*")")) + 
     style_roc(xlab = "")) +
  (plot_roc.a0[[3]] + 
     ggtitle(bquote("Beta"*"("*.(a0_hyper[[3]][1])*", "*.(a0_hyper[[3]][2])*")")) + 
     style_roc(xlab = "")) +
  (plot_roc.a0[[4]] + 
     ggtitle(bquote("Beta"*"("*.(a0_hyper[[4]][1])*", "*.(a0_hyper[[4]][2])*")")) 
  ) +
  plot_layout(ncol = 2)
roc_a0
ggsave("bayesian_subset_selection/actg/results/figures/ppc_a0/ppc_roc_a0.png",
       roc_a0, width = 14, height = 7.2, units = "in", dpi = 300)

roc_all_wip_norm <- 
  (plot_roc.curves.wip +
     ggtitle("Cauchy") +
     style_roc(xlab = "")) +
  (plot_roc.curves.norm[[1]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[1])^2 * I[p^"(m)"]*")")) +
     style_roc(xlab = "")) +
  (plot_roc.curves.norm[[2]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[2])^2 * I[p^"(m)"]*")")) + 
     style_roc(xlab = "")) +
  (plot_roc.curves.norm[[3]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[3])^2 * I[p^"(m)"]*")")) + 
     style_roc()) +
  (plot_roc.curves.norm[[4]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[4])^2 * I[p^"(m)"]*")")) 
  ) + plot_layout(ncol = 2)
roc_all_wip_norm
ggsave("bayesian_subset_selection/actg/results/figures/ppc/ppc_roc_all_wip_norm.png",
       roc_all_wip_norm, width = 10, height = 10.8, units = "in", dpi = 300)

roc_all_a0 <- 
  (plot_roc.curves.a0[[1]] + 
     ggtitle(bquote("Beta"*"("*.(a0_hyper[[1]][1])*", "*.(a0_hyper[[1]][2])*")")) +
     style_roc(xlab = "")) /
  (plot_roc.curves.a0[[2]] + 
     ggtitle(bquote("Beta"*"("*.(a0_hyper[[2]][1])*", "*.(a0_hyper[[2]][2])*")")) + 
     style_roc(xlab = "")) /
  (plot_roc.curves.a0[[3]] + 
     ggtitle(bquote("Beta"*"("*.(a0_hyper[[3]][1])*", "*.(a0_hyper[[3]][2])*")")) + 
     style_roc()) /
  (plot_roc.curves.a0[[4]] + 
     ggtitle(bquote("Beta"*"("*.(a0_hyper[[4]][1])*", "*.(a0_hyper[[4]][2])*")")) 
  ) +
  plot_layout(ncol = 2)
roc_all_a0
ggsave("bayesian_subset_selection/actg/results/figures/ppc_a0/ppc_roc_all_a0.png",
       roc_all_a0, width = 10, height = 7.2, units = "in", dpi = 300)

roc_wip_norm_after_PSM <- 
  (plot_roc.wip_after_PSM +
     ggtitle("Cauchy") +
     style_roc(xlab = "")
  ) +
  (plot_roc.norm_after_PSM[[1]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[1])^2 * I[p^"(m)"]*")")) +
     style_roc(xlab = "")) +
  (plot_roc.norm_after_PSM[[2]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[2])^2 * I[p^"(m)"]*")")) +
     style_roc(xlab = "")) +
  (plot_roc.norm_after_PSM[[3]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[3])^2 * I[p^"(m)"]*")")) +
     style_roc(xlab = "")) +
  (plot_roc.norm_after_PSM[[4]] +
   ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[4])^2 * I[p^"(m)"]*")"))
  ) +
  plot_layout(ncol = 2)
roc_wip_norm_after_PSM
ggsave("bayesian_subset_selection/actg/results/figures/ppc/ppc_roc_wip_norm_after_PSM.png",
       roc_wip_norm_after_PSM, width = 14, height = 10.8, units = "in", dpi = 300)

roc_a0_after_PSM <- 
  (plot_roc.a0_after_PSM[[1]] + 
     ggtitle(bquote("Beta"*"("*.(a0_hyper[[1]][1])*", "*.(a0_hyper[[1]][2])*")")) +
     style_roc(xlab = "")) +
  (plot_roc.a0_after_PSM[[2]] + 
     ggtitle(bquote("Beta"*"("*.(a0_hyper[[2]][1])*", "*.(a0_hyper[[2]][2])*")")) +
     style_roc(xlab = "")) +
  (plot_roc.a0_after_PSM[[3]] + 
     ggtitle(bquote("Beta"*"("*.(a0_hyper[[3]][1])*", "*.(a0_hyper[[3]][2])*")")) +
     style_roc(xlab = "")) +
  (plot_roc.a0_after_PSM[[4]] +
   ggtitle(bquote("Beta"*"("*.(a0_hyper[[4]][1])*", "*.(a0_hyper[[4]][2])*")"))
  ) +
  plot_layout(ncol = 2)
roc_a0_after_PSM
ggsave("bayesian_subset_selection/actg/results/figures/ppc_a0/ppc_roc_a0_after_PSM.png",
       roc_a0_after_PSM, width = 14, height = 7.2, units = "in", dpi = 300)

roc_all_wip_norm_after_PSM <- 
  (plot_roc.curves.wip_after_PSM +
     ggtitle("Cauchy") +
     style_roc(xlab = "")
  ) +
  (plot_roc.curves.norm_after_PSM[[1]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[1])^2 * I[p^"(m)"]*")")) +
     style_roc(xlab = "")) +
  (plot_roc.curves.norm_after_PSM[[2]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[2])^2 * I[p^"(m)"]*")")) +
     style_roc(xlab = "")) +
  (plot_roc.curves.norm_after_PSM[[3]] + 
     ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[3])^2 * I[p^"(m)"]*")")) +
     style_roc()) +
  (plot_roc.curves.norm_after_PSM[[4]] +
   ggtitle(bquote("Normal" * "(" * 0 * "," * .(c0[4])^2 * I[p^"(m)"]*")"))
  ) +
  plot_layout(ncol = 2)
roc_all_wip_norm_after_PSM
ggsave("bayesian_subset_selection/actg/results/figures/ppc/ppc_roc_all_wip_norm_after_PSM.png",
       roc_all_wip_norm_after_PSM, width = 10, height = 10.8, units = "in", dpi = 300)

roc_all_a0_after_PSM <- 
  (plot_roc.curves.a0_after_PSM[[1]] + 
     ggtitle(bquote("Beta"*"("*.(a0_hyper[[1]][1])*", "*.(a0_hyper[[1]][2])*")")) +
     style_roc(xlab = "")) +
  (plot_roc.curves.a0_after_PSM[[2]] + 
     ggtitle(bquote("Beta"*"("*.(a0_hyper[[2]][1])*", "*.(a0_hyper[[2]][2])*")")) +
     style_roc(xlab = "")) +
  (plot_roc.curves.a0_after_PSM[[3]] + 
     ggtitle(bquote("Beta"*"("*.(a0_hyper[[3]][1])*", "*.(a0_hyper[[3]][2])*")")) +
     style_roc()) +
  (plot_roc.curves.a0_after_PSM[[4]] +
   ggtitle(bquote("Beta"*"("*.(a0_hyper[[4]][1])*", "*.(a0_hyper[[4]][2])*")"))
  )
roc_all_a0_after_PSM
ggsave("bayesian_subset_selection/actg/results/figures/ppc_a0/ppc_roc_all_a0_after_PSM.png",
       roc_all_a0_after_PSM, width = 10, height = 7.2, units = "in", dpi = 300)


#-------------------------------------------------------------------------------
# Compare distributions before and after PSM
#-------------------------------------------------------------------------------


load("bayesian_subset_selection/actg/results/figures/ppc/plots_pnew.RData")
load("bayesian_subset_selection/actg/results/figures/ppc/plots_pnew_after_PSM.RData")
load("bayesian_subset_selection/actg/results/figures/ppc_a0/plots_pnew_a0.RData")
load("bayesian_subset_selection/actg/results/figures/ppc_a0/plots_pnew_a0_after_PSM.RData")

pnew_1_before_PSM <- plots_pnew[[1]] +
  xlim(c(0,1)) +
  ggtitle("Before PSM") +
  theme(legend.position = "none")
pnew_1_before_PSM$layers <- pnew_1_before_PSM$layers[-length(pnew_1_before_PSM$layers)]
pnew_1_before_after_PSM <- 
  (pnew_1_before_PSM) | 
  (plots_pnew_after_PSM[[1]])
pnew_1_before_after_PSM
ggsave("bayesian_subset_selection/actg/results/figures/ppc/ppc_pnew_1_before_after_PSM.png",
       pnew_1_before_after_PSM, width = 12, height = 4.5, units = "in", dpi = 300)

pnew_13_before_PSM <- plots_pnew[[13]] +
  xlim(c(0,1)) +
  ggtitle("Before PSM") +
  theme(legend.position = "none")
pnew_13_before_PSM$layers <- 
  pnew_13_before_PSM$layers[-length(pnew_13_before_PSM$layers)]
pnew_13_before_after_PSM <- (pnew_13_before_PSM) | 
  (plots_pnew_after_PSM[[13]] + 
     xlim(c(0,1)) +
     ggtitle("After PSM")
   )
pnew_13_before_after_PSM
ggsave("bayesian_subset_selection/actg/results/figures/ppc/ppc_pnew_13_before_after_PSM.png",
       pnew_13_before_after_PSM, width = 12, height = 4.5, units = "in", dpi = 300)

pnew_a0_1_before_PSM <- plots_pnew_a0[[1]] +
  xlim(c(0,1)) +
  ggtitle("Before PSM") +
  theme(legend.position = "none")
pnew_a0_1_before_PSM$layers <- 
  pnew_a0_1_before_PSM$layers[-length(pnew_a0_1_before_PSM$layers)]
pnew_a0_1_before_after_PSM <- (pnew_a0_1_before_PSM) | 
  (plots_pnew_a0_after_PSM[[1]] + 
     xlim(c(0,1)) +
     ggtitle("After PSM")
   )
pnew_a0_1_before_after_PSM
ggsave("bayesian_subset_selection/actg/results/figures/ppc_a0/ppc_pnew_1_before_after_PSM.png",
       pnew_a0_1_before_after_PSM, width = 12, height = 4.5, units = "in", dpi = 300)

pnew_a0_13_before_PSM <- plots_pnew_a0[[13]] +
  xlim(c(0,1)) +
  ggtitle("Before PSM") +
  theme(legend.position = "none")
pnew_a0_13_before_PSM$layers <- 
  pnew_a0_13_before_PSM$layers[-length(pnew_a0_13_before_PSM$layers)]
pnew_a0_13_before_after_PSM <- (pnew_a0_13_before_PSM) | 
  (plots_pnew_a0_after_PSM[[13]] + 
     xlim(c(0,1)) +
     ggtitle("After PSM")
   )
pnew_a0_13_before_after_PSM
ggsave("bayesian_subset_selection/actg/results/figures/ppc_a0/ppc_pnew_13_before_after_PSM.png",
       pnew_a0_13_before_after_PSM, width = 12, height = 4.5, units = "in", dpi = 300)

roc_wip_norm_before_after_PSM <- 
  (roc_wip_norm + plot_layout(ncol = 1)) | 
  (roc_wip_norm_after_PSM + plot_layout(ncol = 1))
roc_wip_norm_before_after_PSM
ggsave("bayesian_subset_selection/actg/results/figures/ppc/ppc_roc_wip_norm_before_after_PSM.png",
       roc_wip_norm_before_after_PSM, width = 12, height = 20, units = "in", dpi = 300)

roc_all_wip_norm_before_after_PSM <-
  (roc_all_wip_norm + plot_layout(ncol = 1)) | 
  (roc_all_wip_norm_after_PSM + plot_layout(ncol = 1))
roc_all_wip_norm_before_after_PSM
ggsave("bayesian_subset_selection/actg/results/figures/ppc/ppc_roc_all_wip_norm_before_after_PSM.png",
       roc_all_wip_norm_before_after_PSM, width = 12, height = 20, units = "in", dpi = 300)

roc_a0_before_after_PSM <-
  (roc_a0 + plot_layout(ncol = 1)) | 
  (roc_a0_after_PSM + plot_layout(ncol = 1))
roc_a0_before_after_PSM
ggsave("bayesian_subset_selection/actg/results/figures/ppc_a0/ppc_roc_a0_before_after_PSM.png",
       roc_a0_before_after_PSM, width = 12, height = 16, units = "in", dpi = 300)

roc_all_a0_before_after_PSM <-
  (roc_all_a0 + plot_layout(ncol = 1)) | 
  (roc_all_a0_after_PSM + plot_layout(ncol = 1))
roc_all_a0_before_after_PSM
ggsave("bayesian_subset_selection/actg/results/figures/ppc_a0/ppc_roc_all_a0_before_after_PSM.png",
       roc_all_a0_before_after_PSM, width = 12, height = 16, units = "in", dpi = 300)
