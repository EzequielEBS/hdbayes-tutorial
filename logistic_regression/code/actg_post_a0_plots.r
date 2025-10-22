#-------------------------------------------------------------------------------
# This script generates posterior plots for a0 before and after PSM
#-------------------------------------------------------------------------------

# load posterior samples
load("logistic_regression/samples/post_samples_a0_wip.RData")
load("logistic_regression/samples/post_samples_a0_wip_after_PSM.RData")
load("logistic_regression/samples/post_samples_a0_norm.RData")
load("logistic_regression/samples/post_samples_a0_norm_after_PSM.RData")

# load libraries
library(ggplot2)
library(RColorBrewer)
library(grid)
library(wesanderson)
library(patchwork)

plot_wip <- ggplot() +
  geom_density(data = data.frame(y = draws.a0.post.wip[,1]),
               aes(x = y, color = "a01"),
               linewidth = 1) +
  geom_density(data = data.frame(y = draws.a0.post.wip[,2]),
               aes(x = y, color = "a02"),
               linewidth = 1) +
  geom_density(data = data.frame(y = draws.a0.post.wip[,3]),
               aes(x = y, color = "a03"),
               linewidth = 1) +
  geom_density(data = data.frame(y = draws.a0.post.wip[,4]), 
               aes(x = y, color = "a04"),
               linewidth = 1) + 
  scale_color_manual(
    name = NULL,
    values = c(
      "a01" = "black",
      "a02" = brewer.pal(9, "Set1")[7],
      "a03" = brewer.pal(9, "Set1")[8],
      "a04" = brewer.pal(9, "Set1")[9]
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
  theme_bw() +
  theme(legend.position = c(0.81, 0.73),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
  ) +
  xlim(c(0,1)) +
  theme(text = element_text(size = 16),        # Base text size
        axis.title = element_text(size = 18),  # Axis titles
        axis.text = element_text(size = 16),   # Axis tick labels
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  labs(
    x = expression(a[0]),
    y = "",
  )

plot_wip_after_PSM <- ggplot() +
  geom_density(data = data.frame(y = draws.a0.post.wip_after_PSM[,1]),
               aes(x = y, color = "a01"),
               linewidth = 1) +
  geom_density(data = data.frame(y = draws.a0.post.wip_after_PSM[,2]),
               aes(x = y, color = "a02"),
               linewidth = 1) +
  geom_density(data = data.frame(y = draws.a0.post.wip_after_PSM[,3]),
               aes(x = y, color = "a03"),
               linewidth = 1) +
  geom_density(data = data.frame(y = draws.a0.post.wip_after_PSM[,4]), 
               aes(x = y, color = "a04"),
               linewidth = 1) + 
  scale_color_manual(
    name = NULL,
    values = c(
      "a01" = "black",
      "a02" = brewer.pal(9, "Set1")[7],
      "a03" = brewer.pal(9, "Set1")[8],
      "a04" = brewer.pal(9, "Set1")[9]
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
  theme_bw() +
  theme(legend.position = c(0.81, 0.73),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
  ) +
  xlim(c(0,1)) +
  theme(text = element_text(size = 16),        # Base text size
        axis.title = element_text(size = 18),  # Axis titles
        axis.text = element_text(size = 16),   # Axis tick labels
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  labs(
    x = expression(a[0]),
    y = "",
  )


plots_norm <- lapply(seq_along(c0), function(i) {
  ggplot() +
    geom_density(data = data.frame(y = draws.a0.post.norm[[i]][,1]),
                 aes(x = y, color = "a01"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = draws.a0.post.norm[[i]][,2]),
                 aes(x = y, color = "a02"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = draws.a0.post.norm[[i]][,3]),
                 aes(x = y, color = "a03"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = draws.a0.post.norm[[i]][,4]), 
                 aes(x = y, color = "a04"),
                 linewidth = 1) + 
    scale_color_manual(
      name = NULL,
      values = c(
        "a01" = "black",
        "a02" = brewer.pal(9, "Set1")[7],
        "a03" = brewer.pal(9, "Set1")[8],
        "a04" = brewer.pal(9, "Set1")[9]
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
    theme_bw() +
    theme(legend.position = 'none') +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
    ) +
    xlim(c(0,1)) +
    theme(text = element_text(size = 16),        # Base text size
          axis.title = element_text(size = 18),  # Axis titles
          axis.text = element_text(size = 16),   # Axis tick labels
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16)) +
    labs(
      x = expression(a[0]),
      y = ""
    )
})

plots_norm_after_PSM <- lapply(seq_along(c0), function(i) {
  ggplot() +
    geom_density(data = data.frame(y = draws.a0.post.norm_after_PSM[[i]][,1]),
                 aes(x = y, color = "a01"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = draws.a0.post.norm_after_PSM[[i]][,2]),
                 aes(x = y, color = "a02"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = draws.a0.post.norm_after_PSM[[i]][,3]),
                 aes(x = y, color = "a03"),
                 linewidth = 1) +
    geom_density(data = data.frame(y = draws.a0.post.norm_after_PSM[[i]][,4]), 
                 aes(x = y, color = "a04"),
                 linewidth = 1) + 
    scale_color_manual(
      name = NULL,
      values = c(
        "a01" = "black",
        "a02" = brewer.pal(9, "Set1")[7],
        "a03" = brewer.pal(9, "Set1")[8],
        "a04" = brewer.pal(9, "Set1")[9]
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
    theme_bw() +
    theme(legend.position = 'none') +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
    ) +
    xlim(c(0,1)) +
    theme(text = element_text(size = 16),        # Base text size
          axis.title = element_text(size = 18),  # Axis titles
          axis.text = element_text(size = 16),   # Axis tick labels
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16)) +
    labs(
      x = expression(a[0]),
      y = ""
    )
})

plot_before_after_PSM <- 
  ((plot_wip + 
      ggtitle("Before PSM", subtitle = "Cauchy") + 
      xlab("")) +
     (plots_norm[[1]] + 
        ggtitle(label = NULL,
                subtitle = bquote("Normal" * "(" * 0 * "," * .(c0[1])^2 * I[p^"(m)"]*")")) +
        xlab("")) +
     (plots_norm[[2]] + 
        ggtitle(label = NULL,
                subtitle = bquote("Normal" * "(" * 0 * "," * .(c0[2])^2 * I[p^"(m)"]*")")) + 
        xlab("")) +
     (plots_norm[[3]] + 
        ggtitle(label = NULL,
                subtitle = bquote("Normal" * "(" * 0 * "," * .(c0[3])^2 * I[p^"(m)"]*")")) + 
        xlab("") 
     ) +
     (plots_norm[[4]] + 
        ggtitle(label = NULL,
                subtitle = bquote("Normal" * "(" * 0 * "," * .(c0[4])^2 * I[p^"(m)"]*")"))) +
     plot_layout(ncol = 1) & 
     theme(legend.position = 'none') 
  ) |
  ((((plot_wip_after_PSM + 
        ggtitle("After PSM", subtitle = "Cauchy") + 
        xlab("")) +
       (plots_norm_after_PSM[[1]] + 
          ggtitle(label = NULL,
                  subtitle = bquote("Normal" * "(" * 0 * "," * .(c0[1])^2 * I[p^"(m)"]*")")) +
          xlab("")) +
       (plots_norm_after_PSM[[2]] + 
          ggtitle(label = NULL,
                  subtitle = bquote("Normal" * "(" * 0 * "," * .(c0[2])^2 * I[p^"(m)"]*")")) +
          xlab("")) +
       (plots_norm_after_PSM[[3]] + 
          ggtitle(label = NULL,
                  subtitle = bquote("Normal" * "(" * 0 * "," * .(c0[3])^2 * I[p^"(m)"]*")")) +
          xlab("")
       )) & 
      theme(legend.position = 'none') ) +
     (plots_norm_after_PSM[[4]] + 
        ggtitle(label = NULL,
                subtitle = bquote("Normal" * "(" * 0 * "," * .(c0[4])^2 * I[p^"(m)"]*")")) + 
        theme(legend.position = c(0.5, 0.7),
              legend.background = element_rect(fill = "white", color = "black"))) +
     plot_layout(ncol = 1)
  )

plot_before_after_PSM

ggsave("logistic_regression/figures/post_a0_before_after_PSM.png",
       plot_before_after_PSM, 
       width = 14, height = 20, units = "in", dpi = 300)
