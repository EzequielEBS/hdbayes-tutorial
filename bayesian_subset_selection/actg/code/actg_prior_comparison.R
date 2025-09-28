load("bayesian_subset_selection/actg/samples/post_samples_d0.RData")

# Plot using ggplot2
ate_cauchy <- plot_ate(post_samples_cauchy)

ate_cauchy + ggtitle("Cauchy prior on coefficients")
ggsave("bayesian_subset_selection/actg/results/figures/ate_cauchy.png",
       ate_cauchy, width = 10, height = 7, units = "in", dpi = 300)

ate_d0 <- lapply(post_samples_d0, function(x) {
  plot_ate(x)
})

all_ate_d0 <- (ate_d0[[1]] + ggtitle(expression(psi[0] == 0.25))+ xlim(-0.1, 0.05) + xlab("")) /
  (ate_d0[[2]] + ggtitle(expression(psi[0] == 0.5)) + xlim(-0.1, 0.05) + xlab("")) /
  (ate_d0[[3]] + ggtitle(expression(psi[0] == 1)) + xlim(-0.1, 0.05) +xlab("")) /
  (ate_d0[[4]] + ggtitle(expression(psi[0] == 2)) + xlim(-0.1, 0.05)) +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')
all_ate_d0

ggsave("bayesian_subset_selection/actg/results/figures/ate_psi0.png",
       all_ate_d0, width = 10, height = 25, units = "in", dpi = 300)
