# load("bayesian_subset_selection/actg/samples/post_samples_d0.RData") # psi0
# load("bayesian_subset_selection/actg/samples/post_samples_c0.RData") # tau0
load("bayesian_subset_selection/actg/samples/post_samples_c0d0.RData") 
# load("bayesian_subset_selection/actg/samples/post_samples_m0v0.RData")
load("bayesian_subset_selection/actg/samples/post_samples_wip.RData")
load("bayesian_subset_selection/actg/samples/post_samples.RData")
# compare model probabilities

map_models <- post_samples$df_post_ord %>% 
  mutate(model_id = paste0("M", row_number())) %>%
  select(model_id, model)

prob_wip <- post_samples_wip$df_post %>%
  left_join(map_models, by = "model") %>%
  select(model_id, prior_model, post_model) %>%
  pivot_longer(
    cols = c(prior_model, post_model),
    names_to = "type",
    values_to = "value"
  ) %>%
  mutate(model_id = factor(model_id, levels = map_models$model_id))

prob_wip$type <- recode(prob_wip$type,
                        prior_model = "Prior",
                        post_model = "Posterior")
prob_wip$type <- factor(prob_wip$type, levels = c("Prior", "Posterior"))

models_wip <- ggplot(prob_wip, aes(x = model_id, y = value, fill = type)) +
  geom_col(position = "dodge") + # side-by-side
  labs(x = "Model", y = "Value", fill = "") +
  theme_minimal() +
  scale_fill_manual(values = c("Prior" = "#66A8D0", "Posterior" = "#D08E66"))
models_wip

c0 <- c(0.25, 0.5, 1, 2)
models_norm <- lapply(seq_along(c0)
                  , function(i) {
  post <- post_samples_c0d0[[i]]
  c0 <- c0[i]
  df <- post$df_post %>%
    left_join(map_models, by = "model") %>%
    select(model_id, prior_model, post_model) %>%
    pivot_longer(
      cols = c(prior_model, post_model),
      names_to = "type",
      values_to = "value"
    ) %>%
    mutate(model_id = factor(model_id, levels = map_models$model_id)) %>%
    mutate(c0 = c0)
  df$type <- recode(df$type,
                  prior_model = "Prior",
                  post_model = "Posterior")
  df$type <- factor(df$type, levels = c("Prior", "Posterior"))
  ggplot(df, aes(x = model_id, y = value, fill = type)) +
    geom_col(position = "dodge") + # side-by-side
    labs(x = "Model", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = c("Prior" = "#66A8D0", "Posterior" = "#D08E66"))
})

models_wip_norm <- (models_wip + 
                      ggtitle("Cauchy") +
                      labs(x = "",
                           y = "") 
                    ) /
  (models_norm[[1]] + 
     ggtitle(
       bquote("Normal" * "(" * 0 * "," * .(c0[1])^2 * I[p]*")")
     ) +
     labs(x = "",
          y = "")) +
  (models_norm[[2]] + 
     ggtitle(
       bquote("Normal" * "(" * 0 * "," * .(c0[2])^2 * I[p]*")")
     ) +
     labs(x = "",
          y = "")) +
  (models_norm[[3]] + 
     ggtitle(
       bquote("Normal" * "(" * 0 * "," * .(c0[3])^2 * I[p]*")")
     ) +
     labs(x = "",
          y = "")) +
  (models_norm[[4]] + 
     ggtitle(
       bquote("Normal" * "(" * 0 * "," * .(c0[4])^2 * I[p]*")")
     ) +
     labs(y = "")) +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')
models_wip_norm

ggsave("bayesian_subset_selection/actg/results/figures/models_wip_norm.png",
       models_wip_norm, width = 10, height = 15, units = "in", dpi = 300)

df_norm <- post_samples$df_post_ord %>%
  left_join(map_models, by = "model") %>%
  mutate(logml = log(ml)) %>%
  select(model_id, model, prior_model, logml, post_model)

print(xtable::xtable(df_norm, digits = c(3,3,3,3,2,3)), include.rownames = FALSE) 

df_models <- post_samples$df_post_ord %>%
  left_join(map_models, by = "model") %>%
  left_join(post_samples_wip$df_post, by = "model") %>%
  select(model_id, prior_model.x, post_model.x, prior_model.y, post_model.y) %>%
  rename(prior_norm = prior_model.x,
         post_norm = post_model.x,
         prior_wip = prior_model.y,
         post_wip = post_model.y)
print(xtable::xtable(df_models, digits = c(3,3,3,3,3,3)), include.rownames = FALSE)
