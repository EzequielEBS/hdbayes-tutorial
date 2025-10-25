################################################################################
# Checking ate calculations based on 
# https://solomonkurz.netlify.app/blog/2023-04-24-causal-inference-with-logistic-regression/
################################################################################


fit <- glm(outcome ~ ., 
           data = current_data[, colnames(current_data) != "treatment"], 
           family = binomial(link = "logit"))

pred_ctrl <- predict(fit, newdata = current_data_ctrl, type = "response")
pred_trt <- predict(fit, newdata = current_data_trt, type = "response")

mean_glm_ctrl <- mean(pred_ctrl)
mean_glm_trt <- mean(pred_trt)
ate_glm <- mean_glm_trt - mean_glm_ctrl

mean_glm_ctrl
mean_glm_trt
ate_glm


################################################################################

fit_ctrl <- glm(outcome ~ ., 
                data = current_data_ctrl, 
                family = binomial(link = "logit"))
fit_trt <- glm(outcome ~ ., 
               data = current_data_trt, 
               family = binomial(link = "logit"))

pred_ctrl <- predict(fit_ctrl, newdata = current_data, type = "response")
pred_trt <- predict(fit_trt, newdata = current_data, type = "response")

mean_glm_ctrl <- mean(pred_ctrl)
mean_glm_trt <- mean(pred_trt)
ate_glm <- mean_glm_trt - mean_glm_ctrl

mean_glm_ctrl
mean_glm_trt
ate_glm

##############################################################################


glm1 <- glm(
  data = current_data,
  family = binomial,
  outcome ~ treatment
)

# ANCOVA-type model
glm2 <- glm(
  data = current_data,
  family = binomial,
  outcome ~ .
)

# Confidence interval for odds-ratio
tidy(glm1, conf.int = T) %>% 
  filter(term == "treatment") %>% 
  select(estimate, starts_with("conf.")) %>% 
  mutate_all(exp)

# ATE
plogis(coef(glm1)[1] + coef(glm1)[2]) - plogis(coef(glm1)[1])

# Raw ATE
current_data %>% 
  group_by(treatment) %>% 
  summarise(p = mean(outcome == 1)) %>% 
  pivot_wider(names_from = treatment, values_from = p) %>% 
  mutate(sate = `1` - `0`)

# Predict
nd <- tibble(treatment = 0:1)

# log odds metric
predict(glm1, 
        newdata = nd,
        se.fit = TRUE) %>% 
  data.frame() %>% 
  bind_cols(nd)

# probability metric
predict(glm1, 
        newdata = nd,
        se.fit = TRUE,
        type = "response") %>% 
  data.frame() %>% 
  bind_cols(nd)

# sample statistics
current_data %>% 
  group_by(treatment) %>% 
  summarise(p = mean(outcome == 1))

# redefine the data grid
current_data$id <- 1:nrow(current_data)
nd <- current_data %>% 
  select(id) %>% 
  expand_grid(treatment = 0:1)

predict(glm1, 
        newdata = nd,
        se.fit = TRUE,
        # request the probability metric
        type = "response") %>% 
  data.frame() %>% 
  bind_cols(nd) %>% 
  # look at the first 6 rows
  head()

predict(glm1, 
        newdata = nd,
        se.fit = TRUE,
        type = "response") %>% 
  data.frame() %>% 
  bind_cols(nd) %>% 
  select(id, treatment, fit) %>% 
  pivot_wider(names_from = treatment, values_from = fit) %>% 
  summarise(ate = mean(`1` - `0`))


################################################################################

bind_rows(tidy(glm1), tidy(glm2)) %>% 
  filter(term == "treatment") %>% 
  mutate(fit = c("glm1", "glm2"),
         model_type = c("ANOVA", "ANCOVA")) %>%
  rename(`beta[1]` = estimate) %>% 
  select(fit, model_type, `beta[1]`, std.error)

get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

nd <- current_data %>% 
  select(age, race, cd4) %>% 
  expand_grid(treatment = 0:1)

bind_rows(
  avg_comparisons(glm1, newdata = nd, variables = "treatment", conf_level = 0.95),
  avg_comparisons(glm2, newdata = nd, variables = "treatment", conf_level = 0.95)
) %>% 
  data.frame() %>% 
  mutate(fit = c("glm1", "glm2"),
         model_type = c("ANOVA", "ANCOVA")) %>%
  rename(`tau[ATE]` = estimate) %>% 
  select(fit, model_type, `tau[ATE]`, std.error, conf.low, conf.high)