library(hdbayes)
library(posterior)
library(dplyr)
library(parallel)
library(ggplot2)

##################################
### Aux functions
## function to pull out the summaries in a convenient form
get_summaries <- function(fit, pars.interest, digits = 3) {
  fit %>%
    select(all_of(pars.interest)) %>%
    summarise_draws(mean, sd, ~quantile(.x, probs = c(0.025, 0.5, 0.975)),
                    .num_args = list(digits = digits))
}
base.pars      <- c("(Intercept)", "age",
                    "race", "treatment",
                    "cd4")
base.pars.hist <- paste(base.pars, "hist", "1", sep = "_")
##################################

data(actg019)
data(actg036)

hist <-  actg019 %>% mutate(age = scale(age),
                            cd4 = scale(cd4))

curr <- actg036 %>% mutate(age = scale(age),
                           cd4 = scale(cd4))

stacked <- rbind(curr, hist)

formula <- outcome ~  age + race + treatment + cd4
p       <- length(attr(terms(formula), "term.labels")) ## number of predictors
family  <- binomial('logit')

fit.mle.curr  <- glm(formula, family, curr)
fit.mle.hist <- glm(formula, family, hist)
fit.mle.stacked <- glm(formula, family, stacked)

all.data <- list(curr, hist)

round(confint(fit.mle.hist), 3)
round(confint(fit.mle.curr), 3)
round(confint(fit.mle.stacked), 3)

ncores        <- 10
nchains       <- 4 ## number of Markov chains
iter_warmup   <- 1000 ## warmup per chain for MCMC sampling
iter_sampling <- 2500 ## number of samples post warmup per chain

## fit individual models

bayes.fit.curr  <- glm.pp(
  formula = formula, family = family,
  data.list = all.data,
  a0 = 0, ## should zero the likellidhood of hist data
  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
  chains = nchains, parallel_chains = ncores,
  refresh = 0
)

bayes.fit.stacked  <- glm.pp(
  formula = formula, family = family,
  data.list = list(stacked, curr),
  a0 = 0, ## should zero the likellidhood of hist data
  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
  chains = nchains, parallel_chains = ncores,
  refresh = 0
)

get_summaries(bayes.fit.curr,
              pars.interest = base.pars)
confint(fit.mle.curr)

get_summaries(bayes.fit.stacked,
              pars.interest = base.pars)
confint(fit.mle.stacked)

## fit BHM
fit.bhm <- glm.bhm(
  formula, family, all.data,
  meta.mean.mean = 0, meta.mean.sd = 10,
  meta.sd.mean = 0, meta.sd.sd = 0.5,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
  chains = nchains, parallel_chains = ncores,
  refresh = 0
)

get_summaries(fit = fit.bhm,
              pars.interest = c(base.pars, base.pars.hist))


## fit commensurate prior
fit.commensurate <- glm.commensurate(
  formula = formula, family = family, data.list = all.data,
  p.spike = 0.1, spike.mean = 200, spike.sd = 0.1,
  slab.mean = 0, slab.sd = 5,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
  chains = nchains, parallel_chains = ncores,
  refresh = 0
)
base.pars.hist <- paste(base.pars, "hist", sep = "_")
get_summaries(fit = fit.commensurate,
              pars.interest = c(base.pars, base.pars.hist))


## fit RMAP
res.rmap <- glm.rmap(
  formula = formula, family = family, data.list = all.data,
  w = 0.1,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
  chains = nchains, parallel_chains = ncores,
  refresh = 0
)
fit.rmap <- res.rmap[["post.samples"]]
get_summaries(fit.rmap, pars.interest = base.pars)


## fit PP
n0      <- nrow(actg019)
n       <- nrow(actg036)
a0.star <- (n/n0) * 1/2
fit.pp  <- glm.pp(
  formula = formula, family = family, data.list = all.data,
  a0 = a0.star, ## discounting parameter
  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
  chains = nchains, parallel_chains = ncores,
  refresh = 0
)
get_summaries(fit.pp,
              pars.interest = base.pars)


## generate plot of prior density on a_{01} (Figure 1)
library(latex2exp)
devtools::source_url("https://raw.githubusercontent.com/maxbiostat/logPoolR/main/R/beta_elicitator.r")
beta.pars <- elicit_beta_mean_cv(m0 = a0.star, cv = 1)
curve(dbeta(x, shape1 = beta.pars$a, shape2 = beta.pars$b),
      lwd = 3,
      main = TeX("Prior on $a_{0 1}$"),
      ylab = "Density",
      xlab = TeX("$a_{0 1}$"))
abline(v = a0.star, lwd = 2, lty = 3)
legend(x = "topright",
       legend = c(TeX("$\\pi(a_{0 1})$"),
                  TeX("$a_{0 1} = (1/2) (n/n_0)$")),
       lwd = 2, lty = c(1, 3),
       bty = 'n')

a0       <- seq(0, 1, length.out = 21)
histdata <- all.data[[2]]
## wrapper to obtain log normalizing constant in parallel package
logncfun <- function(a0, ...){
  hdbayes::glm.npp.lognc(
    formula = formula, family = family, histdata = histdata, a0 = a0, ...
  )
}
cl <- makeCluster(10)
clusterSetRNGStream(cl, 123)
clusterExport(cl, varlist = c('formula', 'family', 'histdata'))
## call created function
a0.lognc <- parLapply(
  cl = cl, X = a0, fun = logncfun, iter_warmup = 2*iter_warmup,
  iter_sampling = 2*iter_sampling, chains = nchains, refresh = 0
)
stopCluster(cl)

a0.lognc <- data.frame( do.call(rbind, a0.lognc) )
head(a0.lognc) %>%
    mutate(across(where(is.numeric), round, 3))

## generate plot of log normalizing constant versus a_{01} (Figure 2)
ggplot(data = a0.lognc, aes(x = a0, y = lognc)) +
  geom_point(alpha = 0.75) +
  labs(x = TeX("$a_{0 1}$"),
       y = "log normalizing constant") +
  theme_bw(base_size = 15)


## fit NPP
fit.npp <- glm.npp(
  formula = formula, family = family, data.list = all.data,
  a0.lognc = a0.lognc$a0,
  lognc = matrix(a0.lognc$lognc, ncol = 1),
  a0.shape1 = beta.pars$a, a0.shape2 = beta.pars$b, ## beta prior on a_{01}
  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
  chains = nchains, parallel_chains = ncores,
  refresh = 0
)
get_summaries(fit = fit.npp,
              pars.interest = c(base.pars, "a0_hist_1"))


## fit NAPP
fit.napp <- glm.napp(
  formula = formula, family = family, data.list = all.data,
  a0.shape1 = beta.pars$a, a0.shape2 = beta.pars$b,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
  chains = nchains, parallel_chains = ncores,
  refresh = 0
)
get_summaries(fit = fit.napp,
              pars.interest = c(base.pars, "a0_hist_1"))


## fit LEAP
fit.leap <- glm.leap(
  formula = formula, family = family, data.list = all.data,
  K = 2, prob.conc = rep(1, 2),
  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
  chains = nchains, parallel_chains = ncores,
  refresh = 0
)
get_summaries(fit = fit.leap,
              pars.interest = c(base.pars, "gamma"))


## generate plot of coefficient estimates for each of the models (Figure 3)
fit.list <- list('BHM' = fit.bhm,
                 'Commensurate' = fit.commensurate,
                 'RMAP' = fit.rmap,
                 'NAPP' = fit.napp,
                 'NPP' = fit.npp,
                 'PP' = fit.pp,
                 'LEAP' = fit.leap,
                 'Current' = bayes.fit.curr,
                 'Pooled' = bayes.fit.stacked)

## function to get current data regression estimates (posterior mean, sd, and 95% CI) for each model
get_quants <- function(x, name){
  pars <- names(coef(fit.mle.curr))
  summs <- get_summaries(x, pars.interest = pars)
  out <- tibble::tibble(
    parameter = pars,
    mean = as.numeric(unlist(summs[, "mean"])),
    sd = as.numeric(unlist(summs[, "sd"])),
    lwr = as.numeric(unlist(summs[, "2.5%"])),
    upr = as.numeric(unlist(summs[, "97.5%"]))
  )
  out$model <- name
  return(out)
}
model.names  <- names(fit.list)
results.list <- lapply(seq_along(fit.list), function(i){
  get_quants(x = fit.list[[i]], name = model.names[i])
})
results      <- do.call(rbind, results.list)

mle.estimates.past <- data.frame(
  parameter = names(coef(fit.mle.hist)),
  value = coef(fit.mle.hist),
  MLE = "Historical"
)
mle.estimates.curr <- data.frame(
  parameter = names(coef(fit.mle.curr)),
  value = coef(fit.mle.curr),
  MLE = "Current"
)

mle.estimates.stacked <- data.frame(
  parameter = names(coef(fit.mle.stacked)),
  value = coef(fit.mle.stacked),
  MLE = "Stacked"
)
results$model <- as.factor(results$model)
LL <- unique(results$model)


library(scales)
library(ggthemes)

mle.estimates <- rbind(mle.estimates.past,
                       mle.estimates.curr)
ggplot() +
  geom_pointrange(data = results,
                  mapping = aes(x = model,
                                y = mean,
                                ymin = lwr, ymax = upr,
                                colour = model)) +
  scale_color_colorblind() +
  geom_hline(data = mle.estimates,
             aes(yintercept = value, linetype = MLE)) +
  scale_linetype_manual(values = c("dotted", "dashed")) +
  scale_y_continuous("") +
  scale_x_discrete("") +
  facet_wrap(parameter~., scales = "free_y") +
  ggtitle("") +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

results <- results %>% mutate(OR = exp(mean),
                              OR_lwr = exp(lwr),
                              OR_upr = exp(upr))

mle.estimates <- mle.estimates %>% mutate(OR = exp(value))

forPlotting <- subset(results,
                      parameter != "(Intercept)" &
                        model %in% c("BHM", "RMAP",
                                     "NPP", "LEAP",
                                     "Current", "Pooled"))

ggplot() +
  geom_pointrange(data = forPlotting,
                  mapping = aes(x = model,
                                y = OR,
                                ymin = OR_lwr, ymax = OR_upr,
                                colour = model)) +
  scale_color_colorblind() +
  geom_hline(data = subset(mle.estimates, parameter != "(Intercept)"),
             aes(yintercept = exp(value), linetype = MLE)) +
  scale_linetype_manual(values = c("dotted", "dashed")) +
  scale_y_continuous("Odds ratio") +
  scale_x_discrete("") +
  facet_wrap(parameter~., scales = "free_y") +
  ggtitle("") +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggplot() +
  geom_pointrange(data = subset(forPlotting,
                                ! model %in% c("Current") ),
                  mapping = aes(x = model,
                                y = OR,
                                ymin = OR_lwr, ymax = OR_upr,
                                colour = model)) +
  scale_color_colorblind() +
  geom_hline(data = subset(mle.estimates, parameter != "(Intercept)"),
             aes(yintercept = exp(value), linetype = MLE)) +
  scale_linetype_manual(values = c("dotted", "dashed")) +
  scale_y_continuous("Odds ratio") +
  scale_x_discrete("") +
  facet_wrap(parameter~., scales = "free_y") +
  ggtitle("") +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
