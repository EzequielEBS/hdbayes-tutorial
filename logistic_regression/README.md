# Bayesian logistic regression

This directory contains scripts and outputs related to the **logistic regression** analysis performed in the sections 3.2-3.4 and 4.

---

## 📂 Contents

```
logistic_regression/
├── code/                               # R scripts used to run the analysis
|   ├── aux_scripts/                    # Auxiliary R scripts and stan models to perform the analysis (including modifications of functions from hdbayes package)
│   │
│   ├── actg_ate_after_PSM.r            # Generates plots for the posterior distribution of the ATE and odds ratio after PSM
│   ├── actg_ate_bma.r                  # Generates plots for the posterior distribution of the ATE and odds ratio before PSM
│   ├── actg_ate_prior_comp.r           # Generates ATE plots under different priors
│   ├── actg_cauchy_prior.r             # Samples with Cauchy prior
│   ├── actg_com_ate_est.r              # Compares frequentist ATE estimators
│   ├── actg_model_probabilities.r      # Compares model probabilities under different priors
│   ├── actg_norm_sens.r                # Sensitivity analyses for hyperparameters
│   ├── actg_post_a0_plots.r            # Generates posterior plots for a0 before and after PSM under different priors
│   ├── actg_ppc.r                      # Generates prior predictive distributions 
│   ├── actg_ppc_plots.r                # Generates prior predictive checks plots
│   ├── actg_sample_npp_post.r          # Samples from the posterior distribution of a0 under different priors
│   ├── actg_sample_npp_prior.r         # Samples from the NPP prior under different settings
│   └── actg_samples.r                  # Generates posterior model probabilities and ATE useful quantities
│
├── data/                               # Store data after Propensity Score Matching (PSM)
│
├── figures/                            # Store generated figures 
│
├── samples/                            # Store generated samples
|
├── vignettes/                          # Vignettes to explain main analyses
│
└── README.md                           
```

---

## 📘 Related Vignettes

| Vignette | Description |
|-----------|--------------|
| [`vignettes/var_sel.Rmd`](vignettes/var_sel.Rmd) | Explains how to compute the quantities from the variable selection approach (section 3.3.1). |
| [`vignettes/ate.Rmd`](vignettes/ate.Rmd) | Explains how to compute posterior samples from ATE by using model probabilities and BMA (section 3.3.2). |
| [`vignettes/ppc.Rmd`](vignettes/ppc.Rmd) | Explains how to perform prior predicitive checks (PPC) (section 3.3.3). |

