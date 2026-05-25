# Bayesian Logistic Regression

This directory contains the ACTG binary-outcome analysis used in Sections 3.2-3.4 and 4 of the manuscript. The analysis illustrates Bayesian dynamic borrowing for logistic regression, with ACTG036 treated as the current study and ACTG019 treated as the historical study.

The workflow includes Bayesian subset selection, posterior model probabilities, Bayesian model averaging, treatment-effect estimation, prior sensitivity analyses, propensity score matching, and prior predictive checks.

## Directory Contents

```text
logistic_regression/
|-- code/
|   |-- aux_scripts/                  # Shared R helpers plus Stan model files
|   |-- actg_PSM.R                    # Propensity score matching for ACTG019
|   |-- actg_samples.R                # Baseline NPP posterior model probabilities and BMA inputs
|   |-- actg_ate_bma.R                # ATE and odds-ratio plots before PSM
|   |-- actg_ate_after_PSM.r          # ATE and odds-ratio plots after PSM
|   |-- actg_ate_prior_comp.R         # ATE comparison under alternative priors
|   |-- actg_cauchy_prior.R           # Cauchy-prior posterior sampling
|   |-- actg_comp_ate_est.R           # Frequentist ATE estimator comparison
|   |-- actg_model_probabilities.r    # Prior/posterior model probability comparisons
|   |-- actg_norm_sens.R              # Normal-prior hyperparameter sensitivity analysis
|   |-- actg_post_a0_plots.r          # Posterior plots for the borrowing parameter a0
|   |-- actg_ppc.r                    # Prior predictive sampling
|   |-- actg_ppc_plots.R              # Prior predictive check figures
|   |-- actg_sample_npp_post.r        # Posterior samples of a0 under NPP settings
|   `-- actg_sample_npp_prior.r       # Prior samples under NPP settings
|
|-- data/                             # Matched historical data sets
|-- figures/                          # Generated plots used for diagnostics and manuscript figures
|-- samples/                          # Saved posterior, prior, BMA, PPC, ROC, and AUC objects
|-- vignettes/                        # R Markdown walkthroughs
`-- README.md
```

## What Has Been Done

### 1. ACTG Data Setup

The scripts use `actg036` as the current trial and `actg019` as the historical trial from `hdbayes`. Age and CD4 are standardized before fitting models. Several analyses are repeated after propensity score matching, using matched historical data saved in `data/`.

### 2. Bayesian Subset Selection

The subset-selection workflow enumerates candidate logistic-regression models over available covariates, computes historical-data prior model weights, fits NPP posterior models, estimates marginal likelihoods, and converts these into posterior model probabilities.

The conceptual walkthrough is in [vignettes/var_sel.Rmd](vignettes/var_sel.Rmd). The production script is [code/actg_samples.R](code/actg_samples.R), with reusable functions in [code/aux_scripts/functions.R](code/aux_scripts/functions.R).

### 3. Bayesian Model Averaging for Treatment Effects

Posterior draws from each selected model are converted into marginal predicted means under treatment and control. These model-specific estimates are combined using posterior model probabilities to obtain BMA samples for:

- Mean event probabilities by treatment arm.
- Average treatment effect (ATE).
- Odds ratio (OR).

The walkthrough is [vignettes/ate.Rmd](vignettes/ate.Rmd). Main plotting scripts include [code/actg_ate_bma.R](code/actg_ate_bma.R) and [code/actg_ate_after_PSM.r](code/actg_ate_after_PSM.r).

### 4. Propensity Score Matching

[code/actg_PSM.R](code/actg_PSM.R) creates matched versions of the historical ACTG019 data. Downstream scripts repeat the subset-selection, BMA, ATE, odds-ratio, and PPC analyses before and after PSM to evaluate sensitivity to historical-current covariate overlap.

### 5. Prior Sensitivity

Several scripts compare model probabilities and treatment-effect estimates under different prior choices:

- Normal-prior scale sensitivity in [code/actg_norm_sens.R](code/actg_norm_sens.R).
- Cauchy-prior analyses in [code/actg_cauchy_prior.R](code/actg_cauchy_prior.R).
- ATE prior comparisons in [code/actg_ate_prior_comp.R](code/actg_ate_prior_comp.R).
- Model-probability plots in [code/actg_model_probabilities.r](code/actg_model_probabilities.r).

Generated figures are saved under `figures/`.

### 6. Borrowing Parameter Analysis

The analysis studies the normalized power prior borrowing parameter `a0`. Scripts sample from the NPP prior and posterior, then visualize how the distribution changes before and after PSM and under alternative hyperparameters.

Relevant scripts include [code/actg_sample_npp_prior.r](code/actg_sample_npp_prior.r), [code/actg_sample_npp_post.r](code/actg_sample_npp_post.r), and [code/actg_post_a0_plots.r](code/actg_post_a0_plots.r).

### 7. Prior Predictive Checks

The PPC workflow samples prior predictive event probabilities under several Beta hyperparameters for `a0`, then produces:

- Observation-level prior predictive density plots.
- ROC and AUC summaries.
- Before/after PSM comparisons.

The walkthrough is [vignettes/ppc.Rmd](vignettes/ppc.Rmd). Main scripts are [code/actg_ppc.r](code/actg_ppc.r) and [code/actg_ppc_plots.R](code/actg_ppc_plots.R).

## Suggested Run Order

For a full rerun from the repository root:

1. Run `code/actg_PSM.R` to regenerate matched ACTG019 data.
2. Run `code/actg_samples.R` for baseline NPP posterior model probabilities and BMA inputs.
3. Run prior-sensitivity scripts as needed: `actg_cauchy_prior.R`, `actg_norm_sens.R`, and `actg_ate_prior_comp.R`.
4. Run `code/actg_sample_npp_prior.r` and `code/actg_sample_npp_post.r` for `a0` prior/posterior analyses.
5. Run `code/actg_ate_bma.R`, `code/actg_ate_after_PSM.r`, and `code/actg_model_probabilities.r` for final ATE, OR, and model-probability figures.
6. Run `code/actg_ppc.r` and `code/actg_ppc_plots.R` for PPC outputs.
7. Render the vignettes in `vignettes/` for narrative explanations.

## Outputs

- `samples/` contains saved `.RData` objects for posterior samples, BMA draws, prior/posterior `a0` draws, PPC objects, ROC curves, and AUC intervals.
- `figures/` contains generated figures for posterior distributions, ATE/OR estimates, model probabilities, prior predictive checks, and before/after PSM comparisons.
- `figures/ppc/` and `figures/ppc_a0/` contain individual-level PPC plots for multiple observations and PSM settings.

## Notes

Many scripts are computationally intensive. Several use MCMC with `iter_warmup = 1000`, `iter_sampling = 2500`, and parallel execution. Adjust the number of cores in scripts before running on a laptop or shared machine.
