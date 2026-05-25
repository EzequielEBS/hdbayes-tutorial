# Survival Analysis of the E1684 and E1690 Trials

This directory contains the code, figures, and report for the relapse-free survival analysis described in Section 3.1 and Appendix A of the manuscript.

The analysis treats E1690 as the current study and E1684 as the external study. It illustrates Bayesian dynamic borrowing for survival outcomes, focusing on two-year relapse-free survival (RFS) probabilities and treatment effects.

## Directory Contents

```text
survival_analysis/
|-- code/
|   |-- 01_grid_setup.R                 # Defines model/prior/J/arm grid
|   |-- 02_fit_PWE_CurePWE.R            # Fits PWEPH and CurePWEPH models
|   |-- 03_compile_analysis_results.R   # Compiles fit outputs and summaries
|   |-- grid.rds                        # Saved analysis grid
|   `-- wrapper/
|       |-- get_strata_data.R           # PSIPP stratification helper
|       |-- pwe_get_trt_effect.R        # Treatment-effect helper for PWEPH
|       `-- curepwe_get_trt_effect.R    # Treatment-effect helper for CurePWEPH
|
|-- figures/                            # KM, ELPD, ESS, and posterior-density plots
|-- report/
|   |-- survival_analysis_E1684_E1690.Rmd
|   `-- survival_analysis_E1684_E1690.html
`-- README.md
```

## What Has Been Done

### 1. Data Preparation and Descriptive Survival Analysis

The report loads E1684 and E1690 data from `hdbayes`, treats E1684 as external and E1690 as current, adjusts zero-day failure times, centers and scales age, and summarizes baseline covariates by study.

It also includes Kaplan-Meier curves for treatment and control arms in both studies. The rendered figure is saved as `figures/KM_plots.png`.

### 2. Cox Proportional Hazards Benchmarks

The report fits Cox proportional hazards models separately to E1690 and E1684 using:

- Treatment indicator.
- Sex.
- Standardized age.
- Indicator for more than one cancerous lymph node.

These models provide conventional benchmarks before the Bayesian dynamic borrowing analysis.

### 3. Bayesian Survival Outcome Models

Two outcome-model families are evaluated:

- PWEPH: proportional hazards model with a piecewise constant baseline hazard.
- CurePWEPH: mixture cure-rate extension of PWEPH, where a fraction of patients is modeled as cured.

The piecewise baseline hazard uses `J` intervals. The analysis evaluates `J = 2, ..., 9`.

### 4. Dynamic Borrowing Priors

For each model family and interval count, the analysis compares:

- Vague/non-informative prior (`ref`).
- Power prior (`pp`).
- Propensity score-integrated power prior (`psipp`).
- Bayesian hierarchical model (`bhm`).
- Commensurate prior (`cp`).
- Latent exchangeability prior (`leap`).
- Normalized power prior (`npp`).

Non-PSIPP models adjust directly for baseline covariates. PSIPP first estimates study-membership propensity scores, stratifies subjects by propensity-score quantiles, and applies stratum-specific borrowing.

### 5. Grid-Based Model Fitting

[code/01_grid_setup.R](code/01_grid_setup.R) creates a full grid over model family, prior, interval count, and treatment arm, then saves it to `code/grid.rds`.

[code/02_fit_PWE_CurePWE.R](code/02_fit_PWE_CurePWE.R) fits the grid of PWEPH and CurePWEPH models. This is the main computational step and may be intended for a computing cluster.

[code/03_compile_analysis_results.R](code/03_compile_analysis_results.R) combines model outputs, computes ELPD summaries, and prepares posterior summaries of two-year RFS probabilities.

### 6. Model Selection and Treatment Effects

The report selects interval counts using expected log predictive density (ELPD). It then summarizes posterior two-year RFS probabilities for the treatment and control arms and computes the treatment effect:

```text
Delta = S_treatment(2 years) - S_control(2 years)
```

The ELPD figure is saved as `figures/elpd_plots.png`.

### 7. Effective Sample Size and Borrowing

The analysis computes effective sample size (ESS) for each model-prior combination relative to the vague-prior posterior variance. It reports the effective number of external patients borrowed as:

```text
ESS - n_current
```

Figures include:

- `figures/ess_plot_heatmap.png`
- `figures/ess_plot.png`

### 8. Final Posterior Density Comparison

The report compares posterior densities of the treatment effect under an ELPD-selected optimal model-prior pair and a vague-prior reference. The figure is saved as `figures/post_dens_diff.png`.

## Suggested Run Order

From the `survival_analysis/` directory or from the repository root with paths adjusted as needed:

1. Run `code/01_grid_setup.R`.
2. Run `code/02_fit_PWE_CurePWE.R` to fit all model-prior-grid combinations.
3. Run `code/03_compile_analysis_results.R` to combine model outputs.
4. Render `report/survival_analysis_E1684_E1690.Rmd`.

## Reproducibility Notes

The repository includes the grid, report, rendered HTML, and figures. Individual model fit objects and some compiled results are not stored because they are large. They can be regenerated from the scripts above.

The report refers to compiled results under `results/compiled_results/compiled_analysis_results.rds`; regenerate that file before rendering the complete report if it is not present locally.
