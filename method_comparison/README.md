# Method Comparison

This directory contains a compact ACTG analysis script comparing several Bayesian dynamic borrowing approaches for logistic regression.

The script [ACTG_analysis_four_methods.r](ACTG_analysis_four_methods.r) uses ACTG036 as the current study and ACTG019 as the historical study, with the binary outcome modeled by logistic regression.

## What Has Been Done

The script:

- Loads `actg019` and `actg036` from `hdbayes`.
- Standardizes age and CD4.
- Fits classical logistic regressions to the current, historical, and stacked data.
- Fits current-study-only and pooled Bayesian logistic models.
- Fits several dynamic borrowing models from `hdbayes`.
- Summarizes posterior means, standard deviations, and credible intervals for regression coefficients.
- Converts coefficient summaries into odds-ratio summaries.
- Creates comparison plots of coefficient and odds-ratio estimates across models.

## Methods Compared

The script includes:

- Current: Bayesian model using only the current study.
- Pooled: Bayesian model using pooled current and historical data.
- BHM: Bayesian hierarchical model.
- Commensurate prior.
- RMAP: robust meta-analytic predictive prior.
- PP: power prior.
- NPP: normalized power prior.
- NAPP: normalized asymptotic power prior.
- LEAP: latent exchangeability prior.

It also computes frequentist maximum-likelihood estimates for current, historical, and stacked data as reference lines in the plots.

## Main Objects and Outputs

Key fitted objects include:

- `fit.mle.curr`, `fit.mle.hist`, and `fit.mle.stacked`
- `bayes.fit.curr` and `bayes.fit.stacked`
- `fit.bhm`
- `fit.commensurate`
- `fit.rmap`
- `fit.pp`
- `fit.npp`
- `fit.napp`
- `fit.leap`

The script builds a `fit.list`, extracts posterior summaries with `get_summaries()`, and creates plotting data frames for coefficient and odds-ratio comparisons.

## Notes

This file is written as an exploratory comparison script rather than a fully modular pipeline. It includes a remote `devtools::source_url()` call to load a beta-prior elicitation helper, so rerunning it requires internet access unless that helper is copied locally.
