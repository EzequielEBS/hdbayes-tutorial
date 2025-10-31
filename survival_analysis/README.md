# Survival analysis using piecewise exponential and mixture cure rate models

This directory contains the **code, figures, and report** for the survival analysis described in **Section 3.1**
of the manuscript. The analysis illustrates Bayesian dynamic borrowing methods applied to the E1690 melanoma 
trial (current study) while incorporating information from the E1684 trial (external study). The accompanying 
report also includes additional analyses summarized in the **Appendix A** of the paper.

---

## 📂 Contents

```
survival_analysis/
│
├── code/                                   # Scripts and helper functions to perform the analysis
│   ├── 01_grid_setup.R                     # Step 1: Set up the grid of model-prior combinations for the analysis
│   ├── 02_fit_PWE_CurePWE.R                # Step 2: Fit PWEPH and CurePWEPH models under various priors
│   ├── 03_compile_analysis_results.R       # Step 3: Compile and summarize results across all fits
│   ├── grid.rds                            # Saved grid of analysis scenarios
│   └── wrapper/                            # Helper functions used in the main scripts
│       ├── pwe_get_trt_effect.R            # Functino to compute treatment effect for PWEPH models 
│       ├── curepwe_get_trt_effect.R        # Functino to compute treatment effect for CurePWEPH models
│       └── get_strata_data.R               # Helper function for implementing PSIPP
│
├── figures/                                # Plots generated from the compiled analysis
│
└── report/
│   ├── survival_analysis_E1684_E1690.Rmd   # Main R Markdown report (Section 3.1 + Appendix)
│   ├── survival_analysis_E1684_E1690.html  # Rendered HTML output
│
└── README.md
```

## 📘️ Reproducibility Workflow

The survival analysis is organized into **three main steps**, followed by report generation. Each step can be
reproduced using the provided R scripts.

1. **`01_grid_setup.R`**  
   Defines the grid of model, prior, treatment arm, and number of interval combinations (`J`), and saves it to 
   `grid.rds` for use in the subsequent steps.

2. **`02_fit_PWE_CurePWE.R`**  
   Fits two outcome models to the E1684 and E1690 data sets under multiple prior specifications:
   - **PWEPH model** (Proportional Hazards with Piecewise Constant Baseline Hazard)  
   - **CurePWEPH model** (Mixture Cure Rate Extension of PWEPH) 
   
   Supported priors include: *vague, power prior (PP), propensity score–integrated power prior (PSIPP), Bayesian hierarchical model (BHM), commensurate prior (CP), latent exchangeability prior (LEAP), normalized power prior (NPP).*  
   
3. **`03_compile_analysis_results.R`**  
   Combines outputs from all model–prior combinations, computes the **expected log predictive density (elpd)** 
   for model comparison, and summarizes posterior estimates of two-year relapse-free survival (RFS) probabilities.

4. **`report/survival_analysis_E1684_E1690.Rmd`**  
   Generates the final analysis report, including all figures and tables presented in the manuscript’s **Section 3.1** 
   and **Appendix A**.

---

**Note:**  

Individual model fit objects and compiled analysis results are **not stored in this repository** due to their large 
file sizes. However, all results can be **reproduced in full** by running the provided scripts sequentially.

---
