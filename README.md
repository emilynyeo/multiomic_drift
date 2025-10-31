# Multi-omic Predictions Workflow

This repository contains the complete workflow for processing multi-omic data, training prediction models (GLMMLASSO and MERF), and generating publication-ready figures.

## Overview

The workflow consists of four main stages:
1. **Processing Raw Data** → Prepare longitudinal and delta datasets
2. **GLMMLASSO Modeling** → Train GLMMLASSO models and extract feature importances
3. **MERF Modeling** → Train Mixed Effects Random Forest models with SHAP analysis
4. **Statistics and Plotting** → Generate ANOVA tables, feature importance plots, and combined figures

---

## 1. Processing Raw Data

### Input Data
Raw input data should be placed in the appropriate data directory before running processing scripts.

### Processing Scripts
**Location:** `play_scripts/1.processing_omic_inputs/`

#### Key Files:
- **`processing_utils.R`** - Utility functions for data processing
- **`all_processing_long_april2.R`** - Main processing script

#### Outputs:
- **`data/april_processing/long_april29.csv`** - Longitudinal formatted dataset
- **`data/april_processing/all_delta_april29.csv`** - Delta (change) formatted dataset

These processed datasets serve as inputs for all downstream modeling steps.

---

## 2. GLMMLASSO Modeling

GLMMLASSO models are built per-omic type, producing train/test splits and predictions with feature importance rankings for each omic subset.

### Long Models

**Script:** `play_scripts/2.models/glmmlasso/march30_long/march30_glmmlasso_long.R`

**Input:**
- Primary dataset: `data/april_processing/long_april29.csv`

**Outputs (Location: `play_scripts/2.models/glmmlasso/may_basic_plus/long/`):**
- Predictions: `june_lasso_long_predictions_df.csv`
- Feature importance files (per-omic subset):
  - `basic_gl_long_top_features.csv`
  - `meta_gl_long_top_features.csv`
  - `grs_gl_long_top_features.csv`
  - `taxa_gl_long_top_features.csv`
  - `micom_gl_long_top_features.csv`
  - `pathway_gl_long_top_features.csv`
  - `metabo_gl_long_top_features.csv`
  - `all_gl_long_top_features.csv`

### Delta Models

**Script:** `play_scripts/2.models/glmmlasso/march30_delta/glmmlasso_deltas.R`

**Input:**
- Primary dataset: `data/april_processing/all_delta_april29.csv`

**Outputs (Location: `play_scripts/2.models/glmmlasso/may_basic_plus/delta/`):**
- Predictions: `gl_delta_predictions_df.csv`
- Feature importance files (per-omic subset):
  - `basic_gl_delta_top_features.csv`
  - `meta_gl_delta_top_features.csv`
  - `grs_gl_delta_top_features.csv`
  - `taxa_gl_delta_top_features.csv`
  - `micom_gl_delta_top_features.csv`
  - `pathway_gl_delta_top_features.csv`
  - `metabo_gl_delta_top_features.csv`
  - `all_gl_delta_top_features.csv`

---

## 3. MERF (Mixed Effects Random Forest) Modeling

MERF models are trained using Python notebooks with hyperparameter tuning. The workflow produces predictions, R² plots, feature importance plots, and SHAP analysis for multiple model variants (MSE, Prev, PTEV, OOB).

### Supporting Scripts
- **`play_scripts/2.models/merf_python/em_utils.py`** - Core utilities and functions (called by both notebooks)
  - Defines `run_merf_analysis2()` function that writes prediction CSVs

### Input Data
Both MERF notebooks require:
- **`data/april_processing/long_april29.csv`** (for long models)
- **`data/april_processing/all_delta_april29.csv`** (for delta models)

These are produced by: `play_scripts/1.processing_omic_inputs/all_processing_long_april2.R`

### Long MERF Models

**Notebook:** `play_scripts/2.models/merf_python/april/april_long_tune_merf.ipynb`

**Output Location:** `play_scripts/2.models/merf_python/may_basic_plus/long/`

**For each omic subset** (`basic`, `meta_keep`, `only_grs`, `only_taxa`, `only_pathway`, `only_micom`, `only_metabo`, `only_all`):

- **Predictions:** `{key}_long.csv`
- **R² Plots:** 
  - `{key}_long_r2.pdf`
  - `{key}_long_r2_adj.pdf`
- **Feature Importance:** `{key}_long_ft_imp.pdf`
- **SHAP Analysis:**
  - Beeswarm plots: `shap_beeswarm_{key}_{ModelName}.png`, `shap_beeswarm_{ModelName}.png`
  - SHAP values (CSV): `shap_long_{key}_{ModelName}.csv`, `shap_long_{ModelName}.csv`
  - Log file: `shap_log.txt`

**Note:** HTML ANOVA tables and `merf_long_predictions_df.csv` are produced by the downstream R script, not this notebook.

### Delta MERF Models

**Notebook:** `play_scripts/2.models/merf_python/april/april_delta_tune_merf.ipynb`

**Output Location:** `play_scripts/2.models/merf_python/may_basic_plus/delta/`

**For each omic subset** (`basic`, `meta_keep`, `only_grs`, `only_taxa`, `only_pathway`, `only_micom`, `only_metabo`, `only_all`):

- **Predictions:** `{key}_delta.csv`
- **R² Plots:**
  - `{key}_delta_r2.pdf`
  - `{key}_delta_r2_adj.pdf`
- **Feature Importance:** `{key}_delta_ft_imp.pdf`
- **SHAP Analysis:**
  - Beeswarm plots: `shap_beeswarm_{key}_{ModelName}.png`, `shap_beeswarm_{ModelName}.png`
  - SHAP values (CSV): `shap_long_{key}_{ModelName}.csv`, `shap_long_{ModelName}.csv`
  - Log file: `shap_log.txt`

**Note:** HTML ANOVA tables and `merf_delta_predictions_df.csv` are produced by the downstream R script, not this notebook.

### Post-MERF Processing (R Scripts)

After MERF models generate individual predictions, R scripts combine them into lmer scores for statistical analysis.

#### Long Models

**Script:** `play_scripts/2.models/merf_python/april/april_for_long_merf.R`

**Inputs:**
- From `play_scripts/2.models/merf_python/may_basic_plus/long/`:
  - `basic_long.csv`
  - `meta_keep_long.csv`
  - `only_grs_long.csv`
  - `only_taxa_long.csv`
  - `only_pathway_long.csv`
  - `only_micom_long.csv`
  - `only_metabo_long.csv`
  - `only_all_long.csv`
- Test-set outcomes from: `data/april_processing/long_april29.csv`

**Output:**
- `play_scripts/2.models/merf_python/may_basic_plus/long/merf_long_predictions_df.csv`

#### Delta Models

**Script:** `play_scripts/2.models/merf_python/april/april_for_delta_merf.R`

**Inputs:**
- From `play_scripts/2.models/merf_python/may_basic_plus/`:
  - `basic_delta_april29.csv`
  - `meta_keep_delta_april29.csv`
  - `only_grs_delta_april29.csv`
  - `only_taxa_delta_april29.csv`
  - `only_pathway_delta_april29.csv`
  - `only_micom_delta_april29.csv`
  - `only_metabo_delta_april29.csv`
  - `only_all_delta_april29.csv`
- Test-set outcomes from: `data/april_processing/all_delta_april29.csv`

**Output:**
- `play_scripts/2.models/merf_python/may_basic_plus/delta/merf_delta_predictions_df.csv`

---

## 4. Statistics and Plotting

Publication-ready figures, ANOVA tables, and combined visualizations are generated from the model outputs.

### Main Scripts

**Location:** `shiny_scripts/`

- **`rendered_word_plots.Rmd`** - Main R Markdown document for generating plots
- **`plot_input_data.R`** - Loads and processes input data (sourced by Rmd)
- **`plotting_utils.R`** - Utility functions for plotting (sourced by `plot_input_data.R`)

### Input Data Structure

#### MERF Predictions and SHAP

**Long Models (`play_scripts/2.models/merf_python/may_basic_plus/long/`):**

- **Predictions:** `merf_long_predictions_df.csv`
- **Model Inputs:**
  - `basic_long.csv`
  - `meta_keep_long.csv`
  - `only_grs_long.csv`
  - `only_taxa_long.csv`
  - `only_pathway_long.csv`
  - `only_micom_long.csv`
  - `only_metabo_long.csv`
  - `only_all_long.csv`
- **SHAP Files (MSE Model):**
  - `shap_long_basic_BEST_MSE_Model.csv`
  - `shap_long_meta_keep_MSE_Model.csv`
  - `shap_long_only_grs_MSE_Model.csv`
  - `shap_long_only_taxa_MSE_Model.csv`
  - `shap_long_only_pathway_MSE_Model.csv`
  - `shap_long_only_micom_MSE_Model.csv`
  - `shap_long_only_metabo_MSE_Model.csv`
  - `shap_long_only_all_MSE_Model.csv`

**Delta Models (`play_scripts/2.models/merf_python/may_basic_plus/delta/`):**

- **Predictions:** `merf_delta_predictions_df.csv`
- **Model Inputs:**
  - `basic_delta_april29.csv`
  - `meta_keep_delta_april29.csv`
  - `only_grs_delta_april29.csv`
  - `only_taxa_delta_april29.csv`
  - `only_pathway_delta_april29.csv`
  - `only_micom_delta_april29.csv`
  - `only_metabo_delta_april29.csv`
  - `only_all_delta_april29.csv`
- **SHAP Files (MSE Model):**
  - `shap_long_basic_MSE_Model.csv`
  - `shap_long_meta_keep_MSE_Model.csv`
  - `shap_long_only_grs_MSE_Model.csv`
  - `shap_long_only_taxa_MSE_Model.csv`
  - `shap_long_only_pathway_MSE_Model.csv`
  - `shap_long_only_micom_MSE_Model.csv`
  - `shap_long_only_metabo_MSE_Model.csv`
  - `shap_long_only_all_MSE_Model.csv`

#### GLMMLASSO Predictions and Feature Importance

**Long Models (`play_scripts/2.models/glmmlasso/may_basic_plus/long/`):**

- **Predictions:** `june_lasso_long_predictions_df.csv`
- **Feature Importance:**
  - `basic_gl_long_top_features.csv`
  - `meta_gl_long_top_features.csv`
  - `grs_gl_long_top_features.csv`
  - `taxa_gl_long_top_features.csv`
  - `micom_gl_long_top_features.csv`
  - `pathway_gl_long_top_features.csv`
  - `metabo_gl_long_top_features.csv`
  - `all_gl_long_top_features.csv`

**Delta Models (`play_scripts/2.models/glmmlasso/may_basic_plus/delta/`):**

- **Predictions:** `gl_delta_predictions_df.csv`
- **Feature Importance:**
  - `basic_gl_delta_top_features.csv`
  - `meta_gl_delta_top_features.csv`
  - `grs_gl_delta_top_features.csv`
  - `taxa_gl_delta_top_features.csv`
  - `micom_gl_delta_top_features.csv`
  - `pathway_gl_delta_top_features.csv`
  - `metabo_gl_delta_top_features.csv`
  - `all_gl_delta_top_features.csv`

### Output Structure

**Base Output Directory:** `paper_plots/`

#### Generated Outputs:

- **`paper_plots/shap/`**
  - `all_shap_plots_combined_abs.pdf`
  - `all_delta_shap_plots_combined_abs.pdf`

- **`paper_plots/anova_tables/`**
  - `anova_initial_table_with_grs.png`
  - `anova_table_clin_p.png`
  - Word documents:
    - `anova_initial_combined.docx`
    - `anova_clinical_combined.docx`
    - `anova_clinical_combined_with_increase.docx`

- **`paper_plots/aov_feats_combined/`**
  - Combined patchwork PNGs and SHAP variants

- **`paper_plots/venns/`**
  - `venn_4set_*.png`
  - `venn_grid_overall_title.png`

- **`paper_plots/lass_ft_imp/`**
  - GLMMLASSO feature importance PNGs

- **`paper_plots/ft_imp_MERF_*/`**
  - MERF feature importance PNGs (subdirectories by omic type)

- **Heatmaps:**
  - `final_heat_map_scores.png`
  - `score_means_grid.png`

---

## Workflow Summary

```
Raw Data
    ↓
[Processing Scripts] → long_april29.csv, all_delta_april29.csv
    ↓
    ├─→ [GLMMLASSO Models] → Predictions + Feature Importance
    │
    └─→ [MERF Models] → Predictions + SHAP + Feature Importance
            ↓
        [Post-MERF R Scripts] → Combined prediction files
            ↓
        [Plotting Scripts] → Publication figures
```

---

## Quick Start

1. **Process raw data:**
   ```r
   source("play_scripts/1.processing_omic_inputs/all_processing_long_april2.R")
   ```

2. **Run GLMMLASSO models:**
   - Long: `play_scripts/2.models/glmmlasso/march30_long/march30_glmmlasso_long.R`
   - Delta: `play_scripts/2.models/glmmlasso/march30_delta/glmmlasso_deltas.R`

3. **Run MERF models:**
   - Long: `play_scripts/2.models/merf_python/april/april_long_tune_merf.ipynb`
   - Delta: `play_scripts/2.models/merf_python/april/april_delta_tune_merf.ipynb`

4. **Generate combined predictions:**
   - Long: `play_scripts/2.models/merf_python/april/april_for_long_merf.R`
   - Delta: `play_scripts/2.models/merf_python/april/april_for_delta_merf.R`

5. **Generate publication plots:**
   ```r
   rmarkdown::render("shiny_scripts/rendered_word_plots.Rmd")
   ```

---

## Notes

- Ensure all upstream dependencies are completed before running downstream scripts
- Model outputs are organized by omic subset and model type for easy reference
- SHAP analysis uses MSE Model by default for consistency
- Feature importance rankings are available for both GLMMLASSO and MERF models
