# Unused Functions in em_utils.py

## Analysis Date
2025-10-30

## Functions Currently NOT Used by Main MERF Scripts

The main MERF scripts (`april_long_tune_merf.ipynb` and `april_delta_tune_merf.ipynb`, now archived) primarily use:
- `read_data()` - Data loading
- `run_merf_analysis2()` - Main MERF analysis function

### Completely Unused Functions

These functions are **not called anywhere** (neither externally nor internally):

1. **`make_long(wide_data)`** (lines 24-62)
   - Converts wide-format DataFrame to long format
   - Appears to be for data transformation but not used

2. **`create_t_column(df)`** (lines 65-73)
   - Maps timepoints to numeric values
   - Not used in current workflow

3. **`plot_predicted_vs_actual_old(...)`** (lines 76-109)
   - Old version of predicted vs actual plotting
   - Replaced by `plot_predicted_vs_actual()` but not used either

4. **`plot_predicted_vs_actual(...)`** (lines 111-152)
   - Plots predicted vs actual values with metrics
   - Not called by main scripts

5. **`calculate_metrics(Y_true, Y_pred)`** (lines 156-162)
   - Calculates RMSE and R-squared
   - Standalone function not used

6. **`plot_feature_importances(...)`** (lines 165-175)
   - Plots feature importances (all features)
   - Not used (feature plotting is done within `run_merf_analysis2`)

7. **`plot_top_20_feature_importances(...)`** (lines 178-189)
   - Plots top 20 feature importances with customizable color
   - Not used

8. **`create_parameter_grids(df)`** (lines 214-290)
   - Creates parameter grids from tuning results
   - Logic duplicated inline in notebooks; function not called

9. **`run_merf_analysis_old2(...)`** (lines 499-662)
   - Old version of MERF analysis function
   - Replaced by `run_merf_analysis2()` but kept for reference

10. **`compare_r2_values1(...)`** (lines 977-1008)
    - Compares R² values across models (version 1)
    - Not used

11. **`compare_r2_values2(...)`** (lines 1012-1036)
    - Compares R² values across models (version 2)
    - Not used

12. **`compare_r2_values3(...)`** (lines 1040-1067)
    - Compares R² values across models (version 3)
    - Not used

13. **`run_merf(...)`** (lines 1103-1145)
    - Simplified MERF runner function
    - Not used by main scripts

### Functions with Minimal/Legacy Usage

These functions may be called but are not part of the main workflow:

14. **`run_merf_analysis(...)`** (lines 301-494)
    - Earlier version of MERF analysis
    - Replaced by `run_merf_analysis2()` but may still be in some notebooks

15. **`plot_top_feature_importances_comparative(...)`** (lines 1073-1096, 1290-1328)
    - Two versions of this function (duplicate!)
    - One takes feature_importance_dfs, one takes results_list
    - Not used by main scripts

16. **`plot_r2_comparison(results_list)`** (lines 1148-1187)
    - Plots R² comparison (legacy version)
    - Not used by main scripts

17. **`plot_feature_importances_comparison(results_list)`** (lines 1190-1228)
    - Compares feature importances across models
    - Not used by main scripts

18. **`plot_r2_values(results_list, model_names)`** (lines 1238-1283)
    - Plots R² values with color coding
    - Not used by main scripts

---

## Functions Currently USED

These functions are actively used:

1. **`read_data(directory, filename)`** - Used in april notebooks to load CSV files
2. **`manual_r2_score(y_true, y_pred)`** - Used internally by `run_merf_analysis()` and `run_merf_analysis2()`
3. **`plot_shap_beeswarm(...)`** - Used internally by `run_merf_analysis2()` for SHAP visualization
4. **`run_merf_analysis2(...)`** - **Main function** used by april notebooks for MERF analysis

---

## Recommendations

1. **Safe to Remove:** Functions #1-13 (completely unused)
2. **Consider Removing:** Functions #14-18 (legacy/alternative implementations)
3. **Keep:** Functions that are actively used (#1, #2, #3, #4 above)

### Note on Duplicate Function
`plot_top_feature_importances_comparative()` is defined twice (lines 1073 and 1290) with different signatures. Only one version should be kept if either is needed.

---

## Function Statistics

- **Total functions:** 23
- **Actively used:** 4 (`read_data`, `manual_r2_score`, `plot_shap_beeswarm`, `run_merf_analysis2`)
- **Unused:** 13-18 (depending on how strictly you define "unused")
- **Potential cleanup:** ~60-75% of functions could be removed if not needed elsewhere

