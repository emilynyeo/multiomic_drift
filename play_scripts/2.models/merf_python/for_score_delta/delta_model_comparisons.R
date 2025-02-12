library(dplyr)
library(tidyr)
library(tools)
library(broom)
library(stargazer)
library(sjPlot)

# Step 1: Set the directory path (adjust this to your directory)
delta_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/zachs_rerun/drift_fs/csv/all_omic_processed_data/deltas/"

test_all <- read.csv(file.path(delta_dir, 
                               'jan30_all_delta_test_imp_varcheck.csv'))
test <- test_all %>%
  select(subject_id, BMI, range)

predict_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/merf_plots/delta_combined/"

meta <- read.csv(file.path(predict_dir, "merf_results_delta_only_meta_feb12.csv")) %>% 
  select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% 
  dplyr::rename(y_new_meta = y_hat_new)

taxa <- read.csv(file.path(predict_dir, "merf_results_delta_only_meta_taxa_feb12.csv")) %>% 
  select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% 
  dplyr::rename(y_new_taxa = y_hat_new)

pathway <- read.csv(file.path(predict_dir, 
                              "merf_results_delta_meta_taxa_functional_feb12.csv")) %>% 
  select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% 
  dplyr::rename(y_new_pathway = y_hat_new)

micom <- read.csv(file.path(predict_dir, "merf_results_delta_only_micom_feb12.csv")) %>% 
  select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% 
  dplyr::rename(y_new_micom = y_hat_new)

# Merge the result with the next data frame (micom)
merged_meta_micom <- merge(meta, micom, 
                               by = c("Model", "Cluster", "Time"))

# Merge the result with the next data frame (pathway)
merged_meta_micom_pathway <- merge(merged_meta_micom, pathway, 
                                       by = c("Model", "Cluster", "Time"))

# Merge the result with the last data frame (taxa)
merged_all <- merge(merged_meta_micom_pathway, taxa, 
                    by = c("Model", "Cluster", "Time"))

# Merge merged_all with test by matching 'Cluster' with 'all_samples' and 'Time' with 'time'
merged <- merge(merged_all, test, 
                by.x = c("Cluster", "Time"), 
                by.y = c("subject_id", "range")) %>% 
  dplyr::rename(bmi = BMI)

rm(merged_df, merged_df_small, test_all, merged_grs_meta,
   merged_grs_meta_micom, merged_grs_meta_micom_pathway)

### Filter MSE
df_mse <- merged %>%
  dplyr::filter(Model == "MSE Model")
df_prev <- merged %>%
  dplyr::filter(Model == "Prev Model")
df_ptev <- merged %>%
  dplyr::filter(Model == "PTEV Model")
df_oob <- merged %>%
  dplyr::filter(Model == "OOB Model")

# Remove duplicates 
df_mse <- distinct(df_mse)
df_prev <- distinct(df_prev)
df_ptev <- distinct(df_ptev)
df_oob <- distinct(df_oob)

### Make linear models MSE
lm_mse_meta <- lm(bmi ~ y_new_meta, data = df_mse)
lm_mse_grs <- lm(bmi ~ y_new_grs, data = df_mse)
lm_mse_taxa <- lm(bmi ~ y_new_taxa, data = df_mse)
lm_mse_micom <- lm(bmi ~ y_new_micom, data = df_mse)
lm_mse_pathway <- lm(bmi ~ y_new_micom, data = df_mse)
### Make combined MSE models 
lm_mse_meta_grs <- lm(bmi ~ y_new_meta + y_new_grs, data = df_mse)
lm_mse_meta_grs_taxa <- lm(bmi ~ y_new_meta + y_new_grs + y_new_taxa, data = df_mse)
lm_mse_meta_grs_taxa_path <- lm(bmi ~ y_new_meta + y_new_grs + y_new_taxa +
                                  y_new_pathway, data = df_mse)
lm_mse_all_omic <- lm(bmi ~ y_new_meta + y_new_grs + y_new_taxa + 
                        y_new_pathway + y_new_micom, data = df_mse)

# Create a side-by-side comparison table of the models
out_dir <- '/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/for_score_long/'

# Single MSE models 
stargazer(lm_mse_meta, lm_mse_grs, lm_mse_taxa, lm_mse_micom, lm_mse_pathway, 
          type = "html", 
          ci=TRUE, ci.level=0.95,
          title = "Comparison of Single Omics MSE Linear Models",
          out = paste0(out_dir, "MSE_single_model_comparison.html"))

# Combined MSE models 
sjPlot::tab_model(lm_mse_meta, lm_mse_meta_grs, 
                  lm_mse_meta_grs_taxa, lm_mse_meta_grs_taxa_path, 
                  lm_mse_all_omic,
                  title = "Comparing combined models with MSE parameters",
                  string.pred = "Predictors",
                  string.est = "Estimate",
                  string.std = "std. Beta",
                  string.ci = "CI",
                  string.se = "std. Error",
                  p.style = c("numeric_stars"), 
                  p.threshold = c(0.05))

### Make linear models OOB
lm_oob_meta <- lm(bmi ~ y_new_meta, data = df_oob)
lm_oob_grs <- lm(bmi ~ y_new_grs, data = df_oob)
lm_oob_taxa <- lm(bmi ~ y_new_taxa, data = df_oob)
lm_oob_micom <- lm(bmi ~ y_new_micom, data = df_oob)
lm_oob_pathway <- lm(bmi ~ y_new_micom, data = df_oob)

### Make combined OOB models
lm_oob_meta_grs <- lm(bmi ~ y_new_meta + y_new_grs, data = df_oob)
lm_oob_meta_grs_taxa <- lm(bmi ~ y_new_meta + y_new_grs + y_new_taxa, data = df_oob)
lm_oob_meta_grs_taxa_path <- lm(bmi ~ y_new_meta + y_new_grs + y_new_taxa + y_new_pathway, data = df_oob)
lm_oob_all_omic <- lm(bmi ~ y_new_meta + y_new_grs + y_new_taxa + y_new_pathway + y_new_micom, data = df_oob)

# Single OOB models
stargazer(lm_oob_meta, lm_oob_grs, lm_oob_taxa, lm_oob_micom, lm_oob_pathway,
          type = "html",
          ci = TRUE, ci.level = 0.95,
          title = "Comparison of Single Omics OOB Linear Models",
          out = paste0(out_dir, "OOB_single_model_comparison.html"))

# Combined OOB models
sjPlot::tab_model(lm_oob_meta, lm_oob_meta_grs, 
                  lm_oob_meta_grs_taxa, lm_oob_meta_grs_taxa_path, 
                  lm_oob_all_omic,
                  title = "Comparing combined models with OOB parameters", 
                  string.pred = "Predictors",
                  string.est = "Estimate",
                  string.std = "std. Beta",
                  string.ci = "CI",
                  string.se = "std. Error",
                  p.style = c("numeric_stars"),
                  p.threshold = c(0.05))

### PREV PARAMETER COMPARISONS

### Make linear models PREV
lm_prev_meta <- lm(bmi ~ y_new_meta, data = df_prev)
lm_prev_grs <- lm(bmi ~ y_new_grs, data = df_prev)
lm_prev_taxa <- lm(bmi ~ y_new_taxa, data = df_prev)
lm_prev_micom <- lm(bmi ~ y_new_micom, data = df_prev)
lm_prev_pathway <- lm(bmi ~ y_new_micom, data = df_prev)

### Make combined PREV models
lm_prev_meta_grs <- lm(bmi ~ y_new_meta + y_new_grs, data = df_prev)
lm_prev_meta_grs_taxa <- lm(bmi ~ y_new_meta + y_new_grs + y_new_taxa, data = df_prev)
lm_prev_meta_grs_taxa_path <- lm(bmi ~ y_new_meta + y_new_grs + y_new_taxa + y_new_pathway, data = df_prev)
lm_prev_all_omic <- lm(bmi ~ y_new_meta + y_new_grs + y_new_taxa + y_new_pathway + y_new_micom, data = df_prev)

# Single PREV models
stargazer(lm_prev_meta, lm_prev_grs, lm_prev_taxa, lm_prev_micom, lm_prev_pathway,
          type = "html",
          ci = TRUE, ci.level = 0.95,
          title = "Comparison of Single Omics PREV Linear Models",
          out = paste0(out_dir, "PREV_single_model_comparison.html"))

# Combined PREV models
sjPlot::tab_model(lm_prev_meta, lm_prev_meta_grs, 
                  lm_prev_meta_grs_taxa, lm_prev_meta_grs_taxa_path, 
                  lm_prev_all_omic,
                  title = "Comparing combined models with PREV parameters",
                  string.pred = "Predictors",
                  string.est = "Estimate",
                  string.std = "std. Beta",
                  string.ci = "CI",
                  string.se = "std. Error",
                  p.style = c("numeric_stars"),
                  p.threshold = c(0.05))

### PTEV
### Make linear models PTEV
lm_ptev_meta <- lm(bmi ~ y_new_meta, data = df_ptev)
lm_ptev_grs <- lm(bmi ~ y_new_grs, data = df_ptev)
lm_ptev_taxa <- lm(bmi ~ y_new_taxa, data = df_ptev)
lm_ptev_micom <- lm(bmi ~ y_new_micom, data = df_ptev)
lm_ptev_pathway <- lm(bmi ~ y_new_micom, data = df_ptev)

### Make combined PTEV models
lm_ptev_meta_grs <- lm(bmi ~ y_new_meta + y_new_grs, data = df_ptev)
lm_ptev_meta_grs_taxa <- lm(bmi ~ y_new_meta + y_new_grs + y_new_taxa, data = df_ptev)
lm_ptev_meta_grs_taxa_path <- lm(bmi ~ y_new_meta + y_new_grs + y_new_taxa + y_new_pathway, data = df_ptev)
lm_ptev_all_omic <- lm(bmi ~ y_new_meta + y_new_grs + y_new_taxa + y_new_pathway + y_new_micom, data = df_ptev)

# Single PTEV models
stargazer(lm_ptev_meta, lm_ptev_grs, lm_ptev_taxa, lm_ptev_micom, lm_ptev_pathway,
          type = "html",
          ci = TRUE, ci.level = 0.95,
          title = "Comparison of Single Omics PTEV Linear Models",
          out = paste0(out_dir, "PTEV_single_model_comparison.html"))

# Combined PTEV models
sjPlot::tab_model(lm_ptev_meta, lm_ptev_meta_grs, 
                  lm_ptev_meta_grs_taxa, lm_ptev_meta_grs_taxa_path, 
                  lm_ptev_all_omic,
                  title = "Comparing combined models with PTEV parameters",
                  string.pred = "Predictors",
                  string.est = "Estimate",
                  string.std = "std. Beta",
                  string.ci = "CI",
                  string.se = "std. Error",
                  p.style = c("numeric_stars"),
                  p.threshold = c(0.05))


anova(lm_ptev_meta, lm_ptev_all_omic)
anova(lm_prev_meta, lm_prev_all_omic)
anova(lm_oob_meta, lm_oob_all_omic)
anova(lm_mse_meta, lm_mse_all_omic)
