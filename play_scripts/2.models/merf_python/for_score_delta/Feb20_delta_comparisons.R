library(dplyr)
library(tidyr)
library(tools)
library(broom)
library(stargazer)
library(sjPlot)

# Step 1: Set the directory path (adjust this to your directory)
long_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/zachs_rerun/drift_fs/csv/all_omic_processed_data/deltas/"

test_all <- read.csv(file.path(long_dir, 
                               'feb20_all_delta_test.csv'))
test <- test_all %>%
  dplyr::select(-c(Weight, CRP, cholesterol, ghrelin, Weight, CRP, cholesterol, ghrelin, HDL, LDL, HbA1C, 
                   insulin, leptin, peptide_yy, tgcyd, homo_ir, bmi_prs, age.y, outcome_wt_fnl, 
                   outcome_BMI_fnl)) %>% rename(age = age.x)

predict_dir_only <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/merf_plots/delta_combined/feb20"

predict_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/merf_plots/delta_combined/feb20/"

# Function to get best MERF model
get_best_model <- function(df) {
  # Step 1: Calculate average R-squared for each model
  avg_r_squared <- df %>%
    dplyr::group_by(Model) %>%
    dplyr::summarise(avg_R_squared = mean(R_squared, na.rm = TRUE))
  
  # Step 2: Find the model with the highest average R-squared
  best_model <- avg_r_squared %>%
    dplyr::filter(avg_R_squared == max(pull(avg_r_squared, avg_R_squared))) %>%
    pull(Model)
  
  # Step 3: Filter the original data for rows corresponding to the best model
  df_best_model <- df %>%
    dplyr::filter(Model %in% best_model)
  
  return(df_best_model)
}

basic <- read.csv(file.path(predict_dir, 
                            "merf_results_delta_only_basic_feb20.csv")) %>% 
  get_best_model() %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% dplyr::rename(y_new_basic = y_hat_new) %>% unique()

meta <- read.csv(file.path(predict_dir, 
                           "merf_results_delta_only_meta_feb20.csv")) %>% 
  get_best_model() %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% dplyr::rename(y_new_meta = y_hat_new) %>% unique()

# Taxa merf results 
only_taxa <- read.csv(file.path(predict_dir_only, 
                                "merf_results_delta_only_meta_taxa_feb20.csv")) %>% 
  get_best_model() %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% dplyr::rename(y_new_taxa_only = y_hat_new) %>% unique()

# Pathway merf results 
only_pathway <- read.csv(file.path(predict_dir_only, 
                                   "merf_results_delta_only_pathway_feb20.csv")) %>% 
  get_best_model() %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% dplyr::rename(y_new_path_only = y_hat_new) %>% unique()

# Micom merf results 
only_micom <- read.csv(file.path(predict_dir_only, 
                                 "merf_results_delta_only_micom_feb20.csv")) %>% 
  get_best_model() %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% dplyr::rename(y_new_micom_only = y_hat_new) %>% unique()

# Merge meta and basic 
meta_basic <- merge (basic, meta, by = c("Cluster", "Time")) %>% 
  dplyr::select(-c("Model.x", "Model.y"))

merged_meta_micom <- merge(meta_basic, only_micom, 
                               by = c("Cluster", "Time")) %>% 
  dplyr::select(-c("Model"))

merged_meta_micom_pathway <- merge(merged_meta_micom, only_pathway, 
                                   by = c("Cluster", "Time")) %>% 
  dplyr::select(-c("Model"))

merged_all <- merge(merged_meta_micom_pathway, only_taxa, 
                    by = c("Cluster", "Time")) %>% 
  dplyr::select(-c("Model"))

merged_all$Time <- as.factor(merged_all$Time)
test$range <- as.factor(test$range)
merged_only <- merge(merged_all, test, 
                by.x = c("Cluster", "Time"), 
                by.y = c("subject_id", "range"))

### REMOVE 
mod_dat = merged_only 
mod_dat$Time <- as.numeric(mod_dat$Time)

### Make linear models~ # ?lme4 ?nlme ?glm
lmer_basic <- lmer(BMI ~ y_new_basic + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_meta <- lmer(BMI ~ y_new_meta + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_micom <- lmer(BMI ~ y_new_micom_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_path <- lmer(BMI ~ y_new_path_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_tax <- lmer(BMI ~ y_new_taxa_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_time <- lmer(BMI ~ Time+ (1|Cluster), data = mod_dat, REML = FALSE)

### Single plus omic 
lmer_basic <- lmer(BMI ~ y_new_basic + Time + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_meta_b <- lmer(BMI ~ y_new_basic + y_new_meta + Time + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_micom_b <- lmer(BMI ~ y_new_basic + y_new_micom_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_path_b <- lmer(BMI ~ y_new_basic + y_new_path_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_tax_b <- lmer(BMI ~ y_new_basic + y_new_taxa_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_time_b <- lmer(BMI ~ y_new_basic + Time+ (1|Cluster), data = mod_dat, REML = FALSE)

### Combined models 
basic <- lmer(BMI ~ y_new_basic + Time+ (1|Cluster), data = mod_dat)
meta_basic <- lmer(BMI ~ y_new_basic + y_new_meta + Time + (1|Cluster), 
                   data = mod_dat)
meta_basic_tax <- lmer(BMI ~ y_new_basic + y_new_meta + 
                             y_new_taxa_only + Time + 
                             (1|Cluster), data = mod_dat)
meta_basic_tax_path <- lmer(BMI ~ y_new_basic + y_new_meta + 
                                  y_new_taxa_only + 
                                  y_new_path_only + Time + (1|Cluster), 
                                data = mod_dat)
meta_basic_tax_path_micom <- lmer(BMI ~ y_new_basic + y_new_meta + 
                                        y_new_taxa_only + 
                                        y_new_path_only + y_new_micom_only + 
                                        Time+ (1|Cluster), data = mod_dat)

# Combined models sequential 
sjPlot::tab_model(basic, meta_basic,  
                  meta_basic_tax, meta_basic_tax_path, 
                  meta_basic_tax_path_micom,
                  title = "MERF DELTA sequential lmer models",
                  string.pred = "Predictors",
                  string.est = "Estimate",
                  string.std = "std. Beta",
                  string.ci = "95% CI",
                  string.se = "std. Error",
                  p.style = c("numeric"), 
                  p.threshold = c(0.05),
                  #dv.labels = model_titles,
                  auto.label = FALSE)

# https://joshuawiley.com/MonashHonoursStatistics/LMM_Comparison.html#nested-models-in-r
# Basic vs Basic + meta
nobs(lmer_basic)
nobs(lmer_meta_b)
logLik(lmer_basic)
logLik(lmer_meta_b)
anova(lmer_basic, lmer_meta_b, test = "LRT")
anova(lmer_basic, lmer_micom_b, test = "LRT")
anova(lmer_basic, lmer_tax_b, test = "LRT")
anova(lmer_basic, lmer_path_b, test = "LRT")

merf_lmer_models <- list(
  c("lmer_basic", "lmer_meta_b"),
  c("lmer_basic", "lmer_tax_b"),
  c("lmer_basic", "lmer_micom_b"),
  c("lmer_basic", "lmer_path_b"))

library(lme4)  # Make sure the lme4 package is loaded

anova_results <- list()  # empty list to store ANOVA results

for (model_pair in merf_lmer_models) {
  model_1 <- get(model_pair[1])
  model_2 <- get(model_pair[2])
  
  # Perform ANOVA and check if both models are of class lmerMod
  if (inherits(model_1, "lmerMod") && inherits(model_2, "lmerMod")) {
    # Perform ANOVA with mixed models
    anova_result <- anova(model_1, model_2)  # This works with lmerMod objects
    tidied_result <- tidy(anova_result)  # Tidy the ANOVA result using broom::tidy()
    tidied_result$model_comparison <- paste(model_pair[1], "vs", model_pair[2])
    anova_results[[length(anova_results) + 1]] <- tidied_result
  } else {
    message("One of the models in the pair is not a valid lmerMod object.")
  }
}

# Combine all results into one dataframe
anova_table <- do.call(rbind, anova_results)
anova_table_clean <- anova_table %>%
  filter(!is.na(statistic) & !is.na(p.value)) %>%
  mutate(across(where(is.numeric), round, 3)) 

# View the combined table
print(anova_table_clean)
# Create an HTML table from the cleaned anova table
html_table <- kable(anova_table_clean, format = "html", table.attr = "class='table table-striped'")

# Save the table as an HTML file
writeLines(html_table, "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/merf_plots/feb_20_bmi_ranges/MERF_anova_table.html")
