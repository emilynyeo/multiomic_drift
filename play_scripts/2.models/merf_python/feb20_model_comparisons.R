library(dplyr)
library(tidyr)
library(tools)
library(broom)
library(stargazer)
library(sjPlot)

# Step 1: Set the directory path (adjust this to your directory)
long_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/merf_dfs/5.combined/"

test_all <- read.csv(file.path(long_dir, 
                               'feb20_test_merged_all_omics_raw_meta.csv'))
test <- test_all %>%
  dplyr::select(all_samples, outcome_BMI_fnl, time)

predict_dir_only <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/for_score_long/feb20"

predict_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/merf_plots/long_combined/feb20/"

basic <- read.csv(file.path(predict_dir, 
                           "merf_basic_only_bmi_long_feb20.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% dplyr::rename(y_new_basic = y_hat_new)

meta <- read.csv(file.path(predict_dir, 
                           "merf_meta_only_bmi_long_feb20.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% dplyr::rename(y_new_meta = y_hat_new)

only_meta <- read.csv(file.path(predict_dir_only, 
                                "merf_results_long_only_meta_feb20.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% dplyr::rename(y_new_meta_only = y_hat_new)

# Genetic merf resulst 
basic_grs <- read.csv(file.path(predict_dir, 
                                "merf_basic_grs_bmi_long_feb20.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% dplyr::rename(y_new_grs = y_hat_new)

only_grs <- read.csv(file.path(predict_dir_only, 
                                "merf_results_long_grs_feb20.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% dplyr::rename(y_new_grs_only = y_hat_new)

# Taxa merf results 
basic_taxa <- read.csv(file.path(predict_dir, 
                                 "merf_basic_taxa_bmi_long_feb20.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% dplyr::rename(y_new_taxa = y_hat_new)

only_taxa <- read.csv(file.path(predict_dir_only, 
                                   "merf_results_long_only_taxa_feb20.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% dplyr::rename(y_new_taxa_only = y_hat_new)

# Pathway merf results 
basic_pathway <- read.csv(file.path(predict_dir, 
                                    "merf_basic_pathway_bmi_long_feb20.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% dplyr::rename(y_new_pathway = y_hat_new)

only_pathway <- read.csv(file.path(predict_dir_only, 
                                 "merf_results_long_only_pathway_feb20.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% dplyr::rename(y_new_path_only = y_hat_new)

# Micom merf results 
basic_micom <- read.csv(file.path(predict_dir, 
                                  "merf_basic_micom_bmi_long_feb20.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% dplyr::rename(y_new_micom = y_hat_new)

only_micom <- read.csv(file.path(predict_dir_only, 
                                  "merf_results_long_only_micom_feb20.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% dplyr::rename(y_new_micom_only = y_hat_new)

# Merge meta and basic 
meta_basic <- merge (basic, meta, by = c("Model", "Cluster", "Time"))
merged_grs_meta <- merge(basic_grs, meta_basic, by = c("Model", "Cluster", "Time"))

merged_grs_meta_micom <- merge(merged_grs_meta, basic_micom, 
                               by = c("Model", "Cluster", "Time"))

merged_grs_meta_micom_pathway <- merge(merged_grs_meta_micom, basic_pathway, 
                                       by = c("Model", "Cluster", "Time"))

merged_all <- merge(merged_grs_meta_micom_pathway, basic_taxa, 
                    by = c("Model", "Cluster", "Time"))

merged <- merge(merged_all, test, 
                by.x = c("Cluster", "Time"), 
                by.y = c("all_samples", "time")) %>% 
  dplyr::rename(bmi = outcome_BMI_fnl)

# Merge meta and only 
meta_basic_only <- merge(basic, only_meta, by = c("Model", "Cluster", "Time"))
merged_grs_meta_only <- merge(only_grs, meta_basic_only, by = c("Model", "Cluster", "Time"))

merged_grs_meta_micom_only <- merge(merged_grs_meta_only, only_micom, 
                               by = c("Model", "Cluster", "Time"))

merged_grs_meta_micom_pathway_only <- merge(merged_grs_meta_micom_only, only_pathway, 
                                       by = c("Model", "Cluster", "Time"))

merged_only <- merge(merged_grs_meta_micom_pathway_only, only_taxa, 
                    by = c("Model", "Cluster", "Time"))

merged_only <- merge(merged_only, test, 
                by.x = c("Cluster", "Time"), 
                by.y = c("all_samples", "time")) %>% 
  dplyr::rename(bmi = outcome_BMI_fnl)

### REMOVE 
rm(merged_df, test_all, merged_grs_meta, meta_basic, 
   merged_grs_meta_micom, merged_grs_meta_micom_pathway,
   merged_grs_meta_micom_pathway_only, meta_basic_only,
   merged_grs_meta_micom_only, merged_grs_meta_only, merged_all)
rm(basic, basic_grs, basic_micom, basic_pathway, basic_taxa, only_grs,
   only_meta, only_micom, only_pathway, only_taxa, meta)

### Filter MSE
df_mse <- merged_only %>%
  dplyr::filter(Model == "MSE Model")
df_prev <- merged_only %>%
  dplyr::filter(Model == "Prev Model")
df_ptev <- merged_only %>%
  dplyr::filter(Model == "PTEV Model")
df_oob <- merged_only %>%
  dplyr::filter(Model == "OOB Model")

# Remove duplicates 
df_ptev <- distinct(df_ptev)
df_mse <- distinct(df_mse)
df_prev <- distinct(df_prev)
df_oob <- distinct(df_oob)

mod_dat = df_oob
mod_dat$Time <- as.numeric(mod_dat$Time)

### PTEV 
### Make linear models PTEV ~ # ?lme4 ?nlme ?glm
lmer_basic <- lmer(bmi ~ y_new_basic + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_meta <- lmer(bmi ~ y_new_meta_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_grs <- lmer(bmi ~ y_new_grs_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_micom <- lmer(bmi ~ y_new_micom_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_path <- lmer(bmi ~ y_new_path_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_tax <- lmer(bmi ~ y_new_taxa_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_time <- lmer(bmi ~ Time+ (1|Cluster), data = mod_dat, REML = FALSE)

### Single plus omic 
lmer_basic <- lmer(bmi ~ y_new_basic + Time + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_meta_b <- lmer(bmi ~ y_new_basic + y_new_meta_only + Time + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_grs_b <- lmer(bmi ~ y_new_basic + y_new_grs_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_micom_b <- lmer(bmi ~ y_new_basic + y_new_micom_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_path_b <- lmer(bmi ~ y_new_basic + y_new_path_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_tax_b <- lmer(bmi ~ y_new_basic + y_new_taxa_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_time_b <- lmer(bmi ~ y_new_basic + Time+ (1|Cluster), data = mod_dat, REML = FALSE)

### Combined PTEV models 
basic <- lmer(bmi ~ y_new_basic + Time+ (1|Cluster), data = mod_dat)
meta_basic <- lmer(bmi ~ y_new_basic + y_new_meta_only + Time + (1|Cluster), 
                   data = mod_dat)
meta_basic_grs <- lmer(bmi ~ y_new_basic + y_new_meta_only + y_new_grs_only + 
                           Time + (1|Cluster), data = mod_dat)
meta_basic_grs_tax <- lmer(bmi ~ y_new_basic + y_new_meta_only + 
                               y_new_grs_only + y_new_taxa_only + Time + 
                             (1|Cluster), data = mod_dat)
meta_basic_grs_tax_path <- lmer(bmi ~ y_new_basic + y_new_meta_only + 
                                    y_new_grs_only + y_new_taxa_only + 
                                    y_new_path_only + Time + (1|Cluster), 
                                    data = mod_dat)
meta_basic_grs_tax_path_micom <- lmer(bmi ~ y_new_basic + y_new_meta_only + 
                                          y_new_grs_only + y_new_taxa_only + 
                                          y_new_path_only + y_new_micom_only + 
                                          Time+ (1|Cluster), data = mod_dat)

# Combined PTEV models sequential 
sjPlot::tab_model(basic, meta_basic, meta_basic_grs, 
                  meta_basic_grs_tax, meta_basic_grs_tax_path, 
                  meta_basic_grs_tax_path_micom,
                  title = "MERF long sequential lmer models",
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
anova(lmer_basic, lmer_grs_b, test = "LRT")
anova(lmer_basic, lmer_micom_b, test = "LRT")
anova(lmer_basic, lmer_tax_b, test = "LRT")
anova(lmer_basic, lmer_path_b, test = "LRT")

