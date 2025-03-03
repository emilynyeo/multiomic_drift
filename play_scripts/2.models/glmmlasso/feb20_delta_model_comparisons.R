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
  dplyr::select(subject_id, BMI, range)

predict_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/merf_plots/delta_combined/"
predict_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/merf_plots/delta_combined/feb20"

basic <- read.csv(file.path(predict_dir, "merf_results_delta_only_basic_feb20.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% 
  dplyr::rename(y_new_basic = y_hat_new)

meta <- read.csv(file.path(predict_dir, "merf_results_delta_only_meta_feb20.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% 
  dplyr::rename(y_new_meta = y_hat_new)

taxa <- read.csv(file.path(predict_dir, "merf_results_delta_only_meta_taxa_feb20.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% 
  dplyr::rename(y_new_taxa = y_hat_new)

pathway <- read.csv(file.path(predict_dir, 
                              "merf_results_delta_only_pathway_feb20.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% 
  dplyr::rename(y_new_pathway = y_hat_new)

micom <- read.csv(file.path(predict_dir, "merf_results_delta_only_micom_feb20.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% 
  dplyr::rename(y_new_micom = y_hat_new)

# Merge meta and only 
meta_basic_only <- merge(basic, meta, by = c("Model", "Cluster", "Time"))

merged_grs_meta_micom_only <- merge(meta_basic_only, micom, 
                                    by = c("Model", "Cluster", "Time"))

merged_grs_meta_micom_pathway_only <- merge(merged_grs_meta_micom_only, pathway, 
                                            by = c("Model", "Cluster", "Time"))

merged_only <- merge(merged_grs_meta_micom_pathway_only, taxa, 
                     by = c("Model", "Cluster", "Time"))

merged_only$Time <- ifelse(merged_only$Time == 1.9832633040858, 0, 1)
test$range <- ifelse(test$range == -0.991902700743342, 0, 1)
table(merged_only$Time)
table(test$range)

merged_only <- merge(merged_only, test, 
                     by.x = c("Cluster", "Time"), 
                     by.y = c("subject_id", "range")) #%>% 
  #dplyr::rename(bmi = BMI)

rm(merged_df, merged_all, test_all, merged_grs_meta,
   merged_meta_micom, merged_meta_micom_pathway, micom, meta, taxa, pathway)

merged <- merged_only
#merged$Time <- as.factor(merged$Time)

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

mod_dat = df_oob

### Make linear models PTEV ~ # ?lme4 ?nlme ?glm
lmer_basic <- lmer(BMI ~ y_new_basic + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_meta <- lmer(BMI ~ y_new_meta + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_micom <- lmer(BMI ~ y_new_micom + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_path <- lmer(BMI ~ y_new_pathway + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_tax <- lmer(BMI ~ y_new_taxa + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_time <- lmer(BMI ~ Time+ (1|Cluster), data = mod_dat, REML = FALSE)

### Single plus omic 
lmer_basic <- lmer(BMI ~ y_new_basic + Time + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_meta_b <- lmer(BMI ~ y_new_basic + y_new_meta + Time + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_micom_b <- lmer(BMI ~ y_new_basic + y_new_micom + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_path_b <- lmer(BMI ~ y_new_basic + y_new_pathway + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_tax_b <- lmer(BMI ~ y_new_basic + y_new_taxa + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_time_b <- lmer(BMI ~ y_new_basic + Time+ (1|Cluster), data = mod_dat, REML = FALSE)

### Combined PTEV models 
basic <- lmer(BMI ~ y_new_basic + Time+ (1|Cluster), data = mod_dat)
meta_basic <- lmer(BMI ~ y_new_basic + y_new_meta + Time + (1|Cluster), 
                   data = mod_dat)
meta_basic_tax <- lmer(BMI ~ y_new_basic + y_new_meta + y_new_taxa + Time + 
                             (1|Cluster), data = mod_dat)
meta_basic_tax_path <- lmer(BMI ~ y_new_basic + y_new_meta + 
                                  y_new_taxa + 
                                  y_new_pathway + Time + (1|Cluster), 
                                data = mod_dat)
meta_basic_tax_path_micom <- lmer(BMI ~ y_new_basic + y_new_meta + 
                                        y_new_taxa + 
                                        y_new_pathway + y_new_micom + 
                                        Time+ (1|Cluster), data = mod_dat)

# Combined PTEV models sequential 
sjPlot::tab_model(basic, meta_basic, meta_basic, 
                  meta_basic_tax, meta_basic_tax_path, 
                  meta_basic_tax_path_micom,
                  title = "MERF delta sequential lmer models",
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

##########################################################################################
### Make linear models MSE
lmer_mse_meta <- lme(bmi ~ y_new_meta + Time, 
                     random = ~1|Cluster, data = df_mse)
lmer_mse_taxa <- lme(bmi ~ y_new_taxa + Time, 
                     random = ~1|Cluster, data = df_mse)
lmer_mse_micom <- lme(bmi ~ y_new_micom + Time, 
                      random = ~1|Cluster, data = df_mse)
lmer_mse_pathway <- lme(bmi ~ y_new_pathway + Time, 
                        random = ~1|Cluster, data = df_mse)

### Combined models MSE 
mse_meta_taxa <- lme(bmi ~ y_new_meta + y_new_taxa + Time, 
                     random = ~1|Cluster, data = df_mse)
mse_meta_taxa_path <- lme(bmi ~ y_new_meta + y_new_taxa + y_new_pathway + Time, 
                          random = ~1|Cluster, data = df_mse)
mse_meta_taxa_path_micom <- lme(bmi ~ y_new_meta + y_new_taxa + y_new_pathway + 
                                  y_new_micom + Time, random = ~1|Cluster, data = df_mse)
# Combined MSE models 
single_mods <- c("Model 1: Meta", 
                 "Model 2: Taxa", 
                 "Model 3: Pathway", 
                 "Model 4: Micom")
model_titles <- c("Model 1: Meta", 
                  "Model 2: Meta + Taxa", 
                  "Model 3: Meta + Taxa + Pathway", 
                  "Model 4: Meta + Taxa + Pathway + Micom")

sjPlot::tab_model(lmer_mse_meta, lmer_mse_taxa, 
                  lmer_mse_pathway, lmer_mse_micom, 
                  title = "MERF delta single models with MSE parameters",
                  string.pred = "Predictors",
                  string.est = "Estimate",
                  string.std = "std. Beta",
                  string.ci = "95% CI",
                  string.se = "std. Error",
                  p.style = c("numeric"), 
                  p.threshold = c(0.05),
                  dv.labels = single_mods,
                  auto.label = FALSE)

sjPlot::tab_model(lmer_mse_meta, mse_meta_taxa, 
                  mse_meta_taxa_path, mse_meta_taxa_path_micom, 
                  title = "MERF delta combined models with MSE parameters",
                  string.pred = "Predictors",
                  string.est = "Estimate",
                  string.std = "std. Beta",
                  string.ci = "95% CI",
                  string.se = "std. Error",
                  p.style = c("numeric"), 
                  p.threshold = c(0.05),
                  dv.labels = model_titles,
                  auto.label = FALSE)

### Make linear models OOB
lmer_oob_meta <- lme(bmi ~ y_new_meta + Time,  
                     random = ~1|Cluster, data = df_oob)
lmer_oob_taxa <- lme(bmi ~ y_new_taxa + Time, 
                     random = ~1|Cluster, data = df_oob)
lmer_oob_micom <- lme(bmi ~ y_new_micom + Time, 
                      random = ~1|Cluster, data = df_oob)
lmer_oob_pathway <- lme(bmi ~ y_new_pathway + Time, 
                        random = ~1|Cluster, data = df_oob)

### Combined models OOB 
oob_meta_taxa <- lme(bmi ~ y_new_meta + y_new_taxa + Time, 
                     random = ~1|Cluster, data = df_oob)
oob_meta_taxa_path <- lme(bmi ~ y_new_meta + y_new_taxa + y_new_pathway + Time, 
                          random = ~1|Cluster, data = df_oob)
oob_meta_taxa_path_micom <- lme(bmi ~ y_new_meta + y_new_taxa + y_new_pathway +  
                                  y_new_micom + Time, random = ~1|Cluster, data = df_oob)

sjPlot::tab_model(lmer_oob_meta, lmer_oob_taxa, 
                  lmer_oob_pathway, lmer_oob_micom, 
                  title = "MERF delta single models with OOB parameters",
                  string.pred = "Predictors",
                  string.est = "Estimate",
                  string.std = "std. Beta",
                  string.ci = "95% CI",
                  string.se = "std. Error",
                  p.style = c("numeric"), 
                  p.threshold = c(0.05),
                  dv.labels = single_mods,
                  auto.label = FALSE)

# Combined OOB models 
sjPlot::tab_model(lmer_oob_meta, oob_meta_taxa, 
                  oob_meta_taxa_path, oob_meta_taxa_path_micom, 
                  title = "MERF delta combined models with OOB parameters",
                  string.pred = "Predictors",
                  string.est = "Estimate",
                  string.std = "std. Beta",
                  string.ci = "95% CI",
                  string.se = "std. Error",
                  p.style = c("numeric"), 
                  p.threshold = c(0.05),
                  dv.labels = model_titles,
                  auto.label = FALSE)

### Make linear models PREV
lmer_prev_meta <- lme(bmi ~ y_new_meta + Time,  
                      random = ~1|Cluster, data = df_prev)
lmer_prev_taxa <- lme(bmi ~ y_new_taxa + Time, 
                      random = ~1|Cluster, data = df_prev)
lmer_prev_micom <- lme(bmi ~ y_new_micom + Time, 
                       random = ~1|Cluster, data = df_prev)
lmer_prev_pathway <- lme(bmi ~ y_new_pathway + Time, 
                         random = ~1|Cluster, data = df_prev)

### Combined models PREV 
prev_meta_taxa <- lme(bmi ~ y_new_meta + y_new_taxa + Time, 
                      random = ~1|Cluster, data = df_prev)
prev_meta_taxa_path <- lme(bmi ~ y_new_meta + y_new_taxa + y_new_pathway + Time, 
                           random = ~1|Cluster, data = df_prev)
prev_meta_taxa_path_micom <- lme(bmi ~ y_new_meta + y_new_taxa + y_new_pathway +  
                                   y_new_micom + Time, random = ~1|Cluster, data = df_prev)

sjPlot::tab_model(lmer_prev_meta, lmer_prev_taxa, 
                  lmer_prev_pathway, lmer_prev_micom, 
                  title = "MERF delta single models with PREV parameters",
                  string.pred = "Predictors",
                  string.est = "Estimate",
                  string.std = "std. Beta",
                  string.ci = "95% CI",
                  string.se = "std. Error",
                  p.style = c("numeric"), 
                  p.threshold = c(0.05),
                  dv.labels = single_mods,
                  auto.label = FALSE)
# Combined PREV models 
sjPlot::tab_model(lmer_prev_meta, prev_meta_taxa, 
                  prev_meta_taxa_path, prev_meta_taxa_path_micom, 
                  title = "MERF delta combined models with PREV parameters",
                  string.pred = "Predictors",
                  string.est = "Estimate",
                  string.std = "std. Beta",
                  string.ci = "95% CI",
                  string.se = "std. Error",
                  p.style = c("numeric"), 
                  p.threshold = c(0.05),
                  dv.labels = model_titles,
                  auto.label = FALSE)


### Make linear models PTEV
lmer_ptev_meta <- lme(bmi ~ y_new_meta + Time,  
                      random = ~1|Cluster, data = df_ptev)
lmer_ptev_taxa <- lme(bmi ~ y_new_taxa + Time, 
                      random = ~1|Cluster, data = df_ptev)
lmer_ptev_micom <- lme(bmi ~ y_new_micom + Time, 
                       random = ~1|Cluster, data = df_ptev)
lmer_ptev_pathway <- lme(bmi ~ y_new_pathway + Time, 
                         random = ~1|Cluster, data = df_ptev)

### Combined models PTEV 
ptev_meta_taxa <- lme(bmi ~ y_new_meta + y_new_taxa + Time, 
                      random = ~1|Cluster, data = df_ptev)
ptev_meta_taxa_path <- lme(bmi ~ y_new_meta + y_new_taxa + y_new_pathway + Time, 
                           random = ~1|Cluster, data = df_ptev)
ptev_meta_taxa_path_micom <- lme(bmi ~ y_new_meta + y_new_taxa + y_new_pathway +  
                                   y_new_micom + Time, random = ~1|Cluster, data = df_ptev)

sjPlot::tab_model(lmer_ptev_meta, lmer_ptev_taxa, 
                  lmer_ptev_pathway, lmer_ptev_micom, 
                  title = "MERF delta single models with PEV parameters",
                  string.pred = "Predictors",
                  string.est = "Estimate",
                  string.std = "std. Beta",
                  string.ci = "95% CI",
                  string.se = "std. Error",
                  p.style = c("numeric"), 
                  p.threshold = c(0.05),
                  dv.labels = single_mods,
                  auto.label = FALSE)
# Combined PTEV models 
sjPlot::tab_model(lmer_ptev_meta, ptev_meta_taxa, 
                  ptev_meta_taxa_path, ptev_meta_taxa_path_micom, 
                  title = "MERF delta combined models with PTEV parameters",
                  string.pred = "Predictors",
                  string.est = "Estimate",
                  string.std = "std. Beta",
                  string.ci = "95% CI",
                  string.se = "std. Error",
                  p.style = c("numeric"), 
                  p.threshold = c(0.05),
                  dv.labels = model_titles,
                  auto.label = FALSE)

anova(lmer_ptev_meta, ptev_meta_taxa_path_micom)
anova(lmer_ptev_meta, ptev_meta_taxa)
anova(lmer_prev_meta, prev_meta_taxa_path_micom)
anova(lmer_prev_meta, prev_meta_taxa)



