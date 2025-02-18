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

taxa <- read.csv(file.path(predict_dir, "merf_results_delta_only_taxa_feb12.csv")) %>% 
  select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% 
  dplyr::rename(y_new_taxa = y_hat_new)

pathway <- read.csv(file.path(predict_dir, 
                              "merf_results_delta_only_pathway_feb12.csv")) %>% 
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

rm(merged_df, merged_all, test_all, merged_grs_meta,
   merged_meta_micom, merged_meta_micom_pathway, micom, meta, taxa, pathway)

merged$Time <- as.factor(merged$Time)
merged$Time <- ifelse(merged$Time == -0.991902700743342, 0, 1)

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



