## Compare 2 range GLM LASSO

# Trying out long lasso
rm(list = ls())
#source("zc_functions.R") 
#install.packages("rsq")  
library(rsq)
library(pacman)
p_load(tools, reticulate, viridis, tidyplots, patchwork, jsonlite, maps, ggvenn, 
       caret, caretEnsemble, glmnet, xgboost, ggplot2, glmmLasso, corrplot,
       readr, plyr, dplyr, tidyr, purrr, tibble, stringr, psych, randomForest,  
       reshape2, scales, gridExtra, plotly, sf, tidyverse, naniar, VIM, lme4)
library(nlme)
library(gridExtra)
library(sjPlot)
library(htmltools)
library(officer)
library(flextable)
library(webshot)
library(apaTables)
library(MuMIn)
'%ni%' <- Negate('%in%')
r2_general <-function(preds,actual){ 
  return(1- sum((preds - actual) ^ 2)/sum((actual - mean(actual))^2))
}

out_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/feb20/"
data_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/zachs_rerun/drift_fs/csv/all_omic_processed_data/deltas/"
train <- read_csv(paste0(data_dir, "feb20_all_delta_train.csv"))
test <- read_csv(paste0(data_dir, "feb20_all_delta_test.csv"))

train_slim <- train %>% dplyr::select(-c("leptin", "peptide_yy","Weight",
                                         "cholesterol","ghrelin","CRP","tgcyd",
                                         "LDL_Calculated", "HDL_Total_Direct_lipid", 
                                         "Triglyceride_lipid", "outcome_wt_fnl","outcome_BMI_fnl",
                                         "Insulin_endo", "Glucose", "HOMA_IR")) %>% rename(age = age.x)
train_slim$range <- as.factor(train_slim$range)

test_slim <- test %>% dplyr::select(-c("leptin", "peptide_yy","Weight",
                                       "cholesterol","ghrelin","CRP","tgcyd",
                                       "LDL_Calculated", "HDL_Total_Direct_lipid", 
                                       "Triglyceride_lipid", "outcome_wt_fnl","outcome_BMI_fnl",
                                       "Insulin_endo", "Glucose", "HOMA_IR")) %>% rename(age = age.x)
test_slim$range <- as.factor(test_slim$range)

# Define the column names based on your lists
basic <- c('subject_id','BMI', 'range','age', 'sex')
meta_keep <- c('subject_id','BMI', 'range', 'randomized_group', 'sex', 'race', 
               'age', 'HbA1C', 'HDL', 'homo_ir', 'insulin', 'LDL')
only_taxa <- c('subject_id','BMI', 'range', 
               grep("^g__", names(train_slim), value = TRUE))
proton_column <- which(names(train_slim) == "proton")
carbon_dioxide_column <- which(names(train_slim) == "Carbon.dioxide")
only_micom <- c('subject_id','BMI', 'range', 
                names(train_slim)[proton_column:carbon_dioxide_column])
exclude_columns <- unique(c(meta_keep, only_taxa, only_micom,
                            "leptin", "CRP", "ghrelin", "peptide_yy",
                            "cholesterol","tgcyd","Weight"))
only_pathway <- c('subject_id', 'BMI', 'range', 
                  setdiff(names(train_slim), exclude_columns))

# Create data frames based on the columns defined
train_slim$subject_id <- as.factor(train_slim$subject_id)
train_basic <- train_slim[, basic, drop = FALSE]
train_meta <- train_slim[, meta_keep, drop = FALSE]
train_taxa <- train_slim[, only_taxa, drop = FALSE]
train_micom <- train_slim[, only_micom, drop = FALSE]
train_pathway <- train_slim[, only_pathway, drop = FALSE]

test_basic <- test_slim[, basic, drop = FALSE]
test_meta <- test_slim[, meta_keep, drop = FALSE]
test_taxa <- test_slim[, only_taxa, drop = FALSE]
test_micom <- test_slim[, only_micom, drop = FALSE]
test_pathway <- test_slim[, only_pathway, drop = FALSE]

### READ IN MODELS
file_path <- file.path(out_dir, "basic_predictions_feb20.csv")
pred_df_basic <- read.csv(file_path)

file_path <- file.path(out_dir, "meta_predictions_feb20.csv")
pred_df_meta <- read.csv(file_path)

file_path <- file.path(out_dir, "taxa_predictions_feb20.csv")
pred_df_taxa <- read.csv(file_path)

file_path <- file.path(out_dir, "micom_predictions_feb20.csv")
pred_df_micom <- read.csv(file_path)

file_path <- file.path(out_dir, "path_predictions_feb20.csv")
pred_df_path <- read.csv(file_path)

##### Compare models

# Convert the last 3 columns to numeric
pred_df_basic[, tail(names(pred_df_basic), 3)] <- 
  lapply(pred_df_basic[, tail(names(pred_df_basic), 3)], as.numeric)
pred_df_meta[, tail(names(pred_df_meta), 3)] <- 
  lapply(pred_df_meta[, tail(names(pred_df_meta), 3)], as.numeric)
pred_df_taxa[, tail(names(pred_df_taxa), 3)] <- 
  lapply(pred_df_taxa[, tail(names(pred_df_taxa), 3)], as.numeric)
pred_df_micom[, tail(names(pred_df_micom), 3)] <- 
  lapply(pred_df_micom[, tail(names(pred_df_micom), 3)], as.numeric)
pred_df_path[, tail(names(pred_df_path), 3)] <- 
  lapply(pred_df_path[, tail(names(pred_df_path), 3)], as.numeric)

met_basic <- merge(pred_df_meta, pred_df_basic, 
                   by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

met_tax <- merge(met_basic, pred_df_taxa, 
                 by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

met_tax_micom <- merge(met_tax, pred_df_micom, 
                       by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

met_tax_micom_path <- merge(met_tax_micom, pred_df_path, 
                            by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

all_omic <- unique(met_tax_micom_path)
all_omic$predicted_taxa <- as.numeric(all_omic$predicted_taxa)
all_omic$predicted_micom <- as.numeric(all_omic$predicted_micom)
all_omic$predicted_pathway <- as.numeric(all_omic$predicted_pathway)
all_omic[, 3:8] <- scale(all_omic[, 3:8])

########################################################################################
mod_dat = all_omic %>% rename(bmi = actual, 
                              Time = time,
                              Cluster = subject_id,
                              y_new_meta_only = predicted_meta,
                              y_new_basic = predicted_basic,
                              y_new_micom_only = predicted_micom,
                              y_new_path_only = predicted_pathway,
                              y_new_tax_only = predicted_taxa)

mod_dat$Time <- as.numeric(mod_dat$Time)
lmer_basic <- lmer(bmi ~ y_new_basic + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_meta <- lmer(bmi ~ y_new_meta_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_micom <- lmer(bmi ~ y_new_micom_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_path <- lmer(bmi ~ y_new_path_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_tax <- lmer(bmi ~ y_new_tax_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_time <- lmer(bmi ~ Time+ (1|Cluster), data = mod_dat, REML = FALSE)

### Single plus omic 
lmer_basic <- lmer(bmi ~ y_new_basic + Time + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_meta_b <- lmer(bmi ~ y_new_basic + y_new_meta_only + Time + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_micom_b <- lmer(bmi ~ y_new_basic + y_new_micom_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_path_b <- lmer(bmi ~ y_new_basic + y_new_path_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_tax_b <- lmer(bmi ~ y_new_basic + y_new_tax_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_time_b <- lmer(bmi ~ y_new_basic + Time+ (1|Cluster), data = mod_dat, REML = FALSE)

### Combined PTEV models 
basic <- lmer(bmi ~ y_new_basic + Time+ (1|Cluster), data = mod_dat)
meta_basic <- lmer(bmi ~ y_new_basic + y_new_meta_only + Time + (1|Cluster), data = mod_dat)
meta_basic_tax <- lmer(bmi ~ y_new_basic + y_new_meta_only + 
                         y_new_tax_only + Time + (1|Cluster), data = mod_dat)
meta_basic_tax_path <- lmer(bmi ~ y_new_basic + y_new_meta_only + y_new_tax_only + 
                              y_new_path_only + Time + (1|Cluster), data = mod_dat)
meta_basic_tax_path_micom <- lmer(bmi ~ y_new_basic + y_new_meta_only + 
                                    y_new_tax_only + y_new_path_only + y_new_micom_only + 
                                    Time+ (1|Cluster), data = mod_dat)

# Combined PTEV models sequential 
sjPlot::tab_model(basic, meta_basic, 
                  meta_basic_tax, meta_basic_tax_path, 
                  meta_basic_tax_path_micom,
                  title = "glmlasso delta sequential lmer models",
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

#################################################### LONG GLM LASSO 

out_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/feb20_long/"

long_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/merf_dfs/5.combined/"

test <- read.csv(file.path(long_dir, 'feb20_test_merged_all_omics_raw_meta.csv'))
train <- read.csv(file.path(long_dir, 'feb20_training_merged_all_omics_raw_meta.csv'))

train_slim <- train %>% dplyr::select(-c("Unnamed..0_merged_data", "sample_id", 
                                         "subject_id", "record_id"))
train_slim$time <- as.factor(train_slim$time)
train_slim$all_samples <- as.factor(train_slim$all_samples)
train_slim$randomized_group <- as.numeric(train_slim$randomized_group)
train_slim$sex <- as.numeric(train_slim$sex)
train_slim$race <- as.numeric(train_slim$race)
test_slim <- test %>% 
  dplyr::select(-c("Unnamed..0_merged_data", "sample_id", "subject_id",
                   "record_id"))
test_slim$time <- as.factor(test_slim$time)
test_slim$all_samples <- as.factor(test_slim$all_samples)
test_slim$randomized_group <- as.numeric(test_slim$randomized_group)
test_slim$sex <- as.numeric(test_slim$sex)
test_slim$race <- as.numeric(test_slim$race)
train_slim <- train_slim %>% rename(subject_id = all_samples,
                                    BMI = outcome_BMI_fnl,
                                    range = time,
                                    homo_ir = HOMA_IR,
                                    insulin = Insulin_endo,
                                    LDL = LDL_Calculated,
                                    HDL = HDL_Total_Direct_lipid,
                                    HbA1C = Glucose)
test_slim <- test_slim %>% rename(subject_id = all_samples,
                                  BMI = outcome_BMI_fnl,
                                  range = time,
                                  homo_ir = HOMA_IR,
                                  insulin = Insulin_endo,
                                  LDL = LDL_Calculated,
                                  HDL = HDL_Total_Direct_lipid,
                                  HbA1C = Glucose)

# Define the column names based on your lists
basic <- c('subject_id','BMI', 'range','age', 'sex')
meta_keep <- c('subject_id','BMI', 'range', 'randomized_group', 'sex', 'race', 
               'age', 'HbA1C', 'HDL', 'homo_ir', 'insulin', 'LDL')
only_taxa <- c('subject_id','BMI', 'range', 
               grep("^g__", names(train_slim), value = TRUE))
proton_column <- which(names(train_slim) == "proton")
carbon_dioxide_column <- which(names(train_slim) == "Carbon.dioxide")
only_micom <- c('subject_id','BMI', 'range', 
                names(train_slim)[proton_column:carbon_dioxide_column])
exclude_columns <- unique(c(meta_keep, only_taxa, only_micom,
                            "leptin", "CRP", "ghrelin", "peptide_yy",
                            "cholesterol","tgcyd","Weight"))
only_pathway <- c('subject_id', 'BMI', 'range', 
                  setdiff(names(train_slim), exclude_columns))

# Create data frames based on the columns defined
train_basic <- train_slim[, basic, drop = FALSE]
train_meta <- train_slim[, meta_keep, drop = FALSE]
train_taxa <- train_slim[, only_taxa, drop = FALSE]
train_micom <- train_slim[, only_micom, drop = FALSE]
train_pathway <- train_slim[, only_pathway, drop = FALSE]

test_basic <- test_slim[, basic, drop = FALSE]
test_meta <- test_slim[, meta_keep, drop = FALSE]
test_taxa <- test_slim[, only_taxa, drop = FALSE]
test_micom <- test_slim[, only_micom, drop = FALSE]
test_pathway <- test_slim[, only_pathway, drop = FALSE]

train_basic <- unique(train_basic)
train_meta <- unique(train_meta)
train_micom <- unique(train_micom)
train_pathway <- unique(train_pathway)
train_taxa <- unique(train_taxa)

test_basic <- unique(test_basic)
test_meta <- unique(test_meta)
test_micom <- unique(test_micom)
test_pathway <- unique(test_pathway)
test_taxa <- unique(test_taxa)

subject_id_count <- train_meta %>%
  dplyr::filter(range %in% c(0, 6, 12)) %>%
  dplyr::group_by(subject_id) %>%
  dplyr::summarize(range_count = n_distinct(range))  # Count the distinct range values

# Filter `subject_id`s with fewer than 3 unique range values (0, 6, 12)
missing_subjects <- subject_id_count %>%
  dplyr::filter(range_count < 3) %>%
  pull(subject_id)

# Print the result
missing_subjects
# Remove rows from each dataframe where `subject_id` is in `missing_subjects`
train_basic <- train_basic %>%
  dplyr::filter(!subject_id %in% missing_subjects)

train_meta <- train_meta %>%
  dplyr::filter(!subject_id %in% missing_subjects)

train_micom <- train_micom %>%
  dplyr::filter(!subject_id %in% missing_subjects)

train_pathway <- train_pathway %>%
  dplyr::filter(!subject_id %in% missing_subjects)

train_taxa <- train_taxa %>%
  dplyr::filter(!subject_id %in% missing_subjects)

## Read in long lasso models
file_path <- file.path(out_dir, "long_basic_predictions_feb20.csv")
pred_df_basic <- read.csv(file_path)
file_path <- file.path(out_dir, "long_meta_predictions_feb20.csv")
pred_df_meta <- read.csv(file_path)
file_path <- file.path(out_dir, "long_taxa_predictions_feb20.csv")
pred_df_taxa <- read.csv(file_path)
file_path <- file.path(out_dir, "long_micom_predictions_feb20.csv")
pred_df_micom <- read.csv(file_path)
file_path <- file.path(out_dir, "long_path_predictions_feb20.csv")
pred_df_path <- read.csv(file_path)

##### Compare models
# Convert the first 2 columns to factors
pred_df_basic[, head(names(pred_df_basic), 2)] <- 
  lapply(pred_df_basic[, head(names(pred_df_basic), 2)], as.factor)
pred_df_meta[, head(names(pred_df_meta), 2)] <- 
  lapply(pred_df_meta[, head(names(pred_df_meta), 2)], as.factor)
pred_df_taxa[, head(names(pred_df_taxa), 2)] <- 
  lapply(pred_df_taxa[, head(names(pred_df_taxa), 2)], as.factor)
pred_df_micom[, head(names(pred_df_micom), 2)] <- 
  lapply(pred_df_micom[, head(names(pred_df_micom), 2)], as.factor)
pred_df_path[, head(names(pred_df_path), 2)] <- 
  lapply(pred_df_path[, head(names(pred_df_path), 2)], as.factor)

met_basic <- merge(pred_df_basic, pred_df_meta, by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

met_tax <- merge(met_basic, pred_df_taxa, 
                 by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

met_tax_micom <- merge(met_tax, pred_df_micom, 
                       by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

met_tax_micom_path <- merge(met_tax_micom, pred_df_path, 
                            by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

all_omic <- unique(met_tax_micom_path)

# Center and scale the last 5 columns of the dataframe
all_omic[, 3:8] <- scale(all_omic[, 3:8])
head(all_omic)

########################################################################################
mod_dat = all_omic %>% rename(bmi = actual, 
                              Time = time,
                              Cluster = subject_id,
                              y_new_meta_only = predicted_meta,
                              y_new_basic = predicted_basic,
                              y_new_micom_only = predicted_micom,
                              y_new_path_only = predicted_pathway,
                              y_new_tax_only = predicted_taxa)
### Change time to contimuous

lmer_basic <- lmer(bmi ~ y_new_basic + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_meta <- lmer(bmi ~ y_new_meta_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_micom <- lmer(bmi ~ y_new_micom_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_path <- lmer(bmi ~ y_new_path_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_tax <- lmer(bmi ~ y_new_tax_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_time <- lmer(bmi ~ Time+ (1|Cluster), data = mod_dat, REML = FALSE)

### Single plus omic 
lmer_basic <- lmer(bmi ~ y_new_basic + Time + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_meta_b <- lmer(bmi ~ y_new_basic + y_new_meta_only + Time + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_micom_b <- lmer(bmi ~ y_new_basic + y_new_micom_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_path_b <- lmer(bmi ~ y_new_basic + y_new_path_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_tax_b <- lmer(bmi ~ y_new_basic + y_new_tax_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_time_b <- lmer(bmi ~ y_new_basic + Time+ (1|Cluster), data = mod_dat, REML = FALSE)

### Combined PTEV models 
basic <- lmer(bmi ~ y_new_basic + Time+ (1|Cluster), data = mod_dat)
meta_basic <- lmer(bmi ~ y_new_basic + y_new_meta_only + Time + (1|Cluster), data = mod_dat)
meta_basic_tax <- lmer(bmi ~ y_new_basic + y_new_meta_only + 
                         y_new_tax_only + Time + (1|Cluster), data = mod_dat)
meta_basic_tax_path <- lmer(bmi ~ y_new_basic + y_new_meta_only + y_new_tax_only + 
                              y_new_path_only + Time + (1|Cluster), data = mod_dat)
meta_basic_tax_path_micom <- lmer(bmi ~ y_new_basic + y_new_meta_only + 
                                    y_new_tax_only + y_new_path_only + y_new_micom_only + 
                                    Time+ (1|Cluster), data = mod_dat)

# Combined PTEV models sequential 
sjPlot::tab_model(basic, meta_basic, 
                  meta_basic_tax, meta_basic_tax_path, 
                  meta_basic_tax_path_micom,
                  title = "glmlasso long sequential lmer models",
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
########################################################################################

# Make linear models MSE
# Combined models 
lm_meta_tax_time <- lme(actual ~ predicted_meta + predicted_taxa + time, 
                        random = ~1|subject_id, data = all_omic)
summary(lm_meta_tax_time)

### Checks before models 
library(car)
vif(lm(actual ~ predicted_meta + predicted_taxa + predicted_micom + time, data = all_omic))
summary(all_omic)
table(all_omic$subject_id)
cor(all_omic[, c("predicted_meta", "predicted_taxa", "predicted_micom")])
hist(all_omic$actual)
anyDuplicated(all_omic)

# Fit a simpler model first
lm_path <- lme(actual ~ predicted_pathway, random = ~1|subject_id, data = all_omic)
lm_micom <- lme(actual ~ predicted_micom, random = ~1|subject_id, data = all_omic)
lm_tax <- lme(actual ~ predicted_taxa, random = ~1|subject_id, data = all_omic)
lm_meta <- lme(actual ~ predicted_meta, random = ~1|subject_id, data = all_omic)

anova(lm_meta, lm_path)
anova(lm_meta, lm_micom)

model_titles <- c("Model 1: Meta", 
                  "Model 2: Taxa", 
                  "Model 3: Micom", 
                  "Model 4: Pathway")

sjPlot::tab_model(lm_meta, lm_tax, 
                  lm_micom, lm_path, 
                  title = "Comparing single omic glmLASSO models 6m, BL, 12m",
                  string.pred = "Predictors",
                  string.est = "Estimate",
                  string.std = "std. Beta",
                  string.ci = "95% CI",
                  string.se = "std. Error",
                  p.style = c("numeric"), 
                  p.threshold = c(0.05),
                  dv.labels = model_titles,
                  auto.label = FALSE)
### Combined 
lm_meta_tax <- lme(actual ~ predicted_meta + predicted_taxa, 
                   random = ~1|subject_id, data = all_omic)
summary(lm_meta_tax)
lm_meta_tax_micom <- lme(actual ~ predicted_meta + predicted_taxa + 
                           predicted_micom, random = ~1|subject_id, data = all_omic)
summary(lm_meta_tax_micom)
lm_meta_tax_micom_path <- lme(actual ~ predicted_meta + predicted_taxa + 
                                predicted_micom + predicted_pathway, 
                              random = ~1|subject_id, data = all_omic)
lm_no_meta <- lme(actual ~ predicted_taxa + 
                    predicted_micom + predicted_pathway, 
                  random = ~1|subject_id, data = all_omic)
summary(lm_meta_tax_micom_path)
anova(lm_meta, lm_meta_tax)
anova(lm_meta_tax, lm_meta_tax_micom)
anova(lm_meta_tax_micom, lm_meta_tax_micom_path)
anova(lm_meta_tax, lm_meta_tax_micom_path)
anova(lm_meta, lm_meta_tax_micom_path)

model_titles <- c("Model 1: Meta", 
                  "Model 2: Meta + Taxa", 
                  "Model 3: Meta + Taxa + Micom", 
                  "Model 4: Taxa + Micom + Pathway", 
                  "Model 5: Meta + Taxa + Micom + Pathway")

sjPlot::tab_model(lm_meta, lm_meta_tax, 
                  lm_meta_tax_micom, lm_no_meta, lm_meta_tax_micom_path, 
                  title = "Comparing combined omic models 6m, BL, 12m",
                  string.pred = "Predictors",
                  string.est = "Estimate",
                  string.std = "std. Beta",
                  string.ci = "95% CI",
                  string.se = "std. Error",
                  p.style = c("numeric"), 
                  p.threshold = c(0.05),
                  dv.labels = model_titles,
                  auto.label = FALSE)

#### PLOT AIC values
# Extract AIC and BIC for each model
aic_values <- c(AIC(lm_meta), AIC(lm_meta_tax), 
                AIC(lm_meta_tax_micom), AIC(lm_meta_tax_micom_path))
bic_values <- c(BIC(lm_meta), BIC(lm_meta_tax), 
                BIC(lm_meta_tax_micom), BIC(lm_meta_tax_micom_path))

# Create a data frame for plotting
model_names <- c("lm_meta","lm_meta_tax", "lm_meta_tax_micom", "lm_meta_tax_micom_path")
df <- data.frame(Model = model_names, AIC = aic_values, BIC = bic_values)

# Plot AIC and BIC side by side
ggplot(df, aes(x = Model)) +
  geom_bar(aes(y = AIC), stat = "identity", position = "dodge", fill = "blue", alpha = 0.7) +
  geom_bar(aes(y = BIC), stat = "identity", position = "dodge", fill = "red", alpha = 0.7) +
  ylab("Value") +
  ggtitle("Comparison of AIC and BIC for Models") +
  theme_minimal() +
  scale_y_continuous(sec.axis = sec_axis(~., name = "BIC")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))


### Check Correlations
cor_matrix <- cor((all_omic)[3:7], 
                  use = "pairwise.complete.obs", method = "pearson")
melted_cor_matrix <- melt(cor_matrix, na.rm = TRUE)
ggplot(melted_cor_matrix, 
       aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", 
                       high = "red") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 3),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 3)) + 
  labs(title = "Correlations btwn single omics & actual BMI", fill = "Correlation")
