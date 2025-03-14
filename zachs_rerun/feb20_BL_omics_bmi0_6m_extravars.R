#' @author Emily Yeo
#' @email emily.yeo@colorado.edu
#' @purpose BL omic predictors of BL to 6m BMI change 
#' @lab Stanislawski Lab
#' @affiliation University of Colorado Denver - Anschutz Medical, Dept of Biomedical Informatics & Personalized Medicine

###############################
###     Reading R data files    
###############################

# In[1]: Imports ----
rm(list = ls())
source("zc_functions.R") 
library(pacman)
p_load(tools, reticulate, viridis, tidyplots, patchwork, jsonlite, maps, ggvenn, 
       caret, caretEnsemble, readr, plyr, dplyr, tidyr, purrr, tibble, stringr, 
       psych, randomForest, glmnet, xgboost, ggplot2, reshape2, scales, gridExtra, 
       plotly, sf, tidyverse)

###############################
###     Caret Analysis
###############################
# Define peptide variable at the top of the script
extra_var <- "Peptide_YY"  # Change this if needed later

# In[2] Load Datasets ----
data_dir <- "drift_fs/csv/all_omic_processed_data/"

long_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/merf_dfs/5.combined/"

#test <- read.csv(file.path(long_dir, 'feb20_test_merged_all_omics_raw_meta.csv'))
#train <- read.csv(file.path(long_dir, 'feb20_training_merged_all_omics_raw_meta.csv'))
test <- read.csv(file.path(long_dir, 'feb20_test_merged_all_omics_extra_meta.csv'))
train <- read.csv(file.path(long_dir, 'feb20_training_merged_all_omics_extra_meta.csv'))

# FIRST JUST WITH OUTER JOINED
omic_g_ra <- train

### Make BL, 6m and 12m dfs
BL <- omic_g_ra %>%
  dplyr::filter(grepl("0$", time)) %>%
  dplyr::select(-matches("12$|6$"))

m6 <- omic_g_ra %>%
  dplyr::filter(grepl("6$", time)) %>%
  dplyr::select(-matches("0$|12$"))

m12 <- omic_g_ra %>%
  dplyr::filter(grepl("12$", time)) %>%
  dplyr::select(-matches("0$|6$"))

m6_slim <- m6 %>% dplyr::select(c("subject_id", "outcome_BMI_fnl", extra_var)) %>%
  dplyr::rename(outcome_BMI_fnl_6m = outcome_BMI_fnl)

m12_slim <- m12 %>% dplyr::select(c("subject_id", "outcome_BMI_fnl", extra_var)) %>%
  dplyr::rename(outcome_BMI_fnl_12m = outcome_BMI_fnl)

m12_6m <- merge(m12_slim, m6_slim, by = "subject_id", all = TRUE)
all <- merge(m12_6m, BL, by = "subject_id", all = TRUE)
all$bmi_bL_6m <- all$outcome_BMI_fnl_6m - all$outcome_BMI_fnl
all <- all[!is.na(all$bmi_bL_6m), ] # Remove rows with missing values in bmi_bL_6m (30)

### Latent variables
latent_variables_BL <- c(
  "randomized_group", "score_std", "cohort_number", "sex", "race", "age", 
  "Glucose_BL", "HOMA_IR_BL", "Insulin_endo_BL", "HDL_Total_Direct_lipid_BL",             
  "LDL_Calculated_BL", "Triglyceride_lipid_BL", "outcome_BMI_fnl_BL")

# Get all column names in the dataframe "all"
all_columns <- colnames(all)
columns_to_exclude <- c("subject_id", "sample_id", "all_samples", "record_id") # Exclude list
all_latent <- setdiff(all_columns, columns_to_exclude)

### Process DFs
# imputed <- preprocess_data(all, all_latent, "medianImpute")
# remove outcome_BMI_fnl_BL from the genus and species dataframes
imputed <- remove_columns(all, 
                          columns_to_remove = c("sample_id", "all_samples_y",
                                                "record_id", "all_samples",  
                                                "time_y", "outcome_BMI_fnl_12m", 
                                                "outcome_BMI_fnl_6m", "Unnamed..0_merged_data"))
## Set single omic DFs

# Find the column indices for "proton" and "Carbon.dioxide" in train_set
proton_column <- which(colnames(imputed) == "proton")
carbon_dioxide_column <- which(colnames(imputed) == "Carbon.dioxide")
n10_col <- which(colnames(imputed) == "N10.formyl.tetrahydrofolate_biosynthesis")
lval_col <- which(colnames(imputed) == "L.valine_biosynthesis")

# Columns to KEEP for only meta 
only_basic <- c('subject_id', 'bmi_bL_6m', 'outcome_BMI_fnl', 'time', 'age',  
                'sex', extra_var) # Use extra_var here
meta_keep <- c('subject_id', 'bmi_bL_6m',  'outcome_BMI_fnl', 'time',  
               'randomized_group', 'cohort_number', 'sex', 'race', 'age',  
               'Glucose', 'HDL_Total_Direct_lipid', 'HOMA_IR',  
               'Insulin_endo', 'LDL_Calculated', 'Triglyceride_lipid', 
               extra_var) # Use extra_var here

only_grs <- c('subject_id', 'bmi_bL_6m',  'outcome_BMI_fnl', 'bmi_prs', 'time', 
              extra_var) # Use extra_var here
only_taxa <- c('subject_id', 'bmi_bL_6m',  'outcome_BMI_fnl', 'time', 
               extra_var, 
               grep("^g__", colnames(imputed), value = TRUE)) # Use extra_var here
only_micom <- c('subject_id', 'bmi_bL_6m', 'outcome_BMI_fnl', 'time', 
                extra_var, 
                colnames(imputed)[proton_column:carbon_dioxide_column]) # Use extra_var here
only_path <- c('subject_id', 'bmi_bL_6m', 'outcome_BMI_fnl', 'time', 
               extra_var, 
               colnames(imputed)[n10_col:lval_col]) # Use extra_var here

basic <- imputed %>% select(all_of(only_basic))
meta <- imputed %>% select(all_of(meta_keep))
grs <- imputed %>% select(all_of(only_grs))
taxa <- imputed %>% select(all_of(only_taxa))
path <- imputed %>% select(all_of(only_path))
micom <- imputed %>% select(all_of(only_micom))

rm(test_BL, test_m12, test_m12_6m, test_m12_slim, test_m6, 
   test_m6_slim, test_omic_g_ra, test, omic_g_ra, m6_slim, m6, m12, m12_slim, 
   m12_6m, BL, all)

############ RUN on Single Omics ##############################################
datasets <- list(basic, meta, grs, taxa, path, micom)
target_vars <- c(extra_var,extra_var,extra_var,extra_var,extra_var,extra_var)
result_prefixes <- c("basic_BL_bmi0_6m", "meta_BL_bmi0_6m", "grs_BL_bmi0_6m", 
                     "taxa_BL_bmi0_6m", "path_BL_bmi0_6m", "micom_BL_bmi0_6m")
# Append peptide_variable to each result prefix
result_prefixes_with_peptide <- paste0(result_prefixes, "_", extra_var)
set.seed(123)
train_control <- trainControl(method = "cv", number = 5, search = "grid")
train_and_save_multiple_models(datasets, target_vars, train_control, result_prefixes)

# Define base path
base_path <- "drift_fs/csv/results/feb20"

# Define file paths using paste0
file_paths <- list(
  # basic
  basic_beta = paste0("basic_BL_bmi0_6m_", extra_var, "_beta.csv"),
  basic_ft_imp = paste0("basic_BL_bmi0_6m_", extra_var, "_feature_importance.csv"),
  basic_metrics = paste0("basic_BL_bmi0_6m_", extra_var, "_metrics.csv"),
  # grs
  grs_beta = paste0("grs_BL_bmi0_6m_", extra_var, "_beta.csv"),
  grs_ft_imp = paste0("grs_BL_bmi0_6m_", extra_var, "_feature_importance.csv"),
  grs_metrics = paste0("grs_BL_bmi0_6m_", extra_var, "_metrics.csv"),
  # meta
  meta_beta = paste0("meta_BL_bmi0_6m_", extra_var, "_beta.csv"),
  meta_ft_imp = paste0("meta_BL_bmi0_6m_", extra_var, "_feature_importance.csv"),
  meta_metrics = paste0("meta_BL_bmi0_6m_", extra_var, "_metrics.csv"),
  # taxa
  taxa_beta = paste0("taxa_BL_bmi0_6m_", extra_var, "_beta.csv"),
  taxa_ft_imp = paste0("taxa_BL_bmi0_6m_", extra_var, "_feature_importance.csv"),
  taxa_metrics = paste0("taxa_BL_bmi0_6m_", extra_var, "_metrics.csv"),
  # micom
  micom_beta = paste0("micom_BL_bmi0_6m_", extra_var, "_beta.csv"),
  micom_ft_imp = paste0("micom_BL_bmi0_6m_", extra_var, "_feature_importance.csv"),
  micom_metrics = paste0("micom_BL_bmi0_6m_", extra_var, "_metrics.csv"),
  # pathway
  path_beta = paste0("path_BL_bmi0_6m_", extra_var, "_beta.csv"),
  path_ft_imp = paste0("path_BL_bmi0_6m_", extra_var, "_feature_importance.csv"),
  path_metrics = paste0("path_BL_bmi0_6m_", extra_var, "_metrics.csv")
)


# Read all data into a named list using lapply
data_list <- lapply(file_paths, function(path) read.csv(file.path(base_path, path)))
names(data_list) <- names(file_paths) # Assign names to the data list
# Process and plot for all datasets ----
datasets <- list(
  "basic" = list(
    beta = data_list$basic_beta,
    ft_imp = data_list$basic_ft_imp,
    metrics = data_list$basic_metrics
  ),
  "grs" = list(
    beta = data_list$grs_beta,
    ft_imp = data_list$grs_ft_imp,
    metrics = data_list$grs_metrics
  ),
  "meta" = list(
    beta = data_list$meta_beta,
    ft_imp = data_list$meta_ft_imp,
    metrics = data_list$meta_metrics
  ),
  "taxa" = list(
    beta = data_list$taxa_beta,
    ft_imp = data_list$taxa_ft_imp,
    metrics = data_list$taxa_metrics
  ),
  "micom" = list(
    beta = data_list$micom_beta,
    ft_imp = data_list$micom_ft_imp,
    metrics = data_list$micom_metrics
  ),
  "path" = list(
    beta = data_list$path_beta,
    ft_imp = data_list$path_ft_imp,
    metrics = data_list$path_metrics
  )
)

######## PLOT R2 OF MODELS #####################################################
r2_data <- list()
for (dataset_name in names(datasets)) { # Iterate through each dataset
  metrics_df <- datasets[[dataset_name]]$metrics # Get current dataset metrics
  # Filter for Train data and extract model names and R² values
  train_r2 <- metrics_df %>%
    dplyr::filter(DataType == "Train") %>%
    dplyr::select(Model, R2) %>%
    mutate(Dataset = dataset_name) # Add dataset name as a new column
  r2_data[[dataset_name]] <- train_r2 # Append the results to the r2_data list
}
r2_combined <- bind_rows(r2_data) # Combine R² single data frame

# If you only want RF and LASSO models 
r2_combined_rf_lasso <- r2_combined %>%
  dplyr::filter(str_detect(Model, "lasso_model|rf_model"))

custom_colors <- c("basic" = "#85c787",  # Example hex color for lasso_model
                   "grs" = "#ff7f5e",
                   "meta" = "#7a8ec4",
                   "micom" = "#d4ba66",
                   "path" = "#66a4d4", 
                   "taxa" = "#d4668e")
# Plot R² values using ggplot2
plot <- ggplot(r2_combined_rf_lasso, aes(x = Model, y = R2, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(
    title = "R² Values for Train Models Across Datasets",
    x = "Model",
    y = "R²",
    fill = "Dataset") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = custom_colors)
print(plot)

ggsave("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/zachs_rerun/drift_fs/csv/results/feb20/single_omic_BL_bmi0_6m_r2_plot_rf_lass.png", 
       plot = plot, width = 10, height = 8)

######## PLOT FT. IMPORTANCES #################################################

# Create a list to store the feature importance data
ft_imp_data <- list()

for (dataset_name in names(datasets)) { # Iterate through each dataset
  ft_imp_df <- datasets[[dataset_name]]$ft_imp # Get feature importance data for the current dataset
  
  # Reshape the feature importance data
  ft_imp_long <- ft_imp_df %>%
    gather(Model, Importance, -Variable) %>%  # Convert from wide to long format
    mutate(Dataset = dataset_name)  # Add dataset name as a new column
  
  # Sort by importance and select the top 10 features for each dataset and model
  ft_imp_top10 <- ft_imp_long %>%
    group_by(Dataset, Model) %>%
    top_n(10, Importance) %>%  # Select the top 10 features by importance
    ungroup() %>%
    arrange(Dataset, Model, desc(Importance))  # Arrange by importance
  
  ft_imp_data[[dataset_name]] <- ft_imp_top10 # Append the results to the list
}

# Combine the feature importance data for all datasets
ft_imp_combined <- bind_rows(ft_imp_data)

# Shorten the feature names to a maximum of 8 characters
ft_imp_combined$Variable <- substr(ft_imp_combined$Variable, 1, 8)

# Keep the top 10 features for each model and dataset
ft_imp_combined_top10 <- ft_imp_combined %>%
  dplyr::filter(Importance > 0) %>%  # Remove features with importance of 0
  dplyr::group_by(Dataset, Model) %>%  # Group by dataset and model
  top_n(10, Importance) %>%  # Select the top 8 features for each model
  ungroup()  # Ungroup to return to a normal dataframe

# View the updated dataframe
head(ft_imp_combined_top10)
ft_imp_combined_rf_lasso <- ft_imp_combined_top10 %>%
  dplyr::filter(str_detect(Model, "Lasso_Importance|RF_Importance"))

# Define custom colors for each model
custom_colors <- c("basic" = "#85c787",  # Example hex color for lasso_model
                   "grs" = "#ff7f5e",
                   "meta" = "#7a8ec4",
                   "micom" = "#d4ba66",
                   "path" = "#66a4d4", 
                   "tax" = "#d4668e")   # Example hex color for rf_model

# Plot feature importances using ggplot2 with independent y-axes for each facet
plot_ft_imp <- ggplot(ft_imp_combined_rf_lasso, 
                      aes(x = reorder(Variable, Importance), y = Importance, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Dataset, scales = "free_y") +  # Facet by dataset with independent y-axes
  theme_minimal() +
  labs(
    title = "Top 10 Feature Importances for Each Model Across Datasets",
    x = "Feature",
    y = "Importance",
    fill = "Model"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate feature names
  coord_flip() +  # Flip the coordinates for better readability
  scale_fill_manual(values = custom_colors) # Apply custom colors correctly

# Print the plot
print(plot_ft_imp)

# Save plot to a specific directory with a custom filename
ggsave("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/zachs_rerun/drift_fs/csv/results/feb20/single_omic_BL_bmi0_6m_ft_imp_rf_lasso.png", 
       plot = plot_ft_imp, width = 10, height = 8)

###############################################################################
datasets <- list(basic, meta, grs, taxa, path, micom)
# target variables
target_vars <- c("bmi_bL_6m", "bmi_bL_6m", "bmi_bL_6m", 
                 "bmi_bL_6m", "bmi_bL_6m", "bmi_bL_6m")
# result prefixes
result_prefixes <- c("basic_BL_bmi0_6m", "meta_BL_bmi0_6m", "grs_BL_bmi0_6m", 
                     "taxa_BL_bmi0_6m", "path_BL_bmi0_6m", "micom_BL_bmi0_6m")

results <- train_multiple_models(datasets, target_vars, train_control, result_prefixes)

## LASSO results and PREDICTIONS
# basic results
lasso_mod_basic <- results[["basic_BL_bmi0_6m"]]$lasso_model
lasso_pred_basic <- best_lasso_predictions(lasso_mod_basic, 
                                           basic, "bmi_bL_6m", test_basic) %>% 
  rename(predicted_basic = predicted) %>% unique()
# meta
lasso_mod_meta <- results[["meta_BL_bmi0_6m"]]$lasso_model
lasso_pred_meta <- best_lasso_predictions(lasso_mod_meta, 
                                          meta, "bmi_bL_6m", test_meta) %>% 
  rename(predicted_meta = predicted) %>% unique()
# grs
lasso_mod_grs <- results[["grs_BL_bmi0_6m"]]$lasso_model
lasso_pred_grs <- best_lasso_predictions(lasso_mod_grs, 
                                         grs, "bmi_bL_6m", test_grs) %>% 
  rename(predicted_grs = predicted) %>% unique()
# taxa
lasso_mod_taxa <- results[["taxa_BL_bmi0_6m"]]$lasso_model
lasso_pred_taxa <- best_lasso_predictions(lasso_mod_taxa, 
                                          taxa, "bmi_bL_6m", test_taxa) %>% 
  rename(predicted_taxa = predicted) %>% unique()
# micom 
lasso_mod_micom <- results[["micom_BL_bmi0_6m"]]$lasso_model
lasso_pred_micom <- best_lasso_predictions(lasso_mod_micom, 
                                           micom, "bmi_bL_6m", test_micom) %>% 
  rename(predicted_micom = predicted) %>% unique()
# path 
lasso_mod_path <- results[["path_BL_bmi0_6m"]]$lasso_model
lasso_pred_path <- best_lasso_predictions(lasso_mod_path, 
                                          path, "bmi_bL_6m", 
                                          test_path, s = 0.07) %>% 
  rename(predicted_path = predicted) %>% unique()

merged_lasso_preds <- lasso_pred_basic %>%
  left_join(lasso_pred_meta, by = "subject_id") %>%
  left_join(lasso_pred_grs, by = "subject_id") %>%
  left_join(lasso_pred_taxa, by = "subject_id") %>%
  left_join(lasso_pred_micom, by = "subject_id") %>%
  left_join(lasso_pred_path, by = "subject_id") %>% 
  dplyr::select(-c("actual.y.y.y", "time.y.y.y", "actual.x.x.x", "time.x.x.x",
                   "actual.y.y", "time.y.y", "actual.x.x", "time.x.x", "time.y",
                   "actual.y")) %>% rename(time = time.x, actual = actual.x) %>%
  mutate_at(vars(2:9), as.numeric)

### single lasso LM models 
basic_lm <- lm(actual ~ predicted_basic, merged_lasso_preds)
meta_lm <- lm(actual~ predicted_meta, merged_lasso_preds)
grs_lm <- lm(actual ~ predicted_grs, merged_lasso_preds)
taxa_lm <- lm(actual ~ predicted_taxa, merged_lasso_preds)
micom_lm <- lm(actual ~ predicted_micom, merged_lasso_preds)
path_lm <- lm(actual ~ predicted_path, merged_lasso_preds)
anova(basic_lm, meta_lm)
anova(basic_lm, grs_lm)
anova(basic_lm, taxa_lm)
anova(basic_lm, micom_lm)
anova(basic_lm, path_lm)

### combined lasso LM models
basic_lm <- lm(actual ~ predicted_basic, merged_lasso_preds)
basic_meta_lm <- lm(actual~ predicted_basic + predicted_meta, merged_lasso_preds)
basic_grs_lm <- lm(actual ~ predicted_basic + predicted_grs, merged_lasso_preds)
basic_taxa_lm <- lm(actual ~ predicted_basic + predicted_taxa, merged_lasso_preds)
basic_micom_lm <- lm(actual ~ predicted_basic + predicted_micom, merged_lasso_preds)
basic_path_lm <- lm(actual ~ predicted_basic + predicted_path, merged_lasso_preds)
anova(basic_lm, basic_meta_lm)
anova(basic_lm, basic_grs_lm)
anova(basic_lm, basic_taxa_lm)
anova(basic_lm, basic_micom_lm)
anova(basic_lm, basic_path_lm)

### RF results and PREDICTIONS
rf_mod_basic <- results[["basic_BL_bmi0_6m"]]$rf_model
pred_basic_rf <- best_rf_predictions(rf_mod_basic, basic, "bmi_bL_6m", test_basic) %>% 
  rename(pred_basic_rf = predicted) %>% unique()

rf_mod_meta <- results[["meta_BL_bmi0_6m"]]$rf_model
pred_meta_rf <- best_rf_predictions(rf_mod_meta, meta, "bmi_bL_6m", test_meta) %>% 
  rename(pred_meta_rf = predicted) %>% unique()

rf_mod_grs <- results[["grs_BL_bmi0_6m"]]$rf_model
pred_grs_rf <- best_rf_predictions(rf_mod_grs, grs, "bmi_bL_6m", test_grs) %>% 
  rename(pred_grs_rf = predicted) %>% unique()

rf_mod_taxa <- results[["taxa_BL_bmi0_6m"]]$rf_model
pred_taxa_rf <- best_rf_predictions(rf_mod_taxa, taxa, "bmi_bL_6m", test_taxa) %>% 
  rename(pred_taxa_rf = predicted) %>% unique()

rf_mod_micom <- results[["micom_BL_bmi0_6m"]]$rf_model
pred_micom_rf <- best_rf_predictions(rf_mod_micom, micom, "bmi_bL_6m", test_micom) %>% 
  rename(pred_micom_rf = predicted) %>% unique()

rf_mod_path <- results[["path_BL_bmi0_6m"]]$rf_model
pred_path_rf <- best_rf_predictions(rf_mod_path, path, "bmi_bL_6m", test_path) %>% 
  rename(pred_path_rf = predicted) %>% unique()

## Combine predictions
merged_rf_preds <- pred_basic_rf %>%
  left_join(pred_meta_rf, by = "subject_id") %>%
  left_join(pred_grs_rf, by = "subject_id") %>%
  left_join(pred_taxa_rf, by = "subject_id") %>%
  left_join(pred_micom_rf, by = "subject_id") %>%
  left_join(pred_path_rf, by = "subject_id") %>% 
  dplyr::select(-c("actual.y.y.y", "time.y.y.y", "actual.x.x.x", "time.x.x.x",
                   "actual.y.y", "time.y.y", "actual.x.x", "time.x.x", "time.y",
                   "actual.y")) %>% rename(time = time.x, actual = actual.x) %>%
  mutate_at(vars(2:9), as.numeric)

### Get ft imp
models <- list(rf_mod_basic, rf_mod_meta, rf_mod_grs, 
               rf_mod_taxa, rf_mod_micom, rf_mod_path)
labels <- c("rf_basic", "rf_meta", "rf_grs", "rf_taxa", "rf_micom", "rf_path") 
rf_importances <- combine_importances(models, labels) # Combine importances 

### single rf LM models 
basic_lm_rf <- lm(actual ~ pred_basic_rf, merged_rf_preds)
meta_lm_rf <- lm(actual~ pred_meta_rf, merged_rf_preds)
grs_lm_rf <- lm(actual ~ pred_grs_rf, merged_rf_preds)
taxa_lm_rf <- lm(actual ~ pred_grs_rf, merged_rf_preds)
micom_lm_rf <- lm(actual ~ pred_micom_rf, merged_rf_preds)
path_lm_rf <- lm(actual ~ pred_path_rf, merged_rf_preds)
anova(basic_lm_rf, meta_lm_rf)
anova(basic_lm_rf, grs_lm_rf)
anova(basic_lm_rf, taxa_lm_rf)
anova(basic_lm_rf, micom_lm_rf)
anova(basic_lm_rf, path_lm_rf)

### combined lasso LM models
basic_lm_rf <- lm(actual ~ pred_basic_rf, merged_rf_preds)
basic_meta_lm_rf <- lm(actual~ pred_basic_rf + pred_meta_rf, merged_rf_preds)
basic_grs_lm_rf <- lm(actual ~ pred_basic_rf + pred_grs_rf, merged_rf_preds)
basic_taxa_lm_rf <- lm(actual ~ pred_basic_rf + pred_grs_rf, merged_rf_preds)
basic_micom_lm_rf <- lm(actual ~ pred_basic_rf + pred_micom_rf, merged_rf_preds)
basic_path_lm_rf <- lm(actual ~ pred_basic_rf + pred_path_rf, merged_rf_preds)
anova(basic_lm_rf, basic_meta_lm_rf)
anova(basic_lm_rf, basic_grs_lm_rf)
anova(basic_lm_rf, basic_taxa_lm_rf)
anova(basic_lm_rf, basic_micom_lm_rf)
anova(basic_lm_rf, basic_path_lm_rf)

# Combine LASSO and RF
lasso_rf <- left_join(merged_lasso_preds, merged_rf_preds, by = "subject_id") %>% 
  rename(time = time.x, actual = actual.x) %>% dplyr::select(-c(actual.y, time.y))

summary(basic_lm)$adj.r.squared  
summary(basic_lm_rf)$adj.r.squared   


##### IGNORE BELOW ####################################
# remove target var & 'subject_id' 
test_x <- as.data.frame(test_imputed[, -which(names(test_imputed) %in% c("bmi_bL_6m", "subject_id"))])

# Make predictions on the test dataset
predictions <- predict(lasso_model, newdata = test_x)

# If the target variable is continuous (regression), you can now evaluate the predictions:
# Example for regression evaluation (using RMSE):
rmse <- sqrt(mean((predictions - test_data[[target]])^2))
print(paste("Root Mean Squared Error (RMSE):", rmse))

# If it's a classification task, you can check accuracy, confusion matrix, etc.
# Example for classification evaluation (using accuracy):
accuracy <- mean(predictions == test_data[[target]])
print(paste("Accuracy:", accuracy))

###### OLD SCHOOL WAY BELOW - IGNORE ###########################################

set.seed(123)
train_control <- trainControl(method = "cv", number = 5, search = "grid")

# In[5] Regression Models ----
m6_results <- train_and_save_models(
  imputed,
  "bmi_bL_6m",
  train_control,
  "feb20_bL_bmi_6m_BL_train_omic_g_ra_regression")

# describe the data
genus_ra_stats <- describe(omic_g_ra)
imputed_genus_ra_stats <- describe(omic_g_ra)

redundant_columns_genus <- names(omic_g_ra)[
  sapply(omic_g_ra, function(col) mean(col == 0, na.rm = TRUE) > 0.8) |
    genus_ra_stats$mean == 0]

columns_to_remove <- intersect(redundant_columns_genus, colnames(imputed))

# remove all of the redundant columns that are in genus_df_imputed and species_df_imputed
genus_df_imputed_minus_redundant <- remove_columns(imputed, 
                                                   columns_to_remove = columns_to_remove)
# retrain the models
genus_results <- train_and_save_models(
  genus_df_imputed_minus_redundant,
  "bmi_bL_6m",
  train_control,
  "feb20_bL_bmi_6m_BL_train_omic_g_regression_no_redundant")

###############################
###     Figure Analysis
###############################

# In[3]: Define base path and file paths ----
base_path <- "drift_fs/csv/results"

# Define file paths in a structured list
file_paths <- list(
  # genus
  bmi_bL_6m_BL_train_omic_g_ra_regression_beta = "feb20_bL_bmi_6m_BL_train_omic_g_ra_regression_beta.csv",
  bmi_bL_6m_BL_train_omic_g_ra_regression_feature_importance = "feb20_bL_bmi_6m_BL_train_omic_g_ra_regression_feature_importance.csv",
  bmi_bL_6m_BL_train_omic_g_ra_regression_metrics = "feb20_bL_bmi_6m_BL_train_omic_g_ra_regression_metrics.csv",
  # genus no redundant
  bmi_bL_6m_BL_train_omic_g_ra_regression_no_redundant_beta = "feb20_bL_bmi_6m_BL_train_omic_g_regression_no_redundant_beta.csv",
  bmi_bL_6m_BL_train_omic_g_ra_regression_no_redundant_feature_importance = "feb20_bL_bmi_6m_BL_train_omic_g_regression_no_redundant_feature_importance.csv",
  bmi_bL_6m_BL_train_omic_g_ra_regression_no_redundant_metrics = "feb20_bL_bmi_6m_BL_train_omic_g_regression_no_redundant_metrics.csv"
)

# Read all data into a named list using lapply
data_list <- lapply(file_paths, 
                    function(path) read.csv(file.path(base_path, path)))

# Assign names to the data list based on the file paths
names(data_list) <- names(file_paths)

# In[4]: Process and plot for all datasets ----
datasets <- list(
  "bmi_bL_6m_BL_train_omics_Genus_ra" = list(
    beta = data_list$bmi_bL_6m_BL_train_omic_g_ra_regression_beta,
    feature_importance = data_list$bmi_bL_6m_BL_train_omic_g_ra_regression_feature_importance,
    metrics = data_list$bmi_bL_6m_BL_train_omic_g_ra_regression_metrics
  ),
  "bmi_bL_6m_BL_train_omics_Genus_ra_No_Redundant)" = list(
    beta = data_list$bmi_bL_6m_BL_train_omic_g_ra_regression_no_redundant_beta,
    feature_importance = data_list$bmi_bL_6m_BL_train_omic_g_ra_regression_no_redundant_feature_importance,
    metrics = data_list$bmi_bL_6m_BL_train_omic_g_ra_regression_no_redundant_metrics
  )
)

# In[5]: Now from the features, lets do a heatmap for the top 20 features ----

# Get top 20 features for one to test
top_20_features <- get_top_n_features_all_models(data_list$bmi_bL_6m_BL_train_omic_g_ra_regression_feature_importance, 20)
print(top_20_features)

top_20_features_no_r <- get_top_n_features_all_models(data_list$bmi_bL_6m_BL_train_omic_g_ra_regression_no_redundant_feature_importance, 20)
print(top_20_features_no_r)

# In[9]: Metric R ^ 2 for presentation ----

# lets get the testing R^2 for all the different datasets
# Function to extract metrics and calculate max RÃ‚Â²
extract_metrics <- function(dataset_name, datasets) {
  metrics <- datasets[[dataset_name]]$metrics %>%
    dplyr::filter(DataType == "Test") %>%
    select(Model, R2)
  max_r2 <- max(metrics$R2, na.rm = TRUE)
  list(metrics = metrics, max_r2 = max_r2)
}

extract_metrics <- function(dataset) {
  if (is.null(dataset$metrics)) {
    stop("Metrics not found in the dataset.")
  }
  metrics <- dataset$metrics %>%
    dplyr::filter(DataType == "Test") %>%
    dplyr::select(Model, R2)
  max_r2 <- max(metrics$R2, na.rm = TRUE)
  list(metrics = metrics, max_r2 = max_r2)
}

bmi_bL_6m_BL_train_omics_Genus_ra = list(
  beta = data_list$bmi_bL_6m_BL_train_omic_g_ra_regression_beta,
  feature_importance = data_list$bmi_bL_6m_BL_train_omic_g_ra_regression_feature_importance,
  metrics = data_list$bmi_bL_6m_BL_train_omic_g_ra_regression_metrics)

bmi_bL_6m_BL_train_omics_Genus_ra_No_Redundant = list(
  beta = data_list$bmi_bL_6m_BL_train_omic_g_ra_regression_no_redundant_beta,
  feature_importance = data_list$bmi_bL_6m_BL_train_omic_g_ra_regression_no_redundant_feature_importance,
  metrics = data_list$bmi_bL_6m_BL_train_omic_g_ra_regression_no_redundant_metrics)

# Extract metrics and max RÃ‚Â² for all datasets
results_all <- extract_metrics(bmi_bL_6m_BL_train_omics_Genus_ra)
results_allno_re <- extract_metrics(bmi_bL_6m_BL_train_omics_Genus_ra_No_Redundant)

# Calculate overall max RÃ‚Â² for each main category
max_r2_genus <- results_all$max_r2
max_r2_genus_no_re <- results_allno_re$max_r2

# Calculate global max RÃ‚Â² across all categories
max_r2 <- max(max_r2_genus, max_r2_genus_no_re)
max_r2 <- 0.5

# Prepare data and titles for genus
all_g_data_list <- list(results_all$metrics)
all_g_data_list_no_re <- list(results_allno_re$metrics)

# rename the models
all_g_data_list[[1]]$Model <- c("Lasso", "Ridge", "Elastic Net", "Random Forest", "XGBoost")
all_g_data_list_no_re[[1]]$Model <- c("Lasso", "Ridge", "Elastic Net", "Random Forest", "XGBoost")

genus_titles <- c("All Omic Variables Feb 20")
genus_titles_no_redundant <- c("All Omic Variables Non Redundant Feb 20")

# Generate combined genus plot
all_omic_plot_genus_ra <- create_plots(all_g_data_list, max_r2, genus_titles)
pdf("drift_fs/figures/BL_omic_predicting_BL6m_BMI/feb20_BLomic_BL6mBMI_train_g_data_list.pdf", width = 7, height = 7)
print(all_omic_plot_genus_ra)
dev.off()

all_omic_plot_genus_ra_no_redundant <- create_plots(all_g_data_list_no_re, max_r2, genus_titles_no_redundant)
pdf("drift_fs/figures/BL_omic_predicting_BL6m_BMI/feb20_BLomic_BL6mBMI_train_g_data_list_no_re.pdf", width = 7, height = 7)
print(all_omic_plot_genus_ra_no_redundant)
dev.off()

# In[10]: Plotting the top 5-10 features ----

# Extract top features from each dataset
all_omic_genus_no_rendundant_features <- extract_top_features(datasets$`bmi_bL_6m_BL_train_omics_Genus_ra_No_Redundant)`)
all_omic_genus_features <- extract_top_features(datasets$bmi_bL_6m_BL_train_omics_Genus_ra)


# rename the column names
colnames(all_omic_genus_features) <- c("Variable", "Random Forest", "Lasso", 
                                       "Ridge", "Elastic Net",  "XGBoost")

colnames(all_omic_genus_no_rendundant_features) <- c("Variable", "Random Forest", 
                                                     "Lasso", "Ridge", 
                                                     "Elastic Net",  "XGBoost")

all_omic_genus_features <- all_omic_genus_features %>%
  mutate(
    Variable = case_when(
      Variable == "Leptin_BL" ~ "Leptin",
      Variable == "score_std" ~ "Genetic BMI risk score",
      Variable == "HOMA_IR_BL" ~ "Homeostasis Model Assessment",
      Variable == "avg_systolic_BL" ~ "Average Systolic Blood Pressure",
      Variable == "Insulin_endo_BL" ~ "Insulin",
      TRUE ~ Variable))

all_omic_genus_no_rendundant_features <- all_omic_genus_no_rendundant_features %>%
  mutate(
    Variable = case_when(
      Variable == "Leptin_BL" ~ "Leptin",
      Variable == "score_std" ~ "Genetic BMI risk score",
      Variable == "HOMA_IR_BL" ~ "Homeostasis Model Assessment",
      Variable == "avg_systolic_BL" ~ "Average Systolic Blood Pressure",
      Variable == "Insulin_endo_BL" ~ "Insulin",
      TRUE ~ Variable))

# Create and save plots for each dataset
create_feature_plot(
  all_omic_genus_features,
  "Top FI's - BL omics predicting BL-6m BMI change",
  "drift_fs/figures/BL_omic_predicting_BL6m_BMI/feb20_BLomic_BL6mBMI_train_g_feature_plot.pdf")

create_feature_plot(
  all_omic_genus_no_rendundant_features,
  "Top FI's - BL omics predicting BL-6m BMI change (non redundant)",
  "drift_fs/figures/BL_omic_predicting_BL6m_BMI/feb20_BLomic_BL6mBMI_train_g_no_rendundant_feature_plot.pdf")

# In[12]: Plotting the venn diagrams of the top features ----

# Extract top models for each dataset
datasets_names <- c("BL omic Genus ra", "BL omic Genus (No Redundant)")

top_models_list_all_omic_m12 <- bmi_bL_6m_BL_train_omics_Genus_ra$metrics %>%
  filter(DataType == "Test") %>%
  arrange(desc(R2)) %>%
  slice_head(n = 3) %>%
  pull(Model)

top_models_list_all_omic_m12_no_redundant <- bmi_bL_6m_BL_train_omics_Genus_ra_No_Redundant$metrics %>%
  filter(DataType == "Test") %>%
  arrange(desc(R2)) %>%
  slice_head(n = 3) %>%
  pull(Model)

# Print the top models for each data set
print(top_models_list_all_omic_m12)
print(top_models_list_all_omic_m12_no_redundant)