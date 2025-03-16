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
extra_var <- "Leptin"  # Change this to Ghrelin, Leptin, Hemoglobin_A1C, Peptide_YY, bmi_bL_6m

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

new_col_name_6m <- paste0(extra_var, "_6m")
m6_slim <- m6 %>%
  dplyr::select(c("subject_id", "outcome_BMI_fnl", extra_var)) %>%
  dplyr::rename(outcome_BMI_fnl_6m = outcome_BMI_fnl,
                !!new_col_name_6m := extra_var)

new_col_name_12m <- paste0(extra_var, "_12m")
m12_slim <- m12 %>% dplyr::select(c("subject_id", "outcome_BMI_fnl", extra_var)) %>%
  dplyr::rename(outcome_BMI_fnl_12m = outcome_BMI_fnl,
                !!new_col_name_12m := extra_var)

m12_6m <- merge(m12_slim, m6_slim, by = "subject_id", all = TRUE)
all <- merge(m12_6m, BL, by = "subject_id", all = TRUE)
all$bmi_bL_6m <- all$outcome_BMI_fnl_6m - all$outcome_BMI_fnl
all$extra_delta <- all[[new_col_name_6m]] - all[[new_col_name_12m]]
all <- all[!is.na(all$extra_delta), ] # Remove rows with missing values 

# Get all column names in the dataframe "all"
all_columns <- colnames(all)
columns_to_exclude <- c("subject_id", "sample_id", "all_samples", "record_id") # Exclude list
all_latent <- setdiff(all_columns, columns_to_exclude)

### Process DFs
preProc <- preProcess(all, method = 'medianImpute')
imputed_data <- predict(preProc, newdata = all)
imputed <- remove_columns(all, 
                          columns_to_remove = c("sample_id", "all_samples_y",
                                                "record_id", "all_samples",  
                                                "time_y", "outcome_BMI_fnl_12m", 
                                                "outcome_BMI_fnl_6m", "Unnamed..0_merged_data",
                                                new_col_name_6m, new_col_name_12m))
## Set single omic DFs

# Find the column indices for "proton" and "Carbon.dioxide" in train_set
proton_column <- which(colnames(imputed) == "proton")
carbon_dioxide_column <- which(colnames(imputed) == "Carbon.dioxide")
n10_col <- which(colnames(imputed) == "N10.formyl.tetrahydrofolate_biosynthesis")
lval_col <- which(colnames(imputed) == "L.valine_biosynthesis")
extra_delta = "extra_delta"
# Columns to KEEP for only meta 
only_basic <- c('subject_id', 'outcome_BMI_fnl', 'time', 'age',  
                'sex', extra_delta) # Use extra_delta here
meta_keep <- c('subject_id', 'outcome_BMI_fnl', 'time',  
               'randomized_group', 'cohort_number', 'sex', 'race', 'age',  
               'Glucose', 'HDL_Total_Direct_lipid', 'HOMA_IR',  
               'Insulin_endo', 'LDL_Calculated', 'Triglyceride_lipid', 
               extra_delta) # Use extra_delta here

only_grs <- c('subject_id', 'outcome_BMI_fnl', 'bmi_prs', 'time', 
              extra_delta) # Use extra_delta here
only_taxa <- c('subject_id', 'outcome_BMI_fnl', 'time', 
               extra_delta, 
               grep("^g__", colnames(imputed), value = TRUE)) # Use extra_delta here
only_micom <- c('subject_id', 'outcome_BMI_fnl', 'time', 
                extra_delta, 
                colnames(imputed)[proton_column:carbon_dioxide_column]) # Use extra_delta here
only_path <- c('subject_id', 'outcome_BMI_fnl', 'time', 
               extra_delta, 
               colnames(imputed)[n10_col:lval_col]) # Use extra_delta here

basic <- imputed %>% dplyr::select(all_of(only_basic))
meta <- imputed %>% dplyr::select(all_of(meta_keep))
grs <- imputed %>% dplyr::select(all_of(only_grs))
taxa <- imputed %>% dplyr::select(all_of(only_taxa))
path <- imputed %>% dplyr::select(all_of(only_path))
micom <- imputed %>% dplyr::select(all_of(only_micom))

## Repeat on test set
# Filter the test data just like the train data
BL_test <- test %>%
  dplyr::filter(grepl("0$", time)) %>%
  dplyr::select(-matches("12$|6$"))

m6_test <- test %>%
  dplyr::filter(grepl("6$", time)) %>%
  dplyr::select(-matches("0$|12$"))

m12_test <- test %>%
  dplyr::filter(grepl("12$", time)) %>%
  dplyr::select(-matches("0$|6$"))

m6_slim_test <- m6_test %>%
  dplyr::select(c("subject_id", "outcome_BMI_fnl", extra_var)) %>%
  dplyr::rename(outcome_BMI_fnl_6m = outcome_BMI_fnl,
                !!new_col_name_6m := extra_var)

m12_slim_test <- m12_test %>%
  dplyr::select(c("subject_id", "outcome_BMI_fnl", extra_var)) %>%
  dplyr::rename(outcome_BMI_fnl_12m = outcome_BMI_fnl,
                !!new_col_name_12m := extra_var)

# Merge test data for 12m and 6m
m12_6m_test <- merge(m12_slim_test, m6_slim_test, by = "subject_id", all = TRUE)

# Merge with BL data
all_test <- merge(m12_6m_test, BL_test, by = "subject_id", all = TRUE)

# Calculate bmi_bL_6m for the test data
all_test$bmi_bL_6m <- all_test$outcome_BMI_fnl_6m - all_test$outcome_BMI_fnl
all_test$extra_delta <- all_test[[new_col_name_6m]] - all_test[[new_col_name_12m]]
all_test <- all_test[!is.na(all_test$extra_delta), ] # Remove rows with missing values in bmi_bL_6m

# Get all column names in the test data "all_test"
all_columns_test <- colnames(all_test)
columns_to_exclude_test <- c("subject_id", "sample_id", "all_samples", "record_id") # Exclude list

# Select the relevant columns for the test set
all_latent_test <- setdiff(all_columns_test, columns_to_exclude_test)

preProc <- preProcess(all_test, method = 'medianImpute')
imputed_data <- predict(preProc, newdata = all_test)
# Process the test data (remove unwanted columns)
imputed_test <- remove_columns(imputed_data, 
                               columns_to_remove = c("sample_id", "all_samples_y", 
                                                     "record_id", "all_samples",  
                                                     "time_y", "outcome_BMI_fnl_12m",  
                                                     "outcome_BMI_fnl_6m", "Unnamed..0_merged_data",
                                                     new_col_name_6m, new_col_name_12m))

# Repeat the column selections and filter the relevant columns for test data
# For each of the omics data frames:
only_basic_test <- c('subject_id', 'outcome_BMI_fnl', 'time', 'age',  
                     'sex', extra_delta) # Use extra_delta here
meta_keep_test <- c('subject_id', 'outcome_BMI_fnl', 'time',  
                    'randomized_group', 'cohort_number', 'sex', 'race', 'age',  
                    'Glucose', 'HDL_Total_Direct_lipid', 'HOMA_IR',  
                    'Insulin_endo', 'LDL_Calculated', 'Triglyceride_lipid', 
                    extra_delta) # Use extra_delta here

only_grs_test <- c('subject_id', 'outcome_BMI_fnl', 'bmi_prs', 'time', 
                   extra_delta) # Use extra_delta here
only_taxa_test <- c('subject_id', 'outcome_BMI_fnl', 'time', 
                    extra_delta,  
                    grep("^g__", colnames(imputed_test), value = TRUE)) # Use extra_delta here
only_micom_test <- c('subject_id', 'outcome_BMI_fnl', 'time', 
                     extra_delta,  
                     colnames(imputed_test)[proton_column:carbon_dioxide_column]) # Use extra_delta here
only_path_test <- c('subject_id', 'outcome_BMI_fnl', 'time', 
                    extra_delta,  
                    colnames(imputed_test)[n10_col:lval_col]) # Use extra_delta here

# Select the final datasets for the test set
basic_test <- imputed_test %>% dplyr::select(all_of(only_basic_test))
meta_test <- imputed_test %>% dplyr::select(all_of(meta_keep_test))
grs_test <- imputed_test %>% dplyr::select(all_of(only_grs_test))
taxa_test <- imputed_test %>% dplyr::select(all_of(only_taxa_test))
path_test <- imputed_test %>% dplyr::select(all_of(only_path_test))
micom_test <- imputed_test %>% dplyr::select(all_of(only_micom_test))

# Optionally, remove intermediate variables to clear memory
rm(test_BL, test_m12, test_m12_6m, test_m12_slim, test_m6, 
   test_m6_slim, test_omic_g_ra, test, omic_g_ra, m6_slim, m6, m12, m12_slim,  
   m12_6m, BL, all_test)


############ RUN on Single Omics ##############################################
datasets <- list(basic, meta, grs, taxa, path, micom)
target_vars <- c(extra_delta,extra_delta,extra_delta,extra_delta,extra_delta,extra_delta)
result_prefixes <- c("basic_BL_bmi0_6m", "meta_BL_bmi0_6m", "grs_BL_bmi0_6m", 
                     "taxa_BL_bmi0_6m", "path_BL_bmi0_6m", "micom_BL_bmi0_6m")
# Append peptide_variable to each result prefix
result_prefixes_with_extra <- paste0(result_prefixes, "_", extra_delta)
set.seed(123)
train_control <- trainControl(method = "cv", number = 5, search = "grid")
train_and_save_multiple_models(datasets, target_vars, train_control, result_prefixes_with_extra)

# Define base path
base_path <- "drift_fs/csv/results/feb20"

# Define file paths using paste0
file_paths <- list(
  # basic
  basic_beta = paste0("basic_BL_bmi0_6m_", extra_delta, "_beta.csv"),
  basic_ft_imp = paste0("basic_BL_bmi0_6m_", extra_delta, "_feature_importance.csv"),
  basic_metrics = paste0("basic_BL_bmi0_6m_", extra_delta, "_metrics.csv"),
  # grs
  grs_beta = paste0("grs_BL_bmi0_6m_", extra_delta, "_beta.csv"),
  grs_ft_imp = paste0("grs_BL_bmi0_6m_", extra_delta, "_feature_importance.csv"),
  grs_metrics = paste0("grs_BL_bmi0_6m_", extra_delta, "_metrics.csv"),
  # meta
  meta_beta = paste0("meta_BL_bmi0_6m_", extra_delta, "_beta.csv"),
  meta_ft_imp = paste0("meta_BL_bmi0_6m_", extra_delta, "_feature_importance.csv"),
  meta_metrics = paste0("meta_BL_bmi0_6m_", extra_delta, "_metrics.csv"),
  # taxa
  taxa_beta = paste0("taxa_BL_bmi0_6m_", extra_delta, "_beta.csv"),
  taxa_ft_imp = paste0("taxa_BL_bmi0_6m_", extra_delta, "_feature_importance.csv"),
  taxa_metrics = paste0("taxa_BL_bmi0_6m_", extra_delta, "_metrics.csv"),
  # micom
  micom_beta = paste0("micom_BL_bmi0_6m_", extra_delta, "_beta.csv"),
  micom_ft_imp = paste0("micom_BL_bmi0_6m_", extra_delta, "_feature_importance.csv"),
  micom_metrics = paste0("micom_BL_bmi0_6m_", extra_delta, "_metrics.csv"),
  # pathway
  path_beta = paste0("path_BL_bmi0_6m_", extra_delta, "_beta.csv"),
  path_ft_imp = paste0("path_BL_bmi0_6m_", extra_delta, "_feature_importance.csv"),
  path_metrics = paste0("path_BL_bmi0_6m_", extra_delta, "_metrics.csv")
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

ggsave(paste0("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/zachs_rerun/drift_fs/csv/results/feb20/single_omic_BL_bmi0_6m_r2_plot_rf_lass_", extra_delta, ".png"), 
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
ggsave(paste0("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/zachs_rerun/drift_fs/csv/results/feb20/single_omic_BL_bmi0_6m_ft_imp_rf_lasso_", extra_delta,".png"), 
       plot = plot_ft_imp, width = 10, height = 8)

###############################################################################
datasets <- list(basic, meta, grs, taxa, path, micom)
# target variables
target_vars <- c(extra_delta,extra_delta,extra_delta,extra_delta,extra_delta,extra_delta)
result_prefixes <- c("basic_BL_bmi0_6m", "meta_BL_bmi0_6m", "grs_BL_bmi0_6m", 
                     "taxa_BL_bmi0_6m", "path_BL_bmi0_6m", "micom_BL_bmi0_6m")
result_prefixes_with_extra <- paste0(result_prefixes, "_", extra_delta)
results <- train_multiple_models(datasets, target_vars, train_control, result_prefixes_with_extra)

## LASSO results and PREDICTIONS
# basic results
lasso_mod_basic <- results[[paste0("basic_BL_bmi0_6m_",extra_delta)]]$lasso_model
lasso_pred_basic <- best_lasso_predextra(lasso_mod_basic, 
                                           basic, extra_delta, basic_test, extra_delta) %>% 
  rename(predicted_basic = predicted) %>% unique()

# meta
lasso_mod_meta <- results[[paste0("meta_BL_bmi0_6m_", extra_delta)]]$lasso_model
lasso_pred_meta <- best_lasso_predextra(lasso_mod_meta, 
                                          meta, extra_delta, meta_test, extra_delta) %>% 
  rename(predicted_meta = predicted) %>% unique()

# grs
lasso_mod_grs <- results[[paste0("grs_BL_bmi0_6m_", extra_delta)]]$lasso_model
lasso_pred_grs <- best_lasso_predextra(lasso_mod_grs, 
                                       grs, extra_delta, grs_test, extra_delta) %>% 
  rename(predicted_grs = predicted) %>% unique()

# taxa
lasso_mod_taxa <- results[[paste0("taxa_BL_bmi0_6m_", extra_delta)]]$lasso_model
lasso_pred_taxa <- best_lasso_predextra(lasso_mod_taxa, 
                                        taxa, extra_delta, taxa_test, extra_delta) %>% 
  rename(predicted_taxa = predicted) %>% unique()

# micom 
lasso_mod_micom <- results[[paste0("micom_BL_bmi0_6m_", extra_delta)]]$lasso_model
lasso_pred_micom <- best_lasso_predextra(lasso_mod_micom, 
                                         micom, extra_delta, micom_test, extra_delta) %>% 
  rename(predicted_micom = predicted) %>% unique()

# path 
lasso_mod_path <- results[[paste0("path_BL_bmi0_6m_", extra_delta)]]$lasso_model
lasso_pred_path <- best_lasso_predextra(lasso_mod_path, 
                                          path, extra_delta, 
                                          path_test, extra_delta, s = NULL) %>% 
  rename(predicted_path = predicted) %>% unique()

merged_lasso_preds <- lasso_pred_basic %>%
  left_join(lasso_pred_meta, by = "subject_id") %>%
  left_join(lasso_pred_grs, by = "subject_id") %>%
  left_join(lasso_pred_taxa, by = "subject_id") %>%
  left_join(lasso_pred_micom, by = "subject_id") %>%
  left_join(lasso_pred_path, by = "subject_id") %>% 
  dplyr::select(-c("actual.y.y", "time.y.y", "actual.x.x", "time.x.x",
                   "time.y", "actual.y", "time.y.y.y", "actual.y.y.y",
                   "time.x.x.x", "actual.x.x.x" )) %>% 
  rename(time = time.x, actual = actual.x) %>%
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
anova(meta_lm, taxa_lm)
anova(meta_lm, micom_lm)
anova(meta_lm, path_lm)

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
rf_mod_basic <- results[[paste0("basic_BL_bmi0_6m_", extra_delta)]]$rf_model
pred_basic_rf <- best_rf_predictions(rf_mod_basic, basic, extra_delta, basic_test) %>% 
  rename(pred_basic_rf = predicted) %>% unique()

rf_mod_meta <- results[[paste0("meta_BL_bmi0_6m_",extra_delta)]]$rf_model
pred_meta_rf <- best_rf_predictions(rf_mod_meta, meta, extra_delta, meta_test) %>% 
  rename(pred_meta_rf = predicted) %>% unique()

rf_mod_grs <- results[[paste0("grs_BL_bmi0_6m_", extra_delta)]]$rf_model
pred_grs_rf <- best_rf_predictions(rf_mod_grs, grs, extra_delta, grs_test) %>% 
  rename(pred_grs_rf = predicted) %>% unique()

rf_mod_taxa <- results[[paste0("taxa_BL_bmi0_6m_", extra_delta)]]$rf_model
pred_taxa_rf <- best_rf_predictions(rf_mod_taxa, taxa, extra_delta, taxa_test) %>% 
  rename(pred_taxa_rf = predicted) %>% unique()

rf_mod_micom <- results[[paste0("micom_BL_bmi0_6m_", extra_delta)]]$rf_model
pred_micom_rf <- best_rf_predictions(rf_mod_micom, micom, extra_delta, micom_test) %>% 
  rename(pred_micom_rf = predicted) %>% unique()

rf_mod_path <- results[[paste0("path_BL_bmi0_6m_", extra_delta)]]$rf_model
pred_path_rf <- best_rf_predictions(rf_mod_path, path, extra_delta, path_test) %>% 
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