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
# In[2] Load Datasets ----
data_dir <- "drift_fs/csv/all_omic_processed_data/"

long_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/merf_dfs/5.combined/"

test <- read.csv(file.path(long_dir, 'feb20_test_merged_all_omics_raw_meta.csv'))
train <- read.csv(file.path(long_dir, 'feb20_training_merged_all_omics_raw_meta.csv'))

# FIRST JUST WITH OUTER JOINED 
omic_g_ra <- train

### Make BL , 6m and 12m dfs 
BL <- omic_g_ra %>%
  filter(grepl("0$", time)) %>%
  select(-matches("12$|6$"))

m6 <- omic_g_ra %>%
  filter(grepl("6$", time)) %>%
  select(-matches("0$|12$")) 

m12 <- omic_g_ra %>%
  filter(grepl("12$", time)) %>%
  select(-matches("0$|6$"))

m6_slim <- m6 %>% select(c("subject_id", "outcome_BMI_fnl")) %>% 
  rename(outcome_BMI_fnl_6m = outcome_BMI_fnl)
m12_slim <- m12 %>% select(c("subject_id", "outcome_BMI_fnl")) %>% 
  rename(outcome_BMI_fnl_12m = outcome_BMI_fnl)

m12_6m <- merge(m12_slim, m6_slim, by = "subject_id", all = TRUE)
all <- merge(m12_6m, BL, by = "subject_id", all = TRUE)
all$bmi_bL_6m <- all$outcome_BMI_fnl_6m - all$outcome_BMI_fnl
all <- all[!is.na(all$bmi_bL_6m), ] # Remove rows with missing values in bmi_bL_6m (30)

# rm(BL, m12, m6, m12_6m, m12_slim, m6_slim, omic_g_ra, omic_g_ra_inner, omic_g_ra_outer)

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
                          columns_to_remove = c("subject_id", "sample_id", "all_samples_y",
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
only_basic <- c('outcome_BMI_fnl', 'time', 'age', 'sex')
meta_keep <- c('outcome_BMI_fnl', 'time', 'randomized_group', 'cohort_number', 
               'sex', 'race', 'age', 'Glucose', 'HDL_Total_Direct_lipid', 'HOMA_IR', 
               'Insulin_endo', 'LDL_Calculated', 'Triglyceride_lipid')
only_grs <- c('outcome_BMI_fnl', 'bmi_prs', 'time')
only_taxa <- c('outcome_BMI_fnl', 'time', 
               grep("^g__", colnames(imputed), value = TRUE))
only_micom <- c('outcome_BMI_fnl', 'time', 
                colnames(imputed)[proton_column:carbon_dioxide_column])
only_path <- c('outcome_BMI_fnl', 'time',
               colnames(imputed)[n10_col:lval_col])

basic <- imputed %>% select(all_of(only_basic))
meta <- imputed %>% select(all_of(meta_keep))
grs <- imputed %>% select(all_of(only_grs))
taxa <- imputed %>% select(all_of(only_taxa))
path <- imputed %>% select(all_of(only_path))
micom <- imputed %>% select(all_of(only_micom))

############ RUN on Single Omics ##############################################
# Set dfs 
datasets <- list(basic, meta, grs, taxa, path, micom)

# target variables
target_vars <- c("outcome_BMI_fnl", "outcome_BMI_fnl", "outcome_BMI_fnl", 
                 "outcome_BMI_fnl", "outcome_BMI_fnl", "outcome_BMI_fnl")

# result prefixes
result_prefixes <- c("basic_BL_bmi0_6m", "meta_BL_bmi0_6m", "grs_BL_bmi0_6m", 
                     "taxa_BL_bmi0_6m", "path_BL_bmi0_6m", "micom_BL_bmi0_6m")

# Call the function
set.seed(123)
train_control <- trainControl(method = "cv", number = 5, search = "grid")
#train_and_save_multiple_models(datasets, target_vars, train_control, result_prefixes)

base_path <- "drift_fs/csv/results/feb20"

# Define file paths in a structured list
file_paths <- list(
  # basic
  basic_beta = "basic_BL_bmi0_6m_beta.csv",
  basic_ft_imp = "basic_BL_bmi0_6m_feature_importance.csv",
  basic_metrics = "basic_BL_bmi0_6m_metrics.csv",
  # grs
  grs_beta = "grs_BL_bmi0_6m_beta.csv",
  grs_ft_imp = "grs_BL_bmi0_6m_feature_importance.csv",
  grs_metrics = "grs_BL_bmi0_6m_metrics.csv",
  # meta
  meta_beta = "meta_BL_bmi0_6m_beta.csv",
  meta_ft_imp = "meta_BL_bmi0_6m_feature_importance.csv",
  meta_metrics = "meta_BL_bmi0_6m_metrics.csv",
  # taxa
  taxa_beta = "taxa_BL_bmi0_6m_beta.csv",
  taxa_ft_imp = "taxa_BL_bmi0_6m_feature_importance.csv",
  taxa_metrics = "taxa_BL_bmi0_6m_metrics.csv",
  # micom
  micom_beta = "micom_BL_bmi0_6m_beta.csv",
  micom_ft_imp = "micom_BL_bmi0_6m_feature_importance.csv",
  micom_metrics = "micom_BL_bmi0_6m_metrics.csv",
  # pathway
  path_beta = "path_BL_bmi0_6m_beta.csv",
  path_ft_imp = "path_BL_bmi0_6m_feature_importance.csv",
  path_metrics = "path_BL_bmi0_6m_metrics.csv"
)

# Read all data into a named list using lapply
data_list <- lapply(file_paths, function(path) read.csv(file.path(base_path, path)))
names(data_list) <- names(file_paths) # Assign names to the data list

# In[4]: Process and plot for all datasets ----
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

# Create a list to store the R� values
r2_data <- list()
for (dataset_name in names(datasets)) { # Iterate through each dataset
  metrics_df <- datasets[[dataset_name]]$metrics # Get current dataset metrics
  # Filter for Train data and extract model names and R� values
  train_r2 <- metrics_df %>%
    dplyr::filter(DataType == "Train") %>%
    dplyr::select(Model, R2) %>%
    mutate(Dataset = dataset_name) # Add dataset name as a new column
  r2_data[[dataset_name]] <- train_r2 # Append the results to the r2_data list
}
r2_combined <- bind_rows(r2_data) # Combine R� single data frame

# Plot R� values using ggplot2
plot <- ggplot(r2_combined, aes(x = Model, y = R2, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(
    title = "R� Values for Train Models Across Datasets",
    x = "Model",
    y = "R�",
    fill = "Dataset"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plot)

ggsave("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/zachs_rerun/drift_fs/csv/results/feb20/single_omic_BL_bmi0_6m_r2_plot.png", 
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

# Plot feature importances using ggplot2 with independent y-axes for each facet
plot_ft_imp <- ggplot(ft_imp_combined_top10, 
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
  coord_flip()  # Flip the coordinates for better readability

print(plot_ft_imp)
# Save plot to a specific directory with a custom filename
ggsave("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/zachs_rerun/drift_fs/csv/results/feb20/single_omic_BL_bmi0_6m_ft_imp_plot.png", 
       plot = plot_ft_imp, width = 10, height = 8)

###############################################################################
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
# Function to extract metrics and calculate max RÂ²
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

# Extract metrics and max RÂ² for all datasets
results_all <- extract_metrics(bmi_bL_6m_BL_train_omics_Genus_ra)
results_allno_re <- extract_metrics(bmi_bL_6m_BL_train_omics_Genus_ra_No_Redundant)

# Calculate overall max RÂ² for each main category
max_r2_genus <- results_all$max_r2
max_r2_genus_no_re <- results_allno_re$max_r2

# Calculate global max RÂ² across all categories
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