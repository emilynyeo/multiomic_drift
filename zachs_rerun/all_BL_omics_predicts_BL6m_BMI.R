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
omic_g_ra_outer <- read_csv(paste0(data_dir, "jan18_genus_ra_all_omics_outer.csv"))
omic_g_ra_inner <- read_csv(paste0(data_dir, "jan18_genus_ra_all_omics_inner.csv"))

long_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/merf_dfs/5.combined/"

test <- read.csv(file.path(long_dir, 'test_merged_all_omics_raw_meta.csv'))
train <- read.csv(file.path(long_dir, 'training_merged_all_omics_raw_meta.csv'))

# Create the omic training set
omic_train <- omic_g_ra_outer[omic_g_ra_outer$subject_id %in% train$subject_id, ]
omic_test <- omic_g_ra_outer[omic_g_ra_outer$subject_id %in% test$subject_id, ]

# View the first few rows of each set to confirm
head(omic_train)
head(omic_test)

# FIRST JUST WITH OUTER JOINED 
omic_g_ra <- omic_train

### Make BL , 6m and 12m dfs 
BL <- omic_g_ra %>%
  filter(grepl("BL$", SampleID)) %>%
  select(-matches("_12m$|_6m$"))

m6 <- omic_g_ra %>%
  filter(grepl("6m$", SampleID)) %>%
  select(-matches("_BL$|_12m$")) 

m12 <- omic_g_ra %>%
  filter(grepl("12m$", SampleID)) %>%
  select(-matches("_BL$|_6m$"))

m6_slim <- m6 %>% select(c("subject_id", "outcome_BMI_fnl_6m"))
m12_slim <- m12 %>% select(c("subject_id", "outcome_BMI_fnl_12m"))

m12_6m <- merge(m12_slim, m6_slim, by = "subject_id")
all <- merge(m12_6m, BL, by = "subject_id")
all$bmi_bL_6m <- all$outcome_BMI_fnl_6m - all$outcome_BMI_fnl_BL
all <- all[!is.na(all$bmi_bL_6m), ] # Remove rows with missing values in bmi_bL_6m

# rm(BL, m12, m6, m12_6m, m12_slim, m6_slim, omic_g_ra, omic_g_ra_inner, omic_g_ra_outer)

### Latent variables 
latent_variables_BL <- c(
  "randomized_group", "score_std", "cohort_number", "sex", "race", "age", 
  "Glucose_BL", "HOMA_IR_BL", "Insulin_endo_BL", "HDL_Total_Direct_lipid_BL",             
  "LDL_Calculated_BL", "Triglyceride_lipid_BL", "outcome_BMI_fnl_BL")

# Get all column names in the dataframe "all"
all_columns <- colnames(all)
columns_to_exclude <- c("subject_id", "SampleID") # Exclude list
all_latent <- setdiff(all_columns, columns_to_exclude)

### Process DFs 
imputed <- preprocess_data(all, 
                           all_latent, 
                           "medianImpute")

# remove outcome_BMI_fnl_BL from the genus and species dataframes
tail(colnames(imputed), 300)
imputed <- remove_columns(imputed, 
                          columns_to_remove = c("subject_id", "SampleID", 
                                                "outcome_BMI_fnl_12m", "outcome_BMI_fnl_6m"))

set.seed(123)
train_control <- trainControl(method = "cv", number = 5, search = "grid")

# In[5] Regression Models ----
m6_results <- train_and_save_models(
  imputed,
  "bmi_bL_6m",
  train_control,
  "bmi_bL_6m_BL_train_omic_g_ra_regression")

# describe the data
genus_ra_stats <- describe(omic_g_ra)
imputed_genus_ra_stats <- describe(omic_g_ra)

hist(genus_ra_stats$mean)
hist(imputed_genus_ra_stats$mean)

summary(genus_ra_stats$mean)
summary(imputed_genus_ra_stats$mean)

redundant_columns_genus <- names(omic_g_ra)[
  sapply(omic_g_ra, function(col) mean(col == 0, na.rm = TRUE) > 0.8) |
    genus_ra_stats$mean == 0
]

# remove all of the redundant columns that are in genus_df_imputed and species_df_imputed
genus_df_imputed_minus_redundant <- remove_columns(imputed, 
                                                   columns_to_remove = redundant_columns_genus)


# retrain the models
genus_results <- train_and_save_models(
  genus_df_imputed_minus_redundant,
  "bmi_bL_6m",
  train_control,
  "bmi_bL_6m_BL_train_omic_g_regression_no_redundant")

###############################
###     Figure Analysis
###############################


# In[3]: Define base path and file paths ----
base_path <- "drift_fs/csv/results"

# Define file paths in a structured list
file_paths <- list(
  # genus
  bmi_bL_6m_BL_train_omic_g_ra_regression_beta = "bmi_bL_6m_BL_train_omic_g_ra_regression_beta.csv",
  bmi_bL_6m_BL_train_omic_g_ra_regression_feature_importance = "bmi_bL_6m_BL_train_omic_g_ra_regression_feature_importance.csv",
  bmi_bL_6m_BL_train_omic_g_ra_regression_metrics = "bmi_bL_6m_BL_train_omic_g_ra_regression_metrics.csv",
  # genus no redundant
  bmi_bL_6m_BL_train_omic_g_ra_regression_no_redundant_beta = "bmi_bL_6m_BL_train_omic_g_regression_no_redundant_beta.csv",
  bmi_bL_6m_BL_train_omic_g_ra_regression_no_redundant_feature_importance = "bmi_bL_6m_BL_train_omic_g_regression_no_redundant_feature_importance.csv",
  bmi_bL_6m_BL_train_omic_g_ra_regression_no_redundant_metrics = "bmi_bL_6m_BL_train_omic_g_regression_no_redundant_metrics.csv"
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

genus_titles <- c("All Omic Variables (genus ra) - Model Testing RÂ²")
genus_titles_no_redundant <- c("All Omic Variables Non Redundant (genus ra) - Model Testing RÂ²")

# Generate combined genus plot
all_omic_plot_genus_ra <- create_plots(all_g_data_list, max_r2, genus_titles)
pdf("drift_fs/figures/BL_omic_predicting_BL6m_BMI/feb19_BLomic_BL6mBMI_train_g_data_list.pdf", width = 7, height = 7)
print(all_omic_plot_genus_ra)
dev.off()

all_omic_plot_genus_ra_no_redundant <- create_plots(all_g_data_list_no_re, max_r2, genus_titles_no_redundant)
pdf("drift_fs/figures/BL_omic_predicting_BL6m_BMI/feb19_BLomic_BL6mBMI_train_g_data_list_no_re.pdf", width = 7, height = 7)
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
  "drift_fs/figures/BL_omic_predicting_BL6m_BMI/feb19_BLomic_BL6mBMI_train_g_feature_plot.pdf")

create_feature_plot(
  all_omic_genus_no_rendundant_features,
  "Top FI's - BL omics predicting BL-6m BMI change (non redundant)",
  "drift_fs/figures/BL_omic_predicting_BL6m_BMI/feb19_BLomic_BL6mBMI_train_g_no_rendundant_feature_plot.pdf")

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



