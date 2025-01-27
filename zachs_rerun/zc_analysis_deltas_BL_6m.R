#' @author Emily Yeo
#' @email emily.yeo@colorado.edu
#' @purpose Analysis for the Stanislawski Labm
#' @lab Stanislawski Lab
#' @affiliation University of Colorado Denver - Anschutz Medical, Dept of Biomedical Informatics & Personalized Medicine

###############################
###     Reading R data files    
###############################

# In[1]: Imports ----
rm(list = ls())
source("zc_functions.R") 
library(pacman)
p_load(tools, reticulate, viridis, tidyplots, patchwork, jsonlite, maps, ggvenn, caret, caretEnsemble, 
       readr, plyr, dplyr, tidyr, purrr, tibble, stringr, psych, randomForest, glmnet, xgboost, ggplot2, 
       reshape2, scales, gridExtra, plotly, sf, tidyverse)

###############################
###     Caret Analysis
###############################
# In[2] Load Data sets ----
data_dir <- "drift_fs/csv/all_omic_processed_data/deltas/"
omic_g_ra_outer <- read_csv(paste0(data_dir, "jan20_g_ra_all_omics_deltas_outer.csv"))
omic_g_ra_inner <- read_csv(paste0(data_dir, "jan20_g_ra_all_omics_deltas_inner.csv"))

# In[4] Main Analysis ----

# FIRST JUST WITH OUTER JOINED 
omic_g_ra <- omic_g_ra_outer

### Make BL , 6m and 12m dfs 
BL_6m <- omic_g_ra %>%
  filter(grepl("BL$", SampleID)) %>%
  select(-matches("_12m$|_6m$|_BL$")) %>%
  select(-matches(" (12m-BL)$"))
# Remove columns that end with " (12m-BL)"
BL_6m <- BL_6m[, !grepl(" \\(12m-BL\\)$", colnames(BL_6m))]


### Latent variables 
latent_variables_BL_6m <- c(
  "randomized_group", "score_std", "cohort_number", "sex", "race", "age", 
  "Glucose (6m-BL)", "HOMA-IR (6m-BL)", "Insulin (6m-BL)", 
  "HDL (6m-BL)", "LDL (6m-BL)", "Triglyceride lipid (6m-BL)")

### Process DFs 
imputed_BL_6m <- preprocess_data(BL_6m, 
                              latent_variables_BL_6m, 
                              "medianImpute")

# remove outcome_BMI_fnl_BL from the genus and species dataframes
tail(colnames(imputed_BL_6m), 300)
head(colnames(BL_6m), 30)

imputed <- remove_columns(
  imputed_BL_6m,
  columns_to_remove = c("subject_id", "BMI (6m-BL)", "Weight (6m-BL)", 
                        "outcome_wt_fnl_BL","outcome_wt_fnl_6m", "outcome_wt_fnl_12m",
                        "outcome_BMI_fnl_BL", "outcome_BMI_fnl_6m", "outcome_BMI_fnl_12m",          
                        "outcome_wt_fnl_BL_6m", "outcome_BMI_fnl_BL_6m"))

set.seed(123)
train_control <- trainControl(method = "cv", number = 5, search = "grid")

m6_results <- train_and_save_models(
  imputed,
  "BMI (6m-BL)",
  train_control,
  "ALL_deltas_BMI_BL_6m_g_ra_regr")

# describe the data
genus_ra_stats <- describe(omic_g_ra)

hist(genus_ra_stats$mean)
summary(genus_ra_stats$mean)

redundant_columns_genus <- names(omic_g_ra)[
  sapply(omic_g_ra, function(col) mean(col == 0, na.rm = TRUE) > 0.8) |
    genus_ra_stats$mean == 0]

# remove all of the redundant columns that are in genus_df_imputed and species_df_imputed
genus_df_imputed_minus_redundant <- remove_columns(imputed, 
                                                   columns_to_remove = redundant_columns_genus)

# retrain the models
genus_results <- train_and_save_models(
  genus_df_imputed_minus_redundant,
  "BMI (6m-BL)",
  train_control,
  "deltas_BMI_BL_6m_g_ra_reg_no_redundant")

###############################
###     Figure Analysis
###############################

# In[3]: Define base path and file paths ----
base_path <- "drift_fs/csv/results"

# Define file paths in a structured list
file_paths <- list(
  # genus
  deltas_BMI_BL_6m_g_ra_regr_beta = "deltas_BMI_BL_6m_g_ra_regr_beta.csv",
  deltas_BMI_BL_6m_g_ra_regr_feature_importance = "deltas_BMI_BL_6m_g_ra_regr_feature_importance.csv",
  deltas_BMI_BL_6m_g_ra_regr_metrics = "deltas_BMI_BL_6m_g_ra_regr_metrics.csv",
  # genus no redundant
  deltas_BMI_BL_6m_g_ra_reg_no_redundant_beta = "deltas_BMI_BL_6m_g_ra_reg_no_redundant_beta.csv",
  deltas_BMI_BL_6m_g_ra_reg_no_redundant_feature_importance = "deltas_BMI_BL_6m_g_ra_reg_no_redundant_feature_importance.csv",
  deltas_BMI_BL_6m_g_ra_reg_no_redundant_metrics = "deltas_BMI_BL_6m_g_ra_reg_no_redundant_metrics.csv"
)

# Read all data into a named list using lapply
data_list <- lapply(file_paths, 
                    function(path) read.csv(file.path(base_path, path)))

# Assign names to the data list based on the file paths
names(data_list) <- names(file_paths)

# In[4]: Process and plot for all datasets ----
datasets <- list(
  "BL_6m_deltas_genus_ra" = list(
    beta = data_list$deltas_BMI_BL_6m_g_ra_regr_beta,
    feature_importance = data_list$deltas_BMI_BL_6m_g_ra_regr_feature_importance,
    metrics = data_list$deltas_BMI_BL_6m_g_ra_regr_metrics
  ),
  "BL_6m_deltas_genus_ra_no_Redundant" = list(
    beta = data_list$deltas_BMI_BL_6m_g_ra_reg_no_redundant_beta,
    feature_importance = data_list$deltas_BMI_BL_6m_g_ra_reg_no_redundant_feature_importance,
    metrics = data_list$deltas_BMI_BL_6m_g_ra_reg_no_redundant_metrics
  )
)

# In[5]: Now from the features, lets do a heatmap for the top 20 features ----

# Get top 20 features for one to test
top_20_features <- get_top_n_features_all_models(data_list$m6_all_omic_g_ra_regression_feature_importance, 20)
print(top_20_features)

top_20_features_no_r <- get_top_n_features_all_models(data_list$m6_all_omic_g_ra_regression_no_redundant_feature_importance, 20)
print(top_20_features_no_r)

# lets get the testing R^2 for all the different datasets
# Function to extract metrics and calculate max RÂ²
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

deltas_BMI_BL_6m_g_ra_regr = list(
  beta = data_list$deltas_BMI_BL_6m_g_ra_regr_beta,
  feature_importance = data_list$deltas_BMI_BL_6m_g_ra_regr_feature_importance,
  metrics = data_list$deltas_BMI_BL_6m_g_ra_regr_metrics)

deltas_BMI_BL_6m_g_ra_reg_no_redundant = list(
  beta = data_list$deltas_BMI_BL_6m_g_ra_reg_no_redundant_beta,
  feature_importance = data_list$deltas_BMI_BL_6m_g_ra_reg_no_redundant_feature_importance,
  metrics = data_list$deltas_BMI_BL_6m_g_ra_reg_no_redundant_metrics)

# Extract metrics and max RÂ² for all datasets
results_all <- extract_metrics(deltas_BMI_BL_6m_g_ra_regr)
results_allno_re <- extract_metrics(deltas_BMI_BL_6m_g_ra_reg_no_redundant)

# Calculate overall max RÂ² for each main category
max_r2_genus <- results_all$max_r2
max_r2_genus_no_re <- results_allno_re$max_r2

# Calculate global max RÂ² across all categories
max_r2 <- max(max_r2_genus, max_r2_genus_no_re)
max_r2 <- 0.5

# Prepare data and titles for genus
all_g_data_list <- list(results_all$metrics)
all_g_data_list_no_re <- list(results_allno_re$metrics)

# Extract top features from each dataset
all_omic_genus_features <- extract_top_features(datasets$BL_6m_deltas_genus_ra)
all_omic_genus_no_rendundant_features <- extract_top_features(datasets$BL_6m_deltas_genus_ra_no_Redundant) 

all_omic_genus_features <- all_omic_genus_features %>%
  mutate(
    Variable = case_when(
      Variable == "Leptin_6m" ~ "Leptin",
      Variable == "score_std" ~ "Genetic BMI risk score",
      Variable == "HOMA_IR_6m" ~ "Homeostasis Model Assessment",
      Variable == "avg_systolic_6m" ~ "Average Systolic Blood Pressure",
      Variable == "Insulin_endo_6m" ~ "Insulin",
      TRUE ~ Variable))

all_omic_genus_no_rendundant_features <- all_omic_genus_no_rendundant_features %>%
  mutate(
    Variable = case_when(
      Variable == "Leptin_6m" ~ "Leptin",
      Variable == "score_std" ~ "Genetic BMI risk score",
      Variable == "HOMA_IR_6m" ~ "Homeostasis Model Assessment",
      Variable == "avg_systolic_6m" ~ "Average Systolic Blood Pressure",
      Variable == "Insulin_endo_6m" ~ "Insulin",
      TRUE ~ Variable))

# Create and save plots for each dataset
create_feature_plot(
  all_omic_genus_features,
  "Top 10 Features - Non Redundant Genus + All omics",
  "drift_fs/figures/deltas/jan22_deltas_BL_6m_all_g_feature_plot.pdf")

create_feature_plot(
  all_omic_genus_no_rendundant_features,
  "Top 10 Features - Non Redundant Genus + All omics",
  "drift_fs/figures/deltas/jan22_deltas_BL_6m_all_g_no_re_feature_plot.pdf")

# In[12]: Plotting the venn diagrams of the top features ----
# rename the models
all_g_data_list[[1]]$Model <- c("Lasso", "Ridge", "Elastic Net", "Random Forest", "XGBoost")
all_g_data_list_no_re[[1]]$Model <- c("Lasso", "Ridge", "Elastic Net", "Random Forest", "XGBoost")

genus_titles <- c("Delta BL-6m Variables (genus ra) - Model Testing RÂ²")
genus_titles_no_redundant <- c("Delta BL-6m Variables Non Redundant (genus ra) - Model Testing RÂ²")

# Generate combined genus plot
all_omic_plot_genus_ra <- create_plots(all_g_data_list, max_r2, genus_titles)
pdf("drift_fs/figures/deltas/jan22_deltas_BL_6m_all_g_6m.pdf", width = 7, height = 7)
print(all_omic_plot_genus_ra)
dev.off()

all_omic_plot_genus_ra_no_redundant <- create_plots(all_g_data_list_no_re, max_r2, genus_titles_no_redundant)
pdf("drift_fs/figures/deltas/jan22_deltas_BL_6m_all_g_no_re_6m.pdf", width = 7, height = 7)
print(all_omic_plot_genus_ra_no_redundant)
dev.off()

# In[10]: Plotting the top 5-10 features ----

# Extract top features from each dataset
all_omic_genus_features <- extract_top_features(datasets$BL_6m_deltas_genus_ra)
all_omic_genus_no_rendundant_features <- extract_top_features(datasets$BL_6m_deltas_genus_ra_no_Redundant) 

# rename the column names
colnames(all_omic_genus_features) <- c("Variable", "Random Forest", "Lasso", 
                                       "Ridge", "Elastic Net",  "XGBoost")

colnames(all_omic_genus_no_rendundant_features) <- c("Variable", "Random Forest", 
                                                     "Lasso", "Ridge", 
                                                     "Elastic Net",  "XGBoost")

all_omic_genus_features <- all_omic_genus_features %>%
  mutate(
    Variable = case_when(
      Variable == "Leptin_6m" ~ "Leptin",
      Variable == "score_std" ~ "Genetic BMI risk score",
      Variable == "HOMA_IR_6m" ~ "Homeostasis Model Assessment",
      Variable == "avg_systolic_6m" ~ "Average Systolic Blood Pressure",
      Variable == "Insulin_endo_6m" ~ "Insulin",
      TRUE ~ Variable))

all_omic_genus_no_rendundant_features <- all_omic_genus_no_rendundant_features %>%
  mutate(
    Variable = case_when(
      Variable == "Leptin_6m" ~ "Leptin",
      Variable == "score_std" ~ "Genetic BMI risk score",
      Variable == "HOMA_IR_6m" ~ "Homeostasis Model Assessment",
      Variable == "avg_systolic_6m" ~ "Average Systolic Blood Pressure",
      Variable == "Insulin_endo_6m" ~ "Insulin",
      TRUE ~ Variable))

# Create and save plots for each dataset
create_feature_plot(
  all_omic_genus_features,
  "Top 10 Features - Non Redundant Genus + All omics",
  "drift_fs/figures/deltas/jan22_deltas_BL_6m_all_g_feature_plot.pdf")

create_feature_plot(
  all_omic_genus_no_rendundant_features,
  "Top 10 Features - Non Redundant Genus + All omics",
  "drift_fs/figures/deltas/jan22_deltas_BL_6m_all_g_no_re_feature_plot.pdf")

# Extract top models for each dataset

top_models_list_all_omic_m6 <- deltas_BMI_BL_6m_g_ra_regr$metrics %>%
  filter(DataType == "Test") %>%
  arrange(desc(R2)) %>%
  slice_head(n = 3) %>%
  pull(Model)

top_models_list_all_omic_m6_no_redundant <- deltas_BMI_BL_6m_g_ra_reg_no_redundant$metrics %>%
  filter(DataType == "Test") %>%
  arrange(desc(R2)) %>%
  slice_head(n = 3) %>%
  pull(Model)

# Print the top models for each data set
print(top_models_list_all_omic_m6)
print(top_models_list_all_omic_m6_no_redundant)



