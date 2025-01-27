#' @author Emily Yeo
#' @email emily.yeo@colorado.edu
#' @purpose Analysis for the Stanislawski Labm
#' @lab Stanislawski Lab
#' @affiliation University of Colorado Denver - Anschutz Medical, Dept of Biomedical Informatics & Personalized Medicine

###############################
###     Reading R data files    
###############################
rm(list = ls())
source("zc_functions.R") 
library(pacman)
p_load(tools, reticulate, viridis, tidyplots, patchwork, jsonlite, maps, ggvenn, caret, caretEnsemble, 
       readr, plyr, dplyr, tidyr, purrr, tibble, stringr, psych, randomForest, glmnet, xgboost, ggplot2, 
       reshape2, scales, gridExtra, plotly, sf, tidyverse)

###############################
###     Data Preprocessing
###############################

# In[2]: Data Imports ----

# Define the base path for the data files
output_dir <- "drift_fs/csv/unprocessed_data"
zc_pl_dir <- "unprocessed_input/"
local_path <- "drift_fs/csv/unprocessed_data/"
data_dir <- "drift_fs/csv/processed_data/"
df_dir = "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/merf_dfs/5.combined/"
func_dir = "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/micom/aim2/"
m1_dir = "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/clinical/transformed/aim2/"

# Genetic info
updated_analysis <- read_csv(paste0(zc_pl_dir, "grs.diff_110324.csv"))

# Taxa info
genus_clr_data <- read_csv(paste0(local_path, "genus.clr.csv"))
species_clr_data <- read_csv(paste0(local_path, "sp.clr.csv"))

genus_ra_df <- read_csv(paste0(data_dir, "../unprocessed_data/genus.ra.csv"))
sp_ra_df <- read_csv(paste0(data_dir, "../unprocessed_data/sp.ra.csv"))

genus_ra_df <- make_new_columns(genus_ra_df, genus_ra_df$SampleID)
#genus_ra_df <- filter_data(genus_ra_df, genus_ra_df$TIMEPOINT, "BL")

sp_ra_df <- make_new_columns(sp_ra_df, sp_ra_df$SampleID)
#sp_ra_df <- filter_data(sp_ra_df, sp_ra_df$TIMEPOINT, "BL")

# rename the columns
genus_ra_df <- rename_columns_species_to_domain(genus_ra_df)
sp_ra_df <- rename_columns_species_to_domain(sp_ra_df)

# Meta data 
merge_metadata <- read_csv(paste0(zc_pl_dir, "merge_meta_methyl.csv"))
metadata <- read_csv(paste0(zc_pl_dir, "DRIFT_working_dataset_meta_deltas_filtered_05.21.2024.csv"))
full_raw = read_csv(paste0(m1_dir, "a2_meta_not_Transformed_standard_clinical.csv"))

# Combined omics DF
test_all = read_csv(paste0(df_dir, 'test_merged_all_omics_raw_meta.csv'))  # Read in test_micom_no_na
train_all = read_csv(paste0(df_dir, 'training_merged_all_omics_raw_meta.csv'))  # Read in train_micom_no_na

# MICOM data
micom_test = read_csv(paste0(func_dir, "flux_all_clr_testing.csv"))
micom_train = read_csv(paste0(func_dir, "flux_all_clr_training.csv"))
# Add a column to indicate the source
micom_train <- micom_train %>%
  mutate(source = "train")

micom_test <- micom_test %>%
  mutate(source = "test")

# Combine the dataframes
micom_all <- bind_rows(micom_train, micom_test)

# Pathway 
pathway_df <- read_tsv(paste0(zc_pl_dir, "path_abun_unstrat.tsv"))
# Transpose the dataset and ensure it's a DataFrame
pathway_df <- pathway_df %>%
  column_to_rownames("pathway") %>%
  t() %>%
  as.data.frame() %>% 
  rownames_to_column("SampleID")

# In[4]: Data Preprocessing ----

# In[4.1]: Process genus and species clr data ----
genus_clr_data <- make_new_columns(genus_clr_data, genus_clr_data$SampleID)
#genus_clr_data <- filter_data(genus_clr_data, genus_clr_data$TIMEPOINT, "BL")

species_clr_data <- make_new_columns(species_clr_data, species_clr_data$SampleID)
#species_clr_data <- filter_data(species_clr_data, species_clr_data$TIMEPOINT, "BL")

# In[4.2]: Merge the updated_analysis and metadata ----

# Ensure both datasets have 'record_id' for the join
meta_data <- merge_data(updated_analysis,
                        metadata %>% select(-subject_id), 
                        inner_join, 
                        "record_id")

# list of the columns we want from the metadata
columns_to_extract_from_metadata <- c(
  "subject_id",
  "predicted_BL_BMI",
  "differences_BL_BMI",
  "diff_BMI_quartile",
  "diff_BMI_std"
)

# extract the columns from the metadata dataframe
merge_meta_data <- extract_columns(merge_metadata, 
                                   columns_to_extract = columns_to_extract_from_metadata)

# append the columns of merge_meta_data to meta_data
meta_data_df <- cbind(meta_data, 
                      merge_meta_data %>% select(-subject_id))

# only keep the consented samples
meta_data_df <- filter_data(meta_data_df, 
                            meta_data_df$consent, "yes")


### START MERGING GENUS RA with MICOM

# Inner join: Keeps only rows with matches in both dataframes
g_ra_micom_inner <- genus_ra_df %>%
  inner_join(micom_all, by = c("SampleID" = "sample_id"))

# Full join: Keeps all rows, with NA for non-matching rows
g_ra_micom_outer <- genus_ra_df %>%
  full_join(micom_all, by = c("SampleID" = "sample_id"))

### START MERGING GENUS RA & MICOM with Pathway data

# Inner join: Keeps only rows with matches in both dataframes
path_g_ra_micom_inner <- pathway_df %>%
  inner_join(g_ra_micom_inner, by = c("SampleID" = "SampleID"))

# Full join: Keeps all rows, with NA for non-matching rows
path_g_ra_micom_outer <- pathway_df %>%
  full_join(g_ra_micom_outer, by = c("SampleID" = "SampleID"))

### START MERGING GENUS RA & MICOM & Pathway with META data

# Inner join: Keeps only rows with matches in both dataframes
meta_slim <- meta_data_df %>% 
  dplyr::select(c("subject_id", "randomized_group", "score_std", 
                  "cohort_number", "sex", "race", "age", 
                  "outcome_wt_fnl_12m", "Glucose_BL", "HOMA_IR_BL",
                  "Insulin_endo_BL", "HDL_Total_Direct_lipid_BL", 
                  "LDL_Calculated_BL","Triglyceride_lipid_BL", 
                  "Glucose_6m", "HOMA_IR_6m", 
                  "Insulin_endo_6m", "HDL_Total_Direct_lipid_6m", 
                  "LDL_Calculated_6m", "Triglyceride_lipid_6m", 
                  "outcome_BMI_fnl_12m","Glucose_12m", "HOMA_IR_12m", 
                  "Insulin_endo_12m", "HDL_Total_Direct_lipid_12m",
                  "LDL_Calculated_12m", "Triglyceride_lipid_12m"))

meta_path_g_ra_micom_inner <- meta_slim %>%
  inner_join(path_g_ra_micom_inner, by = c("subject_id" = "subject_id"))

# Full join: Keeps all rows, with NA for non-matching rows
meta_path_g_ra_micom_outer <- meta_slim %>%
  full_join(path_g_ra_micom_outer, by = c("subject_id" = "subject_id"))

#rm(meta_data_df, meta_data, merge_meta_data, metadata, merge_metadata)

# In[4.4]: Remove the columns that are not needed ----
colnames(meta_path_g_ra_micom_outer)
tail(colnames(meta_path_g_ra_micom_outer), 300)

# Define columns and pattern to remove
columns_to_remove <- c(
  "all_samples",
  "...1.y",
  "TIMEPOINT",
  "source",
  "...1.x"
)
pattern_to_remove <- "3m|6m|12m|18m"
BL_pattern <- "3m|6m|12m|18m"
pattern_3m <- "BL|6m|12m|18m"
pattern_6m <- "3m|BL|12m|18m"
pattern_12m <- "3m|6m|BL|18m"

# Remove columns from genus dataset
g_ra_all_BL <- remove_columns(meta_path_g_ra_micom_outer, 
                              columns_to_remove = columns_to_remove, 
                              pattern = BL_pattern)

g_ra_all_outer <- meta_path_g_ra_micom_outer %>% 
  dplyr::select(-c("all_samples", "...1.y", "TIMEPOINT",
                   "source", "...1.x"))

g_ra_all_inner <- meta_path_g_ra_micom_inner %>% 
  dplyr::select(-c("all_samples", "...1.y", "TIMEPOINT",
                   "source", "...1.x"))

# save these dataframes
save_dir <- "drift_fs/csv/all_omic_processed_data/"

# check if the directory exists
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

write.csv(g_ra_all_outer, 
          paste0(save_dir, 
                 "jan22_wt_12m_genus_ra_all_omics_outer.csv"), 
          row.names = FALSE)

write.csv(g_ra_all_inner, 
          paste0(save_dir, 
                 "jan22_wt_12m_genus_ra_all_omics_inner.csv"), 
          row.names = FALSE)

###############################
###     Caret Analysis
###############################
rm(list = ls())
source("zc_functions.R") 
# In[2] Load Datasets ----
data_dir <- "drift_fs/csv/all_omic_processed_data/"
omic_g_ra_outer <- read_csv(paste0(data_dir, "jan22_wt_12m_genus_ra_all_omics_outer.csv"))
omic_g_ra_inner <- read_csv(paste0(data_dir, "jan22_wt_12m_genus_ra_all_omics_inner.csv"))

# In[4] Main Analysis ----

# FIRST JUST WITH OUTER JOINED 
omic_g_ra <- omic_g_ra_outer

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

### Latent variables 
latent_variables_BL <- c(
  "randomized_group", "score_std", "cohort_number", "sex", "race", "age", 
  "Glucose_BL", "HOMA_IR_BL", "Insulin_endo_BL", "HDL_Total_Direct_lipid_BL",             
  "LDL_Calculated_BL", "Triglyceride_lipid_BL", "outcome_wt_fnl_BL")

latent_variables_6m <- c(
  "randomized_group", "score_std", "cohort_number", "sex", "race", "age", 
  "Glucose_6m", "HOMA_IR_6m", "Insulin_endo_6m", "HDL_Total_Direct_lipid_6m",             
  "LDL_Calculated_6m", "Triglyceride_lipid_6m", "outcome_wt_fnl_6m")

latent_variables_12m <- c(
  "randomized_group", "score_std", "cohort_number", "sex", "race", "age", 
  "Glucose_12m", "HOMA_IR_12m", "Insulin_endo_12m", "HDL_Total_Direct_lipid_12m",             
  "LDL_Calculated_12m", "Triglyceride_lipid_12m", "outcome_wt_fnl_12m")

### Process DFs 
imputed_12m <- preprocess_data(m12, 
                              latent_variables_12m, 
                              "medianImpute")

# remove outcome_BMI_fnl_BL from the genus and species dataframes
tail(colnames(imputed_12m), 300)
imputed <- remove_columns(imputed_12m, 
                          columns_to_remove = c("subject_id", "SampleID", "outcome_BMI_fnl_12m"))

set.seed(123)
train_control <- trainControl(method = "cv", number = 5, search = "grid")

# In[5] Regression Models ----
m12_results <- train_and_save_models(
  imputed,
  "outcome_wt_fnl_12m",
  train_control,
  "m12_wt_all_omic_g_ra_regression")

# describe the data
genus_ra_stats <- describe(omic_g_ra)

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
  "outcome_wt_fnl_12m",
  train_control,
  "m12_wt_all_omic_g_regression_no_redundant"
)

###############################
###
###     Figure Analysis
###
###############################

# In[3]: Define base path and file paths ----
base_path <- "drift_fs/csv/results"

# Define file paths in a structured list
file_paths <- list(
  # genus
  m12_all_omic_g_ra_regression_beta = "m12_wt_all_omic_g_ra_regression_beta.csv",
  m12_all_omic_g_ra_regression_feature_importance = "m12_wt_all_omic_g_ra_regression_feature_importance.csv",
  m12_all_omic_g_ra_regression_metrics = "m12_wt_all_omic_g_ra_regression_metrics.csv",
  
  # genus no redundant
  m12_all_omic_g_ra_regression_no_redundant_beta = "m12_wt_all_omic_g_regression_no_redundant_beta.csv",
  m12_all_omic_g_ra_regression_no_redundant_feature_importance = "m12_wt_all_omic_g_regression_no_redundant_feature_importance.csv",
  m12_all_omic_g_ra_regression_no_redundant_metrics = "m12_wt_all_omic_g_regression_no_redundant_metrics.csv"
)

# Read all data into a named list using lapply
data_list <- lapply(file_paths, 
                    function(path) read.csv(file.path(base_path, path)))

# Assign names to the data list based on the file paths
names(data_list) <- names(file_paths)

# In[4]: Process and plot for all datasets ----
datasets <- list(
  "m12_all_omics_Genus_ra" = list(
    beta = data_list$m12_all_omic_g_ra_regression_beta,
    feature_importance = data_list$m12_all_omic_g_ra_regression_feature_importance,
    metrics = data_list$m12_all_omic_g_ra_regression_metrics
  ),
  "m12_all_omics_Genus_ra_No_Redundant)" = list(
    beta = data_list$m12_all_omic_g_ra_regression_no_redundant_beta,
    feature_importance = data_list$m12_all_omic_g_ra_regression_no_redundant_feature_importance,
    metrics = data_list$m12_all_omic_g_ra_regression_no_redundant_metrics
  )
)

# In[5]: Now from the features, lets do a heatmap for the top 20 features ----

# Get top 20 features for one to test
top_20_features <- get_top_n_features_all_models(data_list$m12_all_omic_g_ra_regression_feature_importance, 20)
print(top_20_features)

top_20_features_no_r <- get_top_n_features_all_models(data_list$m12_all_omic_g_ra_regression_no_redundant_feature_importance, 20)
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

m12_all_omics_Genus_ra = list(
  beta = data_list$m12_all_omic_g_ra_regression_beta,
  feature_importance = data_list$m12_all_omic_g_ra_regression_feature_importance,
  metrics = data_list$m12_all_omic_g_ra_regression_metrics)

m12_all_omics_Genus_ra_No_Redundant = list(
  beta = data_list$m12_all_omic_g_ra_regression_no_redundant_beta,
  feature_importance = data_list$m12_all_omic_g_ra_regression_no_redundant_feature_importance,
  metrics = data_list$m12_all_omic_g_ra_regression_no_redundant_metrics)

# Extract metrics and max RÂ² for all datasets
results_all <- extract_metrics(m12_all_omics_Genus_ra)
results_allno_re <- extract_metrics(m12_all_omics_Genus_ra_No_Redundant)

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
pdf("drift_fs/figures/all_omics_genus_ra/jan22_WT_12m_all_g_data_list.pdf", width = 7, height = 7)
print(all_omic_plot_genus_ra)
dev.off()

all_omic_plot_genus_ra_no_redundant <- create_plots(all_g_data_list_no_re, max_r2, genus_titles_no_redundant)
pdf("drift_fs/figures/all_omics_genus_ra/jan22_WT_12m_all_g_data_list_no_re.pdf", width = 7, height = 7)
print(all_omic_plot_genus_ra_no_redundant)
dev.off()

# In[10]: Plotting the top 5-10 features ----

# Extract top features from each dataset
all_omic_genus_features <- extract_top_features(datasets$m12_all_omics_Genus_ra)
all_omic_genus_no_rendundant_features <- extract_top_features(datasets$`m12_all_omics_Genus_ra_No_Redundant`)

# rename the column names
colnames(all_omic_genus_features) <- c("Variable", "Random Forest", "Lasso", "Ridge", "Elastic Net",  "XGBoost")
colnames(all_omic_genus_no_rendundant_features) <- c("Variable", "Random Forest", "Lasso", "Ridge", "Elastic Net",  "XGBoost")

all_omic_genus_features <- all_omic_genus_features %>%
  mutate(
    Variable = case_when(
      Variable == "Leptin_m12" ~ "Leptin",
      Variable == "score_std" ~ "Genetic BMI risk score",
      Variable == "HOMA_IR_m12" ~ "Homeostasis Model Assessment",
      Variable == "avg_systolic_m12" ~ "Average Systolic Blood Pressure",
      Variable == "Insulin_endo_m12" ~ "Insulin",
      TRUE ~ Variable))

all_omic_genus_no_rendundant_features <- all_omic_genus_no_rendundant_features %>%
  mutate(
    Variable = case_when(
      Variable == "Leptin_m12" ~ "Leptin",
      Variable == "score_std" ~ "Genetic BMI risk score",
      Variable == "HOMA_IR_m12" ~ "Homeostasis Model Assessment",
      Variable == "avg_systolic_m12" ~ "Average Systolic Blood Pressure",
      Variable == "Insulin_endo_m12" ~ "Insulin",
      TRUE ~ Variable))

# Create and save plots for each dataset
create_feature_plot(
  all_omic_genus_features,
  "Top 10 Features - Non Redundant Genus + All omics",
  "drift_fs/figures/all_omics_genus_ra/jan22_WT_all_m12_genus_feature_plot.pdf")

create_feature_plot(
  all_omic_genus_no_rendundant_features,
  "Top 10 Features - Non Redundant Genus + All omics",
  "drift_fs/figures/all_omics_genus_ra/jan22_WT_all_m12_genus_no_rendundant_feature_plot.pdf")
