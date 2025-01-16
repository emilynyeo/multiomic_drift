#' @author Zachary Caterer
#' @email ztcaterer@colorado.edu
#' @purpose Analysis for the Stanislawski Lab and Interdiciplinary Quantitative Biology Program
#' @lab Stanislawski Lab
#' @affiliation University of Colorado Denver - Anschutz Medical Campus, Department of Biomedical Informatics and Personalized Medicine

###############################
###
###     Reading R data files    
###
###############################

# In[1]: Imports ----
rm(list = ls())
source("zc_functions.R") 
p_load(tools, reticulate, viridis, tidyplots, patchwork, jsonlite, maps, ggvenn, caret, caretEnsemble, 
readr, plyr, dplyr, tidyr, purrr, tibble, stringr, psych, randomForest, glmnet, xgboost, ggplot2, 
reshape2, scales, gridExtra, plotly, sf, tidyverse)

# In[2]: Functions ----

# Function to load RData file into a new environment and return the environment
load_rdata <- py$load_rdata
# Function to automatically save all data frames/matrices in an environment as CSV files
save_env_to_csv <- py$save_env_to_csv

# In[3]: Taxa Data ----
# Load another RData file for taxa data and save to CSV
genus_tables <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/zachs_rerun/unprocessed_input/Genus_Sp_tables.RData"
env_taxa <- load_rdata(genus_tables)

output_dir <- "drift_fs/csv/unprocessed_data"
save_env_to_csv(env_taxa, output_dir)

###############################
###     Data Preprocessing
###############################

# In[2]: Data Imports ----

# Define the base path for the data files
zc_pl_dir <- "unprocessed_input/"
local_path <- "drift_fs/csv/unprocessed_data/"

updated_analysis <- read_csv(paste0(zc_pl_dir, "grs.diff_110324.csv"))
genus_clr_data <- read_csv(paste0(local_path, "genus.clr.csv"))
species_clr_data <- read_csv(paste0(local_path, "sp.clr.csv"))
merge_metadata <- read_csv(paste0(zc_pl_dir, "merge_meta_methyl.csv"))
metadata <- read_csv(paste0(zc_pl_dir, "DRIFT_working_dataset_meta_deltas_filtered_05.21.2024.csv"))

# In[3]: Functions ----
make_new_columns <- py$make_new_columns
filter_data <- py$filter_data
merge_data <- py$merge_data
remove_columns <- py$remove_columns
extract_columns <- py$extract_columns
rename_columns_species_to_domain <- py$rename_columns_species_to_domain

# In[4]: Data Preprocessing ----

# In[4.1]: Process genus and species clr data ----

genus_clr_data <- make_new_columns(genus_clr_data, genus_clr_data$SampleID)
genus_clr_data <- filter_data(genus_clr_data, genus_clr_data$TIMEPOINT, "BL")

species_clr_data <- make_new_columns(species_clr_data, species_clr_data$SampleID)
species_clr_data <- filter_data(species_clr_data, species_clr_data$TIMEPOINT, "BL")

# In[4.2]: Merge the updated_analysis and metadata ----

# Ensure both datasets have 'record_id' for the join
meta_data <- merge_data(updated_analysis,metadata %>% select(-subject_id), inner_join, "record_id")

# list of the columns we want from the metadata
columns_to_extract_from_metadata <- c(
  "subject_id",
  "predicted_BL_BMI",
  "differences_BL_BMI",
  "diff_BMI_quartile",
  "diff_BMI_std"
)

# extract the columns from the metadata dataframe
merge_meta_data <- extract_columns(merge_metadata, columns_to_extract = columns_to_extract_from_metadata)

# append the columns of merge_meta_data to meta_data
meta_data_df <- cbind(meta_data, merge_meta_data %>% select(-subject_id))

# only keep the consented samples
meta_data_df <- filter_data(meta_data_df, meta_data_df$consent, "yes")

# In[4.3]: Merge the genus and species clr data with the metadata ----
genus_clr_data <- merge_data(genus_clr_data, meta_data_df, inner_join, "subject_id")
species_clr_data <- merge_data(species_clr_data, meta_data_df, inner_join, "subject_id")

rm(meta_data_df, meta_data, merge_meta_data, metadata, merge_metadata)

# In[4.4]: Remove the columns that are not needed ----

# Define columns and pattern to remove
columns_to_remove <- c(
  "SampleID",
  "TIMEPOINT",
  "record_id",
  "withdrawal_date_check",
  "start_treatment"
)
pattern_to_remove <- "3m|6m|12m|18m"

# Remove columns from genus dataset
genus_clr_data <- remove_columns(genus_clr_data, columns_to_remove = columns_to_remove, pattern = pattern_to_remove)

# Remove columns from species dataset
species_clr_data <- remove_columns(species_clr_data, columns_to_remove = columns_to_remove, pattern = pattern_to_remove)

# In[4.5]: Finalize Genus and Species datasets ----
print(colnames(genus_clr_data))
# columns we want to use in the analysis
# any column that contains "__"
latent_variables_to_use <- c(
  "subject_id",
  "age",
  "sex",
  "cohort_number",
  "race",
  "ethnicity",
  "education",
  "rmr_kcald_BL",
  "spk_EE_int_kcal_day_BL",
  "avg_systolic_BL",
  "avg_diastolic_BL",
  "C_Reactive_Protein_BL",
  "Cholesterol_lipid_BL",
  "Ghrelin_BL",
  "Glucose_BL",
  "HDL_Total_Direct_lipid_BL",
  "Hemoglobin_A1C_BL",
  "Insulin_endo_BL",
  "LDL_Calculated_BL",
  "Leptin_BL",
  "Peptide_YY_BL",
  "Triglyceride_lipid_BL",
  "HOMA_IR_BL",
  "outcome_BMI_fnl_BL",
  # prediction variables for the regression model are
  # "differences_BL_BMI",
  "diff_std_bmi_score"
)

# ensure all the columns are present in the data
genus_clr_latent_unclean <- extract_columns(
  genus_clr_data,
  columns_to_extract = latent_variables_to_use,
  pattern = "__"
)

species_clr_latent_unclean <- extract_columns(
  species_clr_data,
  columns_to_extract = latent_variables_to_use,
  pattern = "__"
)
print(colnames(genus_clr_latent_unclean))

# Apply the function to your dataframe
genus_clr_latent_clean <- rename_columns_species_to_domain(genus_clr_latent_unclean)
species_clr_latent_clean <- rename_columns_species_to_domain(species_clr_latent_unclean)

# Verify the results
print(colnames(genus_clr_latent_clean))
print(colnames(species_clr_latent_clean))

# save these dataframes
save_dir <- "drift_fs/csv/processed_data/"

# check if the directory exists
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

write.csv(genus_clr_latent_clean, paste0(save_dir, "genus_latent.csv"), row.names = FALSE)
write.csv(species_clr_latent_clean, paste0(save_dir, "species_latent.csv"), row.names = FALSE)

###############################
###     Caret Analysis
###############################

# In[2] Load Datasets ----
data_dir <- "drift_fs/csv/processed_data/"
species_df <- read_csv(paste0(data_dir, "species_latent.csv"))
genus_df <- read_csv(paste0(data_dir, "genus_latent.csv"))

# In[4] Main Analysis ----
latent_variables_to_use <- c(
  "subject_id",
  "age",
  "sex",
  "cohort_number",
  "race",
  "ethnicity",
  "education",
  "rmr_kcald_BL",
  "spk_EE_int_kcal_day_BL",
  "avg_systolic_BL",
  "avg_diastolic_BL",
  "C_Reactive_Protein_BL",
  "Cholesterol_lipid_BL",
  "Ghrelin_BL",
  "Glucose_BL",
  "HDL_Total_Direct_lipid_BL",
  "Hemoglobin_A1C_BL",
  "Insulin_endo_BL",
  "LDL_Calculated_BL",
  "Leptin_BL",
  "Peptide_YY_BL",
  "Triglyceride_lipid_BL",
  "HOMA_IR_BL",
  # prediction variables for the regression model are
  # "differences_BL_BMI",
  # "outcome_BMI_fnl_BL",
  "diff_std_bmi_score"
)

genus_df_imputed <- preprocess_data(genus_df, latent_variables_to_use, "medianImpute")
species_df_imputed <- preprocess_data(species_df, latent_variables_to_use, "medianImpute")

# remove outcome_BMI_fnl_BL from the genus and species dataframes
genus_df_imputed <- remove_columns(genus_df_imputed, columns_to_remove = "outcome_BMI_fnl_BL")
species_df_imputed <- remove_columns(species_df_imputed, columns_to_remove = "outcome_BMI_fnl_BL")

set.seed(123)
train_control <- trainControl(method = "cv", number = 5, search = "grid")

# In[5] Regression Models ----

genus_results <- train_and_save_models(
  genus_df_imputed,
  "diff_std_bmi_score",
  train_control,
  "genus_regression"
)

species_results <- train_and_save_models(
  species_df_imputed,
  "diff_std_bmi_score",
  train_control,
  "species_regression"
)

# In[7] Remove columns that have missing data from clr datasets ----

genus_ra_df <- read_csv(paste0(data_dir, "../unprocessed_data/genus.ra.csv"))
sp_ra_df <- read_csv(paste0(data_dir, "../unprocessed_data/sp.ra.csv"))

genus_ra_df <- make_new_columns(genus_ra_df, genus_ra_df$SampleID)
genus_ra_df <- filter_data(genus_ra_df, genus_ra_df$TIMEPOINT, "BL")

sp_ra_df <- make_new_columns(sp_ra_df, sp_ra_df$SampleID)
sp_ra_df <- filter_data(sp_ra_df, sp_ra_df$TIMEPOINT, "BL")

genus_ra_df <- remove_columns(genus_ra_df, c("SampleID", "TIMEPOINT", "...1"))
sp_ra_df <- remove_columns(sp_ra_df, c("SampleID", "TIMEPOINT", "...1"))

# rename the columns
genus_ra_df <- rename_columns_species_to_domain(genus_ra_df)
sp_ra_df <- rename_columns_species_to_domain(sp_ra_df)

# describe the data
genus_ra_df_stats <- describe(genus_ra_df)
species_ra_df_stats <- describe(sp_ra_df)

redundant_columns_genus <- names(genus_ra_df)[
  sapply(genus_ra_df, function(col) mean(col == 0, na.rm = TRUE) > 0.8) |
    genus_ra_df_stats$mean == 0
]

redundant_columns_sp <- names(sp_ra_df)[
  sapply(sp_ra_df, function(col) mean(col == 0, na.rm = TRUE) > 0.8) |
    species_ra_df_stats$mean == 0
]

# remove all of the redundant columns that are in genus_df_imputed and species_df_imputed
genus_df_imputed_minus_redundant <- remove_columns(genus_df_imputed, columns_to_remove = redundant_columns_genus)
species_df_imputed_minus_redundant <- remove_columns(species_df_imputed, columns_to_remove = redundant_columns_sp)

# retrain the models
genus_results <- train_and_save_models(
  genus_df_imputed_minus_redundant,
  "diff_std_bmi_score",
  train_control,
  "genus_regression_no_redundant"
)

species_results <- train_and_save_models(
  species_df_imputed_minus_redundant,
  "diff_std_bmi_score",
  train_control,
  "species_regression_no_redundant"
)

# remove the latent variables from the genus and species dataframes except for the target variable
genus_df_imputed_minus_redundant_no_latent <- remove_columns(genus_df_imputed_minus_redundant, latent_variables_to_use[-which(latent_variables_to_use %in% c("diff_std_bmi_score"))])
species_df_imputed_minus_redundant_no_latent <- remove_columns(species_df_imputed_minus_redundant, latent_variables_to_use[-which(latent_variables_to_use %in% c("diff_std_bmi_score"))])

genus_results_no_latent <- train_and_save_models(
  genus_df_imputed_minus_redundant_no_latent,
  "diff_std_bmi_score",
  train_control,
  "genus_regression_no_redundant_no_latent"
)

species_results_no_latent <- train_and_save_models(
  species_df_imputed_minus_redundant_no_latent,
  "diff_std_bmi_score",
  train_control,
  "species_regression_no_redundant_no_latent"
)

# remove the latent variables from the genus and species dataframes except for the target variable
genus_df_imputed_no_latent <- remove_columns(genus_df_imputed, latent_variables_to_use[-which(latent_variables_to_use %in% c("diff_std_bmi_score"))])
species_df_imputed_no_latent <- remove_columns(species_df_imputed, latent_variables_to_use[-which(latent_variables_to_use %in% c("diff_std_bmi_score"))])

print(colnames(genus_df_imputed_no_latent))
print(head(genus_df_imputed_no_latent))

genus_results_no_latent <- train_and_save_models(
  genus_df_imputed_no_latent,
  "diff_std_bmi_score",
  train_control,
  "genus_regression_no_latent"
)

species_results_no_latent <- train_and_save_models(
  species_df_imputed_no_latent,
  "diff_std_bmi_score",
  train_control,
  "species_regression_no_latent"
)

# In[8] pathway analysis ----

# Load the pathway data
pathway_df <- read_tsv(paste0(zc_pl_dir, "path_abun_unstrat.tsv"))

# extract only the latent variables from the genus_df_imputed
genus_df_imputed <- genus_df_imputed %>%
  select(all_of(latent_variables_to_use))

# Extract all columns containing "BL" in the name and the pathway column
pathway_df <- pathway_df %>%
  select(matches("BL") | matches("^pathway$")) %>%
  select(-matches("3m|6m|12m|18m"))

# Rename column names by removing ".BL"
colnames(pathway_df) <- gsub("\\.BL", "", colnames(pathway_df))

subject_id <- colnames(pathway_df)[!colnames(pathway_df) %in% c("pathway")]

# Transpose the dataset and ensure it's a DataFrame
pathway_df <- pathway_df %>%
  column_to_rownames("pathway") %>%
  t() %>%
  as.data.frame()

# Add subject_id as a column
pathway_df$subject_id <- subject_id

# innerjoin the pathway_df with the genus_df_imputed  by subject_id
pathway_df <- inner_join(genus_df_imputed, pathway_df, by = "subject_id")
print(colnames(pathway_df))

# run the models on the pathway data
pathway_results <- train_and_save_models(
  pathway_df,
  "diff_std_bmi_score",
  train_control,
  "pathway_regression"
)

# remove the redundant columns from the pathway data
pathway_df_stats <- describe(pathway_df)

# Identify columns with more than 80% zeros or with zero mean
redundant_columns_pathway <- names(pathway_df)[
  sapply(pathway_df, function(col) mean(col == 0, na.rm = TRUE) > 0.8) |
    pathway_df_stats$mean == 0
]

# Remove redundant columns from the pathway data
pathway_df_minus_redundant <- remove_columns(pathway_df, columns_to_remove = redundant_columns_pathway)

# retrain the models
pathway_results_minus_rendundant <- train_and_save_models(
  pathway_df_minus_redundant,
  "diff_std_bmi_score",
  train_control,
  "pathway_regression_no_redundant"
)

# remove the latent variables from the genus and species dataframes except for the target variable
pathway_df_no_latent <- remove_columns(pathway_df, columns_to_remove = latent_variables_to_use[-which(latent_variables_to_use %in% c("diff_std_bmi_score"))])

pathway_df_no_latent_results <- train_and_save_models(
  pathway_df_no_latent,
  "diff_std_bmi_score",
  train_control,
  "pathway_regression_no_latent"
)

pathway_df_minus_redundant_no_latent <- remove_columns(pathway_df_minus_redundant, columns_to_remove = latent_variables_to_use[-which(latent_variables_to_use %in% c("diff_std_bmi_score"))])

pathway_df_minus_redundant_no_latent_results <- train_and_save_models(
  pathway_df_minus_redundant_no_latent,
  "diff_std_bmi_score",
  train_control,
  "pathway_regression_no_redundant_no_latent"
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
  genus_regression_beta = "genus_regression_beta.csv",
  genus_regression_feature_importance = "genus_regression_feature_importance.csv",
  genus_regression_metrics = "genus_regression_metrics.csv",
  
  # genus no redundant
  genus_regression_no_redundant_beta = "genus_regression_no_redundant_beta.csv",
  genus_regression_no_redundant_feature_importance = "genus_regression_no_redundant_feature_importance.csv",
  genus_regression_no_redundant_metrics = "genus_regression_no_redundant_metrics.csv",
  
  # genus no latent
  genus_regression_no_latent_beta = "genus_regression_no_latent_beta.csv",
  genus_regression_no_latent_feature_importance = "genus_regression_no_latent_feature_importance.csv",
  genus_regression_no_latent_metrics = "genus_regression_no_latent_metrics.csv",
  
  # pathway regression
  pathway_regression_beta = "pathway_regression_beta.csv",
  pathway_regression_feature_importance = "pathway_regression_feature_importance.csv",
  pathway_regression_metrics = "pathway_regression_metrics.csv",
  
  # pathway regression no redundant
  pathway_regression_no_redundant_beta = "pathway_regression_no_redundant_beta.csv",
  pathway_regression_no_redundant_feature_importance = "pathway_regression_no_redundant_feature_importance.csv",
  pathway_regression_no_redundant_metrics = "pathway_regression_no_redundant_metrics.csv", 
  
  # pathway regression no latent
  pathway_regression_no_latent_beta = "pathway_regression_no_latent_beta.csv",
  pathway_regression_no_latent_feature_importance = "pathway_regression_no_latent_feature_importance.csv",
  pathway_regression_no_latent_metrics = "pathway_regression_no_latent_metrics.csv",
  
  # pathway subset regression
  # pathway_subset_regression_beta = "pathway_subset_regression_beta.csv",
  # pathway_subset_regression_feature_importance = "pathway_subset_regression_feature_importance.csv",
  # pathway_subset_regression_metrics = "pathway_subset_regression_metrics.csv",
  
  # pathway subset regression no redundant
  # pathway_subset_regression_no_redundant_beta = "pathway_subset_regression_no_redundant_beta.csv",
  # pathway_subset_regression_no_redundant_feature_importance = "pathway_subset_regression_no_redundant_feature_importance.csv",
  # pathway_subset_regression_no_redundant_metrics = "pathway_subset_regression_no_redundant_metrics.csv",
  
  # species regression
  species_regression_beta = "species_regression_beta.csv",
  species_regression_feature_importance = "species_regression_feature_importance.csv",
  species_regression_metrics = "species_regression_metrics.csv",
  
  # species regression no redundant
  species_regression_no_redundant_beta = "species_regression_no_redundant_beta.csv",
  species_regression_no_redundant_feature_importance = "species_regression_no_redundant_feature_importance.csv",
  species_regression_no_redundant_metrics = "species_regression_no_redundant_metrics.csv",
  
  # species regression no latent
  species_regression_no_latent_beta = "species_regression_no_latent_beta.csv",
  species_regression_no_latent_feature_importance = "species_regression_no_latent_feature_importance.csv",
  species_regression_no_latent_metrics = "species_regression_no_latent_metrics.csv"
)

# Read all data into a named list using lapply
data_list <- lapply(file_paths, function(path) read.csv(file.path(base_path, path)))

# Assign names to the data list based on the file paths
names(data_list) <- names(file_paths)

# In[4]: Process and plot for all datasets ----
datasets <- list(
  "Genus" = list(
    beta = data_list$genus_regression_beta,
    feature_importance = data_list$genus_regression_feature_importance,
    metrics = data_list$genus_regression_metrics
  ),
  "Genus (No Redundant)" = list(
    beta = data_list$genus_regression_no_redundant_beta,
    feature_importance = data_list$genus_regression_no_redundant_feature_importance,
    metrics = data_list$genus_regression_no_redundant_metrics
  ),
  "Genus (No Latent)" = list(
    beta = data_list$genus_regression_no_latent_beta,
    feature_importance = data_list$genus_regression_no_latent_feature_importance,
    metrics = data_list$genus_regression_no_latent_metrics
  ),
  "Pathway" = list(
    beta = data_list$pathway_regression_beta,
    feature_importance = data_list$pathway_regression_feature_importance,
    metrics = data_list$pathway_regression_metrics
  ),
  "Pathway (No Redundant)" = list(
    beta = data_list$pathway_regression_no_redundant_beta,
    feature_importance = data_list$pathway_regression_no_redundant_feature_importance,
    metrics = data_list$pathway_regression_no_redundant_metrics
  ),
  "Pathway (No Latent)" = list(
    beta = data_list$pathway_regression_no_latent_beta,
    feature_importance = data_list$pathway_regression_no_latent_feature_importance,
    metrics = data_list$pathway_regression_no_latent_metrics
  ),
  # "Pathway Subset" = list(
  #     beta = data_list$pathway_subset_regression_beta,
  #     feature_importance = data_list$pathway_subset_regression_feature_importance,
  #     metrics = data_list$pathway_subset_regression_metrics
  # ),
  # "Pathway Subset (No Redundant)" = list(
  #     beta = data_list$pathway_subset_regression_no_redundant_beta,
  #     feature_importance = data_list$pathway_subset_regression_no_redundant_feature_importance,
  #     metrics = data_list$pathway_subset_regression_no_redundant_metrics
  # ),
  "Species" = list(
    beta = data_list$species_regression_beta,
    feature_importance = data_list$species_regression_feature_importance,
    metrics = data_list$species_regression_metrics
  ),
  "Species (No Redundant)" = list(
    beta = data_list$species_regression_no_redundant_beta,
    feature_importance = data_list$species_regression_no_redundant_feature_importance,
    metrics = data_list$species_regression_no_redundant_metrics
  ),
  "Species (No Latent)" = list(
    beta = data_list$species_regression_no_latent_beta,
    feature_importance = data_list$species_regression_no_latent_feature_importance,
    metrics = data_list$species_regression_no_latent_metrics
  )
)

# Process and plot data for each dataset
# lapply(names(datasets), function(dataset_name) {
#     process_and_plot_data(datasets[[dataset_name]], dataset_name)
# })

# In[5]: Now from the features, lets do a heatmap for the top 20 features ----

# Get top 20 features for one to test
top_20_features <- get_top_n_features_all_models(data_list$latent_feature_importance, 20)
print(top_20_features)

# In[9]: Metric R ^ 2 for presentation ----

# lets get the testing R^2 for all the different datasets
# Function to extract metrics and calculate max RÂ²
extract_metrics <- function(dataset_name, datasets) {
  metrics <- datasets[[dataset_name]]$metrics %>%
    filter(DataType == "Test") %>%
    select(Model, R2)
  max_r2 <- max(metrics$R2, na.rm = TRUE)
  list(metrics = metrics, max_r2 = max_r2)
}

# List of dataset names
dataset_names <- c(
  "Genus", "Genus (No Latent)", "Genus (No Redundant)",
  "Pathway", "Pathway (No Latent)", "Pathway (No Redundant)",
  "Species", "Species (No Latent)", "Species (No Redundant)"
)

# Extract metrics and max RÂ² for all datasets
results <- lapply(dataset_names, extract_metrics, datasets = datasets)

# Organize results into a list for each main category
genus_results <- results[1:3]
pathway_results <- results[4:6]
species_results <- results[7:9]

# Calculate overall max RÂ² for each main category
max_r2_genus <- max(sapply(genus_results, function(res) res$max_r2))
max_r2_pathway <- max(sapply(pathway_results, function(res) res$max_r2))
max_r2_species <- max(sapply(species_results, function(res) res$max_r2))

# Calculate global max RÂ² across all categories
max_r2 <- max(max_r2_genus, max_r2_pathway, max_r2_species)

# Prepare data and titles for genus
genus_data_list <- list(genus_results[[1]]$metrics, genus_results[[2]]$metrics, genus_results[[3]]$metrics)
print(genus_data_list)
# rename the models
genus_data_list[[1]]$Model <- c("Lasso", "Ridge", "Elastic Net", "Random Forest", "XGBoost")
genus_data_list[[2]]$Model <- c("Lasso", "Ridge", "Elastic Net", "Random Forest", "XGBoost")
genus_data_list[[3]]$Model <- c("Lasso", "Ridge", "Elastic Net", "Random Forest", "XGBoost")

genus_titles <- c(
  "All Genus + Clinical Variables - Model Testing RÂ²",
  "Only Genus - Model Testing RÂ²",
  "Non Redundant Genus + Clinical Variables - Model Testing RÂ²"
)
# Generate combined genus plot
combined_plot_genus <- create_plots(genus_data_list, max_r2, genus_titles)
pdf("drift_fs/figures/genus_combined_plot.pdf", width = 7, height = 7)
print(combined_plot_genus)
dev.off()

# Prepare data and titles for species
species_data_list <- list(species_results[[1]]$metrics, species_results[[2]]$metrics, species_results[[3]]$metrics)
species_titles <- c(
  "All Species + Clinical Variables - Model Testing RÂ²",
  "Only Species - Model Testing RÂ²",
  "Non Redundant Species + Clinical Variables - Model Testing RÂ²"
)

# rename the models
species_data_list[[1]]$Model <- c("Lasso", "Ridge", "Elastic Net", "Random Forest", "XGBoost")
species_data_list[[2]]$Model <- c("Lasso", "Ridge", "Elastic Net", "Random Forest", "XGBoost")
species_data_list[[3]]$Model <- c("Lasso", "Ridge", "Elastic Net", "Random Forest", "XGBoost")

# Generate combined species plot
combined_plot_species <- create_plots(species_data_list, max_r2, species_titles)
pdf("drift_fs/figures/species_combined_plot.pdf", width = 7, height = 7)
print(combined_plot_species)
dev.off()

# Prepare data and titles for pathway
pathway_data_list <- list(pathway_results[[1]]$metrics, pathway_results[[2]]$metrics, pathway_results[[3]]$metrics)
pathway_titles <- c(
  "All Pathway + Clinical Variables - Model Testing RÂ²",
  "Only Pathway - Model Testing RÂ²",
  "Non Redundant Pathway + Clinical Variables - Model Testing RÂ²"
)

# rename the models
pathway_data_list[[1]]$Model <- c("Lasso", "Ridge", "Elastic Net", "Random Forest", "XGBoost")
pathway_data_list[[2]]$Model <- c("Lasso", "Ridge", "Elastic Net", "Random Forest", "XGBoost")
pathway_data_list[[3]]$Model <- c("Lasso", "Ridge", "Elastic Net", "Random Forest", "XGBoost")

# Generate combined pathway plot
combined_plot_pathway <- create_plots(pathway_data_list, max_r2, pathway_titles)
pdf("drift_fs/figures/pathway_combined_plot.pdf", width = 7, height = 7)
print(combined_plot_pathway)
dev.off()

# In[10]: Plotting the top 5-10 features ----

# Extract top features from each dataset
genus_no_rendundant_features <- extract_top_features(datasets[["Genus (No Redundant)"]])
pathway_features <- extract_top_features(datasets[["Pathway"]])
species_no_rendundant_features <- extract_top_features(datasets[["Species (No Redundant)"]])

# rename the column names
colnames(genus_no_rendundant_features) <- c("Variable", "Random Forest", "Lasso", "Ridge", "Elastic Net",  "XGBoost")
colnames(pathway_features) <- c("Variable", "Random Forest", "Lasso", "Ridge", "Elastic Net", "XGBoost")
colnames(species_no_rendundant_features) <- c("Variable", "Random Forest", "Lasso", "Ridge", "Elastic Net", "XGBoost")

genus_no_rendundant_features <- genus_no_rendundant_features %>%
  mutate(
    Variable = case_when(
      Variable == "Leptin_BL" ~ "Leptin",
      Variable == "rmr_kcald_BL" ~ "Resting Metabolic Rate",
      Variable == "spk_EE_int_kcal_day_BL" ~ "Spontaneous Energy Expenditure",
      Variable == "HOMA_IR_BL" ~ "Homeostasis Model Assessment",
      Variable == "avg_systolic_BL" ~ "Average Systolic Blood Pressure",
      Variable == "Insulin_endo_BL" ~ "Insulin",
      TRUE ~ Variable
    )
  )

pathway_features <- pathway_features %>%
  mutate(
    Variable = case_when(
      Variable == "Leptin_BL" ~ "Leptin",
      Variable == "rmr_kcald_BL" ~ "Resting Metabolic Rate",
      Variable == "spk_EE_int_kcal_day_BL" ~ "Spontaneous Energy Expenditure",
      Variable == "HOMA_IR_BL" ~ "Homeostasis Model Assessment",
      Variable == "avg_systolic_BL" ~ "Average Systolic Blood Pressure",
      Variable == "Insulin_endo_BL" ~ "Insulin",
      TRUE ~ Variable
    )
  )

species_no_rendundant_features <- species_no_rendundant_features %>%
  mutate(
    Variable = case_when(
      Variable == "Leptin_BL" ~ "Leptin",
      Variable == "rmr_kcald_BL" ~ "Resting Metabolic Rate",
      Variable == "spk_EE_int_kcal_day_BL" ~ "Spontaneous Energy Expenditure",
      Variable == "HOMA_IR_BL" ~ "Homeostasis Model Assessment",
      Variable == "avg_systolic_BL" ~ "Average Systolic Blood Pressure",
      Variable == "Insulin_endo_BL" ~ "Insulin",
      TRUE ~ Variable
    )
  )

# Create and save plots for each dataset
create_feature_plot(
  genus_no_rendundant_features,
  "Top 10 Features - Non Redundant Genus + Clinical Variables",
  "drift_fs/figures/genus_no_rendundant_feature_plot.pdf"
)

create_feature_plot(
  pathway_features,
  "Top 10 Features - All Pathways + Clinical Variables",
  "drift_fs/figures/pathway_feature_plot.pdf"
)

create_feature_plot(
  species_no_rendundant_features,
  "Top 10 Features - Non Redundant Species + Clinical Variables",
  "drift_fs/figures/species_no_rendundant_feature_plot.pdf"
)

# In[12]: Plotting the venn diagrams of the top features ----

# Extract top models for each dataset
datasets_names <- c("Genus (No Redundant)", "Pathway", "Species (No Redundant)")
top_models_list <- datasets_names %>%
  set_names() %>%
  map(~ get_top_models(datasets[[.x]]))

# Print the top models for each dataset
print(top_models_list)

# Map for model names and feature importance columns
model_names <- c("lasso_model", "ridge_model", "elastic_net_model")
model_names_map_to <- c("Lasso_Importance", "Ridge_Importance", "Enet_Importance")

# Get top features for each dataset
top_features_list <- datasets_names %>%
  set_names() %>%
  map(get_features_for_dataset)

# Print the top features for each dataset
print(top_features_list)

# get a set of the variable names for model in each dataset
genus_lasso <- top_features_list[["Genus (No Redundant)"]]$lasso_model$Variable
genus_enet <- top_features_list[["Genus (No Redundant)"]]$elastic_net_model$Variable
genus_ridge <- top_features_list[["Genus (No Redundant)"]]$ridge_model$Variable

# make a set from these variables
genus_set <- list(Lasso = genus_lasso, "Elastic Net" = genus_enet, Ridge = genus_ridge)
print(genus_set)
ggvenn(
  data = genus_set,
  fill_color = viridis(3),
  show_percentage = FALSE,
  text_size = 10, set_name_size = 15
)

pdf("drift_fs/figures/genus_venn_diagrams.pdf", width = 10, height = 10)

species_lasso <- top_features_list[["Species (No Redundant)"]]$lasso_model$Variable
species_enet <- top_features_list[["Species (No Redundant)"]]$elastic_net_model$Variable
species_ridge <- top_features_list[["Species (No Redundant)"]]$ridge_model$Variable
species_set <- list(Lasso = species_lasso, "Elastic Net" = species_enet, Ridge = species_ridge)

ggvenn(
  data = species_set,
  fill_color = viridis(3),
  show_percentage = FALSE,
  text_size = 10, set_name_size = 15
) 
pdf("drift_fs/figures/species_venn_diagrams.pdf", width = 10, height = 10)

pathway_lasso <- top_features_list[["Pathway"]]$lasso_model$Variable
pathway_enet <- top_features_list[["Pathway"]]$elastic_net_model$Variable
pathway_ridge <- top_features_list[["Pathway"]]$ridge_model$Variable

pathway_set <- list(Lasso = pathway_lasso, "Elastic Net" = pathway_enet, Ridge = pathway_ridge)

ggvenn(
  data = pathway_set,
  fill_color = viridis(3),
  show_percentage = FALSE,
  text_size = 10, set_name_size = 15
)
pdf("drift_fs/figures/pathway_venn_diagrams.pdf", width = 10, height = 10)