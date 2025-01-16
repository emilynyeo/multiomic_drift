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

library(tools)
library(reticulate)
library(viridis)
library(tidyplots)
library(patchwork)
library(jsonlite)
library(maps)
library(ggvenn)
library(caret)
library(caretEnsemble)
library(readr)
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)
library(psych)
library(randomForest)
library(glmnet)
library(xgboost)
library(ggplot2)
library(reshape2)
library(scales)
library(gridExtra)
library(plotly)

# difficult packages
library(sf)
library(tidyverse)

# In[2]: Functions ----

# Function to load RData file into a new environment and return the environment
load_rdata <- function(file_path) {
    if (file_ext(file_path) == "RData") {
        # Create a new environment to load the RData file into
        env <- new.env()
        load(file_path, envir = env) # Load the RData file into the new environment
        message("RData file loaded successfully.")

        # Return the environment for further inspection
        return(env)
    } else {
        stop("The file is not an RData file.")
    }
}

# Function to automatically save all data frames/matrices in an environment as CSV files
save_env_to_csv <- function(env, output_dir) {
    # Ensure the output directory exists
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    # Loop through the objects in the environment
    for (obj_name in ls(envir = env)) {
        obj <- get(obj_name, envir = env)

        # Check if the object is a data frame or matrix
        if (is.data.frame(obj) || is.matrix(obj)) {
            file_path <- file.path(output_dir, paste0(obj_name, ".csv"))
            write.csv(obj, file = file_path)
            message(paste("Saved:", file_path))
        } else {
            message(paste("Skipping:", obj_name, "as it is not a data frame or matrix."))
        }
    }
}

# In[3]: Taxa Data ----
# Load another RData file for taxa data and save to CSV
genus_tables <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/zachs_rerun/unprocessed_input/Genus_Sp_tables.RData"
env_taxa <- load_rdata(genus_tables)

output_dir <- "drift_fs/csv/unprocessed_data"
save_env_to_csv(env_taxa, output_dir)

###############################
###
###     Data Preprocessing
###
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

make_new_columns <- function(data, column_name) {
    data <- data %>%
        mutate(
            subject_id = str_split(column_name, "\\.", simplify = TRUE)[, 1],
            TIMEPOINT = str_split(column_name, "\\.", simplify = TRUE)[, 2]
        )
    return(data)
}

filter_data <- function(data, column, value) {
    data <- data %>%
        filter(column == value)
    return(data)
}

merge_data <- function(data1, data2, join_type, columnname) {
    data <- join_type(data1, data2, by = columnname)

    # Remove the duplicated columns and rename as necessary
    data <- data %>%
        select(-matches(paste0(columnname, "\\.y$"))) %>%
        rename_with(~ gsub("\\.x$", "", .), ends_with(".x"))

    return(data)
}

remove_columns <- function(data, columns_to_remove = NULL, pattern = NULL) {
    # Remove specified columns if provided
    if (!is.null(columns_to_remove)) {
        data <- data %>% select(-all_of(columns_to_remove))
    }

    # Remove columns matching the pattern if provided
    if (!is.null(pattern)) {
        data <- data %>% select(-matches(pattern))
    }

    return(data)
}

extract_columns <- function(data, columns_to_extract = NULL, pattern = NULL) {
    # Case 1: Both columns_to_extract and pattern are specified
    if (!is.null(columns_to_extract) && !is.null(pattern)) {
        data <- data %>% 
                select(all_of(intersect(names(data), 
                                        columns_to_extract)), 
                       matches(pattern))

        # Case 2: Only pattern is specified
    } else if (is.null(columns_to_extract) && !is.null(pattern)) {
        data <- data %>% select(matches(pattern))

        # Case 3: Only columns_to_extract is specified
    } else if (!is.null(columns_to_extract) && is.null(pattern)) {
        data <- data %>% select(all_of(columns_to_extract))
    }

    # Return data with the selected columns
    return(data)
}

# Python function definition
# py_run_string("
# import re
# 
# order = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
# def rename_columns_species_to_domain(dataframe):
#     # Try to get the lowest level of taxonomy
#     # If the lowest level is not present, then get the next lowest level
# 
#     # Get the columns of the dataframe
#     columns = dataframe.columns
# 
#     # Get all the columns containing all levels of taxonomy (starting with 'd__')
#     species_columns = [column for column in columns if column.startswith('d__')]
# 
#     for column in species_columns:
#         # split by _{any single letter}__
#         split_column = re.split(r'_[a-z]__', column)
#         split_column[0] = split_column[0].split('d__')[1]
# 
#         # remove the elements of the list that are empty
#         split_column = [x for x in split_column if x]
# 
#         # the length of the split_column will be the lowest level of taxonomy
#         order_index = len(split_column) - 1
# 
#         #join the split_column to get the lowest level of taxonomy
#         new_column_name = order[order_index] + split_column[-1]
# 
#         # rename the column
#         dataframe.rename(columns={column: new_column_name}, inplace=True)
# 
#     # Return dataframe without any changes for now
#     return dataframe
# ")

# Access the Python function from R
# rename_columns_species_to_domain <- py$rename_columns_species_to_domain

# chatgpt translation of the python function
rename_columns_species_to_domain <- function(dataframe) {
  # Define the order of taxonomic levels
  order <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__")
  
  # Get all column names
  columns <- colnames(dataframe)
  
  # Filter columns starting with 'd__'
  species_columns <- columns[grepl("^d__", columns)]
  
  # Loop through each species column
  for (column in species_columns) {
    # Split by '_{any single letter}__'
    split_column <- unlist(strsplit(column, "_[a-z]__"))
    split_column[1] <- sub("d__", "", split_column[1]) # Remove 'd__' prefix
    
    # Remove empty elements
    split_column <- split_column[split_column != ""]
    
    # Get the lowest level of taxonomy
    order_index <- length(split_column) - 1
    new_column_name <- paste0(order[order_index + 1], split_column[length(split_column)])
    
    # Rename the column in the data frame
    colnames(dataframe)[colnames(dataframe) == column] <- new_column_name
  }
  
  return(dataframe)
}


# In[4]: Data Preprocessing ----

# In[4.1]: Process genus and species clr data ----

genus_clr_data <- make_new_columns(genus_clr_data, genus_clr_data$SampleID)
genus_clr_data <- filter_data(genus_clr_data, genus_clr_data$TIMEPOINT, "BL")

species_clr_data <- make_new_columns(species_clr_data, species_clr_data$SampleID)
species_clr_data <- filter_data(species_clr_data, species_clr_data$TIMEPOINT, "BL")

# In[4.2]: Merge the updated_analysis and metadata ----

# Ensure both datasets have 'record_id' for the join
meta_data <- merge_data(updated_analysis, metadata %>% select(-subject_id), inner_join, "record_id")

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
###
###     Caret Analysis
###
###############################

# In[2] Load Datasets ----

data_dir <- "drift_fs/csv/processed_data/"

species_df <- read_csv(paste0(data_dir, "species_latent.csv"))
genus_df <- read_csv(paste0(data_dir, "genus_latent.csv"))

# In[3] Functions ----

process_data <- function(data, columns_to_remove, columns_to_standardize, impute_method = "medianImpute") {
    data_cleaned <- remove_columns(data, columns_to_remove)
    data_standardized <- preprocess_data(data_cleaned, columns_to_standardize, impute_method)
    return(data_standardized)
}

preprocess_data <- function(data, columns_to_standardize, imputation_method) {
    data_imputed <- predict(preProcess(data, method = c(imputation_method)), data)
    data_standardized <- predict(preProcess(data_imputed[, columns_to_standardize], method = c("center", "scale")), data_imputed)
    return(data_standardized)
}

train_all_models <- function(data, target, train_control) {
    lasso_model <- train(
        x = as.data.frame(data[, -which(names(data) %in% c(target, "subject_id"))]),
        y = as.numeric(data[[target]]),
        method = "glmnet",
        trControl = train_control,
        tuneGrid = expand.grid(alpha = 1, lambda = seq(0.1, 1, 0.1))
    )

    ridge_model <- train(
        x = as.data.frame(data[, -which(names(data) %in% c(target, "subject_id"))]),
        y = as.numeric(data[[target]]),
        method = "glmnet",
        trControl = train_control,
        tuneGrid = expand.grid(alpha = 0, lambda = seq(0.1, 1, 0.1))
    )

    elastic_net_model <- train(
        x = as.data.frame(data[, -which(names(data) %in% c(target, "subject_id"))]),
        y = as.numeric(data[[target]]),
        method = "glmnet",
        trControl = train_control,
        tuneGrid = expand.grid(alpha = seq(0.1, 1, 0.1), lambda = seq(0.1, 1, 0.1))
    )

    rf_model <- train(
        x = as.data.frame(data[, -which(names(data) %in% c(target, "subject_id"))]),
        y = as.numeric(data[[target]]),
        method = "rf",
        trControl = train_control,
        tuneGrid = expand.grid(mtry = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)),
        ntree = 10000
    )

    xgb_model <- train(
        x = as.data.frame(data[, -which(names(data) %in% c(target, "subject_id"))]),
        y = as.numeric(data[[target]]),
        method = "xgbTree",
        trControl = train_control
    )

    caret.list <- caretList(
        x = as.data.frame(data[, -which(names(data) %in% c(target, "subject_id"))]),
        y = as.numeric(data[[target]]),
        trControl = train_control,
        methodList = c("glmnet", "rf", "xgbTree")
    )

    ens <- caretEnsemble(caret.list)

    return(list(
        lasso_model = lasso_model,
        ridge_model = ridge_model,
        elastic_net_model = elastic_net_model,
        rf_model = rf_model,
        xgb_model = xgb_model,
        ens = ens
    ))
}

extract_importance_df <- function(model, label) {
    importance <- varImp(model)$importance %>%
        as.data.frame() %>%
        rownames_to_column("Variable")
    colnames(importance)[2] <- label
    return(importance)
}

combine_importances <- function(model_list, labels) {
    importance_dfs <- map2(model_list, labels, extract_importance_df)
    reduce(importance_dfs, full_join, by = "Variable")
}

extract_best_betas <- function(model_list, labels) {
    beta_dfs <- map2(model_list, labels, function(model, label) {
        best_lambda <- model$bestTune$lambda
        betas <- as.data.frame(as.matrix(coef(model$finalModel, s = best_lambda))) %>%
            rownames_to_column("Variable")
        colnames(betas)[2] <- label
        return(betas)
    })
    beta_combined <- reduce(beta_dfs, full_join, by = "Variable") %>%
        mutate(across(everything(), ~ replace_na(., 0)))
    return(beta_combined)
}

replace_na <- function(x, value) {
    ifelse(is.na(x), value, x)
}

# Helper function to calculate model performance metrics for either training or testing data
calculate_metrics <- function(model, data, target_var, model_name, data_type) {
    predictions <- predict(model, data)
    actuals <- data[[target_var]]

    r2 <- caret::R2(predictions, actuals)
    mae <- caret::MAE(predictions, actuals)
    rmse <- caret::RMSE(predictions, actuals)

    return(data.frame(Model = model_name, DataType = data_type, R2 = r2, MAE = mae, RMSE = rmse))
}

# Update the train_and_save_models function to include training and testing metric calculation
train_and_save_models <- function(data, target_var, train_control, result_prefix, test_size = 0.3) {
    # Split data into training and testing sets
    set.seed(123) # Ensure reproducibility
    train_indices <- sample(seq_len(nrow(data)), size = (1 - test_size) * nrow(data))
    train_data <- data[train_indices, ]
    test_data <- data[-train_indices, ]

    # Train models
    results <- train_all_models(train_data, target_var, train_control)

    # Combine importances
    feature_importance <- combine_importances(
        list(results$rf_model, results$lasso_model, results$ridge_model, results$elastic_net_model, results$xgb_model),
        c("RF_Importance", "Lasso_Importance", "Ridge_Importance", "Enet_Importance", "XGBoost_Importance")
    )
    # ensure output directory exists
    if (!dir.exists("drift_fs/csv/results/")) {
        dir.create("drift_fs/csv/results/", recursive = TRUE)
    }

    write.csv(feature_importance, paste0("drift_fs/csv/results/", result_prefix, "_feature_importance.csv"), row.names = FALSE)

    # Extract best betas
    beta_coefficients <- extract_best_betas(
        list(results$lasso_model, results$ridge_model, results$elastic_net_model),
        c("Lasso_Beta", "Ridge_Beta", "Enet_Beta")
    )
    write.csv(beta_coefficients, paste0("drift_fs/csv/results/", result_prefix, "_beta.csv"), row.names = FALSE)

    # Initialize an empty DataFrame to store performance metrics
    metrics_df <- data.frame(Model = character(), DataType = character(), R2 = numeric(), MAE = numeric(), RMSE = numeric(), stringsAsFactors = FALSE)

    # Calculate and store metrics for each model on both training and testing data
    for (model_name in names(results[-length(results)])) {
        model <- results[[model_name]]

        # Training metrics
        train_metrics <- calculate_metrics(model, train_data, target_var, model_name, "Train")
        metrics_df <- rbind(metrics_df, train_metrics)

        # Testing metrics
        test_metrics <- calculate_metrics(model, test_data, target_var, model_name, "Test")
        metrics_df <- rbind(metrics_df, test_metrics)
    }

    # Save the metrics DataFrame as CSV
    write.csv(metrics_df, paste0("drift_fs/csv/results/", result_prefix, "_metrics.csv"), row.names = FALSE)

    # ensure the model directory exists
    if (!dir.exists("drift_fs/models/")) {
        dir.create("drift_fs/models/", recursive = TRUE)
    }
    # Save model results
    saveRDS(results, paste0("drift_fs/models/", result_prefix, "_results.rds"))

    return(results)
}

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

# remove the latent variables from the genus and species dataframes except for the target variable
# genus_df_imputed <- remove_columns(genus_df_imputed, latent_variables_to_use[-which(latent_variables_to_use %in% c("outcome_BMI_fnl_BL"))])
# species_df_imputed <- remove_columns(species_df_imputed, latent_variables_to_use[-which(latent_variables_to_use %in% c("diff_std_bmi_score"))])

# In[6] Data Exploration ----

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

# In[2]: Functions ----

# Function to get top N important features for a given model
get_top_n_features <- function(feature_importance, model_importance_column, n = 20) {
    feature_importance %>%
        select(Variable, all_of(model_importance_column)) %>%
        arrange(desc(get(model_importance_column))) %>%
        head(n)
}

# Function to get top N important features for all models in the feature importance dataframe
get_top_n_features_all_models <- function(feature_importance, n = 20) {
    # Identify columns ending with '_Importance'
    importance_columns <- names(feature_importance)[grepl("_Importance$", names(feature_importance))]

    # Initialize a list to store results for each model
    top_features_list <- list()

    # Loop through each importance column and get top N features
    for (column in importance_columns) {
        model_name <- gsub("_Importance", "", column) # Extract model name
        top_features_list[[model_name]] <- get_top_n_features(feature_importance, column, n)
    }

    return(top_features_list)
}

# Function to plot Venn diagram for top features
plot_venn_diagram <- function(feature_sets, colors = NULL, figurename = "venn_diagram.png", output_dir = "drift_fs/figures/") {
    if (is.null(colors)) {
        colors <- viridis(length(feature_sets))
    }

    venn.plot <- venn.diagram(
        x = feature_sets,
        category.names = names(feature_sets),
        filename = NULL, # Plot directly to the object
        output = TRUE,
        col = "transparent",
        fill = colors,
        alpha = 0.3,
        cex = 1.5,
        cat.cex = 1.2,
        cat.pos = 0,
        margin = 0.1
    )

    ggsave(file.path(output_dir, figurename), plot = venn.plot, dpi = 600, bg = "white")
}

# Function to plot feature importance or beta values
plot_importance_or_beta <- function(data, value_column, plot_title, y_label, figurename, palette = "viridis") {
    long_format <- data %>%
        pivot_longer(
            cols = ends_with(value_column),
            names_to = "Model",
            values_to = value_column
        )

    ggplot(long_format, aes(y = reorder(Variable, !!sym(value_column)), x = !!sym(value_column), fill = Model)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
        scale_fill_viridis_d(option = palette) +
        labs(
            title = plot_title,
            x = paste(value_column, "Score"),
            y = y_label
        ) +
        theme_minimal() +
        theme(
            axis.text.y = element_text(size = 5), # Adjust text size
            axis.title.y = element_text(vjust = 1)
        ) +
        scale_y_discrete(expand = expansion(mult = c(0.1, 0.2))) # Add padding between features

    ggsave(file.path("drift_fs/figures/", figurename), dpi = 600, bg = "white")
}

# Function to plot model performance metrics
plot_performance_metrics <- function(metrics, dataset_name) {
    metrics_long <- metrics %>%
        pivot_longer(
            cols = c("R2", "MAE", "RMSE"),
            names_to = "Metric",
            values_to = "Value"
        )

    # Reorder factors for better readability
    metrics_long$DataType <- factor(metrics_long$DataType, levels = c("Train", "Test"))
    metrics_long$Metric <- factor(metrics_long$Metric, levels = c("R2", "MAE", "RMSE"))

    ggplot(metrics_long, aes(x = Model, y = Value, fill = DataType)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
        labs(
            title = paste("Model Performance Metrics (R², MAE, RMSE) -", dataset_name),
            x = "Model",
            y = "Metric Value",
            fill = "Dataset"
        ) +
        scale_fill_viridis(discrete = TRUE) +
        facet_wrap(~Metric, scales = "free_y", nrow = 3) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 10),
            strip.text = element_text(size = 12),
            legend.position = "top",
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12)
        )

    ggsave(file.path("drift_fs/figures/", paste0(dataset_name, "_performance_metrics.png")), dpi = 600, bg = "white")
}

# Function to process data and generate plots for a dataset
process_and_plot_data <- function(data_list, dataset_name, n = 20) {
    beta <- data_list$beta
    feature_importance <- data_list$feature_importance
    metrics <- data_list$metrics

    # Get top N features for all models
    top_20_features <- get_top_n_features_all_models(feature_importance, n)

    # Create feature sets for Venn diagram
    feature_sets <- lapply(top_20_features, function(df) df$Variable)

    # Plot Venn diagram
    plot_venn_diagram(feature_sets, figurename = paste0(dataset_name, "_venn_diagram.png"))

    # Extract top features and sort by total importance
    all_top_features <- unique(unlist(feature_sets))
    filtered_feature_importance <- feature_importance %>%
        filter(Variable %in% all_top_features) %>%
        mutate(Total_Importance = rowSums(select(., ends_with("_Importance")), na.rm = TRUE)) %>%
        arrange(desc(Total_Importance))

    # Remove total importance column
    filtered_feature_importance <- select(filtered_feature_importance, -Total_Importance)

    # Plot Feature Importances
    plot_importance_or_beta(filtered_feature_importance, "Importance", "Top 20 Feature Importances", "Features", paste0(dataset_name, "_feature_importance.png"))

    # Plot model performance metrics
    plot_performance_metrics(metrics, dataset_name)
}

# Function to add points and recalculate regression
generate_plot <- function(x, y, modified_x = NULL, modified_y = NULL) {
    # Create the data frame
    plot_data <- data.frame(x = x, y = y)

    # Add modified points if provided
    if (!is.null(modified_x) && !is.null(modified_y)) {
        modified_data <- data.frame(x = modified_x, y = modified_y)
        plot_data <- rbind(plot_data, modified_data)
    }

    # Fit the linear regression model
    model <- lm(y ~ x, data = plot_data)
    r_squared <- summary(model)$r.squared

    # Generate the base plot
    p <- ggplot(plot_data, aes(x = x, y = y)) +
        geom_point(aes(color = ifelse(x %in% modified_x, "Modified", "Original")), size = 5) +
        scale_color_manual(values = c("Original" = "#1C7C54", "Modified" = "red")) +
        theme_minimal() +
        xlim(0, max(plot_data$x)) +
        ylim(min(plot_data$y) - 2, max(plot_data$y) + 2) +
        theme(
            legend.position = "none",
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.ticks.length = unit(0.3, "cm"), # Adjust tick size
            axis.title.x = element_text(size = 22, face = "bold"), # Larger x-axis title
            axis.title.y = element_text(size = 22, face = "bold") # Larger y-axis title
        )

    # Add regression line and R² annotation
    p1 <- p +
        geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray", alpha = 0.5) +
        annotate("text",
            x = max(plot_data$x) - 2, y = mean(plot_data$y),
            label = paste("R² =", round(r_squared, 2)),
            color = "black", size = 8
        )

    # Return both plots in a list
    return(list(base_plot = p, regression_plot = p1))
}

# Define a function to create plots
create_plots <- function(data_list, max_r2, titles) {
    plots <- lapply(seq_along(data_list), function(i) {
        ggplot(data_list[[i]], aes(x = R2, y = Model, fill = Model)) +
            geom_bar(stat = "identity", position = "dodge") +
            labs(
                title = titles[i],
                x = "R²",
                y = "Model",
                fill = "Model"
            ) +
            coord_cartesian(xlim = c(0, max_r2)) +
            theme_minimal() +
            scale_fill_viridis_d() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
                axis.text.y = element_text(size = 15, angle = 45),
                legend.position = "none",
                axis.title.x = element_text(size = 20),
                axis.title.y = element_text(size = 20),
                plot.title = element_text(size = 20)
            )
    })
    # Combine the plots vertically
    Reduce(`/`, plots)
}

# Updated function to create and save plots
create_feature_plot <- function(features, title, save_path) {
    # Prepare data for plotting, excluding the `Model` column from pivoting
    features_long <- features %>%
        pivot_longer(
            cols = c(-Variable), # Exclude `Variable` and `Model`
            names_to = "Model",
            values_to = "Importance"
        )

    # Create the plot
    feature_plot <- ggplot(features_long, aes(x = Importance, y = reorder(Variable, Importance), fill = Model)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(
            title = title,
            x = "Importance",
            y = "Feature",
            fill = "Model"
        ) +
        scale_fill_viridis_d() +
        theme_minimal() +
        theme(
            axis.text.x = element_text(size = 13),
            axis.text.y = element_text(size = 13),
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15),
            title = element_text(size = 15),
            legend.position = "top"
        )

    # Print and save the plot
    print(feature_plot)
    pdf(file = save_path, width = 10, height = 10)
    print(feature_plot)
    dev.off()
}

# Function to extract top features
extract_top_features <- function(dataset, top_n = 10) {
    # Extract columns ending with '_Importance'
    importance_columns <- names(dataset$feature_importance) %>%
        grep("_Importance$", ., value = TRUE)

    # Check if importance columns exist
    if (length(importance_columns) == 0) {
        stop("No columns ending with '_Importance' found in feature_importance.")
    }

    # get a column with the total importance
    dataset$feature_importance <- dataset$feature_importance %>%
        mutate(Total_Importance = rowSums(select(., ends_with("_Importance")), na.rm = TRUE))
    
    # sort by total importance from high to low
    dataset$feature_importance <- dataset$feature_importance %>%
        arrange(desc(Total_Importance))
    
    # get the top n features
    top_features <- dataset$feature_importance %>%
        select(Variable, Total_Importance) %>%
        head(10)
    
    # remove the total importance column
    dataset$feature_importance <- select(dataset$feature_importance, -Total_Importance) %>% head(10)

    return(dataset$feature_importance)
}
# Function to get top 20 features for each model
get_top_features <- function(dataset, model_name, importance_col) {
    dataset$feature_importance %>%
        select(Variable, all_of(importance_col)) %>%
        arrange(desc(get(importance_col))) %>%
        head(20)
}

# Extract top 20 features for each model and dataset
get_features_for_dataset <- function(dataset_name) {
    dataset <- datasets[[dataset_name]]
    map2(
        model_names, model_names_map_to,
        ~ get_top_features(dataset, .x, .y)
    ) %>%
        set_names(model_names)
}

# Function to get top models based on R²
get_top_models <- function(dataset, top_n = 3) {
    dataset$metrics %>%
        filter(DataType == "Test") %>%
        arrange(desc(R2)) %>%
        slice_head(n = top_n) %>%
        pull(Model)
}

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

# In[6]: world map obesity plot ----

# Fetch the data
df <- read.csv("https://ourworldindata.org/grapher/share-of-adults-defined-as-obese.csv?v=1&csvType=full&useColumnShortNames=true")

# Get the most recent year for each country
latest_data <- df %>%
    group_by(Entity) %>%
    filter(Year == max(Year)) %>%
    ungroup()

# Standardize country names to match map_data
latest_data <- latest_data %>%
    mutate(Entity = case_when(
        Entity == "United States" ~ "USA",
        Entity == "United Kingdom" ~ "UK",
        Entity == "Czechia" ~ "Czech Republic", # Example; add more as needed
        TRUE ~ Entity
    ))

# Load the world map
world <- map_data("world")

# Merge data with the map
map_data <- world %>%
    left_join(latest_data, by = c("region" = "Entity"))

# Find the maximum value in the prevalence column
max_prevalence <- max(
    latest_data$prevalence_of_obesity_among_adults__bmi__gt__30__crude_estimate__pct__sex_both_sexes__age_group_18plus__years,
    na.rm = TRUE
)

# Plot the map with dynamic max value
ggplot(map_data, aes(long, lat,
    group = group,
    fill = prevalence_of_obesity_among_adults__bmi__gt__30__crude_estimate__pct__sex_both_sexes__age_group_18plus__years
)) +
    geom_polygon(color = "black") +
    scale_fill_gradient(
        name = "Obesity Prevalence (%)",
        low = "#eaeaea",
        high = "#1C7C54",
        na.value = "grey50",
        limits = c(0, 40) # Set limits dynamically
    ) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.title = element_blank())

ggsave("drift_fs/figures/world_map_obesity.png", height = 6, width = 12, dpi = 600, bg = "white")

# In[8]: Showing how linear regression works ----

# Step 1: Generate initial data
set.seed(123)
x <- seq(1, 20, by = 1) # Independent variable
y <- 3 * x + rnorm(length(x), mean = 0, sd = 10) # Dependent variable with noise

# Step 2: Display the original plot
original_plots <- generate_plot(x, y)

# Print the base plot and regression plot
print(original_plots$base_plot)
print(original_plots$regression_plot)

set.seed(456) # Ensure reproducibility
modified_x <- runif(20, min = 0, max = 20) # 20 random points for x
modified_y <- runif(20, min = 0, max = 70) # 20 random points for y

# Step 4: Display updated plot with modified points
updated_plots <- generate_plot(x, y, modified_x, modified_y)

# Print the updated plots
print(updated_plots$base_plot)
print(updated_plots$regression_plot)

# Save the plots
ggsave("drift_fs/figures/linear_regression_example.png", plot = original_plots$base_plot, dpi = 600, bg = "white")
ggsave("drift_fs/figures/linear_regression_example_modified.png", plot = updated_plots$base_plot, dpi = 600, bg = "white")
ggsave("drift_fs/figures/linear_regression_example_regression.png", plot = original_plots$regression_plot, dpi = 600, bg = "white")
ggsave("drift_fs/figures/linear_regression_example_modified_regression.png", plot = updated_plots$regression_plot, dpi = 600, bg = "white")

# In[9]: Metric R ^ 2 for presentation ----

# lets get the testing R^2 for all the different datasets
# Function to extract metrics and calculate max R²
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

# Extract metrics and max R² for all datasets
results <- lapply(dataset_names, extract_metrics, datasets = datasets)

# Organize results into a list for each main category
genus_results <- results[1:3]
pathway_results <- results[4:6]
species_results <- results[7:9]

# Calculate overall max R² for each main category
max_r2_genus <- max(sapply(genus_results, function(res) res$max_r2))
max_r2_pathway <- max(sapply(pathway_results, function(res) res$max_r2))
max_r2_species <- max(sapply(species_results, function(res) res$max_r2))

# Calculate global max R² across all categories
max_r2 <- max(max_r2_genus, max_r2_pathway, max_r2_species)

# Prepare data and titles for genus
genus_data_list <- list(genus_results[[1]]$metrics, genus_results[[2]]$metrics, genus_results[[3]]$metrics)
print(genus_data_list)
# rename the models
genus_data_list[[1]]$Model <- c("Lasso", "Ridge", "Elastic Net", "Random Forest", "XGBoost")
genus_data_list[[2]]$Model <- c("Lasso", "Ridge", "Elastic Net", "Random Forest", "XGBoost")
genus_data_list[[3]]$Model <- c("Lasso", "Ridge", "Elastic Net", "Random Forest", "XGBoost")

genus_titles <- c(
    "All Genus + Clinical Variables - Model Testing R²",
    "Only Genus - Model Testing R²",
    "Non Redundant Genus + Clinical Variables - Model Testing R²"
)
# Generate combined genus plot
combined_plot_genus <- create_plots(genus_data_list, max_r2, genus_titles)
pdf("drift_fs/figures/genus_combined_plot.pdf", width = 7, height = 7)
print(combined_plot_genus)
dev.off()

# Prepare data and titles for species
species_data_list <- list(species_results[[1]]$metrics, species_results[[2]]$metrics, species_results[[3]]$metrics)
species_titles <- c(
    "All Species + Clinical Variables - Model Testing R²",
    "Only Species - Model Testing R²",
    "Non Redundant Species + Clinical Variables - Model Testing R²"
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
    "All Pathway + Clinical Variables - Model Testing R²",
    "Only Pathway - Model Testing R²",
    "Non Redundant Pathway + Clinical Variables - Model Testing R²"
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

# because of this plot we will hence forth use the
# genus no rendundant
# pathway
# species no rendundant

# In[10]: Plotting the top 5-10 features ----

# Load required libraries

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