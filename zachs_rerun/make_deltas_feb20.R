#' @author Emily Yeo
#' @email emily.yeo@colorado.edu
#' @purpose Analysis for the Stanislawski Lab 
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
###     Data Preprocessing
###############################

# In[2]: Data Imports ----

# Define the base path for the data files
output_dir <- "drift_fs/csv/unprocessed_data"
zc_pl_dir <- "unprocessed_input/"
local_path <- "drift_fs/csv/unprocessed_data/"
data_dir <- "drift_fs/csv/processed_data/"
df_dir = "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/merf_dfs/5.combined/"

gen_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/genetic/"
tax_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/taxa/aim2_transformed/genus/"
path_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/functional/aim2/"
micom_dir = "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/micom/aim2/"
m1_dir = "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/clinical/transformed/aim2/"

# Genetic info
gen_all <- read_csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/genetic/genetic_risk.csv")
gen_test <- read_csv(paste0(gen_dir, "genetic_risk_testing.csv"))
gen_train <- read_csv("~/projects/research/Stanislawski/comps/mutli-omic-predictions/data/genetic/genetic_risk_training.csv")

gen_all <- gen_all %>% dplyr::select(c("subject_id", "bmi_prs"))

# Taxa info
genus_all <- read_csv(paste0(tax_dir, "clr_taxa_all_feb20.csv"))
genus_test <- read_csv(paste0(tax_dir, "aim2_clr_testing_feb20.csv"))
genus_train <- read_csv(paste0(tax_dir, "aim2_clr_training_feb20.csv"))

#genus_all <- make_new_columns(genus_all, genus_all$SampleID)
#genus_test <- make_new_columns(genus_test, genus_test$SampleID)
#genus_train <- make_new_columns(genus_train, genus_train$SampleID)

genus_all <- rename_columns_species_to_domain(genus_all) %>% 
  dplyr::rename("sample_id" = "...1")
genus_test <- rename_columns_species_to_domain(genus_test) # rename the columns
genus_train <- rename_columns_species_to_domain(genus_train)

# Meta data 
merge_metadata <- read_csv(paste0(zc_pl_dir, "merge_meta_methyl.csv"))
metadata <- read_csv(paste0(zc_pl_dir, "DRIFT_working_dataset_meta_deltas_filtered_05.21.2024.csv"))
full_raw = read_csv(paste0(m1_dir, "a2_meta_not_Transformed_standard_clinical.csv"))

meta_all <- read_csv(paste0(m1_dir, "a2_meta_Transformed_standard_clinical_feb20.csv"))
meta_test <- read_csv(paste0(m1_dir, "a2_test_samples_standard_clinical_feb20.csv"))
meta_train <- read_csv(paste0(m1_dir, "a2_train_samples_standard_clinical_feb20.csv"))
  
# MICOM data
micom_all = read_csv(paste0(micom_dir, "flux_all_feb20.csv"))
micom_test = read_csv(paste0(micom_dir, "flux_all_testing_feb20.csv"))
micom_train = read_csv(paste0(micom_dir, "flux_all_training_feb20.csv"))

micom_all$id <- micom_all$sample_id
micom_all <- micom_all %>%
  separate(sample_id, into = c("subject_id", "time"), sep = "\\.")

# Pathway 
path_all <- read_csv(paste0(path_dir, "clr_taxa_all_feb20.csv"))
path_test <- read_csv(paste0(path_dir, "all_clr_testing_feb20.csv"))
path_train <- read_csv(paste0(path_dir, "all_clr_training_feb20.csv"))

path_all$id <- path_all$SampleID
path_all <- path_all %>%
  separate(SampleID, 
           into = c("subject_id", 
                    "time"), sep = "\\.")

rm(micom_test, micom_train, path_test, path_train, meta_test, meta_train, full_raw,
   genus_test, genus_train, gen_test, gen_train)

# In[4.2]: Merge the updated_analysis and metadata ----

# Ensure both datasets have 'record_id' for the join
meta_data <- merge_data(gen_all, meta_all, 
                        inner_join, "subject_id")

meta_full <- gen_all %>% full_join(meta_all, by = "subject_id")

### START MERGING GENUS RA with MICOM

# Full join: Keeps all rows, with NA for non-matching rows
meta_full_path <- meta_full %>% full_join(path_all, by = "subject_id")

g_ra_micom_path <- meta_full_path %>%
  full_join(micom_all, by = "id")

g_ra_micom_path <- g_ra_micom_path %>%
  full_join(genus_all, by = c("id" = "sample_id"))

### START MERGING GENUS RA & MICOM with Pathway data
meta_slim <- metadata %>% 
  dplyr::select("subject_id",
                "Glucose (6m-BL)", "Glucose (12m-BL)",
                "HOMA-IR (6m-BL)", "HOMA-IR (12m-BL)",
                "Insulin (6m-BL)", "Insulin (12m-BL)",
                "HDL (6m-BL)", "HDL (12m-BL)",
                "LDL (6m-BL)", "LDL (12m-BL)",
                "Triglyceride lipid (6m-BL)", "Triglyceride lipid (12m-BL)",
                "BMI (6m-BL)", "BMI (12m-BL)",
                "Weight (6m-BL)", "Weight (12m-6m)", "Weight (12m-BL)" ,
                "outcome_wt_fnl_BL", "outcome_wt_fnl_6m", "outcome_wt_fnl_12m")

# Full join: Keeps all rows, with NA for non-matching rows
meta_path_g_ra_micom <- meta_slim %>%
  full_join(g_ra_micom_path, by = c("subject_id" = "subject_id.x"))

#rm(meta_data_df, meta_data, merge_meta_data, metadata, merge_metadata)

# In[4.4]: Remove the columns that are not needed ----
colnames(meta_path_g_ra_micom)
tail(colnames(meta_path_g_ra_micom), 300)

# Define columns and pattern to remove
columns_to_remove <- c(
  "...1.y", "...1.x", "id",
  "subject_id.y", "time.y"
)

pattern_to_remove <- "3m|6m|12m|18m"
BL_pattern <- "3m|6m|12m|18m"
pattern_3m <- "BL|6m|12m|18m"
pattern_6m <- "3m|BL|12m|18m"
pattern_12m <- "3m|6m|BL|18m"

# Remove columns from genus dataset
g_ra_all_BL <- remove_columns(meta_path_g_ra_micom, 
                              columns_to_remove = columns_to_remove, 
                              pattern = BL_pattern)

meta_path_g_ra_micom <- meta_path_g_ra_micom %>% 
  dplyr::select(-c("subject_id.y", "...1", "time.y",
                   "...1.y", "...1.x",
                   "time.x"))

# Remove rows that have 3 and 18m
meta_path_g_ra_micom_fil <- meta_path_g_ra_micom %>%
  filter(!grepl("\\.(3m|18m)$", id))

# save these dataframes
save_dir <- "drift_fs/csv/all_omic_processed_data/deltas/"

# check if the directory exists
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

write.csv(meta_path_g_ra_micom_fil, 
          paste0(save_dir, 
                 "feb20_g_ra_all_omics_deltas_outer.csv"), 
          row.names = FALSE)

# Look at missingness for innerjoined df 
library(mice)
md.pattern(meta_path_g_ra_micom_fil)

# Percentage of missing data by column
colSums(is.na(meta_path_g_ra_micom_fil)) / nrow(meta_path_g_ra_micom_fil) * 100

# Overall percentage of missing data
sum(is.na(meta_path_g_ra_micom_fil)) / (nrow(meta_path_g_ra_micom_fil) * ncol(meta_path_g_ra_micom_fil)) * 100

# Create a correlation matrix for missingness
missing_corr <- cor(is.na(meta_path_g_ra_micom_fil))
print(meta_path_g_ra_micom_fil)
