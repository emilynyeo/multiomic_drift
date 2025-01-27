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
                  "outcome_BMI_fnl_BL", "Glucose_BL", "HOMA_IR_BL",
                  "Insulin_endo_BL", "HDL_Total_Direct_lipid_BL", 
                  "LDL_Calculated_BL","Triglyceride_lipid_BL", 
                  "outcome_BMI_fnl_6m","Glucose_6m", "HOMA_IR_6m", 
                  "Insulin_endo_6m", "HDL_Total_Direct_lipid_6m", 
                  "LDL_Calculated_6m", "Triglyceride_lipid_6m", 
                  "outcome_BMI_fnl_12m","Glucose_12m", "HOMA_IR_12m", 
                  "Insulin_endo_12m", "HDL_Total_Direct_lipid_12m",
                  "LDL_Calculated_12m", "Triglyceride_lipid_12m",
                  "BMI (6m-BL)", "BMI (12m-BL)",
                  "Weight (6m-BL)", "Weight (12m-6m)", "Weight (12m-BL)" ,
                  "outcome_wt_fnl_BL", "outcome_wt_fnl_6m", "outcome_wt_fnl_12m"))

meta_deltas <- meta_data_df %>% 
  dplyr::select("subject_id", "randomized_group", "score_std", 
                "cohort_number", "sex", "race", "age",
                "Glucose (6m-BL)", "Glucose (12m-BL)",
                "HOMA-IR (6m-BL)", "HOMA-IR (12m-BL)",
                "Insulin (6m-BL)", "Insulin (12m-BL)",
                "HDL (6m-BL)", "HDL (12m-BL)",
                "LDL (6m-BL)", "LDL (12m-BL)",
                "Triglyceride lipid (6m-BL)", "Triglyceride lipid (12m-BL)",
                "BMI (6m-BL)", "BMI (12m-BL)",
                "Weight (6m-BL)", "Weight (12m-6m)", "Weight (12m-BL)" ,
                "outcome_wt_fnl_BL", "outcome_wt_fnl_6m", "outcome_wt_fnl_12m",
                "outcome_BMI_fnl_BL", "outcome_BMI_fnl_6m", "outcome_BMI_fnl_12m")

meta_path_g_ra_micom_inner <- meta_deltas %>%
  inner_join(path_g_ra_micom_inner, by = c("subject_id" = "subject_id"))

# Full join: Keeps all rows, with NA for non-matching rows
meta_path_g_ra_micom_outer <- meta_deltas %>%
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
save_dir <- "drift_fs/csv/all_omic_processed_data/deltas/"

# check if the directory exists
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

write.csv(g_ra_all_outer, 
          paste0(save_dir, 
                 "jan20_g_ra_all_omics_deltas_outer.csv"), 
          row.names = FALSE)

write.csv(g_ra_all_inner, 
          paste0(save_dir, 
                 "jan20_g_ra_all_omics_deltas_inner.csv"), 
          row.names = FALSE)

# Look at missingness for innerjoined df 
library(mice)
md.pattern(g_ra_all_inner)

# Percentage of missing data by column
colSums(is.na(g_ra_all_inner)) / nrow(g_ra_all_inner) * 100

# Overall percentage of missing data
sum(is.na(g_ra_all_inner)) / (nrow(g_ra_all_inner) * ncol(g_ra_all_inner)) * 100

# Create a correlation matrix for missingness
missing_corr <- cor(is.na(g_ra_all_inner))
print(missing_corr)

library(ComplexHeatmap)
# Visualize missingness as a heatmap
heatmap(as.matrix(is.na(g_ra_all_inner)), 
        Rowv = NA, Colv = NA, scale = "none")

# Install and load naniar
install.packages("naniar")
library(naniar)
# Visualize missingness
dev.off() 
gg_miss_var(g_ra_all_inner)  # Missing values by variable
vis_miss(g_ra_all_inner)     # Overall missingness

gg_miss_var(g_ra_all_inner, subject_id, 
            show_pct = TRUE) + ylim(0, 100)

ggsave("drift_fs/csv/all_omic_processed_data/deltas/missing_values_plot_misvis.png", 
       plot = vis_miss(g_ra_all_inner[, 1:30]))

# Look at missingness for outer joined df
ggsave("drift_fs/csv/all_omic_processed_data/deltas/outer_missing_values_plot_misvis.png", 
       plot = vis_miss(g_ra_all_outer[, 1:30]))

ggsave("drift_fs/csv/all_omic_processed_data/deltas/outer_missing_values_plot.png", 
       plot = gg_miss_var(g_ra_all_outer[, 1:30]))

### Do the same for normal 
data_dir <- "drift_fs/csv/all_omic_processed_data/"
omic_g_ra_outer1 <- read_csv(paste0(data_dir, "jan18_genus_ra_all_omics_outer.csv"))
omic_g_ra_inner1 <- read_csv(paste0(data_dir, "jan18_genus_ra_all_omics_inner.csv"))

ggsave("drift_fs/csv/all_omic_processed_data/outer_missing_values_plot_misvis.png", 
       plot = vis_miss(omic_g_ra_outer1[, 1:30]))

ggsave("drift_fs/csv/all_omic_processed_data/outer_missing_values_plot.png", 
       plot = gg_miss_var(omic_g_ra_outer1[, 1:30]))
