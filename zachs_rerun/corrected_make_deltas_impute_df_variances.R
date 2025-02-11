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
p_load(tools, reticulate, viridis, tidyplots, patchwork, jsonlite, maps, ggvenn, 
       caret, caretEnsemble, readr, plyr, dplyr, tidyr, purrr, tibble, stringr, 
       psych, randomForest, glmnet, xgboost, ggplot2, reshape2, scales, 
       gridExtra, plotly, sf, tidyverse)

###############################
###     Caret Analysis
###############################
# In[2] Load Data sets ----
data_dir <- "drift_fs/csv/all_omic_processed_data/deltas/"
omic_g_ra_outer <- read_csv(paste0(data_dir, "jan20_g_ra_all_omics_deltas_outer.csv"))
omic_g_ra_inner <- read_csv(paste0(data_dir, "jan20_g_ra_all_omics_deltas_inner.csv"))

deltas_BL_6m <- read_csv(paste0(data_dir, "jan22_all_delta_BL_6m.csv"))
deltas_6m_12m <- read_csv(paste0(data_dir, "jan22_all_delta_6m_12m.csv"))

metadata <- read_csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/zachs_rerun/unprocessed_input/DRIFT_working_dataset_meta_deltas_filtered_05.21.2024.csv")

# Use gsub to modify column names
colnames(metadata) <- gsub(" \\(6m-BL\\)$", "_BL_6m", colnames(metadata))

# Create a regular expression to match the desired column name suffixes
pattern <- "(.*_6m$|.*_12m$|.*_id$|group$|sex$|race$|age$|consent$)"
meta_filtered <- metadata[, grep(pattern, colnames(metadata))] # Keep only columns whose names match the pattern

meta_filtered$CRP_6m_12m <- meta_filtered$C_Reactive_Protein_12m - meta_filtered$C_Reactive_Protein_6m
meta_filtered$cholesterol_6m_12m <- meta_filtered$Cholesterol_lipid_12m - meta_filtered$Cholesterol_lipid_6m
meta_filtered$ghrelin_6m_12m <- meta_filtered$Ghrelin_12m - meta_filtered$Ghrelin_6m
meta_filtered$HDL_6m_12m <- meta_filtered$HDL_Total_Direct_lipid_12m - meta_filtered$HDL_Total_Direct_lipid_6m
meta_filtered$LDL_6m_12m <- meta_filtered$LDL_Calculated_12m - meta_filtered$LDL_Calculated_6m
meta_filtered$HbA1C_6m_12m <- meta_filtered$Hemoglobin_A1C_12m - meta_filtered$Hemoglobin_A1C_6m
meta_filtered$insulin_6m_12m <- meta_filtered$Insulin_endo_12m - meta_filtered$Insulin_endo_6m
meta_filtered$leptin_6m_12m <- meta_filtered$Leptin_12m - meta_filtered$Leptin_6m
meta_filtered$peptide_yy_6m_12m <- meta_filtered$Peptide_YY_12m - meta_filtered$Peptide_YY_6m
meta_filtered$tgcyd_6m_12m <- meta_filtered$Triglyceride_lipid_12m - meta_filtered$Triglyceride_lipid_6m
meta_filtered$homo_ir_6m_12m <- meta_filtered$HOMA_IR_12m - meta_filtered$HOMA_IR_6m
meta_filtered$BMI_6m_12m <- meta_filtered$outcome_BMI_fnl_12m - meta_filtered$outcome_BMI_fnl_6m
meta_filtered$Weight_6m_12m <- meta_filtered$outcome_wt_fnl_12m - meta_filtered$outcome_wt_fnl_6m

metadata$test_bmi_12m_BL <- metadata$outcome_BMI_fnl_12m - metadata$outcome_BMI_fnl_BL

meta_6m_12m <- meta_filtered %>% 
  select(c("subject_id", 
           "randomized_group", "consent", "age", "sex", "race",
           "CRP_6m_12m", "cholesterol_6m_12m", "ghrelin_6m_12m",
           "HDL_6m_12m", "LDL_6m_12m", "HbA1C_6m_12m",
           "insulin_6m_12m", "leptin_6m_12m", "peptide_yy_6m_12m",
           "tgcyd_6m_12m", "homo_ir_6m_12m",
           "BMI_6m_12m", "Weight_6m_12m"))
meta_BL_6m <- meta_filtered %>% 
  select(c("subject_id", "randomized_group", "consent",
           "age", "sex", "race",
           "BMI_BL_6m", "Weight_BL_6m", 
           "C Reactive Protein_BL_6m", "Cholesterol lipid_BL_6m",
           "Ghrelin_BL_6m","HDL_BL_6m", "LDL_BL_6m", "HbA1C_BL_6m", 
           "Insulin_BL_6m", "Leptin_BL_6m", "Peptide YY_BL_6m", 
           "Triglyceride lipid_BL_6m", "HOMA-IR_BL_6m")) %>% 
  rename(CRP_BL_6m = "C Reactive Protein_BL_6m",
         cholesterol_BL_6m = "Cholesterol lipid_BL_6m",
         ghrelin_BL_6m = "Ghrelin_BL_6m",
         insulin_BL_6m = "Insulin_BL_6m",
         leptin_BL_6m = "Leptin_BL_6m",
         peptide_yy_BL_6m = "Peptide YY_BL_6m",
         tgcyd_BL_6m = "Triglyceride lipid_BL_6m",
         homo_ir_BL_6m = "HOMA-IR_BL_6m")

# In[4] Main Analysis ----

# FIRST JUST WITH OUTER JOINED 
omic_g_ra <- omic_g_ra_outer

# Assuming omic_g_ra is your data frame
df_BL <- omic_g_ra[grep("\\.BL$", omic_g_ra$SampleID), ]
df_6m <- omic_g_ra[grep("\\.6m$", omic_g_ra$SampleID), ]
df_12m <- omic_g_ra[grep("\\.12m$", omic_g_ra$SampleID), ]

# Trim them before Merging
df_6m_trim <- df_6m[, c(1, 31:1057)]
df_12m_trim <- df_12m[, c(1, 31:1057)]

# Add "_6m" to the column names of columns 31:1057
colnames(df_BL)[31:1057] <- paste0(colnames(df_BL)[31:1057], "_BL")
colnames(df_6m_trim)[31:1057] <- paste0(colnames(df_6m_trim)[31:1057], "_6m")
colnames(df_12m_trim)[31:1057] <- paste0(colnames(df_12m_trim)[31:1057], "_12m")

# merge
df_BL_6m <- merge(df_BL, df_6m_trim, by.x = "subject_id", by.y = "subject_id", all = TRUE)
df_BL_6m_12m <- merge(df_BL_6m, df_12m_trim, by.x = "subject_id", by.y = "subject_id", all = TRUE)

colnames(df_BL_6m_12m)
tail(colnames(df_BL_6m_12m), 100)

colnames_with_BL <- grep("_BL$", colnames(df_BL_6m_12m), value = TRUE)
print(colnames_with_BL)

colnames_with_6m <- grep("_6m$", colnames(df_BL_6m_12m), value = TRUE)
print(colnames_with_6m)

colnames_with_12m <- grep("_12m$", colnames(df_BL_6m_12m), value = TRUE)
print(colnames_with_12m)

# Make the deltas columns
# Get column names for each time point
BL_cols <- grep("_BL$", colnames(df_BL_6m_12m), value = TRUE)
sixm_cols <- grep("_6m$", colnames(df_BL_6m_12m), value = TRUE)
twelvem_cols <- grep("_12m$", colnames(df_BL_6m_12m), value = TRUE)

# Extract the base names (prefixes) for alignment
BL_base <- gsub("_BL$", "", BL_cols)
sixm_base <- gsub("_6m$", "", sixm_cols)
twelvem_base <- gsub("_12m$", "", twelvem_cols)

# Create new columns for differences
for (base in BL_base) {
  # BL - 6m
  if (base %in% sixm_base) {
    df_BL_6m_12m[[paste0(base, "_BL_6m")]] <- df_BL_6m_12m[[paste0(base, "_6m")]] - df_BL_6m_12m[[paste0(base, "_BL")]]
  }
  
  # BL - 12m
  if (base %in% twelvem_base) {
    df_BL_6m_12m[[paste0(base, "_BL_12m")]] <- df_BL_6m_12m[[paste0(base, "_12m")]] - df_BL_6m_12m[[paste0(base, "_BL")]]
  }
}

for (base in sixm_base) {
  # 6m - 12m
  if (base %in% twelvem_base) {
    df_BL_6m_12m[[paste0(base, "_6m_12m")]] <- df_BL_6m_12m[[paste0(base, "_12m")]] - df_BL_6m_12m[[paste0(base, "_6m")]]
  }
}
# debugging the last one 
common_bases <- intersect(sixm_base, twelvem_base)
print(common_bases)

# Check if the _6m and _12m columns have data for one base
count(is.na(df_BL_6m_12m[[paste0(common_bases[1], "_6m")]]))
count(is.na(df_BL_6m_12m[[paste0(common_bases[1], "_12m")]]))

# first check if there are NAs, then make delta 
df_BL_6m_12m[[paste0(base, "_6m_12m")]] <- 
  ifelse(is.na(df_BL_6m_12m[[paste0(base, "_6m")]]) | is.na(df_BL_6m_12m[[paste0(base, "_12m")]]),
         NA,
         df_BL_6m_12m[[paste0(base, "_12m")]] - df_BL_6m_12m[[paste0(base, "_6m")]])


# Get the names of the new columns
new_columns <- grep("(_BL_6m|_BL_12m|_6m_12m)$", colnames(df_BL_6m_12m), value = TRUE)
print(df_BL_6m_12m[, new_columns])

BL_6m_cols <- grep("(_BL_6m)$", colnames(df_BL_6m_12m), value = TRUE)
BL_12m_cols <- grep("(_BL_12m)$", colnames(df_BL_6m_12m), value = TRUE)
m6_12m_cols <- grep("(_6m_12m)$", colnames(df_BL_6m_12m), value = TRUE)

# separate new delta dfs
# Select columns from BL_6m_cols and columns 1-30
all_delta_BL_6m <- df_BL_6m_12m[, c(1:30, 
                                    which(colnames(df_BL_6m_12m) %in% 
                                            BL_6m_cols))]

all_delta_BL_12m <- df_BL_6m_12m[, c(1:30, 
                                     which(colnames(df_BL_6m_12m) %in% 
                                             BL_12m_cols))]

all_delta_6m_12m <- df_BL_6m_12m[, c(1:30, 
                                     which(colnames(df_BL_6m_12m) %in% 
                                             m6_12m_cols))]

#vis_miss(all_delta_BL_6m)
#vis_miss(all_delta_BL_12m)
#vis_miss(all_delta_6m_12m)

# In[4] Main Analysis ----
colnames(all_delta_BL_6m[, c(1, 1:35)])

### Make dfs 
BL_6m <- all_delta_BL_6m 
BL_6m <- BL_6m[, !grepl(" \\(12m-BL\\)$", colnames(BL_6m))]
BL_6m <- BL_6m[, !grepl(" \\(12m-6m\\)$", colnames(BL_6m))]

BL_12m <- all_delta_BL_12m 
BL_12m <- BL_12m[, !grepl(" \\(6m-BL\\)$", colnames(BL_12m))]
BL_12m <- BL_12m[, !grepl(" \\(12m-6m\\)$", colnames(BL_12m))]

m6_12m <- all_delta_6m_12m 
m6_12m <- m6_12m[, !grepl(" \\(6m-BL\\)$", colnames(m6_12m))]
m6_12m <- m6_12m[, !grepl(" \\(12m-BL\\)$", colnames(m6_12m))]

#save_dir <- "drift_fs/csv/all_omic_processed_data/deltas/"
#write.csv(all_delta_BL_6m, 
#          paste0(save_dir, "jan22_all_delta_BL_6m.csv"), row.names = FALSE)

#write.csv(all_delta_BL_12m, 
#          paste0(save_dir, "jan22_all_delta_BL_12m.csv"), row.names = FALSE)

#write.csv(all_delta_6m_12m, 
#          paste0(save_dir, "jan22_all_delta_6m_12m.csv"), row.names = FALSE)

# Make a delta long dataframe

# Step 1: Add 'range' column to each data frame
new_BL_6m <- BL_6m 
new_BL_6m$range <- "1" # BL_6m
new_m6_12m <- m6_12m
new_m6_12m$range <- "2" # m6_12m

# Step 2: remove _BL_6 and _6m_12m from each of the data frame columns 

# Keep only certain columns
new_BL_6m_filtered <- new_BL_6m[, c(1, 24:ncol(new_BL_6m))]
new_m6_12m_filtered <- new_m6_12m[, c(1, 17:ncol(new_m6_12m))]

# add meta data
new_BL_6m_meta <- merge(meta_BL_6m, new_BL_6m_filtered, by = "subject_id")
new_m6_12m_meta <- merge(meta_6m_12m, new_m6_12m_filtered, by = "subject_id")

# Assuming m6_12m is your data frame
# Get column names from the 17th column onward
colnames(new_BL_6m_meta)[2:ncol(new_BL_6m_meta)] <- gsub("_BL_6m", "", 
                                                         colnames(new_BL_6m_meta)[2:ncol(new_BL_6m_meta)])
colnames(new_m6_12m_meta)[2:ncol(new_m6_12m_meta)] <- gsub("_6m_12m", "", 
                                                           colnames(new_m6_12m_meta)[2:ncol(new_m6_12m_meta)])

# Step 3: combine dataframes 
setdiff(colnames(new_BL_6m_meta), colnames(new_m6_12m_meta))
combined_data <- rbind(new_BL_6m_meta, new_m6_12m_meta)

# Look at missingness of combined data. 
library(naniar)
library(VIM)
library(ggplot2)

gg_miss_var(combined_data[1:101])
mis_sum <- miss_var_summary(combined_data)
mis_sum

# Heatmap using VIM
aggr(combined_data[1:202], 
     col=c("navyblue", "yellow"), 
     numbers=TRUE, sortVars=TRUE, 
     labels=names(combined_data[1:202]), cex.axis=0.7)

vis_miss(combined_data[1:101])
vis_miss(combined_data[505:605])
vis_miss(combined_data[707:808])

# Visualizing correlations in missingness
matrix_miss <- apply(combined_data[7:19], 2, function(x) is.na(x))
corrplot::corrplot(cor(matrix_miss, 
                       use="pairwise.complete.obs"), 
                   method="circle", 
                   type="upper")

# Look at row missingness
# Calculate the number of missing values per row
missing_per_row <- apply(combined_data, 1, function(x) sum(is.na(x)))
threshold <- ncol(combined_data) * 0.9 # % missing
rows_with_more_than_x_missing <- sum(missing_per_row > threshold)
# Percentage of rows with more than 90% missing data
percentage_missing <- (rows_with_more_than_x_missing / nrow(combined_data)) * 100
percentage_missing

# Filter out rows where more than 90% of the data is missing (29% of rows)
combined_data_cleaned <- combined_data[missing_per_row <= threshold, ]
vis_miss(combined_data_cleaned[1:202])

count(unique(combined_data$subject_id)) # 156
count(unique(combined_data_cleaned$subject_id)) # 114 (-42 humans)

# Separate into test and train sets 
m1_dir = "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/clinical/transformed/aim2/"
test = read_csv(paste0(m1_dir, "a2_test_samples_standard_clinical.csv"))
train = read_csv(paste0(m1_dir, "a2_train_samples_standard_clinical.csv"))

count(intersect(test$subject_id, combined_data$subject_id)) # 24
count(intersect(train$subject_id, combined_data$subject_id)) # 96

combined_train <- combined_data %>% filter(subject_id %in% train$subject_id)
combined_test <- combined_data %>% filter(subject_id %in% test$subject_id)

count(unique(combined_train$subject_id)) # 96 
count(unique(combined_test$subject_id)) # 24

table(combined_train$consent)
table(combined_test$consent)

train <- combined_train %>% 
  select(-c("consent")) %>% 
  mutate(range = as.numeric(range))

test <- combined_test %>% 
  select(-c("consent")) %>% 
  mutate(range = as.numeric(range))

########## Check variances and impute ############

# Look at variance of train set 
numeric_train <- train[, sapply(train, is.numeric)]
variances <- apply(numeric_train, 2, var, na.rm = TRUE)

variance_df <- data.frame(
  variable = names(variances),
  variance = variances)

variance_df_sorted <- variance_df[order(-variance_df$variance), ]
top_50 <- variance_df_sorted[1:50, ]
bottom_100 <- variance_df_sorted[(nrow(variance_df_sorted) - 100):nrow(variance_df_sorted), ]

### Plot the variance of each variable
ggplot(bottom_100, aes(x = reorder(variable, variance), y = variance)) +
  geom_bar(stat = "identity", fill = "lightgreen") +
  coord_flip() +  # Flip the coordinates for better readability
  labs(title = "Variance of Variables",
       x = "Variable",
       y = "Variance") +
  theme_minimal()

### Remove variables with zero variance
non_zero_variance_df <- variance_df[variance_df$variance != 0, ]
train_no_zero_var <- train[, c("subject_id", non_zero_variance_df$variable)]

gg_miss_var(train_no_zero_var[1:101])
mis_sum <- miss_var_summary(train_no_zero_var[1:101])
mis_sum

# impute missing of test set 
stardardizing_cols <- colnames(train_no_zero_var[9:920])
train_imputed <- preprocess_data(train_no_zero_var, 
                                 stardardizing_cols, 
                                 "medianImpute")

# Look at variance of test set 
numeric_test <- test[, sapply(test, is.numeric)]
variances_test <- apply(numeric_test, 2, var, na.rm = TRUE)

variance_df_test <- data.frame(
  variable = names(variances_test),
  variance = variances_test)

variance_df_sorted_test <- variance_df_test[order(-variance_df_test$variance), ]
top_50_test <- variance_df_sorted_test[1:50, ]
bottom_100_test <- variance_df_sorted_test[(nrow(variance_df_sorted_test) - 100):nrow(variance_df_sorted_test), ]

### Plot the variance of each variable
ggplot(bottom_100_test, aes(x = reorder(variable, variance), y = variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Flip the coordinates for better readability
  labs(title = "Variance of bottom 100 Variables in test set",
       x = "Variable",
       y = "Variance") +
  theme_minimal()

### Remove variables with zero variance
non_zero_variance_df_test <- variance_df_test[variance_df_test$variance != 0, ]
test_no_zero_var <- test[, c("subject_id", non_zero_variance_df_test$variable)]

gg_miss_var(test_no_zero_var[1:101])
mis_sum <- miss_var_summary(test_no_zero_var[1:101])
mis_sum

# impute missing of test set 
stardardizing_cols_test <- colnames(test_no_zero_var[8:825])
test_imputed_novar <- preprocess_data(test_no_zero_var, 
                                stardardizing_cols_test, 
                                "medianImpute")

test_imputed <- preprocess_data(test, 
                                stardardizing_cols_test, 
                                "medianImpute")

setdiff(colnames(test_imputed), colnames(train_imputed))
diff_cols <- setdiff(colnames(test_imputed), colnames(train_imputed))
diff_cols_train <- setdiff(colnames(train_imputed), colnames(test_imputed))

dim(train_imputed)
dim(test_imputed)

# Identify the common column names between 'train_imputed' and 'test_imputed'
common_cols <- intersect(colnames(train_imputed), colnames(test_imputed))
train_imputed_trimmed <- train_imputed[, common_cols]
test_imputed_trimmed <- test_imputed[, common_cols]


#save_dir <- "drift_fs/csv/all_omic_processed_data/deltas/"
#write.csv(train_imputed_trimmed, 
#          paste0(save_dir, "jan30_all_delta_train_imp_varcheck.csv"), 
#          row.names = FALSE)

#write.csv(test_imputed_trimmed, 
#          paste0(save_dir, "jan30_all_delta_test_imp_varcheck.csv"), 
#          row.names = FALSE)
