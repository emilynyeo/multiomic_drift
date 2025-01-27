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
latent_variables_BL_6m <- c(
  "randomized_group", "score_std", "cohort_number", "sex", "race", "age", 
  "Glucose (6m-BL)", "HOMA-IR (6m-BL)", "Insulin (6m-BL)", 
  "HDL (6m-BL)", "LDL (6m-BL)", "Triglyceride lipid (6m-BL)",           
  "outcome_BMI_fnl_BL")

latent_variables_BL_12m <- c(
  "randomized_group", "score_std", "cohort_number", "sex", "race", "age", 
  "Glucose (12m-BL)", "HOMA-IR (12m-BL)", "Insulin (12m-BL)", 
  "HDL (12m-BL)", "LDL (12m-BL)", "Triglyceride lipid (12m-BL)",           
  "BMI (6m-BL)", "BMI (12m-BL)", 
  "Weight (6m-BL)", "Weight (12m-BL)")


### Process DFs 
imputed_BL_6m <- preprocess_data(BL_6m, 
                              latent_variables_BL_6m, 
                              "medianImpute")

imputed_BL_12m <- preprocess_data(BL_12m, 
                                 latent_variables_BL_12m, 
                                 "medianImpute")

# remove outcome_BMI_fnl_BL from the genus and species dataframes
tail(colnames(imputed_6m), 300)

imputed <- remove_columns(imputed_BL_6m, 
                          columns_to_remove = "subject_id", 
                          "SampleID")

set.seed(123)
train_control <- trainControl(method = "cv", number = 5, search = "grid")






