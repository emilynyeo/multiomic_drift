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
    df_BL_6m_12m[[paste0(base, "_BL_6m")]] <- df_BL_6m_12m[[paste0(base, "_BL")]] - df_BL_6m_12m[[paste0(base, "_6m")]]
  }
  
  # BL - 12m
  if (base %in% twelvem_base) {
    df_BL_6m_12m[[paste0(base, "_BL_12m")]] <- df_BL_6m_12m[[paste0(base, "_BL")]] - df_BL_6m_12m[[paste0(base, "_12m")]]
  }
}

for (base in sixm_base) {
  # 6m - 12m
  if (base %in% twelvem_base) {
    df_BL_6m_12m[[paste0(base, "_6m_12m")]] <- df_BL_6m_12m[[paste0(base, "_6m")]] - df_BL_6m_12m[[paste0(base, "_12m")]]
  }
}
# debugging the last one 
common_bases <- intersect(sixm_base, twelvem_base)
print(common_bases)

# Check if the _6m and _12m columns have data for one base
count(is.na(df_BL_6m_12m[[paste0(common_bases[1], "_6m")]]))
count(is.na(df_BL_6m_12m[[paste0(common_bases[1], "_12m")]]))

df_BL_6m_12m[[paste0(base, "_6m_12m")]] <- 
  ifelse(is.na(df_BL_6m_12m[[paste0(base, "_6m")]]) | is.na(df_BL_6m_12m[[paste0(base, "_12m")]]),
         NA,
         df_BL_6m_12m[[paste0(base, "_6m")]] - df_BL_6m_12m[[paste0(base, "_12m")]])


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

vis_miss(all_delta_BL_6m)
vis_miss(all_delta_BL_12m)
vis_miss(all_delta_6m_12m)

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

save_dir <- "drift_fs/csv/all_omic_processed_data/deltas/"
write.csv(all_delta_BL_6m, 
          paste0(save_dir, "jan22_all_delta_BL_6m.csv"), row.names = FALSE)

write.csv(all_delta_BL_12m, 
          paste0(save_dir, "jan22_all_delta_BL_12m.csv"), row.names = FALSE)

write.csv(all_delta_6m_12m, 
          paste0(save_dir, "jan22_all_delta_6m_12m.csv"), row.names = FALSE)
