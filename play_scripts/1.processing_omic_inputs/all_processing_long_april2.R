# April 2nd 
# Combined all input processing script 
rm(list = ls())
pacman::p_load(knitr, data.table, dplyr, tidyr, tableone, kableExtra, readxl,
               readr, car, RColorBrewer, gridExtra, mlbench, earth, ggplot2, missForest,
               AppliedPredictiveModeling, caret, reshape2, corrplot, stringr,
               summarytools, grid, mice, plyr, mlmRev, cowplot, ape, e1071,
               jtools, broom, patchwork, phyloseq, microbiome, glmnet, ISLR,
               MicrobiomeStat, ANCOMBC, ape, vegan, zCompositions, janitor, MASS,
               RColorBrewer, DT, ggpubr, microbiomeutilities, compositions, VIM)
library(naniar)
library(tibble)
## Functions
# Function to process the names BL
process_names_bl <- function(names) {
  bl_names <- grep("\\.BL$", names, value = TRUE) # ending with "BL"
  extracted_numbers_BL <- sub(".*-(\\d+)\\.BL$", 
                              "\\1", bl_names) # Numbers after "-", before "."
  cleaned_numbers_bl <- sub("^0+", "", extracted_numbers_BL) # Remove leading zeros
  cleaned_numbers_bl} # Return the cleaned numbers
# Function to process the names ALL
process_names_all <- function(names) {
  cleaned_numbers_all <- character(length(names)) #vector to store cleaned numbers
  all_names <- grep("\\..*$", names, value = TRUE) #names matching the pattern
  for (i in seq_along(names)) { # Extract numbers and remove leading zeros
    if (names[i] %in% all_names) {
      extracted_number <- sub(".*-(\\d+)\\..*$", "\\1", names[i])
      cleaned_number <- sub("^0+", "", extracted_number)
      cleaned_numbers_all[i] <- cleaned_number
    } else {
      cleaned_numbers_all[i] <- names[i] # If name don't match, keep the original
    }
  }
  return(cleaned_numbers_all)  # Return the cleaned numbers
}
# Plot Density
my_plt_density <- function(df, cols, xlim_vals, plot_title) {
  path_long <- pivot_longer(df, cols = cols, names_to = "variable", values_to = "value")
  
  # Create a color palette based on the number of variables
  num_vars <- length(unique(path_long$variable))
  color_palette <- rep("skyblue", num_vars)  # Skyblue for all variables (or change as needed)
  color_palette_lines <- rep("black", num_vars)  # Black for all variables (or change as needed)
  
  ggplot(path_long, aes(x = value, fill = variable, color = variable)) + 
    geom_density(alpha = 0.5) + 
    theme_minimal() + 
    xlim(xlim_vals) + 
    labs(title = plot_title, x = "Value", y = "Density") + 
    scale_fill_manual(values = color_palette) + 
    scale_color_manual(values = color_palette_lines) +
    guides(fill = FALSE, color = FALSE)  # Hide the color and fill legends
}
# Check normality and skewness 
check_normality_and_skewness <- function(df, cols) {
  # Shapiro-Wilk test for normality
  normality_summary <- table(ifelse(apply(df[, cols], 2, 
                                          function(x) shapiro.test(x)$p.value) > 0.05, 
                                    "Normal", "Non-Normal"))
  
  # Skewness calculation
  skew_summary <- data.frame(Variable = names(df)[cols],
                             Skewness = apply(df[, cols], 2, skewness, na.rm = TRUE))
  
  # Summary of skewness (skewed or not)
  skew_summary_count <- table(ifelse(abs(skew_summary$Skewness) > 1, "Skewed", "Not Skewed"))
  
  # Automatically print normality and skewness results
  print("Normality Summary:")
  print(normality_summary)
  
  print("Skewness Summary:")
  print(skew_summary_count)
  
  # Return the results as a list
  return(list(
    normality_summary = normality_summary,
    skew_summary_count = skew_summary_count
  ))
}
# Make unique
make_unique_names <- function(names) {
  name_counts <- table(names)  # Create a table of name counts
  occurrence <- integer(length(names)) # Vector to track occurrences
  
  # Loop through each name and assign a unique suffix if necessary
  for (i in seq_along(names)) {
    if (name_counts[names[i]] > 1) {
      occurrence[i] <- sum(names[1:i] == names[i])  # Count occurrences up to current index
      names[i] <- paste0(names[i], "_", occurrence[i])  # Append the occurrence number
    }
  }
  
  # Ensure unique names
  for (i in seq_along(names)) {
    while (sum(names[i] == names) > 1) {
      occurrence[i] <- occurrence[i] + 1
      names[i] <- sub("_\\d+$", "", names[i])  # Remove existing suffix
      names[i] <- paste0(names[i], "_", occurrence[i])  # Append new occurrence number
    }
  }
  return(names)
}
# fancy_processing
fancy_process <- function(data, min_sample_threshold = 0.1, 
                          var_threshold = 0.1, 
                          correlation_threshold = 0.95, 
                          center = TRUE, scale = TRUE) {
  
  # Check if the data is a data.frame or matrix
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input data must be a data frame or matrix.")
  }
  
  # Initialize the removed_info list
  removed_info <- list()
  
  # Step 0: Determine the most suitable transformation based on data distribution
  cat("\nStep 0: Determining the most suitable transformation...\n")
  
  # Function to check normality with Shapiro-Wilk test
  normality_test <- function(x) {
    if (length(x) < 3) return(TRUE)  # Handle small data cases
    p_value <- shapiro.test(x)$p.value
    return(p_value > 0.05)  # p-value > 0.05 means data is normally distributed
  }
  
  # Isolate numeric columns only
  numeric_data <- data[sapply(data, is.numeric)]
  
  # Step 0.1: Test if the numeric data is normal (Shapiro-Wilk test) and skewness
  normal_distrib <- apply(numeric_data, 2, normality_test)  # Apply to all numeric columns
  skew_vals <- apply(numeric_data, 2, function(x) skewness(x, na.rm = TRUE))
  
  # Check if the majority of the data columns are normal
  if (all(normal_distrib)) {
    transformation <- "zscore"
    cat("Data appears to be normally distributed. Using Z-score normalization.\n")
  } else {
    # If data is not normally distributed, check for skewness
    skew_columns <- which(abs(skew_vals) > 0.5)  # Skew threshold (positive or negative skew)
    
    if (length(skew_columns) > 0) {
      if (any(skew_vals[skew_columns] > 0)) {
        # If positively skewed
        transformation <- "log10"
        cat("Data is positively skewed. Using Log10 transformation.\n")
      } else {
        # If negatively skewed or zero/negative values present
        transformation <- "arcsinh"
        cat("Data is negatively skewed or has zero/negative values. Using Arcsinh transformation.\n")
      }
    } else {
      # If data isn't skewed, use Box-Cox for stabilization
      transformation <- "boxcox"
      cat("Data is not skewed. Using Box-Cox transformation.\n")
    }
  }
  
  # Step 0.2: Detect if the data is relative abundance (compositional data)
  is_relative_abundance <- all(abs(rowSums(numeric_data) - 1) < 0.01)  # Adjust threshold as needed
  if (is_relative_abundance) {
    transformation <- "clr"
    cat("Data appears to be relative abundance data. Using CLR transformation.\n")
  }
  
  # Step 1: Normalize the numeric data using selected transformation
  cat("\nStep 1: Normalizing the data...\n")
  if (transformation == "log10") {
    numeric_data <- log10(numeric_data + 1)
    cat("Applied log10 transformation.\n")
  } else if (transformation == "clr") {
    gm <- apply(numeric_data, 1, function(x) exp(sum(log(x[x > 0])) / sum(x > 0)))  # Geometric mean/row
    numeric_data <- t(apply(numeric_data, 1, function(x) log(x / gm)))
    cat("Applied CLR transformation.\n")
  } else if (transformation == "zscore") {
    numeric_data <- scale(numeric_data)
    cat("Applied Z-score normalization.\n")
  } else if (transformation == "boxcox") {
    numeric_data <- apply(numeric_data, 2, function(x) ifelse(all(x > 0), boxcox(x ~ 1)$y, x))
    cat("Applied Box-Cox transformation.\n")
  } else if (transformation == "arcsinh") {
    numeric_data <- asinh(numeric_data)
    cat("Applied Arcsinh transformation.\n")
  }
  
  removed_info$step1 <- list(removed_rows = 0, removed_columns = 0)
  
  # Step 2: Centering and Scaling (standardize data)
  cat("\nStep 2: Centering and Scaling the data...\n")
  if (center) {
    numeric_data <- scale(numeric_data, center = TRUE, scale = FALSE)  # Centering only
    cat("Applied centering (subtracting mean).\n")
  }
  if (scale) {
    numeric_data <- scale(numeric_data, center = FALSE, scale = TRUE)  # Scaling only
    cat("Applied scaling (dividing by standard deviation).\n")
  }
  
  removed_info$step2 <- list(removed_rows = 0, removed_columns = 0)
  
  # Step 3: Remove variables not present in a minimum % threshold of samples
  cat("\nStep 3: Removing variables not present in a minimum percentage of samples...\n")
  min_samples <- floor(min_sample_threshold * nrow(numeric_data))
  cols_before_step3 <- ncol(numeric_data)
  numeric_data <- numeric_data[, colSums(!is.na(numeric_data)) >= min_samples]
  cols_removed_step3 <- cols_before_step3 - ncol(numeric_data)
  cat(cols_removed_step3, "columns removed based on the minimum sample threshold of", 
      min_sample_threshold * 100, "%.\n")
  removed_info$step3 <- list(removed_rows = 0, removed_columns = cols_removed_step3)
  
  # Step 4: Remove variables with variance below the variance threshold
  cat("\nStep 4: Removing variables with variance below", var_threshold, "...\n")
  vars_before_step4 <- ncol(numeric_data)
  numeric_data <- numeric_data[, apply(numeric_data, 2, var, na.rm = TRUE) >= var_threshold]
  vars_removed_step4 <- vars_before_step4 - ncol(numeric_data)
  cat(vars_removed_step4, "columns removed due to variance below the threshold.\n")
  removed_info$step4 <- list(removed_rows = 0, removed_columns = vars_removed_step4)
  
  # Step 5: Remove highly collinear features (Pearson's > 0.95)
  cat("\nStep 5: Removing highly collinear features (Pearson's > 0.95)...\n")
  non_zero_variance_data <- numeric_data[, apply(numeric_data, 2, var, na.rm = TRUE) > 0]
  if (ncol(non_zero_variance_data) > 1) {
    correlation_matrix <- cor(non_zero_variance_data, use = "pairwise.complete.obs")
    upper_tri <- correlation_matrix[upper.tri(correlation_matrix)]
    highly_correlated <- which(abs(upper_tri) > correlation_threshold, arr.ind = TRUE)
    if (length(highly_correlated) > 0) {
      high_corr_idx <- unique(c(highly_correlated[, 1], highly_correlated[, 2]))
      cols_before_step5 <- ncol(numeric_data)
      numeric_data <- numeric_data[, -high_corr_idx]
      cols_removed_step5 <- cols_before_step5 - ncol(numeric_data)
      cat(cols_removed_step5, "columns removed due to high correlation ( > 0.95).\n")
      removed_info$step5 <- list(removed_rows = 0, removed_columns = cols_removed_step5)
    } else {
      cat("No highly correlated features found (  > 0.95).\n")
      removed_info$step5 <- list(removed_rows = 0, removed_columns = 0)
    }
  } else {
    cat("Not enough columns left for correlation analysis.\n")
    removed_info$step5 <- list(removed_rows = 0, removed_columns = 0)
  }
  
  # Replace numeric data back into the original data frame, keeping non-numeric columns unchanged
  non_numeric_data <- data[sapply(data, Negate(is.numeric))]
  final_data <- cbind(non_numeric_data, numeric_data)
  
  # Return the processed data
  return(final_data)
}

# Define a function to print taxa removal summary
print_taxa_summary <- function(before, after, step) {
  removed <- before - after
  print(paste0("# Taxa removed in step ", step, ": ", removed))
  print(paste0("# Taxa remaining after step ", step, ": ", after))
}

# Process each omic seperately #################################################

# Metabolomics data 

# Meta data

# GRS data 

# Taxa data

# Pathway data 

# Micom data 

########### METABOLOMICS #######################################################
tabo <- read_csv("data/metabolomics/UTF-8SendSafely_NGH19889 results_20250401T080522Z/19889 COinterv-01-Apr-2025-Results.csv")
colnames(tabo)
head(tabo[1:30])
# Remove 2-21 cols which are empty & keep only 0,6,12m rows
tabo <- tabo[, c(1, 23:ncol(tabo))] %>% 
  filter(!grepl("\\.3m$", `Sample ID`),
         !grepl("\\.18m$", `Sample ID`)) %>% 
  dplyr::rename(SampleID = `Sample ID`)
tabo$Glycerol <- as.numeric(tabo$Glycerol)

# heatmap
vis_miss(tabo[50:100])
vis_miss(tabo[150:190])
cor_mat <- cor(tabo[2:ncol(tabo)], use = "pairwise.complete.obs")
heatmap(cor_mat)

# Extract highly correlaeted pairs
high_corr <- which(abs(cor_mat) > 0.99 & abs(cor_mat) < 1, arr.ind = TRUE) %>% 
             as.data.frame() %>% rownames_to_column("cor_vars")
# keep <- high_corr %>% dplyr::filter(!grepl("\\.(1|2|3|4)$", cor_vars))

tabo2 <- tabo %>% dplyr::select(colnames(tabo)[!grepl("^(L_|XL_|S_|M_|XS_|XXL_)", colnames(tabo))])
heatmap(cor(tabo2[2:ncol(tabo2)], use = "pairwise.complete.obs"))

vis_miss(tabo2)
my_plt_density(tabo2, 2:ncol(tabo2), c(0.01, 5), "Metabolomics Conc. Raw")
check_normality_and_skewness(tabo2, 2:ncol(tabo2))

# center scale and remove var corr of meta 
preProcValues <- preProcess(tabo2, 
                            method = c("scale", "center", "nzv", "corr"), # Same for metabolomics 
                            thresh = 0.95, # % variance
                            na.remove = TRUE, verbose = TRUE, 
                            freqCut = 95/5, # ratio of most common val to 2nd most common val.
                            uniqueCut = 10, # % distinct vals / total no. samples
                            cutoff = 0.99) # correlation cut off 
preProcValues
tab_cs <- predict(preProcValues, tabo2) 
cor_mat_2 <- cor(tab_cs[2:ncol(tab_cs)], use = "pairwise.complete.obs")
heatmap(cor_mat_2)

library(tidyr)

tab_cs <- tab_cs %>%
  separate(SampleID, into = c("record_id", "time"), sep = "\\.")

tab_cs <- tab_cs %>%
  mutate(time = case_when(grepl("BL", time) ~ 0, grepl("6m", time) ~ 6, 
                          grepl("12m", time) ~ 12, TRUE ~ NA_real_)) %>%
  filter(!grepl("\\.3m$", time))

########### META ###############################################################

meta <- read_csv("~/projects/research/Stanislawski/BMI_risk_scores/data/correct_meta_files/ashleys_meta/DRIFT_working_dataset_meta_deltas_filtered_05.21.2024.csv")

# Replace spaces with underscores in column names
colnames(meta) <- gsub(" ", "_", colnames(meta))
length(unique(meta$subject_id))
meta <- meta[meta$consent != "no", ]
length(unique(meta$subject_id))
a2_extra <- meta %>% dplyr::select(c(record_id, subject_id, randomized_group, consent,
                                     cohort_number, sex, race, completer, age,
                                     # BL #outcome_BMI_fnl_BL
                                     Glucose_BL, HOMA_IR_BL, 
                                     Insulin_endo_BL, HDL_Total_Direct_lipid_BL, 
                                     LDL_Calculated_BL, Triglyceride_lipid_BL,
                                     Peptide_YY_BL, Ghrelin_BL, Leptin_BL, Hemoglobin_A1C_BL,
                                     # 6m # outcome_BMI_fnl_6m
                                     Glucose_6m, HOMA_IR_6m,
                                     Insulin_endo_6m, HDL_Total_Direct_lipid_6m, 
                                     LDL_Calculated_6m, Triglyceride_lipid_6m,
                                     Peptide_YY_6m, Ghrelin_6m, Leptin_6m, Hemoglobin_A1C_6m,
                                     # 12m # outcome_BMI_fnl_12m
                                     Glucose_12m, HOMA_IR_12m,
                                     Insulin_endo_12m, HDL_Total_Direct_lipid_12m,
                                     LDL_Calculated_12m, Triglyceride_lipid_12m,
                                     Peptide_YY_12m, Ghrelin_12m, Leptin_12m, Hemoglobin_A1C_12m))  %>% 
  mutate(subject_id = as.factor(subject_id),
         record_id = as.factor(record_id), 
         randomized_group = as.factor(randomized_group),
         consent = as.factor(consent),
         sex = as.factor(sex), 
         race = as.factor(race),
         completer = as.factor(completer),
         cohort_number = as.factor(cohort_number))

my_plt_density(a2_extra, 9:39, c(-2, 80), "Meta Count Raw")
check_normality_and_skewness(a2_extra, 9:39)
vis_miss(a2_extra)

# center scale and remove var corr of meta 
preProcValues <- preProcess(a2_extra, 
                            method = c("scale", "center", "nzv", "corr"), # Same for metabolomics 
                            thresh = 0.95, # % variance
                            na.remove = TRUE, verbose = TRUE, 
                            freqCut = 95/5, # ratio of most common val to 2nd most common val.
                            uniqueCut = 10, # % distinct vals / total no. samples
                            cutoff = 1) # correlation cut off 
preProcValues
meta_cs <- predict(preProcValues, a2_extra) 
my_plt_density(meta_cs, 9:39, c(-10, 10), "META CENTER SCALED")
check_normality_and_skewness(meta_cs, 9:39)

### Seperate outcome vars ######################################################
siy <- meta %>% dplyr::select(c(subject_id, outcome_BMI_fnl_BL, 
                                outcome_BMI_fnl_6m, outcome_BMI_fnl_12m)) 
my_plt_density(siy, 2:4, c(-5, 50), "raw Y")

### GRS #########################################################################
grs <- read_csv("/Users/emily/projects/research/Stanislawski/BMI_risk_scores/full_cohort_pulling_snps/bigsnpr/made_scores/merge_meta_methyl.csv") %>%
  dplyr::select(c(subject_id, raw_score)) %>%
  dplyr::rename(bmi_prs = raw_score)
my_plt_density(grs, 2:2, c(-1, 1), "GRS")
check_normality_and_skewness(grs, 2:2)
grs <- grs %>% distinct(subject_id, .keep_all = TRUE)

preProcValues <- preProcess(grs[,2:2], 
                            method = c("nzv", "scale", "center"), 
                            thresh = 0.95, na.remove = TRUE, 
                            numUnique = 15, verbose = TRUE,  freqCut = 95/5, 
                            uniqueCut = 10, cutoff = 0.95)
preProcValues
grs_transformed <- predict(preProcValues, grs[, 1:2]) 
my_plt_density(grs_transformed, 2:2, c(-10, 10), "GRS nzv")

### Taxa #################################################################################################### 
load("~/projects/research/Stanislawski/BMI_risk_scores/microbiome_rs/data/PhyloseqObj.RData")
load("/Users/emily/projects/research/Stanislawski/BMI_risk_scores/microbiome_rs/data/Genus_Sp_tables.RData")
print_ps(drift.phy.count)

# Step 1: Remove taxa not seen more than 3 times in at least 10% of the samples
initial_taxa_count <- ntaxa(drift.phy.count)
GP_count <- filter_taxa(drift.phy.count, function(x) sum(x > 3) > (0.1*length(x)), TRUE)
print_taxa_summary(initial_taxa_count, ntaxa(GP_count), "1")

# Step 2: Apply coefficient of variation cutoff
initial_taxa_count_gp <- ntaxa(GP_count)
gpsf_count <- filter_taxa(GP_count, function(x) sd(x)/mean(x) > 1.0, TRUE)
print_taxa_summary(initial_taxa_count_gp, ntaxa(gpsf_count), "2")

# Get genus --- double check MATTS email 
tax_genus <- tax_glom(gpsf_count, "Genus")
genus_filtered_count <- otu_table(tax_genus) 
genus_filtered_count <- t(genus_filtered_count) %>% as.data.frame()
colnames(genus_filtered_count) <- as.data.frame(tax_table(tax_genus))$Genus
rm(drift.phy.clr,drift.phy.count,GP_count,gpsf_count,sp.clr, tax_genus,sp.count,
   drift.phy.count.r21116, drift.phy.ra, genus.clr, genus.count, genus.ra, sp.ra)

# Function to make column names unique
colnames(genus_filtered_count) <- make_unique_names(colnames(genus_filtered_count)) # rename the columns
genus_count <- genus_filtered_count %>% rownames_to_column(var = "subject_id") 
my_plt_density(genus_count, 2:180, c(-100, 100), "Genus Count")
check_normality_and_skewness(genus_count, 2:180)

# Relative abundance conversion
#Genus_relative_abundance <- genus_filtered_count %>%
#  rownames_to_column(var = "subject_id") %>%  dplyr::rowwise() %>%
#  mutate(total = sum(c_across(g__Parabacteroides_B_862066:`g__Massilistercora`), na.rm = TRUE)) %>%
#  mutate(across(g__Parabacteroides_B_862066:`g__Massilistercora`, ~ .x / total)) %>%
#  dplyr::select(-total) %>% column_to_rownames(var = "subject_id")

##  removes any column where the proportion of zeros is greater than or equal to 80%.
threshold <- 0.80
zero_percentage <- colSums(genus_count == 0, na.rm = TRUE) / nrow(genus_count) # % zeros / column
genus_count_cleaned <- genus_count[, zero_percentage < threshold] # only colz < 20% zeros
dim(genus_count) - dim(genus_count_cleaned) #

# CLR Transformation -  ON THE COUNT DATA 
genus_clr_transformed <- cbind(genus_count_cleaned[, 1, drop = FALSE], 
                               apply(genus_count_cleaned[, -1], 2, clr) %>% as.data.frame())


my_plt_density(genus_clr_transformed, 2:133, c(-1, 1), "Genus RA CLR")

# remove low var and high corr
preProcValues <- preProcess(genus_clr_transformed[,1:133], 
                            method = c( "nzv", "corr"), 
                            thresh = 0.95, na.remove = TRUE, fudge = 0.2, 
                            numUnique = 15, verbose = TRUE,  freqCut = 95/5, 
                            uniqueCut = 10, cutoff = 0.95)
preProcValues
cs_transformed <- predict(preProcValues, genus_clr_transformed[, 1:133]) 
heatmap(cor(cs_transformed[, 2:122]))
my_plt_density(cs_transformed, 2:122, c(-1, 1), "Genus CLR Var Corr")
check_normality_and_skewness(cs_transformed, 2:122)
cs_transformed$all_samples <- process_names_all(cs_transformed$subject_id)

### Functional ####################################################################################################
pathways <- fread("/Users/emily/projects/research/Stanislawski/BMI_risk_scores/picrust2/june7/pathways_out/path_abun_unstrat_descrip.tsv") %>% .[, -1] %>% 
  t() %>% as.data.frame() %>% row_to_names(1) %>% rownames_to_column("SampleID") %>% mutate(across(-1, as.numeric))
my_plt_density(pathways, 2:391, c(-50, 50), "Pathways")

# Remove variables with low presence thresholds 
threshold <- 0.80
zero_percentage <- colSums(pathways == 0, na.rm = TRUE) / nrow(pathways) # % zeros / column
path_all_time_cleaned <- pathways[, zero_percentage < threshold] # only colz < 20% zeros
dim(pathways) - dim(path_all_time_cleaned) 
any(path_all_time_cleaned < 0, na.rm = TRUE)
my_plt_density(path_all_time_cleaned, 2:312, c(-5, 5), "Pathways")
check_normality_and_skewness(path_all_time_cleaned, 2:312)

# MAYBE CLR TRANSFORMATION for PATHWAYs
path_clr_transformed <- cbind(path_all_time_cleaned[, 1, drop = FALSE], 
                               apply(path_all_time_cleaned[, -1], 2, clr) %>% as.data.frame())

# remove low varianve and high correlation 
preProcValues <- preProcess(path_clr_transformed[,2:312], 
                            method = c("nzv", "corr"), thresh = 0.95, fudge = 0.2,
                            na.remove = TRUE, numUnique = 15, verbose = TRUE, 
                            freqCut = 95/5, uniqueCut = 10, cutoff = 0.90)
preProcValues
path_all_time_cs <- predict(preProcValues, path_clr_transformed[,1:312])
heatmap(cor(path_all_time_cs[, 2:128]))
my_plt_density(path_all_time_cs, 2:128, c(-5, 5), "Pathways CLR var corr ")
check_normality_and_skewness(path_all_time_cs, 2:128)
# Process samples 
path_all_time_cs$all_samples <- process_names_all(path_all_time_cs$SampleID)
colnames(path_all_time_cs) <- make.names(colnames(path_all_time_cs), unique = TRUE)
# make time column 
path_all_time_cs <- path_all_time_cs %>%
  mutate(time = case_when(grepl("BL", SampleID) ~ 0, grepl("6m", SampleID) ~ 6, 
         grepl("12m", SampleID) ~ 12, TRUE ~ NA_real_)) %>% filter(!grepl("\\.3m$", SampleID)) 

### MICOM ####################################################################################################
load("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/micom/WD_Cplex/products_all.RData")
flux <- products.all %>% 
  dplyr::filter(diet == "Western", timepoint2 %in% c("BL", "6m", "12m")) %>% 
  dplyr::select(sample_id, description, flux) %>% 
  pivot_wider(names_from = description, values_from = flux, values_fill = list(flux = 0))
any(flux < 0, na.rm = TRUE)

# Remove variables with % of zero values is greater than or equal to the threshold
threshold <- 0.80
zero_percentage <- colSums(flux == 0, na.rm = TRUE) / nrow(flux)
flux_cleaned <- flux[, zero_percentage < threshold]
dim(flux) - dim(flux_cleaned) # 0 
my_plt_density(flux_cleaned, 2:104, c(-0.001, 0.001), "MICOM raw")
check_normality_and_skewness(flux_cleaned, 2:104)
heatmap(cor(flux_cleaned[, 2:104]))

# CHECK MICOM DATA TRANSFORMATIONS - SEE WHAT THEY DO - APPLIED USING MICOM ON HUMAN DATA
# https://github.com/Gibbons-Lab/scfa_predictions/blob/main/notebooks/summarized_exvivos.ipynb
# This paper uses z-score normalization for micom

# Z socre tramsformation as per the paper 
flux_cleaned <- as.data.frame(flux_cleaned)
flux_cleaned[] <- lapply(flux_cleaned, function(x) if(is.numeric(x)) scale(x) else x)
my_plt_density(flux_cleaned, 2:121, c(-1, 1), "MICOM z-scored")

flux_cleaned[, 2:121] <- lapply(flux_cleaned[, 2:121], as.numeric)
preProcValues <- preProcess(flux_cleaned[, 2:121], 
                            method = c("nzv", "corr"), thresh = 0.95, fudge = 0.2, 
                            numUnique = 15, verbose = TRUE, freqCut = 95/5, 
                            uniqueCut = 10, cutoff = 0.85, na.remove = TRUE)

preProcValues
flux_all_time_cs <- predict(preProcValues, flux_cleaned)
heatmap(cor(flux_all_time_cs[, 2:ncol(flux_all_time_cs)]))

dim(flux_cleaned) - dim(flux_all_time_cs) # 60 
my_plt_density(flux_all_time_cs, 2:ncol(flux_all_time_cs), 
               c(-2.5, 2.5), "MICOM z-scored var corr")
check_normality_and_skewness(flux_all_time_cs, 2:ncol(flux_all_time_cs))

# Process names
flux_all_time_cs$all_samples <- process_names_all(flux_all_time_cs$sample_id)
colnames(flux_all_time_cs) <- make.names(colnames(flux_all_time_cs), unique = TRUE)

flux_all_time_cs <- flux_all_time_cs %>%
  mutate(time = case_when(grepl("BL", sample_id) ~ 0, grepl("6m", sample_id) ~ 6, 
                          grepl("12m", sample_id) ~ 12, TRUE ~ NA_real_)) %>%
  filter(!grepl("\\.3m$", sample_id))

##### Merging DFs ###########################################################################
siy_meta <- merge(siy, meta_cs, by = "subject_id", all = TRUE)
siy_meta_grs <- merge(siy_meta, grs, by = "subject_id", all = TRUE)

meta_cs_0 <- siy_meta_grs %>% setNames(gsub("_BL$", "", colnames(.))) %>% 
             mutate(time = 0) %>% dplyr::select(-which(grepl("(_6m$|_12m$)", colnames(.))))
meta_cs_6 <- siy_meta_grs %>% setNames(gsub("_6m$", "", colnames(.))) %>% 
             mutate(time = 6) %>% dplyr::select(-which(grepl("(_BL$|_12m$)", colnames(.))))
meta_cs_12 <- siy_meta_grs %>% setNames(gsub("_12m$", "", colnames(.))) %>% 
              mutate(time = 12) %>% dplyr::select(-which(grepl("(_BL$|_6m$)", colnames(.))))
lg_met <- bind_rows(meta_cs_0, meta_cs_6, meta_cs_12) %>% unique()
vis_miss(lg_met)
colnames(lg_met)

# Merge Metabolomics
setdiff(lg_met$record_id, tab_cs$record_id)
setdiff(tab_cs$record_id, lg_met$record_id)
lg_met_met <- merge(lg_met, tab_cs, by = c("record_id", "time"), all.x = TRUE)

# Merge other long data
tx_pth <- merge(cs_transformed, path_all_time_cs, by.x = "subject_id", by.y = "SampleID")
tx_pth_micom <- merge(tx_pth, flux_all_time_cs, by.x = "subject_id", by.y = "sample_id") %>% 
  dplyr::select(-c("all_samples.x", "all_samples.y", "subject_id", "time.y"))

long <- merge(lg_met_met, tx_pth_micom, 
              by.x = c("record_id", "time"), by.y = c("all_samples", "time.x"),
              all = TRUE) %>% arrange(time) %>% unique()

vis_miss(long)
vis_miss(long[1:30])
md.pattern(long[, c(1:30, (ncol(long)-9):ncol(long))], rotate.names = TRUE)

### PROCESS LONG DATA FRAME 

# Summarize 'outcome_BMI_fnl' by 'time' column
long$time <- as.factor(long$time)
table(long$time)
# Drop put subject_id's
missing_6_12m_bmi <- long %>% 
          dplyr::filter(time != 0, is.na(outcome_BMI_fnl)) %>% 
          dplyr::select(subject_id) %>% count() 
count(unique(missing_6_12m_bmi$subject_id)) # 45 

# This DF will be used for the long models :
long_for_models <- long %>% dplyr::filter(!is.na(outcome_BMI_fnl)) %>% 
                            dplyr::filter(rowSums(is.na(.)) / ncol(.) <= 0.8)
length(unique(long_for_models$subject_id))
table(table(long_for_models$subject_id)) #  LBL-047 is repeated? 

# remove LBL-047
print(names(which(table(long_for_models$subject_id) == 5)))
long_for_models2 <- long_for_models[!grepl("LBL-047", long_for_models$subject_id), ] 
print(names(which(table(long_for_models2$subject_id) == 5)))
table(table(long_for_models2$subject_id))

# impute. test
long_imputed <- cbind(long_for_models2[, !sapply(long_for_models2, is.numeric)], 
                      missForest(long_for_models2[, sapply(long_for_models2, is.numeric)])$ximp) %>% 
  unique()
table(table(long_imputed$subject_id))
vis_miss(long_imputed) # This DF will be for longitudinal aims 


# This DF will be used for the DELTA making 
print(names(which(table(long$subject_id) == 5)))
long2 <- long[!grepl("LBL-047", long$subject_id), ]
long_no_drops_0_6m <- long2 %>% 
  dplyr::filter(!(time = 6 & is.na(outcome_BMI_fnl))) %>% arrange(time)
dim(long) - dim(long_no_drops_0_6m)
vis_miss(long_no_drops_0_6m)
long_min_na_0_6m <- long_no_drops_0_6m %>%
  filter(time != 12) %>%  # remove 12 month vars 
  filter(rowSums(is.na(.)) / ncol(.) <= 0.80) # remove >80% missig data
vis_miss(long_min_na_0_6m)
vis_miss(long_min_na_0_6m[1:35])
table(table(long_min_na_0_6m$subject_id.x))

long_no_drops_0_12m <- long2 %>% 
  dplyr::filter(!(time != 0 & is.na(outcome_BMI_fnl))) %>% arrange(time)
dim(long) - dim(long_no_drops_0_12m)
vis_miss(long_no_drops_0_6m)
long_min_na_0_12m <- long_no_drops_0_12m %>%
  filter(rowSums(is.na(.)) / ncol(.) <= 0.80) # remove >80% missig data
vis_miss(long_min_na_0_12m)
vis_miss(long_min_na_0_12m[1:35])
table(table(long_min_na_0_12m$subject_id.x))

# Take out black box line and impute the rest 
print(names(which(table(long_min_na_0_12m$subject_id) == 5)))

###### MAKE DELTAS #############################################################
# Filter the data into two data frames based on 'time' column (0 and 6)
df_time_6 <- long_min_na_0_6m %>% filter(time == 6)
df_time_0 <- long_min_na_0_6m %>% filter(time == 0) 

# Subset df_time_0 and df_time_6 to only include matching subject_id.x values
df_time_0 <- df_time_0 %>% filter(subject_id %in% df_time_6$subject_id)
df_time_6 <- df_time_6 %>% filter(subject_id %in% df_time_0$subject_id)
nrow(df_time_0)
nrow(df_time_6)
  
# Identify numeric columns in both data frames
numeric_cols_0 <- df_time_0 %>% select_if(is.numeric) %>% dplyr::select(-c(age, bmi_prs))
numeric_cols_6 <- df_time_6 %>% select_if(is.numeric) %>% dplyr::select(-c(age, bmi_prs))

# Calculate the change (difference between time 6 and time 0 for numeric columns)
change_cols <- numeric_cols_6 - numeric_cols_0

# Select non-numeric columns
non_numeric_cols <- df_time_6 %>%
  dplyr::select( which(!sapply(., is.numeric)),  # non-numeric columns
    age, bmi_prs)

# Combine the numeric changes with non-numeric columns
change_df <- bind_cols(non_numeric_cols, change_cols)

# View the change dataframe
head(change_df)
change_notx_pth <- change_df %>% dplyr::select(-c(101:348)) # remove taxa and patway data
change_notx_pth_0_6_imputed <- cbind(change_notx_pth[, !sapply(change_notx_pth, is.numeric)], 
                                     missForest(change_notx_pth[, sapply(change_notx_pth, is.numeric)])$ximp) %>% 
  unique()
vis_miss(change_notx_pth_0_6_imputed)

# Load Delta change taxa
#load("~/projects/research/Stanislawski/comps/mutli-omic-predictions/data/taxa_change/Change Data/DRIFT2-K01-GenusChangeCLR_0_6_pldist.RData")
load('~/projects/research/Stanislawski/comps/mutli-omic-predictions/data/taxa_change/Change Data/tax_0_6_pldist_CLR.RData')
tax0_6m <- quant.clr %>% as.data.frame() %>% rownames_to_column("subject_id")
length(unique(tax0_6m$subject_id))
length(setdiff(change_notx_pth_0_6_imputed$subject_id.x, tax0_6m$subject_id))
# only keep the taxa that are also in the long datasets above 
setdiff(colnames(cs_transformed), colnames(tax0_6m))
setdiff(colnames(tax0_6m), colnames(cs_transformed))
tax0_6m_slim <- tax0_6m[, intersect(colnames(tax0_6m), colnames(cs_transformed))]

# Make Pathway change 0-6 #
load('~/projects/research/Stanislawski/comps/mutli-omic-predictions/data/taxa_change/Change Data/pathway_0_6_pldist_CLR.RData')
path_0_6 <- quant.clr %>% as.data.frame() %>% rownames_to_column("subject_id")
length(unique(path_0_6$subject_id))
path_0_6_slim <- path_0_6[, c("subject_id", intersect(colnames(path_0_6)[-1], colnames(path_all_time_cs)))]

# Merge taxa and meta data
met_tax_0_6 <- merge(change_notx_pth_0_6_imputed, tax0_6m_slim, 
                     by.x = "subject_id", by.y = "subject_id")
# And pathway data 
change_all_0_6 <- merge(met_tax_0_6, path_0_6_slim, 
                        by.x = "subject_id", by.y = "subject_id")

vis_miss(change_all_0_6)

#### Repeat for 6-12m ##########################################################
# Filter the data into two data frames based on 'time' column (0 and 6)
df_time_12 <- long_min_na_0_12m %>% filter(time == 12)
df_time_6 <- long_min_na_0_12m %>% filter(time == 6) 

# Subset df_time_6 and df_time_6 to only include matching subject_id.x values
df_time_6 <- df_time_6 %>% filter(subject_id %in% df_time_12$subject_id)
df_time_12 <- df_time_12 %>% filter(subject_id %in% df_time_6$subject_id)
nrow(df_time_6)
nrow(df_time_12)

# Identify numeric columns in both data frames
numeric_cols_6 <- df_time_6 %>% select_if(is.numeric) %>% dplyr::select(-c(age, bmi_prs))
numeric_cols_12 <- df_time_12 %>% select_if(is.numeric) %>% dplyr::select(-c(age, bmi_prs))

# Calculate the change (difference between time 6 and time 0 for numeric columns)
change_cols <- numeric_cols_12 - numeric_cols_6

# Select non-numeric columns
non_numeric_cols <- df_time_12 %>% dplyr::select( which(!sapply(., is.numeric)),  # non-numeric columns
                                                  age, bmi_prs)

# Combine the numeric changes with non-numeric columns
change_df <- bind_cols(non_numeric_cols, change_cols)

# View the change dataframe
head(change_df)
change_notx_pth <- change_df %>% dplyr::select(-c(101:348)) 
change_notx_pth_0_6_imputed <- cbind(change_notx_pth[, !sapply(change_notx_pth, is.numeric)], 
                                     missForest(change_notx_pth[, sapply(change_notx_pth, is.numeric)])$ximp) %>% 
  unique()
vis_miss(change_notx_pth_0_6_imputed)

# Load Delta change taxa
#load("~/projects/research/Stanislawski/comps/mutli-omic-predictions/data/taxa_change/Change Data/DRIFT2-K01-GenusChangeCLR_0_6_pldist.RData")
load('~/projects/research/Stanislawski/comps/mutli-omic-predictions/data/taxa_change/Change Data/tax_6_12_pldist_CLR.RData')
tax0_6m <- quant.clr %>% as.data.frame() %>% rownames_to_column("subject_id")
length(unique(tax0_6m$subject_id))
length(setdiff(change_notx_pth_0_6_imputed$subject_id.x, tax0_6m$subject_id))
# only keep the taxa that are also in the long datasets above 
setdiff(colnames(cs_transformed), colnames(tax0_6m))
setdiff(colnames(tax0_6m), colnames(cs_transformed))
tax0_6m_slim <- tax0_6m[, intersect(colnames(tax0_6m), colnames(cs_transformed))]

# Make Pathway change 6-12 #
load('~/projects/research/Stanislawski/comps/mutli-omic-predictions/data/taxa_change/Change Data/pathway_6_12_pldist_CLR.RData')
path_0_6 <- quant.clr %>% as.data.frame() %>% rownames_to_column("subject_id")
length(unique(path_0_6$subject_id))
path_0_6_slim <- path_0_6[, c("subject_id", intersect(colnames(path_0_6)[-1], colnames(path_all_time_cs)))]

# Merge taxa and meta data
met_tax_0_6 <- merge(change_notx_pth_0_6_imputed, tax0_6m_slim, 
                     by.x = "subject_id", by.y = "subject_id")
# And pathway data 
change_all_6_12 <- merge(met_tax_0_6, path_0_6_slim, 
                        by.x = "subject_id", by.y = "subject_id")

vis_miss(change_all_6_12)

## Merge Delta Files ##
setdiff(colnames(change_all_0_6), colnames(change_all_6_12))
all_delta <- rbind(change_all_0_6, change_all_6_12)

# SAVE FILES ###############################################################################

# long 
#write.csv(long_imputed, file = '/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/april_processing/long.csv')

# deltas
#write.csv(change_all_0_6, file = '/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/april_processing/delta_0_6.csv')
#write.csv(change_all_6_12, file = '/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/april_processing/delta_6_12.csv')
#write.csv(all_delta, file = '/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/april_processing/all_delta.csv')

# Make t_plus1 files fpr long data ###############################################################################

long_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/april_processing/"

long <- read.csv(file.path(long_dir, 'long.csv')) %>% 
  dplyr::select(-c("consent", "record_id", "completer", 
                   "Peptide_YY", "Ghrelin", "Leptin")) %>%
  dplyr::mutate(time = as.factor(time),
                subject_id = as.factor(subject_id),
                randomized_group = as.factor(randomized_group),
                sex = as.numeric(sex),
                randomized_group = as.numeric(randomized_group),
                cohort_number = as.numeric(cohort_number),
                race = as.numeric(race)) %>% dplyr::rename(BMI = outcome_BMI_fnl,
                                                    range = time,
                                                    homo_ir = HOMA_IR,
                                                    insulin = Insulin_endo,
                                                    LDL = LDL_Calculated,
                                                    HDL = HDL_Total_Direct_lipid,
                                                    HbA1c = Hemoglobin_A1C) 

# Split into A and B 
A_df <- long %>%
  dplyr::select(subject_id, range, BMI) %>%
  dplyr::mutate(range = dplyr::recode(as.numeric(as.character(range)), 
                                      `0` = 0, `6` = 1, `12` = 2)) %>% 
  dplyr::filter(!range == 0)

B_df <- long %>%
  dplyr::mutate(range = dplyr::recode(as.numeric(as.character(range)), 
                                      `0` = 1, `6` = 2, `12` = 3)) %>% 
  dplyr::filter(!range == 3)

long_t_plus <- merge(A_df, B_df, by = c("subject_id", "range")) # -4 


# Make t_plus1 files fpr delta data ###############################################################################







