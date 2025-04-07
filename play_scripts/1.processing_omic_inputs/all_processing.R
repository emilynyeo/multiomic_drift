# March 19 
# Combined all input processing script 

pacman::p_load(knitr, data.table, dplyr, tidyr, tableone, kableExtra, readxl,
               readr, car, RColorBrewer, gridExtra, mlbench, earth, ggplot2, missForest,
               AppliedPredictiveModeling, caret, reshape2, corrplot, stringr,
               summarytools, grid, mice, plyr, mlmRev, cowplot, ape, e1071,
               jtools, broom, patchwork, phyloseq, microbiome, glmnet, ISLR,
               MicrobiomeStat, ANCOMBC, ape, vegan, zCompositions, janitor, MASS,
               RColorBrewer, DT, ggpubr, microbiomeutilities, compositions, VIM)

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

# Scaling ensures that all features contribute equally to the analysis, especially when the features have different units or ranges. 
# Centering removes any biases caused by differences in the central values of the features.
# The operations are applied in this order: zero-variance filter, near-zero variance filter, correlation filter, Box-Cox/Yeo-Johnson/exponential transformation, centering, scaling, range, imputation, PCA, ICA then spatial sign
# se Exponential Transformations when the data represents rates, decay, or growth over time.

#### META ##########################################################################################
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

my_plt_density(a2_extra, 9:39, c(-2, 80), "Meta Count")
check_normality_and_skewness(a2_extra, 9:39)
vis_miss(a2_extra)
# LOOk 
preProcValues <- preProcess(a2_extra, 
                            method = c("scale", "center", "nzv", #"YeoJohnson", 
                                       "corr"), # "bagImpute"
                            thresh = 0.95, # % variance
                            na.remove = TRUE, verbose = TRUE, 
                            freqCut = 95/5, # ratio of most common val to 2nd most common val.
                            uniqueCut = 10, # % distinct vals / total no. samples
                            cutoff = 1) # correlation cut off 
                            #rangeBounds = c(0, 1)) # interval for range transformation
preProcValues
meta_cs <- predict(preProcValues, a2_extra) 
heatmap(cor(meta_cs[, 9:39]))
my_plt_density(meta_cs, 9:39, c(-10, 10), "META CENTER SCALED")
check_normality_and_skewness(meta_cs, 9:39)

#### SET Y Vars ####################################################################################
siy <- meta %>% dplyr::select(c(subject_id, outcome_BMI_fnl_BL, 
                                outcome_BMI_fnl_6m, outcome_BMI_fnl_12m)) 
my_plt_density(siy, 2:4, c(-5, 50), "raw Y")

### GRS ####################################################################################################
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

dub <- grs_transformed %>% filter(duplicated(subject_id) | duplicated(subject_id, fromLast = TRUE))

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
Genus_relative_abundance <- genus_filtered_count %>%
  rownames_to_column(var = "subject_id") %>%  
  dplyr::rowwise() %>%
  mutate(total = sum(c_across(g__Parabacteroides_B_862066:`g__Massilistercora`), na.rm = TRUE)) %>%
  mutate(across(g__Parabacteroides_B_862066:`g__Massilistercora`, ~ .x / total)) %>%
  dplyr::select(-total) %>%
  column_to_rownames(var = "subject_id")

# Plot RA distribuitions 
any(genus_count < 0, na.rm = TRUE)
any(Genus_relative_abundance < 0, na.rm = TRUE)
my_plt_density(Genus_relative_abundance, 2:179, c(-0.001, 0.09), "Genus RA")
check_normality_and_skewness(Genus_relative_abundance, 2:179)

##  removes any column where the proportion of zeros is greater than or equal to 80%.
threshold <- 0.80
zero_percentage <- colSums(Genus_relative_abundance == 0, na.rm = TRUE) / nrow(Genus_relative_abundance) # % zeros / column
Genus_relative_abundance_cleaned <- Genus_relative_abundance[, zero_percentage < threshold] # only colz < 20% zeros
dim(Genus_relative_abundance) - dim(Genus_relative_abundance_cleaned) # 47 

# CLR Transformation - CHANGE THIS TO ONLY BE DONE ON THE COUNT
genus_clr_transformed <- apply(Genus_relative_abundance_cleaned, 2, clr) %>% as.data.frame()
my_plt_density(genus_clr_transformed, 1:132, c(-0.005, 0.005), "Genus RA CLR")
my_plt_density(genus_clr_transformed, 1:132, c(-0.1, 0.1), "Genus RA CLR")
check_normality_and_skewness(genus_clr_transformed, 1:132)

#  with more than 95% of its values being the same (or near-zero variance) will be excluded from the dataset.
# CENTER AND SCALE
preProcValues <- preProcess(genus_clr_transformed[,1:132], 
                            method = c( "nzv", "corr"), # , "expoTrans""scale", "center",
                            thresh = 0.95, na.remove = TRUE, 
                            fudge = 0.2, 
                            numUnique = 15, verbose = TRUE,  freqCut = 95/5, 
                            uniqueCut = 10, cutoff = 0.95, rangeBounds = c(0, 1))

preProcValues
cs_transformed <- predict(preProcValues, genus_clr_transformed[, 1:132]) %>% 
                   rownames_to_column(var = "subject_id")
heatmap(cor(cs_transformed[, 2:133]))
my_plt_density(cs_transformed, 2:133, c(-1, 1), "GenusRA CLR")
check_normality_and_skewness(cs_transformed, 2:133)
cs_transformed$all_samples <- process_names_all(cs_transformed$subject_id)

### Functional ####################################################################################################
pathways <- fread("/Users/emily/projects/research/Stanislawski/BMI_risk_scores/picrust2/june7/pathways_out/path_abun_unstrat_descrip.tsv") %>% .[, -1] %>% t() %>% as.data.frame() %>% 
  row_to_names(1) %>% rownames_to_column("SampleID") %>% 
  mutate(across(-1, as.numeric))
my_plt_density(pathways, 2:391, c(-50, 50), "Pathways")

# Remove variables with low presence thresholds 
threshold <- 0.80
zero_percentage <- colSums(pathways == 0, na.rm = TRUE) / nrow(pathways) # % zeros / column
path_all_time_cleaned <- pathways[, zero_percentage < threshold] # only colz < 20% zeros
dim(pathways) - dim(path_all_time_cleaned) #79
any(path_all_time_cleaned < 0, na.rm = TRUE)
my_plt_density(path_all_time_cleaned, 2:312, c(-5, 5), "Pathways")
check_normality_and_skewness(path_all_time_cleaned, 2:312)

# MAYBE CLR TRANSFORMATION

# imputation and scaling 
preProcValues <- preProcess(path_all_time_cleaned[,2:312], 
                            method = c("scale", "center", "nzv", #  "YeoJohnson"
                                       "corr"),
                            thresh = 0.95,  na.remove = TRUE, 
                            fudge = 0.2, numUnique = 15, verbose = TRUE,  
                            freqCut = 95/5, uniqueCut = 10, cutoff = 0.90)
preProcValues
path_all_time_cs <- predict(preProcValues, path_all_time_cleaned[,1:312])
heatmap(cor(path_all_time_cs[, 2:121]))
my_plt_density(path_all_time_cs, 2:121, c(-5, 5), "Pathways Center Scaled")
check_normality_and_skewness(path_all_time_cs, 2:121)
# Process samples 
path_all_time_cs$all_samples <- process_names_all(path_all_time_cs$SampleID)
colnames(path_all_time_cs) <- make.names(colnames(path_all_time_cs), unique = TRUE)
# make time column 
path_all_time_cs <- path_all_time_cs %>%
  mutate(time = case_when(grepl("BL", SampleID) ~ 0, grepl("6m", SampleID) ~ 6, 
                          grepl("12m", SampleID) ~ 12, TRUE ~ NA_real_)) %>%
  filter(!grepl("\\.3m$", SampleID)) # Remove 3 month cols

### MICOM ####################################################################################################
load("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/micom/WD_Cplex/products_all.RData")
flux <- products.all %>% 
  filter(diet == "Western", timepoint2 %in% c("BL", "6m", "12m")) %>% 
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
my_plt_density(flux_cleaned, 2:104, c(-0.001, 0.001), "MICOM z-scored")

# remove low variance 
preProcValues <- preProcess(flux_cleaned[, 2:104], 
                            method = c("nzv", "corr"), thresh = 0.95, fudge = 0.2, 
                            numUnique = 15, verbose = TRUE, freqCut = 95/5, 
                            uniqueCut = 10, cutoff = 0.95, na.remove = TRUE)
preProcValues
flux_all_time_cs <- predict(preProcValues, flux_cleaned[,1:104])
heatmap(cor(flux_all_time_cs[, 2:61]))

dim(flux_cleaned) - dim(flux_all_time_cs) # 60 
my_plt_density(flux_all_time_cs, 2:61, c(-2.5, 2.5), "MICOM z-scored var corr")
check_normality_and_skewness(flux_all_time_cs, 2:61)

# Process names
flux_all_time_cs$all_samples <- process_names_all(flux_all_time_cs$sample_id)
colnames(flux_all_time_cs) <- make.names(colnames(flux_all_time_cs), unique = TRUE)
flux_all_time_cs <- flux_all_time_cs %>%
  mutate(time = case_when(grepl("BL", sample_id) ~ 0, grepl("6m", sample_id) ~ 6, 
                          grepl("12m", sample_id) ~ 12, TRUE ~ NA_real_)) %>%
  filter(!grepl("\\.3m$", sample_id))

##### Merging DFs ###########################################################################
# siy, meta_cs, grs_transformed, 
# cs_transformed, path_all_time_cs, flux_all_time_cs

head(cs_transformed$subject_id)
head(path_all_time_cs$SampleID)
head(flux_all_time_cs$sample_id)
head(siy)
head(meta_cs)
head(grs_transformed)

tx_pth <- merge(cs_transformed, path_all_time_cs, by.x = "subject_id", by.y = "SampleID")
tx_pth_micom <- merge(tx_pth, flux_all_time_cs, by.x = "subject_id", by.y = "sample_id")

# BL 
tx_p_m_BL <- tx_pth_micom %>% filter(time.x == 0.00) 
colnames(tx_p_m_BL)[-1] <- paste0(colnames(tx_p_m_BL)[-1], "_BL")
tx_p_m_BL$subject_id_6digits <- substr(tx_p_m_BL$subject_id, 1, 7)
# 6m
tx_p_m_6m <- tx_pth_micom %>% filter(time.x == 6.000) 
colnames(tx_p_m_6m)[-1] <- paste0(colnames(tx_p_m_6m)[-1], "_6m")
tx_p_m_6m$subject_id_6digits <- substr(tx_p_m_6m$subject_id, 1, 7)
# 12m 
tx_p_m_12m <- tx_pth_micom %>% filter(time.x == 12.000)
colnames(tx_p_m_12m)[-1] <- paste0(colnames(tx_p_m_12m)[-1], "_12m")
tx_p_m_12m$subject_id_6digits <- substr(tx_p_m_12m$subject_id, 1, 7)

y_met_grs <- siy %>%
  left_join(meta_cs, by = "subject_id") %>%
  left_join(grs_transformed, by = "subject_id")

y_met_grs_BL <- y_met_grs %>% dplyr::select(1:12, "bmi_prs" , ends_with("BL"))
y_met_grs_6m <- y_met_grs %>% dplyr::select(1:12, "bmi_prs" , ends_with("6m"))
y_met_grs_12m <- y_met_grs %>% dplyr::select(1:12, "bmi_prs", ends_with("12m"))

BL <- y_met_grs_BL %>% merge(tx_p_m_BL, by.x = "subject_id", 
                             by.y = "subject_id_6digits", all.y = TRUE) %>% 
  unique() %>% dplyr::filter(grepl("BL$", subject_id.y))

m6 <- y_met_grs_6m %>% merge(tx_p_m_6m, by.x = "subject_id", 
                             by.y = "subject_id_6digits", all.y = TRUE) %>% unique()
m12 <- y_met_grs_12m %>% merge(tx_p_m_12m, by.x = "subject_id", 
                             by.y = "subject_id_6digits", all.y = TRUE) %>% unique()

#Now merge the dataframes based on the new 7-character columns
BL_6m <- full_join(BL, m6, by = c("subject_id" = "subject_id")) %>% unique()
m6_m12 <- full_join(m6, m12, by = c("subject_id" = "subject_id")) %>% unique()

length(unique(BL_6m$subject_id))
length(unique(m6_m12$subject_id))
# % of rows with >45% missing data 
vis_miss(BL_6m)
cat("BL_6m missing >45% data:", mean(apply(BL_6m, 1, function(x) sum(is.na(x)) / length(x)) > 0.45) * 100, "%\n")
vis_miss(m6_m12)
cat("m6_m12 missing >45% data:", mean(apply(m6_m12, 1, function(x) sum(is.na(x)) / length(x)) > 0.45) * 100, "%\n")

# DONT REMOVE HERE 
# Remove rows where more than 45% of the data is missing
BL_6m_clean <- BL_6m[apply(BL_6m, 1, function(x) sum(is.na(x)) / length(x)) <= 0.45, ]
m6_m12_clean <- m6_m12[apply(m6_m12, 1, function(x) sum(is.na(x)) / length(x)) <= 0.45, ]
length(unique(BL_6m_clean$subject_id))
length(unique(m6_m12_clean$subject_id))

########### Make Long DF ##################################################################
# Seeing how many matched there are. 
# Remove the suffixes and compare the base names between the dataframes
remove_suffix <- function(cols) sub("(_BL$|_6m$|_12m$)", "", cols)
# Create a summary table comparing the base names between BL, m6, and m12
summary_table <- data.frame(
  Comparison = c("BL vs m6", "BL vs m12", "m6 vs m12"),
  Matches = c(
    sum(remove_suffix(colnames(BL)) %in% remove_suffix(colnames(m6))),
    sum(remove_suffix(colnames(BL)) %in% remove_suffix(colnames(m12))),
    sum(remove_suffix(colnames(m6)) %in% remove_suffix(colnames(m12)))
  ),
  Non_matches = c(
    length(colnames(BL)) - sum(remove_suffix(colnames(BL)) %in% remove_suffix(colnames(m6))),
    length(colnames(BL)) - sum(remove_suffix(colnames(BL)) %in% remove_suffix(colnames(m12))),
    length(colnames(m6)) - sum(remove_suffix(colnames(m6)) %in% remove_suffix(colnames(m12)))
  )
)
print(summary_table)


# Rename columns by removing the "_BL", "_6m", and "_12m" suffixes
BL_clean <- BL %>% dplyr::select(-c("outcome_BMI_fnl_6m","outcome_BMI_fnl_12m")) %>% 
  rename_with(~ sub("(_BL$|_6m$|_12m$)", "", .)) %>% mutate(time = "BL")

m6_clean <- m6 %>% dplyr::select(-c("outcome_BMI_fnl_BL","outcome_BMI_fnl_12m")) %>% 
  rename_with(~ sub("(_BL$|_6m$|_12m$)", "", .)) %>% mutate(time = "6m")

m12_clean <- m12 %>% dplyr::select(-c("outcome_BMI_fnl_6m","outcome_BMI_fnl_BL")) %>% 
  rename_with(~ sub("(_BL$|_6m$|_12m$)", "", .)) %>% mutate(time = "12m")

# Combine the dataframes using rbind, ensuring columns match
long_df <- bind_rows(BL_clean, m6_clean, m12_clean) %>% unique()
length(long_df$subject_id)
length(unique(long_df$subject_id))
vis_miss(long_df[1:50])
vis_miss(long_df)

aggr_plot <- aggr(long_df[1:50], col=c('navyblue','red'), numbers=TRUE, 
                  sortVars=TRUE, labels=names(long_df), cex.axis=.7, 
                  gap=3, ylab=c("Histogram of missing long data","Pattern"))

cols_with_na <- colnames(long_df)[colSums(is.na(long_df)) > 0]
sapply(long_df[cols_with_na], class) # Check the class of these columns
# Remove rows with NA in specific columns
long_df <- long_df[complete.cases(long_df$outcome_BMI_fnl), ]
table(table(long_df$subject_id))

# Impute missing data using knn 
library(DMwR)
long_df_imputed_knn <- cbind(long_df[, !sapply(long_df, is.numeric)], 
                         knnImputation(long_df[, sapply(long_df, is.numeric)])) %>% 
  dplyr::select(-c(all_samples, all_samples.x, all_samples.y, subject_id.y, time.y))
vis_miss(long_df_imputed_knn[1:10])
length(unique(long_df_imputed_knn$subject_id))
table(table(long_df_imputed_knn$subject_id))

# Impute missing values using missForest <- DONT DO IMPUTATION ALL TOGETHER - Mayne do it up above and only for meta BL 
long_df_imputed <- cbind(long_df[, !sapply(long_df, is.numeric)], 
                         missForest(long_df[, sapply(long_df, is.numeric)])$ximp) %>%
  dplyr::select(-c(all_samples, all_samples.x, all_samples.y, subject_id.y, time.y))

# Visualize missing data in the first 10 rows
vis_miss(long_df_imputed[1:100])
length(unique(long_df_imputed$subject_id))
table(table(long_df_imputed$subject_id))
my_plt_density(long_df_imputed, 10:335, c(-5, 40), "long df bag imputed")

### DONT DO BELOW> ALL SHOULD BE DONE ON INDIVIDUAL OMICS 
# scale longitdinally 
long_df_scaled <- long_df_imputed %>% 
  mutate(across(where(is.numeric), scale))

my_plt_density(long_df_scaled, 10:335, c(-25, 25), "long df bag imputed scaled")

# Center and Scale
long_df_center_scaled <- long_df_imputed %>%
  mutate(across(where(is.numeric), ~ scale(.)))  # Apply scale function column-wise

my_plt_density(long_df_center_scaled, 10:335, c(-25, 25), "long df bag imputed scaled")
check_normality_and_skewness(long_df_center_scaled, 10:361)

# Range Bounds 
preprocess_long_df_range <- preProcess(long_df_imputed, method =  "range", rangeBounds = c(0, 1))
preprocess_long_df_range
long_imputed_range <- predict(preprocess_long_df_range, long_df_imputed)

my_plt_density(long_imputed_range, 10:335, c(-0.15, 1.1), "long df bag imputed scaled range 0-1")
heatmap(cor(long_imputed_range[, 10:361]))
check_normality_and_skewness(long_imputed_range, 10:361)
heatmap(cor(long_imputed_range[, 10:361]))

########### Make delta's ##################################################################

# CHECK PETA LIBRARY - for taxa - CHANGE DATA - genus change CLR 
# fine for micom and meta 
#BL to 6m
colnames_BL_6m <- colnames(BL_6m_clean)
base_names <- sub("_BL$|_6m$", "", colnames_BL_6m)
unique_base_names <- unique(base_names)
change_BL_m6 <- data.frame(matrix(ncol = 0, nrow = nrow(BL_6m_clean)))
# Loop
for (base_name in unique(base_names)) {
  bl_col <- paste0(base_name, "_BL")
  m6_col <- paste0(base_name, "_6m")
  # Check if both "_BL" and "_6m" columns exist in the dataframe
  if (bl_col %in% colnames_BL_6m & m6_col %in% colnames_BL_6m) {
    # Extract the values of the columns and ensure they are numeric
    bl_values <- as.numeric(BL_6m_clean[[bl_col]])
    m6_values <- as.numeric(BL_6m_clean[[m6_col]])
    change_BL_m6[[base_name]] <- m6_values - bl_values # change
  }
}
# Now include columns that do not end in "_BL" or "_6m" unchanged in the change_df
for (col in colnames_BL_6m) {
  # Only add columns that don't match the "_BL" or "_6m" pattern
  if (!grepl("_BL$", col) & !grepl("_6m$", col)) {
    change_BL_m6[[col]] <- BL_6m_clean[[col]]
  }
}
change_BL_m6 <- change_BL_m6 %>% dplyr::select(-c(subject_id.y.x, subject_id.y.y.y, all_samples, all_samples.x,
                                                  all_samples.y, record_id.y, outcome_BMI_fnl_6m.y, outcome_BMI_fnl_BL.y,
                                                  outcome_BMI_fnl_12m.y, time.y, consent.y, cohort_number.y, completer.y,
                                                  sex.y, race.y, age.y, randomized_group.y))
# Remove '.x' from column names
colnames(change_BL_m6) <- gsub("\\.x$", "", colnames(change_BL_m6))
# Remove rows with NA in specific columns
change_BL_m6 <- change_BL_m6[complete.cases(change_BL_m6$outcome_BMI_fnl_BL, change_BL_m6$outcome_BMI_fnl_6m), ]
# Make bmi change var
change_BL_m6$BMI0_6m <- change_BL_m6$outcome_BMI_fnl_6m - change_BL_m6$outcome_BMI_fnl_BL

table(table(change_BL_m6$subject_id))
vis_miss(change_BL_m6[1:20])
vis_miss(change_BL_m6[21:364])

# Impute missing data using miss Forest 
colnames(change_BL_m6) <- make.names(colnames(change_BL_m6))
#change_BL_m6_imputed <- mice(change_BL_m6, method = 'pmm', m = 1)
change_BL_m6_imputed <- cbind(change_BL_m6[, !sapply(change_BL_m6, is.numeric)], 
                         missForest(change_BL_m6[, sapply(change_BL_m6, is.numeric)])$ximp)

change_BL_m6_imputed <- complete(change_BL_m6_imputed) # Get the completed dataset
vis_miss(change_BL_m6_imputed)

my_plt_density(change_BL_m6_imputed, 9:364, c(-5, 5), "BL_m6")
  
### PREPROCESS DELTA
# first just remove low presence/ vaiance and highly correlated vars 
preprocess_BL_6m <- preProcess(change_BL_m6_imputed, 
                               method = c("nzv", "corr"),
                               thresh = 0.95, # % variance
                               na.remove = TRUE, verbose = TRUE, 
                               freqCut = 95/5, # ratio of most common val to 2nd most common val.
                               uniqueCut = 10, # % distinct vals / total no. samples
                               cutoff = 1) # correlation cut off 
preprocess_BL_6m
change_BL_m6_imputed_var_corr <- predict(preprocess_BL_6m, change_BL_m6_imputed)
my_plt_density(change_BL_m6_imputed_var_corr, 9:355, c(-5, 5), "BL_m6 var corr checked")
check_normality_and_skewness(change_BL_m6_imputed_var_corr, 9:355)
setdiff(colnames(change_BL_m6_imputed), colnames(change_BL_m6_imputed_var_corr))

# then center and scale 
change_BL_m6_center_scaled <- change_BL_m6_imputed_var_corr %>%
  mutate(across(where(is.numeric), ~ scale(.)))  # Apply scale function column-wise

my_plt_density(change_BL_m6_center_scaled, 10:355, c(-5, 5), "BL-6m imputed var corr center scaled")
check_normality_and_skewness(change_BL_m6_center_scaled, 10:355)

# Maybe convert to a range 
preprocess_BL_6m_range <- preProcess(change_BL_m6_imputed_var_corr, 
                                       method =  c("range"), rangeBounds = c(0, 1))
preprocess_BL_6m_range
BL_6m_range <- predict(preprocess_BL_6m_range, change_BL_m6_imputed_var_corr)

my_plt_density(BL_6m_range, 10:355, c(-0.1, 1.1), "BL-6m imputed var corr range")
check_normality_and_skewness(BL_6m_range, 10:355)
heatmap(cor(BL_6m_range[, 10:355]))

### Repeat for 6m - 12m ########################################################
colnames_m6_m12 <- colnames(m6_m12_clean)
base_names <- sub("_6m$|_12m$", "", colnames_m6_m12)
unique_base_names <- unique(base_names)
change_m6_m12 <- data.frame(matrix(ncol = 0, nrow = nrow(m6_m12_clean)))
# Loop
for (base_name in unique(base_names)) {
  m6_col <- paste0(base_name, "_6m")
  m12_col <- paste0(base_name, "_12m")
  # Check if both "_BL" and "_6m" columns exist in the dataframe
  if (m6_col %in% colnames_m6_m12 & m12_col %in% colnames_m6_m12) {
    # Extract the values of the columns and ensure they are numeric
    m6_values <- as.numeric(m6_m12_clean[[m6_col]])
    m12_values <- as.numeric(m6_m12_clean[[m12_col]])
    change_m6_m12[[base_name]] <- m12_values - m6_values # change
  }
}
# Now include columns that do not end in "_BL" or "_6m" unchanged in the change_df
for (col in colnames_m6_m12) {
  # Only add columns that don't match the "_BL" or "_6m" pattern
  if (!grepl("_6m$", col) & !grepl("_12m$", col)) {
    change_m6_m12[[col]] <- m6_m12_clean[[col]]
  }
}

change_m6_m12 <- change_m6_m12 %>% 
  dplyr::select(-c(subject_id.y.x, subject_id.y.y.y, all_samples, all_samples.x,
                   all_samples.y, record_id.y, outcome_BMI_fnl_6m.y, outcome_BMI_fnl_BL.y,
                   outcome_BMI_fnl_12m.y, time.y, consent.y, cohort_number.y, completer.y,
                   sex.y, race.y, age.y, randomized_group.y))
# Remove '.x' from column names
colnames(change_m6_m12) <- gsub("\\.x$", "", colnames(change_m6_m12))
# Remove rows with NA in specific columns
change_m6_m12 <- change_m6_m12[complete.cases(change_m6_m12$outcome_BMI_fnl_12m, change_m6_m12$outcome_BMI_fnl_6m), ]
# make change var 
change_m6_m12$BMI6_12m <- change_m6_m12$outcome_BMI_fnl_12m - change_m6_m12$outcome_BMI_fnl_6m

# impute final missing with miss Forest 
colnames(change_m6_m12) <- make.names(colnames(change_m6_m12))
#change_m6_m12_imputed <- mice(change_m6_m12, method = 'pmm', m = 1)
change_m6_m12_imputed <- cbind(change_m6_m12[, !sapply(change_m6_m12, is.numeric)], 
                              missForest(change_m6_m12[, sapply(change_m6_m12, is.numeric)])$ximp)

change_m6_m12_imputed <- complete(change_m6_m12_imputed) # Get the completed dataset

table(table(change_m6_m12_imputed$subject_id))
vis_miss(change_m6_m12_imputed[1:10])
colnames(change_m6_m12_imputed)[colSums(is.na(change_m6_m12_imputed)) > 0]
my_plt_density(change_m6_m12_imputed, 9:364, c(-5, 5), "m6_m12")

### PREPROCESS DELTA
# first just remove low presence/ vaiance and highly correlated vars 
preprocessL_6m_12m <- preProcess(change_m6_m12_imputed, 
                               method = c("nzv", "corr"),
                               thresh = 0.95, # % variance
                               na.remove = TRUE, verbose = TRUE, 
                               freqCut = 95/5, # ratio of most common val to 2nd most common val.
                               uniqueCut = 10, # % distinct vals / total no. samples
                               cutoff = 1) # correlation cut off 
preprocessL_6m_12m
change_m6_12m_imputed_var_corr <- predict(preprocessL_6m_12m, change_m6_m12_imputed)
my_plt_density(change_m6_12m_imputed_var_corr, 9:355, c(-5, 5), "m6-m12var corr checked")
check_normality_and_skewness(change_m6_12m_imputed_var_corr, 9:355)
setdiff(colnames(change_m6_m12_imputed), colnames(change_m6_12m_imputed_var_corr))

# then center and scale 
change_m6_m12_center_scaled <- change_m6_12m_imputed_var_corr %>%
  mutate(across(where(is.numeric), ~ scale(.)))  # Apply scale function column-wise

my_plt_density(change_m6_m12_center_scaled, 10:355, c(-5, 5), "6m-12m imputed var corr center scaled")
check_normality_and_skewness(change_m6_m12_center_scaled, 10:355)

# Maybe convert to a range 
preprocess_6m_12m_range <- preProcess(change_m6_12m_imputed_var_corr, 
                                     method =  c("range"), rangeBounds = c(0, 1))
preprocess_6m_12m_range
m6_m12_range <- predict(preprocess_6m_12m_range, change_m6_12m_imputed_var_corr)

my_plt_density(m6_m12_range, 10:355, c(-0.1, 1.1), "6m-12m imputed var corr range")
check_normality_and_skewness(m6_m12_range, 10:355)

heatmap(cor(m6_m12_range[, 10:355]))

# Save files made
# Define the base path
base_path <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/march_20/"
# List of filenames and corresponding data frames
files <- list(
  list(filename = "long_df_imputed.csv", data = long_df_imputed),
  list(filename = "long_df_imputed_cent_scale.csv", data = long_df_center_scaled),
  list(filename = "long_df_imputed_range.csv", data = long_imputed_range),
  # BL-6m
  list(filename = "BL_6m.csv", data = change_BL_m6_imputed_var_corr),
  list(filename = "BL_6m_cent_scale.csv", data = change_BL_m6_center_scaled),
  list(filename = "BL_6m_range.csv", data = BL_6m_range),
  # 6m-12m
  list(filename = "m6_m12.csv", data = change_m6_12m_imputed_var_corr),
  list(filename = "m6_m12_cent_scale.csv", data = change_m6_m12_center_scaled),
  list(filename = "m6_m12_range.csv", data = m6_m12_range)
)
lapply(files, function(x) write.csv(x$data, file.path(base_path, x$filename), row.names = FALSE))



# Save them 
#march_20 <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/march_20/all_processing_long.csv"
#write.csv(long_df_imputed, march_20, row.names = FALSE)
#march_25 <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/march_20/all_processing_0_6m.csv"
#write.csv(change_BL_m6_imputed, march_25, row.names = FALSE)
#march_30 <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/march_20/all_processing_6_12m.csv"
#write.csv(change_m6_m12_imputed, march_30, row.names = FALSE)









