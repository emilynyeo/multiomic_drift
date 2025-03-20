# March 19 
# Combined all input processing script 

pacman::p_load(knitr, data.table, dplyr, tidyr, tableone, kableExtra, readxl,
               readr, car, RColorBrewer, gridExtra, mlbench, earth, ggplot2, 
               AppliedPredictiveModeling, caret, reshape2, corrplot, stringr,
               summarytools, grid, mice, plyr, mlmRev, cowplot, ape, e1071,
               jtools, broom, patchwork, phyloseq, microbiome, glmnet, ISLR,
               MicrobiomeStat, ANCOMBC, ape, vegan, zCompositions, janitor,
               RColorBrewer, DT, ggpubr, microbiomeutilities, compositions)
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

#### META ###############################################################################################
meta <- read_csv("~/projects/research/Stanislawski/BMI_risk_scores/data/correct_meta_files/ashleys_meta/DRIFT_working_dataset_meta_deltas_filtered_05.21.2024.csv")

# Replace spaces with underscores in column names
colnames(meta) <- gsub(" ", "_", colnames(meta))
length(unique(meta$subject_id))
meta <- meta[meta$consent != "no", ]
length(unique(meta$subject_id))
a2_extra <- meta %>% dplyr::select(c(record_id, subject_id, randomized_group, consent,
                                     cohort_number, sex, race, completer, age,
                                     # BL
                                     outcome_BMI_fnl_BL, Glucose_BL, HOMA_IR_BL, 
                                     Insulin_endo_BL, HDL_Total_Direct_lipid_BL, 
                                     LDL_Calculated_BL, Triglyceride_lipid_BL,
                                     Peptide_YY_BL, Ghrelin_BL, Leptin_BL, Hemoglobin_A1C_BL,
                                     # 6m
                                     outcome_BMI_fnl_6m, Glucose_6m, HOMA_IR_6m,
                                     Insulin_endo_6m, HDL_Total_Direct_lipid_6m, 
                                     LDL_Calculated_6m, Triglyceride_lipid_6m,
                                     Peptide_YY_6m, Ghrelin_6m, Leptin_6m, Hemoglobin_A1C_6m,
                                     # 12m
                                     outcome_BMI_fnl_12m, Glucose_12m, HOMA_IR_12m,
                                     Insulin_endo_12m, HDL_Total_Direct_lipid_12m,
                                     LDL_Calculated_12m, Triglyceride_lipid_12m,
                                     Peptide_YY_12m, Ghrelin_12m, Leptin_12m, Hemoglobin_A1C_12m))

my_plt_density(a2_extra, 9:42, c(-2, 150), "Meta Count")
check_normality_and_skewness(genus_count, 9:42)

preProcValues <- preProcess(a2_extra[,9:42], 
                            method = c("scale", #"center", 
                                       "nzv"),
                            thresh = 0.95, pcaComp = NULL, na.remove = TRUE, 
                            k = 5, knnSummary = mean, fudge = 0.2, 
                            numUnique = 15, verbose = TRUE,  freqCut = 95/5, 
                            uniqueCut = 10, cutoff = 0.85, rangeBounds = c(0, 1))

meta_cs <- predict(preProcValues, a2_extra[,1:42]) 
my_plt_density(meta_cs, 9:42, c(-10, 100), "META DATA CENTER SCALED")
check_normality_and_skewness(meta_cs, 9:42)

### GRS ####################################################################################################
grs <- read_csv("/Users/emily/projects/research/Stanislawski/BMI_risk_scores/full_cohort_pulling_snps/bigsnpr/made_scores/merge_meta_methyl.csv") %>%
  dplyr::select(c(subject_id, raw_score)) %>%
  dplyr::rename(bmi_prs = raw_score)

### Taxa #################################################################################################### 
load("~/projects/research/Stanislawski/BMI_risk_scores/microbiome_rs/data/PhyloseqObj.RData")
load("/Users/emily/projects/research/Stanislawski/BMI_risk_scores/microbiome_rs/data/Genus_Sp_tables.RData")
print_ps(drift.phy.count)
# Remove taxa not seen more than 2 times in at least 20% of the samples.
GP_count = filter_taxa(drift.phy.count, 
                       function(x) sum(x > 3) > (0.1*length(x)), TRUE)
print(paste0("#Taxa with >3 counts in > 10% samples if using count: ", ntaxa(GP_count)))

# Coefficient of variation cut off 
gpsf_count = filter_taxa(GP_count, function(x) sd(x)/mean(x) > 1.0, TRUE)
print(paste0("# Taxa with coefficient of variation > 10 if using count: ", ntaxa(gpsf_count)))

# Get gennus 
tax_genus <- tax_glom(gpsf_count, "Genus")
genus_filtered_count <- otu_table(tax_genus) 
genus_filtered_count <- t(genus_filtered_count) %>% as.data.frame()
colnames(genus_filtered_count) <- as.data.frame(tax_table(tax_genus))$Genus
rm(drift.phy.clr,drift.phy.count,GP_count,gpsf_count,sp.clr, tax_genus,sp.count,
   drift.phy.count.r21116, drift.phy.ra, genus.clr, genus.count, genus.ra, sp.ra)

# Function to make column names unique
colnames(genus_filtered_count) <- make_unique_names(colnames(genus_filtered_count)) # rename the columns
genus_count <- genus_filtered_count %>% rownames_to_column(var = "subject_id") 
my_plt_density(genus_count, 2:180, c(-10, 10), "Genus Count")
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
my_plt_density(Genus_relative_abundance, 2:179, c(-0.001, 0.01), "Genus RA")
check_normality_and_skewness(Genus_relative_abundance, 2:179)

# CLR Transformation
genus_clr_transformed <- apply(Genus_relative_abundance, 2, clr) %>% as.data.frame()
my_plt_density(genus_clr_transformed, 1:179, c(-0.005, 0.005), "Genus RA CLR")
check_normality_and_skewness(genus_clr_transformed, 1:179)

# CENTER AND SCALE
clr_transformed_0 <- genus_clr_transformed 
clr_transformed_0$SampleID <- rownames(genus_clr_transformed)
preProcValues <- preProcess(genus_clr_transformed[,1:179], 
                            method = c("scale", #"center", 
                                       "nzv"),
                            thresh = 0.95, pcaComp = NULL, na.remove = TRUE, 
                            k = 5, knnSummary = mean, fudge = 0.2, 
                            numUnique = 15, verbose = TRUE,  freqCut = 95/5, 
                            uniqueCut = 10, cutoff = 0.85, rangeBounds = c(0, 1))

clr_transformed <- predict(preProcValues, genus_clr_transformed[, 1:179]) %>% 
                   rownames_to_column(var = "subject_id")
my_plt_density(clr_transformed, 2:180, c(-0.1, 0.1), "Genus RA CLR CENTER SCALED")
check_normality_and_skewness(clr_transformed, 2:180)
# Process names 
clr_transformed$all_samples <- process_names_all(rownames(clr_transformed))

### Functional ####################################################################################################
pathways <- fread("/Users/emily/projects/research/Stanislawski/BMI_risk_scores/picrust2/june7/pathways_out/path_abun_unstrat_descrip.tsv") %>% .[, -1] %>% t() %>% as.data.frame() %>% 
  row_to_names(1) %>% rownames_to_column("SampleID") %>% 
  mutate(across(-1, as.numeric))

# Remove variables with low presence thresholds 
threshold <- 0.20
zero_percentage <- colSums(pathways == 0, na.rm = TRUE) / nrow(pathways) # % zeros / column
path_all_time_cleaned <- pathways[, zero_percentage < threshold] # only colz < 20% zeros
dim(pathways) - dim(path_all_time_cleaned) 
any(path_all_time_cleaned < 0, na.rm = TRUE)
my_plt_density(path_all_time_cleaned, 2:267, c(-10, 10), "Pathways")
check_normality_and_skewness(path_all_time_cleaned, 2:267)

# imputation and scaling 
preProcValues <- preProcess(path_all_time_cleaned[,2:267], 
                            method = c("scale", #"center", 
                                       "nzv"),
                            thresh = 0.95, pcaComp = NULL, na.remove = TRUE, 
                            k = 5, knnSummary = mean, fudge = 0.2, 
                            numUnique = 15, verbose = TRUE,  freqCut = 95/5, 
                            uniqueCut = 10, cutoff = 0.85, rangeBounds = c(0, 1))

path_all_time_cs <- predict(preProcValues, path_all_time_cleaned[,1:267])
my_plt_density(path_all_time_cs, 2:266, c(-1, 10), "Pathways Center Scaled")
check_normality_and_skewness(path_all_time_cs, 2:266)
# Process samples 
path_all_time_cs$all_samples <- process_names_all(path_all_time_cs$SampleID)
# Make pathway names unique
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

# Remove variables with low presence thresholds 
threshold <- 0.20
zero_percentage <- colSums(flux == 0, na.rm = TRUE) / nrow(flux)
flux_cleaned <- flux[, zero_percentage < threshold]
my_plt_density(flux_cleaned, 2:104, c(-0.1, 0.1), "MICOM raw")
check_normality_and_skewness(flux_cleaned, 2:104)

log_all_micom <- flux_cleaned %>% mutate(across(-1, log1p))
my_plt_density(log_all_micom, 2:104, c(-0.01, 0.01), "MICOM logged")
check_normality_and_skewness(log_all_micom, 2:104)

# Preprocess the data - imputation and Centering
preProcValues <- preProcess(flux_cleaned[,2:104], 
                            method = c("scale", #"center", 
                                       "nzv"),
                            thresh = 0.95, pcaComp = NULL, na.remove = TRUE, 
                            k = 5, knnSummary = mean, fudge = 0.2, 
                            numUnique = 15, verbose = TRUE,  freqCut = 95/5, 
                            uniqueCut = 10, cutoff = 0.85, rangeBounds = c(0, 1))

flux_all_time_cs <- predict(preProcValues, flux_cleaned[,1:104])
my_plt_density(flux_all_time_cs, 2:92, c(-0.1, 0.1), "MICOM Center Scaled")
check_normality_and_skewness(flux_all_time_cs, 2:92)

# Process names
flux_all_time_cs$all_samples <- process_names_all(flux_all_time_cs$sample_id)
colnames(flux_all_time_cs) <- make.names(colnames(flux_all_time_cs), unique = TRUE)
flux_all_time_cs <- flux_all_time_cs %>%
  mutate(time = case_when(grepl("BL", sample_id) ~ 0, grepl("6m", sample_id) ~ 6, 
                          grepl("12m", sample_id) ~ 12, TRUE ~ NA_real_)) %>%
  filter(!grepl("\\.3m$", sample_id))
