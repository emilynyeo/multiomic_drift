# Trying out long lasso
rm(list = ls())
#source("zc_functions.R") 
#install.packages("rsq")  
library(rsq)
library(pacman)
p_load(tools, reticulate, viridis, tidyplots, patchwork, jsonlite, maps, ggvenn, 
       caret, caretEnsemble, glmnet, xgboost, ggplot2, glmmLasso, corrplot,
       readr, plyr, dplyr, tidyr, purrr, tibble, stringr, psych, randomForest,  
       reshape2, scales, gridExtra, plotly, sf, tidyverse, naniar, VIM, kableExtra)
#library(nlme)
library(gridExtra)
library(sjPlot)
library(htmltools)
library(officer)
library(flextable)
library(webshot)
library(apaTables)
library(kneedle)
library(MuMIn)
'%ni%' <- Negate('%in%')
r2_general <-function(preds,actual){ 
  return(1- sum((preds - actual) ^ 2)/sum((actual - mean(actual))^2))
}

#out_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_delta/new_split/"
out_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/delta/"

#data_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/zachs_rerun/drift_fs/csv/all_omic_processed_data/deltas/"
data_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/april_processing/"
all_deltas <- read_csv(paste0(data_dir, "all_delta.csv")) %>% 
                       dplyr::select(-c("...1", "consent", "completer", 
                                        "Peptide_YY", "Ghrelin", "Leptin")) %>%
  dplyr::mutate(time = as.factor(time),
                subject_id = as.factor(subject_id),
                randomized_group = as.factor(randomized_group),
                sex = as.numeric(sex),
                race = as.numeric(race)) %>% rename(BMI = outcome_BMI_fnl,
                                                    range = time,
                                                    homo_ir = HOMA_IR,
                                                    insulin = Insulin_endo,
                                                    LDL = LDL_Calculated,
                                                    HDL = HDL_Total_Direct_lipid,
                                                    HbA1c = Hemoglobin_A1C)

# Define the column names based on your lists
basic <- c('subject_id','BMI', 'range', 'sex', 'age', 'randomized_group') #
meta_keep <- c('subject_id','BMI', 'range', 'randomized_group', 'sex', 'race', 
               'age', 'HbA1c', 'HDL', 'homo_ir', 'insulin', 'LDL', 'Glucose.x') # 
meta_keep_no_age_sex <- c('subject_id','BMI', 'range', 'race', 
                          'HbA1c', 'HDL', 'homo_ir', 'insulin', 'LDL', 'Glucose.x')
only_grs <- c('subject_id','BMI', 'range','bmi_prs')
only_taxa <- c('subject_id','BMI', 'range', grep("^g__", names(all_deltas), value = TRUE))

micom_start <- which(names(all_deltas) == "Diacetyl")
micom_end <- which(names(all_deltas) == "aldehydo.D.xylose")
only_micom <- c('subject_id','BMI', 'range', names(all_deltas)[micom_start:micom_end])

path_start <- which(names(all_deltas) == "arginine..ornithine.and.proline.interconversion")
path_end <- which(names(all_deltas) == "UDP.N.acetyl.D.glucosamine.biosynthesis.I")
only_pathway <- c('subject_id','BMI', 'range', names(all_deltas)[path_start:path_end])

tabo_start <- which(names(all_deltas) == "non_HDL_C")
tabo_end <- which(names(all_deltas) == "IDL_TG_pct")
only_tabo <- c('subject_id','BMI', 'range', names(all_deltas)[tabo_start:tabo_end])

all_col <- c('subject_id','BMI', 'range',
             'randomized_group', 'sex', 'race', 
             'age', 'HbA1c', 'HDL', 'homo_ir', 'insulin', 'LDL', 'Glucose.x',
             grep("^g__", names(all_deltas), value = TRUE),
             names(all_deltas)[micom_start:micom_end],
             names(all_deltas)[path_start:path_end],
             names(all_deltas)[tabo_start:tabo_end])

all_col_no_age_sex <- c('subject_id','BMI', 'range','race', 
                        'HbA1c', 'HDL', 'homo_ir', 'insulin', 'LDL', 'Glucose.x',
                        grep("^g__", names(all_deltas), value = TRUE),
                        names(all_deltas)[micom_start:micom_end],
                        names(all_deltas)[path_start:path_end],
                        names(all_deltas)[tabo_start:tabo_end])

all_col_no_clin <- c('subject_id','BMI', 'range',
                     grep("^g__", names(all_deltas), value = TRUE),
                     names(all_deltas)[micom_start:micom_end],
                     names(all_deltas)[path_start:path_end],
                     names(all_deltas)[tabo_start:tabo_end])

# Create data frames based on the columns defined
basic <- all_deltas[, basic, drop = FALSE] %>% unique()
meta <- all_deltas[, meta_keep, drop = FALSE] %>% unique()
meta_no_age_sex <- all_deltas[, meta_keep_no_age_sex, drop = FALSE] %>% unique()
grs <- all_deltas[, only_grs, drop = FALSE] %>% unique()
taxa <- all_deltas[, only_taxa, drop = FALSE] %>% unique()
micom <- all_deltas[, only_micom, drop = FALSE] %>% unique()
pathway <- all_deltas[, only_pathway, drop = FALSE] %>% unique()
metabo <- all_deltas[, only_tabo, drop = FALSE] %>% unique()
all <- all_deltas[, all_col, drop = FALSE] %>% unique() %>% 
       dplyr::mutate(randomized_group = as.numeric(randomized_group))
all_no_age_sex <- all_deltas[, all_col_no_age_sex, drop = FALSE] %>% unique()
all_no_clin <- all_deltas[, all_col_no_clin, drop = FALSE] %>% unique() 

# Make train and test set
# Test sample names
#test_names <- c("ABR-079", "AGA-071", "AHE-055", "ALI-121", "ALO-163", "AMA-031", "ASO-013", "AWI-167", "BMO-164", "CWA-183", "DSC-024", "EBE-130", "EHI-177", "EJO-092", "GFU-188", "HGI-010", "JCA-109", "JGO-100", "KBU-085", "KCE-034", "KHE-170", "LDO-148", "LST-186", "LZD-142", "MAR-119", "MCA-088", "MJA-153", "MWE-112", "NPO-149", "RAE-114", "SBO-020", "SEG-080", "SKA-195", "SLO-178", "SSH-028", "TDU-086","TFA-016", "VCA-041")

# Taking out the odd one : "AHE-055"

test_names <- c("ASO-013", "NTA-021", "KGI-029", "KPA-042", "AWA-052",  "COW-066", "NBI-069", "CEL-073", "CAL-074", "ABR-079", "SEG-080", "NKA-090", "NEL-094", "LJA-101", "ADA-105", "MLU-106", "MDI-107", "JER-110", "TRO-113", "MFB-118", "ALI-121", "KWA-122", "RAF-125", "EBE-130", "CGA-134", "LZD-142", "NPO-149", "HDE-154", "AMC-155", "SAB-160", "QNG-166", "NCO-171", "BSA-174", "EHI-177", "LST-186", "MBA-187", "BAN-193", "AHE-055")

# Train sample names
#train_names <- c("AAL-144", "ACO-053", "ADA-105", "AKE-009", "AKI-011", "AKO-139", "AMC-155", "AME-128", "AME-157", "ATA-129", "AWA-052", "AWA-083", "BAN-193", "BHO-014", "BIN-201", "BKN-104", "BMI-156", "BSA-174", "CAM-057", "CCO-189", "CED-026", "CEL-073", "CGA-134", "CIS-077", "CKR-078", "CLE-049", "COW-066", "CRO-108", "CWA-161", "EBE-051", "EKA-135", "EKR-045", "ELA-159", "EPO-182", "EVO-184", "FWI-098", "GHA-035", "HDE-154", "IBE-120", "JDI-140", "JER-110", "JFU-027", "JJO-093", "JKN-127", "JPO-022", "JUG-116", "JUT-032", "JVE-126", "KAN-138", "KBR-162", "KEL-185", "KEL-199", "KGI-029", "KHU-196", "KPA-042", "KRI-072", "KVA-038", "KWA-122", "KWA-141", "LBL-047", "LBU-015", "LEL-147", "LFI-003", "LJA-101", "LMC-111", "LPF-198", "LVA-017", "MBA-187", "MCW-065", "MDI-107", "MES-068", "MFB-118", "MGA-076", "MHO-117", "MKE-192", "MMA-036", "MRT-179", "MSH-091", "MST-039", "MWE-143", "MWO-133", "MWY-152", "NAR-099", "NBI-048", "NBI-069", "NCO-171", "NDI-067", "NEL-094", "NKA-090", "NMO-151", "NTA-021", "PBE-123", "QNG-166", "RAF-125", "RAM-050", "RHP-023", "RLA-132", "ROL-006", "SAB-160", "SCA-043", "SCR-061", "SDA-150", "SGA-062", "SKA-087", "SRO-194", "TBU-115", "TFA-172", "TRO-113", "TSH-146", "TSL-056", "WPE-005", "YOR-103", "YSU-097", "ZVU-096")

train_names <- c("SDA-150", "LBU-015", "CIS-077", "ATA-129", "KHU-196", "MWY-152", "AGA-071", "AME-157", "CWA-183", "RHP-023", "MST-025", "SSH-028", "JUG-116", "EJO-092", "VCA-041", "NMO-151", "BHO-014", "KBU-085", "SBO-020", "MWO-133", "KRI-072", "AAL-144", "ALO-163", "AKI-011", "MHO-117", "TSH-146", "RAE-114", "FWI-098", "MAR-119", "JGO-100", "CAM-057", "YOR-103", "HGI-010", "KAN-138", "SGA-062", "CKR-078", "MWE-112", "ROL-006", "MMA-036", "DSC-024", "LDO-148", "MCA-088", "CPU-075", "AKO-139", "LFI-003", "KWA-141", "GFU-188", "BMO-164", "JPO-022", "EVO-184", "LPF-198", "TBU-115", "SRO-194", "KEL-199", "JFU-027", "SKA-195", "IBE-120", "TSL-056", "NDI-067", "AWA-083", "CWA-161", "TDU-086", "JCA-109", "CBO-004", "NAR-099", "MES-068", "AMA-031", "SLO-178", "SCA-043", "AWI-167",  "KBR-162", "TFA-172", "BIN-201", "NBI-048", "KHE-170", "CSH-012", "BMI-156", "MWE-143", "EKA-135", "WPE-005", "AKE-009", "YSU-097", "MCW-065", "EBE-051", "ZVU-096", "JJO-093", "KVA-038", "ACO-053", "RLA-132", "MBA-176", "CED-026", "JDI-140", "CCO-189", "EKR-045", "MJA-153", "CLE-049", "LMC-111", "SKA-087", "JUT-032", "MKE-192", "JVE-126", "KCE-034", "KEL-185", "MRT-179", "JKN-127", "LEL-147", "BKN-104", "AME-128", "MSH-091", "MGA-076", "LVA-017", "EPO-182")

cat("Length of test names:", length(test_names), "\n")
cat("Length of train names:", length(train_names), "\n")

subject_id_count <- meta %>%
  dplyr::filter(range %in% c("0","6", "12")) %>%
  dplyr::group_by(subject_id) %>%
  dplyr::summarize(range_count = n_distinct(range))  # Count the distinct range values
table(subject_id_count$range_count)
missing_subjects  <- subject_id_count %>% dplyr::filter(range_count != 3)

# CHECK ALL CORRELATIONS
preProcValues <- preProcess(all, 
                            method = c("nzv", "corr"), thresh = 0.95, fudge = 0.2, 
                            numUnique = 15, verbose = TRUE, freqCut = 95/5, 
                            uniqueCut = 10, cutoff = 0.70, na.remove = TRUE)
preProcValues
all <- predict(preProcValues, all)
heatmap(cor(all[, c(2, 5:ncol(all))]))

## Check all_no_age_sex correlation
preProcValues_train <- preProcess(all_no_age_sex[, c(2, 4:ncol(all_no_age_sex))], 
                                  method = c("nzv", "corr"), thresh = 0.95, fudge = 0.2, 
                                  numUnique = 15, verbose = TRUE, freqCut = 95/5, 
                                  uniqueCut = 10, cutoff = 0.70, na.remove = TRUE)
preProcValues_train
all_no_age_sex <- predict(preProcValues_train, all_no_age_sex)

## Check all_no_clin correlation
preProcValues_train <- preProcess(all_no_clin[, c(2, 4:ncol(all_no_clin))], 
                                  method = c("nzv", "corr"), thresh = 0.95, fudge = 0.2, 
                                  numUnique = 15, verbose = TRUE, freqCut = 95/5, 
                                  uniqueCut = 10, cutoff = 0.70, na.remove = TRUE)
preProcValues_train
all_no_clin <- predict(preProcValues_train, all_no_clin)

# Make test and tain sets for each omic 
data_frames <- c("basic", "meta", "meta_no_age_sex", "grs", "taxa", "pathway", 
                 "micom", "metabo", "all", "all_no_age_sex", "all_no_clin")
for (df in data_frames) {
  df_data <- get(df)  # Get the data frame by name
  df_data <- df_data %>% dplyr::filter(!df_data$subject_id %in% missing_subjects)  # Filter rows
  # Split into training and test sets
  train_set <- df_data[df_data$subject_id %in% train_names, ]
  test_set <- df_data[df_data$subject_id %in% test_names, ]
  # Assign the new sets back to their respective variables if needed
  assign(paste0(df, "_train"), train_set)
  assign(paste0(df, "_test"), test_set)
}

# CHECK MICOM CORRELATIONS
preProcValues <- preProcess(micom_train[, c(2, 4:ncol(micom_train))], 
                            method = c("nzv", "corr"), thresh = 0.95, fudge = 0.2, 
                            numUnique = 15, verbose = TRUE, freqCut = 95/5, 
                            uniqueCut = 10, cutoff = 0.90, na.remove = TRUE)
preProcValues
micom_train <- predict(preProcValues, micom_train)
heatmap(cor(micom_train[, c(2, 4:ncol(micom_train))]))

# CHECK ALL CORRELATIONS
preProcValues <- preProcess(all_train, 
                            method = c("nzv", "corr"), thresh = 0.95, fudge = 0.2, 
                            numUnique = 15, verbose = TRUE, freqCut = 95/5, 
                            uniqueCut = 10, cutoff = 0.75, na.remove = TRUE)
preProcValues
all_train <- predict(preProcValues, all_train)
heatmap(cor(all[, c(2, 5:ncol(all_train))]))

# CHECK ALL_NO_AGE_SEX CORRELATIONS
preProcValues <- preProcess(all_no_age_sex_train, 
                            method = c("nzv", "corr"), thresh = 0.95, fudge = 0.2, 
                            numUnique = 15, verbose = TRUE, freqCut = 95/5, 
                            uniqueCut = 10, cutoff = 0.75, na.remove = TRUE)
preProcValues
all_no_age_sex_train <- predict(preProcValues, all_no_age_sex_train)

# CHECK ALL_NO_CLINICAL CORRELATIONS
preProcValues <- preProcess(all_no_clin_train, 
                            method = c("nzv", "corr"), thresh = 0.95, fudge = 0.2, 
                            numUnique = 15, verbose = TRUE, freqCut = 95/5, 
                            uniqueCut = 10, cutoff = 0.75, na.remove = TRUE)
preProcValues
all_no_clin_train <- predict(preProcValues, all_no_clin_train)

# STEP 1
#data_frames <- c("all") # "basic",
data_frames <- c("basic", "meta", "meta_no_age_sex", "grs", "taxa", "pathway", 
                 "micom", "metabo", "all", "all_no_age_sex", "all_no_clin") #
for (df_name in data_frames) {
  train_data <- get(paste0(df_name, "_train"))
  test_data <- get(paste0(df_name, "_test"))
  numvar <- c() 
  #lambdavec <- seq(from = 25, to = 100, by = 1)
  
  # Use smaller lambda sweep for basic and grs, default for others
  if (df_name %in% c("grs", "basic")) {
    lambdavec <- seq(from = 0.01, to = 1, by = 0.05)
  } else {
    lambdavec <- seq(from = 25, to = 100, by = 1)
  }
  
  # Loop through each lambda value to perform Lasso regression
  for (lambdy in lambdavec) {
    predictors <- setdiff(names(train_data), c("BMI", "subject_id", "range"))
    predictors_escaped <- paste0("`", predictors, "`", collapse = " + ")
    fix_formula <- as.formula(paste("BMI ~", predictors_escaped))
    
    # Run the Lasso regression with GLMM
    lm1 <- glmmLasso(fix = fix_formula,
                     data = train_data,
                     rnd = list(subject_id = ~ 1, range = ~ 1),
                     lambda = lambdy,
                     family = gaussian(link = "identity"))
    summary(lm1)
    # Get the non-zero lasso features
    lassoFeatures <- names(lm1$coefficients[which(lm1$coefficients != 0)])
    lassoFeatures <- lassoFeatures[lassoFeatures %ni% c("(Intercept)")]
    lassoFeatures <- unique(c(lassoFeatures))
    numvar <- c(numvar, length(lassoFeatures)) # Store the number of variables included for this lambda
  }
  lam <- plot(lambdavec, numvar,  xlab = "lambda", ylab = "numvars", 
              pch=21, col="blue", bg="lightblue", type = "b")
  lam
  #lambda_value <- kneedle(lambdavec, numvar)
  lambda_value <- tryCatch(
    {
      kneedle(lambdavec, numvar)
    },
    error = function(e) {
      message("kneedle() failed: ", e$message)
      25
    }
  )
  lambda_value <- lambda_value[1]
  assign(paste0(df_name, "lambda_value"), lambda_value)
  
  # STEP 2 Create the final model on the training data using the chosen lambda
  best_model <- glmmLasso(fix = fix_formula,
                          data = train_data,
                          rnd = list(subject_id = ~ 1, range = ~ 1),
                          lambda = lambda_value,  # Use the lambda_value set manually
                          family = gaussian(link = "identity"))
  assign(paste0(df_name, "_best_model"), best_model)
  
  # Extract lasso features
  lassoFeatures <- names(best_model$coefficients)
  lassoFeatures <- lassoFeatures[lassoFeatures %ni% c("(Intercept)")]
  lassoFeatures <- unique(c(lassoFeatures))
  lassoFeatures <- lassoFeatures[grep("subject_id|BMI|range", lassoFeatures, invert = TRUE)]  # Exclude specific features
  lassoFeatures <- unique(c(lassoFeatures, "subject_id", "BMI", "range"))
  assign(paste0(df_name, "_lassoFeatures"), lassoFeatures)
  
  # Extract coefficients and plot the top features
  coef_df <- as.data.frame(summary(best_model)$coefficients)
  coef_df$Feature <- rownames(coef_df)
  coef_df <- coef_df %>% arrange(desc(Estimate))  # Sort by coefficients
  assign(paste0(df_name, "_coef_df"), coef_df)
  
  # Filter top 10 features by absolute coefficient magnitude
  top_features <- coef_df %>%
    mutate(abs_estimate = abs(Estimate)) %>%
    dplyr::filter(Feature != "(Intercept)") %>%
    mutate(Feature = str_to_title(Feature),  # Capitalize the first letter of each word
           Feature = str_replace_all(Feature, "[._]", " ")) %>%
    arrange(desc(abs_estimate)) 
  
  write.csv(top_features, file = paste0(out_dir, paste0(df_name, "_gl_delta_top_features.csv")), row.names = FALSE)
  
  plot_feat <- top_features %>%
    dplyr::filter(!str_detect(tolower(Feature), "time")) %>%  # Remove rows with "time" in Feature
    slice_head(n = 10) 
  
  # Plot top features
  features <- ggplot(plot_feat, aes(x = reorder(Feature, Estimate), y = Estimate)) +
    geom_bar(stat = "identity", fill = "#1C4C98") +
    coord_flip() + theme_bw() +
    ggtitle(paste("Top Features & Coefficients from", df_name)) +
    xlab("Feature") + ylab("Coefficient Estimate") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12))
  features 
  
  # Save the plot using ggsave
  ggsave(filename = paste0(out_dir, paste0(df_name, "_top_deltas_features.png")),
         plot = features, width = 8, height = 6, units = "in", dpi = 300)
  
  # STEP 3:  Predict on test data 
  test_data <- test_data[complete.cases(test_data),] %>% as.data.frame() # Filter complete cases
  pred_risk_scores <- predict(best_model, test_data)
  range_long <- test_data$range
  sub_ids <- test_data$subject_id
  actual_bmi <- test_data$BMI
  pred_df <- data.frame(
    subject_id = as.character(sub_ids),
    time = as.numeric(range_long), 
    actual = as.numeric(actual_bmi),
    predicted = as.numeric(pred_risk_scores))
  
  assign(paste0(df_name, "_delta_pred_df"), pred_df) # Dynamically assign the prediction results
  
  # Calculate R-squared value
  actual <- test_data$BMI
  predicted <- pred_risk_scores
  pred_plot <- plot(predicted, actual, xlab = paste0(df_name, " Predicted BMI"), ylab = "Actual BMI",
                    pch = 21, col = "blue", bg = "lightblue") +
    abline(lm(actual ~ predicted), col = "red", lwd = 2)
  pred_plot
  
  sse <- sum((actual - predicted)^2)
  sst <- sum((actual - mean(actual))^2)
  testing_rsq <- 1 - sse / sst
  print(paste(df_name, "R-squared:", round(testing_rsq, 3)))
  
  mean_actual <- mean(actual)
  ss_total <- sum((actual - mean_actual)^2)
  ss_residual <- sum((actual - predicted)^2)
  r_squared <- 1 - (ss_residual / ss_total)
  assign(paste0(df_name, "_r_squared"), r_squared)
  
  # Print R-squared for the current dataset
  print(paste(df_name, "R-squared:", round(r_squared, 3)))
  
  # Save predictions to a CSV file
  file_path <- file.path(out_dir, paste0(df_name, "_delta_predictions.csv"))
  write.csv(pred_df, file_path, row.names = FALSE)
}

# 1. Compute and plot the correlation distribuition
cor_matrix <- cor(all[sapply(all, is.numeric)], use = "pairwise.complete.obs")
cor_values <- cor_matrix[lower.tri(cor_matrix)]
hist(cor_values, 
     main = "Distribution of Pairwise Correlations", 
     xlab = "Correlation Coefficient", 
     breaks = 20, 
     col = "skyblue", 
     border = "white")

################################################################################
meta_delta_pred_df <- meta_delta_pred_df %>% rename(meta_predicted = predicted)
meta_no_age_sex_delta_pred_df <- meta_no_age_sex_delta_pred_df %>% rename(meta_no_age_sex_predicted = predicted)
basic_delta_pred_df <- basic_delta_pred_df %>% rename(basic_predicted = predicted)
grs_delta_pred_df <- grs_delta_pred_df %>% rename(grs_predicted = predicted)
taxa_delta_pred_df <- taxa_delta_pred_df %>% rename(taxa_predicted = predicted)
micom_delta_pred_df <- micom_delta_pred_df %>% rename(micom_predicted = predicted)
pathway_delta_pred_df <- pathway_delta_pred_df %>% rename(pathway_predicted = predicted)
metabo_delta_pred_df <- metabo_delta_pred_df %>% rename(metabo_predicted = predicted)
all_delta_pred_df <- all_delta_pred_df %>% rename(all_predicted = predicted)
all_no_age_sex_delta_pred_df <- all_no_age_sex_delta_pred_df %>% rename(all_no_age_sex_predicted = predicted)
all_no_clin_delta_pred_df <- all_no_clin_delta_pred_df %>% rename(all_no_clin_predicted = predicted)

##### Compare models
met_tax <- merge(meta_delta_pred_df, taxa_delta_pred_df, 
                 by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

met_tax_micom <- merge(met_tax, micom_delta_pred_df, 
                       by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

met_tax_micom_path <- merge(met_tax_micom, pathway_delta_pred_df, 
                            by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

met_tax_micom_path_tab <- merge(met_tax_micom_path, metabo_delta_pred_df, 
                                by = c("subject_id", "time")) %>%
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)
  
met_tax_micom_path_tab_all <- merge(met_tax_micom_path_tab, all_delta_pred_df, 
                  by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x) %>% unique()

all_but_b <- merge(met_tax_micom_path_tab_all, grs_delta_pred_df,
                   by = c("subject_id", "time")) %>% 
             dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

all_but_b_mnoas <- merge(all_but_b, meta_no_age_sex_delta_pred_df,
                         by = c("subject_id", "time")) %>% 
                dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

all_but_b_anoas <- merge(all_but_b_mnoas, all_no_age_sex_delta_pred_df,
                         by = c("subject_id", "time")) %>% 
                   dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

all_but_b_anoc <- merge(all_but_b_anoas, all_no_clin_delta_pred_df,
                         by = c("subject_id", "time")) %>% 
                    dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

all_omic <- merge(basic_delta_pred_df, all_but_b_anoc, 
                  by = c("subject_id", "time")) %>% 
            dplyr::select(-c(actual.y)) %>% rename(actual = actual.x) %>% unique()

#all_omic[, 3:8] <- scale(all_omic[, 3:8])

########################################################################################
mod_dat = all_omic %>% dplyr::rename(bmi = actual, 
                              Time = time,
                              Cluster = subject_id,
                              y_new_basic_only = basic_predicted, 
                              y_new_meta_only = meta_predicted,
                              y_new_meta_noas_only = meta_no_age_sex_predicted,
                              y_new_grs_only = grs_predicted,
                              y_new_micom_only = micom_predicted,
                              y_new_path_only = pathway_predicted,
                              y_new_tax_only = taxa_predicted,
                              y_new_metab_only = metabo_predicted,
                              y_new_all_only = all_predicted,
                              y_new_all_noas_only = all_no_age_sex_predicted,
                              y_new_all_nclin_only = all_no_clin_predicted)

write.csv(mod_dat, file = paste0(out_dir, "gl_delta_predictions_df.csv"))

mod_dat$Time <- as.numeric(mod_dat$Time)

### Single basic model plus each omic 
lmer_basic <- lmer(bmi ~ y_new_basic_only + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_meta_b <- lmer(bmi ~ y_new_basic_only + y_new_meta_only + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_grs_b <- lmer(bmi ~ y_new_basic_only + y_new_grs_only + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_micom_b <- lmer(bmi ~ y_new_basic_only + y_new_micom_only + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_path_b <- lmer(bmi ~ y_new_basic_only + y_new_path_only + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_tax_b <- lmer(bmi ~ y_new_basic_only + y_new_tax_only + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_metabo_b <- lmer(bmi ~ y_new_basic_only + y_new_metab_only + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_all_b <- lmer(bmi ~ y_new_basic_only + y_new_all_only + (1|Cluster), data = mod_dat, REML = FALSE)

anova(lmer_basic, lmer_meta_b)
anova(lmer_basic, lmer_grs_b)
anova(lmer_basic, lmer_micom_b)
anova(lmer_basic, lmer_path_b)
anova(lmer_basic, lmer_tax_b)
anova(lmer_basic, lmer_metabo_b)
anova(lmer_basic, lmer_all_b)

# https://joshuawiley.com/MonashHonoursStatistics/LMM_Comparison.html#nested-models-in-r
# Basic vs Basic + meta
nobs(lmer_basic)
nobs(lmer_meta_b)
logLik(lmer_basic)
logLik(lmer_meta_b)
anova(lmer_basic, lmer_meta_b, test = "LRT")
anova(lmer_basic, lmer_grs_b, test = "LRT")
anova(lmer_basic, lmer_micom_b, test = "LRT")
anova(lmer_basic, lmer_tax_b, test = "LRT")
anova(lmer_basic, lmer_path_b, test = "LRT")
anova(lmer_basic, lmer_metabo_b, test = "LRT")
anova(lmer_basic, lmer_all_b, test = "LRT")
glmmlass_lmer_models <- list(
  c("lmer_basic", "lmer_meta_b"),
  c("lmer_basic", "lmer_grs_b"),
  c("lmer_basic", "lmer_tax_b"),
  c("lmer_basic", "lmer_micom_b"),
  c("lmer_basic", "lmer_path_b"),
  c("lmer_basic", "lmer_metabo_b"),
  c("lmer_basic", "lmer_all_b"))

library(lme4)  # Make sure the lme4 package is loaded

anova_results <- list()  # empty list to store ANOVA results

for (model_pair in glmmlass_lmer_models) {
  model_1 <- get(model_pair[1])
  model_2 <- get(model_pair[2])
  
  # Perform ANOVA and check if both models are of class lmerMod
  if (inherits(model_1, "lmerMod") && inherits(model_2, "lmerMod")) {
    # Perform ANOVA with mixed models
    anova_result <- anova(model_1, model_2)  # This works with lmerMod objects
    tidied_result <- tidy(anova_result)  # Tidy the ANOVA result using broom::tidy()
    tidied_result$model_comparison <- paste(model_pair[1], "vs", model_pair[2])
    anova_results[[length(anova_results) + 1]] <- tidied_result
  } else {
    message("One of the models in the pair is not a valid lmerMod object.")
  }
}

# Combine all results into one dataframe
anova_table <- do.call(rbind, anova_results)
anova_table_clean <- anova_table %>%
  filter(!is.na(statistic) & !is.na(p.value)) %>%
  mutate(across(where(is.numeric), round, 3)) 

# View the combined table
print(anova_table_clean)
# Create an HTML table from the cleaned anova table
html_table <- kable(anova_table_clean, format = "html", table.attr = "class='table table-striped'")

# Save the table as an HTML file
writeLines(html_table, paste0("glmmlasso_delta_anova_table_noT.html"))

########################################################################################
# library(nlme)
# all_omic$subject_id <- as.factor(all_omic$subject_id)
# all_omic$time <- as.factor(all_omic$time)
# 
# # Make linear models 
# basic_model <- lme(actual ~ basic_predicted + time, 
#                    random = ~1|subject_id, data = all_omic)
# 
# # Combined models 
# lm_meta_tax_time <- lme(actual ~ meta_predicted + taxa_predicted + time, 
#                         random = ~1|subject_id, data = all_omic)
# summary(lm_meta_tax_time)
# 
# ### Checks before models 
# vif(lm(actual ~ meta_predicted + taxa_predicted + micom_predicted + time, data = all_omic))
# summary(all_omic)
# table(all_omic$subject_id)
# cor(all_omic[, c("meta_predicted", "taxa_predicted", "micom_predicted")])
# hist(all_omic$actual)
# anyDuplicated(all_omic)
# 
# # Fit a simpler model first
# lm_path <- lme(actual ~ pathway_predicted, random = ~1|subject_id, data = all_omic)
# lm_micom <- lme(actual ~ micom_predicted, random = ~1|subject_id, data = all_omic)
# lm_tax <- lme(actual ~ taxa_predicted, random = ~1|subject_id, data = all_omic)
# lm_meta <- lme(actual ~ meta_predicted, random = ~1|subject_id, data = all_omic)
# lm_basic <- lme(actual ~ basic_predicted, random = ~1|subject_id, data = all_omic)
# 
# model_titles <- c("1: Basic",
#                   "2: Meta", 
#                   "3: Taxa", 
#                   "4: Micom", 
#                   "5: Pathway")
# 
# # Create the table as a grid object
# sjPlot::tab_model(lm_basic, lm_meta, lm_tax, 
#                   lm_micom, lm_path, 
#                   title = "Comparing single omic glmLASSO models deltas",
#                   string.pred = "Predictors",
#                   string.est = "Estimate",
#                   string.std = "std. Beta",
#                   string.ci = "95% CI",
#                   string.se = "std. Error",
#                   p.style = c("numeric"), 
#                   p.threshold = c(0.05),
#                   dv.labels = model_titles,
#                   auto.label = FALSE)
# ### Combined 
# lm_basic <- lme(actual ~ basic_predicted, 
#                 random = ~1|subject_id, data = all_omic)
# lm_meta_basic <- lme(actual ~ basic_predicted + meta_predicted, 
#                      random = ~1|subject_id, data = all_omic)
# lm_basic_tax <- lme(actual ~ basic_predicted + taxa_predicted, 
#                     random = ~1|subject_id, data = all_omic)
# lm_basic_micom <- lme(actual ~ basic_predicted + micom_predicted, 
#                       random = ~1|subject_id, data = all_omic)
# lm_basic_path <- lme(actual ~ basic_predicted + pathway_predicted, 
#                      random = ~1|subject_id, data = all_omic)
# 
# # Define model titles for clarity
# model_titles <- c("Model 0: Basic" ,
#                   "Model 1: Basic + Meta", 
#                   "Model 2: Basic + Taxa", 
#                   "Model 3: Basic + Micom", 
#                   "Model 4: Basic + Pathway")
# 
# sjPlot::tab_model(lm_basic, lm_meta_basic, lm_basic_tax, 
#                   lm_basic_micom, lm_basic_path, 
#                   title = "Comparing combined omic delta models",
#                   string.pred = "Predictors",
#                   string.est = "Estimates",
#                   string.std = "std. Beta",
#                   string.ci = "95% CI",
#                   string.se = "std. Error",
#                   p.style = c("numeric"), 
#                   p.threshold = c(0.05),
#                   dv.labels = model_titles,
#                   auto.label = FALSE)

# Save the plot
#output_file <- file.path(out_dir, "delta_combined_omics_models_feb20.png")
#ggsave(output_file, plot_model, width = 10, height = 6)

# ### Check Correlations
# cor_matrix <- cor((all_omic)[3:7], 
#                   use = "pairwise.complete.obs", method = "pearson")
# melted_cor_matrix <- melt(cor_matrix, na.rm = TRUE)
# ggplot(melted_cor_matrix, 
#        aes(Var1, Var2, fill = value)) +
#   geom_tile() +
#   scale_fill_gradient2(low = "blue", 
#                        high = "red") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, size = 3),
#     axis.text.y = element_text(angle = 0, hjust = 1, size = 3)) + 
#   labs(title = "Correlations btwn HIGH corr pathway2", fill = "Correlation")
# 
# # Run the ANOVA tests
# anova_basic_meta <- anova(lm_basic, lm_meta_basic)
# anova_basic_tax <- anova(lm_basic, lm_basic_tax)
# anova_basic_micom <- anova(lm_basic, lm_basic_micom)
# anova_basic_path <- anova(lm_basic, lm_basic_path)
# anova_basic_all <- anova(lm_basic, lm_basic_all_omic)
# 
# meta_r <- rsq.lmm(lm_meta,adj=TRUE)
# basic_r <- rsq.lmm(lm_basic,adj=TRUE)
# meta_basic_r <- rsq.lmm(lm_meta_basic,adj=TRUE)
# tax_basic_r <- rsq.lmm(lm_basic_tax,adj=TRUE)
# micom_basic_r <- rsq.lmm(lm_basic_micom,adj=TRUE)
# path_basic_r <- rsq.lmm(lm_basic_path,adj=TRUE)
# all_r <- rsq.lmm(lm_basic_all_omic,adj=TRUE)
# 
# ### Plot r-squared 
# # Create a data frame with model names and their corresponding R-squared values
# r_squared_data <- data.frame(
#   Model = rep(c("lm_basic", "lm_basic_meta", "lm_basic_tax", 
#                 "lm_basic_micom", "lm_basic_path", "basic_all_omic"), each = 3),
#   Type = rep(c("Model_Total", "Fixed_Effects", "Random_Effects"), times = 6),
#   R_squared = c(path_basic_r$model, path_basic_r$fixed, path_basic_r$random,
#                 meta_basic_r$model, meta_basic_r$fixed, meta_basic_r$random,
#                 tax_basic_r$model, tax_basic_r$fixed, tax_basic_r$random,
#                 micom_basic_r$model, micom_basic_r$fixed, micom_basic_r$random,
#                 basic_r$model, basic_r$fixed, basic_r$random,
#                 all_r$model, all_r$fixed, all_r$random))
# 
# # Plot the data using ggplot2
# r2_plot <- ggplot(r_squared_data, aes(x = Model, y = R_squared, fill = Type)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   theme_minimal() +
#   labs(x = "Model", y = "R Squared", 
#        title = "R Squared for Different Models (glmlasso delta)") +
#   scale_fill_manual(values = c("lightblue", "lightgreen", "lightcoral")) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   ylim(0, 0.6)
# 
# output_file <- file.path(out_dir, "glm_delta_r2_models_feb20.png")
# ggsave(output_file, r2_plot, width = 10, height = 6)
# 
# ### PLOT ANOVA RESULTS
# anova_results <- data.frame(
#   Model_Comparison = c("basic vs meta + basic", "basic vs basic + tax", 
#                        "basic vs basic + micom", "basic vs basic + path",
#                        "basic vs basic + all omic"),
#   L_Ratio = c(anova_basic_meta$L.Ratio[2], anova_basic_tax$L.Ratio[2], 
#               anova_basic_micom$L.Ratio[2], anova_basic_path$L.Ratio[2],
#               anova_basic_all$L.Ratio[2]),
#   p_value = c(anova_basic_meta$`p-value`[2], anova_basic_tax$`p-value`[2], 
#               anova_basic_micom$`p-value`[2], anova_basic_path$`p-value`[2],
#               anova_basic_all$`p-value`[2]))
# 
# # Plot ANOVA
# anova_plot <- ggplot(anova_results, aes(x = Model_Comparison, y = L_Ratio)) +
#   geom_bar(stat = "identity", fill = "skyblue") +  # Create bars
#   geom_text(aes(label = round(p_value, 3)), vjust = -0.5) +  # Add p_value above the bars
#   theme_minimal() +  # Clean theme
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
#   labs(x = "Model Comparison", y = "L Ratio", 
#        title = "L Ratio for Each Model Comparison (glmlasso delta models)") +
#   ylim(0, 7)
# 
# output_file <- file.path(out_dir, "glm_delta_anova_models_feb20.png")
# ggsave(output_file, anova_plot, width = 10, height = 6)
