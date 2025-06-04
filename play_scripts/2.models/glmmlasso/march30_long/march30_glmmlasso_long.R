# Trying out long lasso
rm(list = ls())
#source("zc_functions.R") 
library(pacman)
p_load(tools, reticulate, viridis, tidyplots, patchwork, jsonlite, maps, ggvenn, 
       caret, caretEnsemble, glmnet, xgboost, ggplot2, glmmLasso, corrplot,
       readr, plyr, dplyr, tidyr, purrr, tibble, stringr, psych, randomForest,  
       reshape2, scales, gridExtra, plotly, sf, tidyverse, naniar, VIM, kableExtra)
#remotes::install_github("thepira/cv.glmmLasso")
library(cv.glmmLasso)
library(kneedle)
'%ni%' <- Negate('%in%')

#out_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/feb20_long/"
out_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_long/new_split_may/"
out_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/long/"

#long_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/merf_dfs/5.combined/"
#long_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/march_20/"
long_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/april_processing/"

long <- read.csv(file.path(long_dir, 'long_april29.csv')) %>% 
                 dplyr::select(-c("consent", "record_id", "completer", 
                                  "Peptide_YY", "Ghrelin", "Leptin")) %>%
  dplyr::mutate(time = as.factor(time),
    subject_id = as.factor(subject_id),
    randomized_group = as.factor(randomized_group),
    sex = as.numeric(sex),
    randomized_group = as.numeric(randomized_group),
    cohort_number = as.numeric(cohort_number),
    race = as.numeric(race)) %>% rename(BMI = outcome_BMI_fnl,
                                    range = time,
                                    homo_ir = HOMA_IR,
                                    insulin = Insulin_endo,
                                    LDL = LDL_Calculated,
                                    HDL = HDL_Total_Direct_lipid,
                                    HbA1c = Hemoglobin_A1C)

# Create the plot with lines between dots of the same subject_id
ggplot(long, aes(x = range, y = bmi_prs, color = factor(subject_id), group = factor(subject_id))) + 
  geom_point() +  # Scatter plot
  geom_line() +   # Add lines between points of the same subject_id
  labs(title = "BMI PRS by Range and Subject ID", 
       x = "Time", 
       y = "BMI PRS") + 
  theme(legend.position = "none") +  # Remove the legend
  scale_color_manual(values = rainbow(length(unique(long$subject_id)))) 

# Define the column names based on your lists
basic <- c('subject_id','BMI', 'range','age', 'sex', 'randomized_group')
meta_keep <- c('subject_id','BMI', 'range', 'randomized_group', 'sex', 'race', 
               'age', 'HbA1c', 'HDL', 'homo_ir', 'insulin', 'LDL', 'Glucose.x')
only_grs <- c('subject_id', 'BMI', 'range',  'bmi_prs')
only_taxa <- c('subject_id','BMI', 'range', grep("^g__", names(long), value = TRUE))

micom_start <- which(names(long) == "Diacetyl")
micom_end <- which(names(long) == "aldehydo.D.xylose")
only_micom <- c('subject_id','BMI', 'range', names(long)[micom_start:micom_end])

path_start <- which(names(long) == "arginine..ornithine.and.proline.interconversion")
path_end <- which(names(long) == "UDP.N.acetyl.D.glucosamine.biosynthesis.I")
only_pathway <- c('subject_id','BMI', 'range', names(long)[path_start:path_end])

tabo_start <- which(names(long) == "non_HDL_C")
tabo_end <- which(names(long) == "IDL_TG_pct")
only_tabo <- c('subject_id','BMI', 'range', names(long)[tabo_start:tabo_end])

all_col <- c('subject_id','BMI', 'range',
             'randomized_group', 'sex', 'race', 
             'age', 'HbA1c', 'HDL', 'homo_ir', 'insulin', 'LDL', 'Glucose.x',
             grep("^g__", names(long), value = TRUE),
             names(long)[micom_start:micom_end],
             names(long)[path_start:path_end],
             names(long)[tabo_start:tabo_end])

# Create data frames based on the columns defined
basic <- long[, basic, drop = FALSE] %>% unique()
meta <- long[, meta_keep, drop = FALSE] %>% unique()
grs <- long[, only_grs, drop = FALSE] %>% unique()
taxa <- long[, only_taxa, drop = FALSE] %>% unique()
micom <- long[, only_micom, drop = FALSE] %>% unique()
pathway <- long[, only_pathway, drop = FALSE] %>% unique()
metabo <- long[, only_tabo, drop = FALSE] %>% unique()
all <- long[, all_col, drop = FALSE] %>% unique()

## Check MICOM correlatioon 
heatmap(cor(micom[, 4:ncol(micom)]))
heatmap(cor(all[, 4:ncol(all)]))

## Check all correlation
preProcValues_train <- preProcess(all[, c(2, 4:ncol(all))], 
                                  method = c("nzv", "corr"), thresh = 0.95, fudge = 0.2, 
                                  numUnique = 15, verbose = TRUE, freqCut = 95/5, 
                                  uniqueCut = 10, cutoff = 0.75, na.remove = TRUE)
preProcValues_train
all <- predict(preProcValues_train, all)

subject_id_count <- meta %>%
  dplyr::filter(range %in% c("0","6", "12")) %>%
  dplyr::group_by(subject_id) %>%
  dplyr::summarize(range_count = n_distinct(range))  # Count the distinct range values
table(subject_id_count$range_count)

# Filter `subject_id`s with fewer than 3 unique range values (0, 6, 12)
missing_subjects  <- subject_id_count %>% dplyr::filter(range_count != 3)
missing_subjects

# Make train and test set
# Test sample names
#test_names <- c("ABR-079", "AGA-071", "AHE-055", "ALI-121", "ALO-163", "AMA-031", "ASO-013", "AWI-167", "BMO-164", "CWA-183", "DSC-024", "EBE-130", "EHI-177", "EJO-092", "GFU-188", "HGI-010", "JCA-109", "JGO-100", "KBU-085", "KCE-034", "KHE-170", "LDO-148", "LST-186", "LZD-142", "MAR-119", "MCA-088", "MJA-153", "MWE-112", "NPO-149", "RAE-114", "SBO-020", "SEG-080", "SKA-195", "SLO-178", "SSH-028", "TDU-086","TFA-016", "VCA-041")

test_names <- c("ASO-013", "NTA-021", "KGI-029", "KPA-042", "AWA-052", "AHE-055", "COW-066", "NBI-069", "CEL-073", "CAL-074", "ABR-079", "SEG-080", "NKA-090", "NEL-094", "LJA-101", "ADA-105", "MLU-106", "MDI-107", "JER-110", "TRO-113", "MFB-118", "ALI-121", "KWA-122", "RAF-125", "EBE-130", "CGA-134", "LZD-142", "NPO-149", "HDE-154", "AMC-155", "SAB-160", "QNG-166", "NCO-171", "BSA-174", "EHI-177", "LST-186", "MBA-187", "BAN-193")
  
# Train sample names
#rain_names <- c("AAL-144", "ACO-053", "ADA-105", "AKE-009", "AKI-011", "AKO-139", "AMC-155", "AME-128", "AME-157", "ATA-129", "AWA-052", "AWA-083", "BAN-193", "BHO-014", "BIN-201", "BKN-104", "BMI-156", "BSA-174", "CAM-057", "CCO-189", "CED-026", "CEL-073", "CGA-134", "CIS-077", "CKR-078", "CLE-049", "COW-066", "CRO-108", "CWA-161", "EBE-051", "EKA-135", "EKR-045", "ELA-159", "EPO-182", "EVO-184", "FWI-098", "GHA-035", "HDE-154", "IBE-120", "JDI-140", "JER-110", "JFU-027", "JJO-093", "JKN-127", "JPO-022", "JUG-116", "JUT-032", "JVE-126", "KAN-138", "KBR-162", "KEL-185", "KEL-199", "KGI-029", "KHU-196", "KPA-042", "KRI-072", "KVA-038", "KWA-122", "KWA-141", "LBL-047", "LBU-015", "LEL-147", "LFI-003", "LJA-101", "LMC-111", "LPF-198", "LVA-017", "MBA-187", "MCW-065", "MDI-107", "MES-068", "MFB-118", "MGA-076", "MHO-117", "MKE-192", "MMA-036", "MRT-179", "MSH-091", "MST-039", "MWE-143", "MWO-133", "MWY-152", "NAR-099", "NBI-048", "NBI-069", "NCO-171", "NDI-067", "NEL-094", "NKA-090", "NMO-151", "NTA-021", "PBE-123", "QNG-166", "RAF-125", "RAM-050", "RHP-023", "RLA-132", "ROL-006", "SAB-160", "SCA-043", "SCR-061", "SDA-150", "SGA-062", "SKA-087", "SRO-194", "TBU-115", "TFA-172", "TRO-113", "TSH-146", "TSL-056", "WPE-005", "YOR-103", "YSU-097", "ZVU-096")

train_names <- c("SDA-150", "LBU-015", "CIS-077", "ATA-129", "KHU-196", "MWY-152", "AGA-071", "AME-157", "CWA-183", "RHP-023", "MST-025", "SSH-028", "JUG-116", "EJO-092", "VCA-041", "NMO-151", "BHO-014", "KBU-085", "SBO-020", "MWO-133", "KRI-072", "AAL-144", "ALO-163", "AKI-011", "MHO-117", "TSH-146", "RAE-114", "FWI-098", "MAR-119", "JGO-100", "CAM-057", "YOR-103", "HGI-010", "KAN-138", "SGA-062", "CKR-078", "MWE-112", "ROL-006", "MMA-036", "DSC-024", "LDO-148", "MCA-088", "CPU-075", "AKO-139", "LFI-003", "KWA-141", "GFU-188", "BMO-164", "JPO-022", "EVO-184", "LPF-198", "TBU-115", "SRO-194", "KEL-199", "JFU-027", "SKA-195", "IBE-120", "TSL-056", "NDI-067", "AWA-083", "CWA-161", "TDU-086", "JCA-109", "CBO-004", "NAR-099", "MES-068", "AMA-031", "SLO-178", "SCA-043", "AWI-167",  "KBR-162", "TFA-172", "BIN-201", "NBI-048", "KHE-170", "CSH-012", "BMI-156", "MWE-143", "EKA-135", "WPE-005", "AKE-009", "YSU-097", "MCW-065", "EBE-051", "ZVU-096", "JJO-093", "KVA-038", "ACO-053", "RLA-132", "MBA-176", "CED-026", "JDI-140", "CCO-189", "EKR-045", "MJA-153", "CLE-049", "LMC-111", "SKA-087", "JUT-032", "MKE-192", "JVE-126", "KCE-034", "KEL-185", "MRT-179", "JKN-127", "LEL-147", "BKN-104", "AME-128", "MSH-091", "MGA-076", "LVA-017", "EPO-182")

cat("Length of test names:", length(test_names), "\n")
cat("Length of train names:", length(train_names), "\n")

# Make test and tain sets for each omic 
data_frames <- c("basic", "meta", "grs","micom", "pathway", "taxa", "metabo", "all")
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

# CHECK CORRELATIONS
#preProcValues_train <- preProcess(all_train[, c(2, 4:ncol(all_train))], 
#                            method = c("nzv", "corr"), thresh = 0.95, fudge = 0.2, 
#                            numUnique = 15, verbose = TRUE, freqCut = 95/5, 
#                            uniqueCut = 10, cutoff = 0.75, na.remove = TRUE)
#preProcValues_train
#all_train <- predict(preProcValues_train, all_train)
#heatmap(cor(all_train[, c(2, 4:ncol(all_train))]))

#preProcValues_test <- preProcess(all_test[, c(2, 4:ncol(all_test))], 
#                            method = c("nzv", "corr"), thresh = 0.95, fudge = 0.2, 
#                            numUnique = 15, verbose = TRUE, freqCut = 95/5, 
#                            uniqueCut = 10, cutoff = 0.75, na.remove = TRUE)
#preProcValues_test
#all_test <- predict(preProcValues_test, all_test)
#### Test CV LASSO

# Step 1: run glmmlasso through a grid of lambdas for the best one
# Step 2: re-run glmmlasso using the best lambdas
# Step 4: Use that model in step 3 to predict BMI in the test set
# Step 5: That prediction in step 4 becomes the "risk score" for that omic

# STEP 1
#data_frames <- c("grs")
#df_name <- "metabo"
data_frames <- c("basic", "meta", "grs", "taxa", "pathway", "micom", "metabo", "all")
for (df_name in data_frames) {
  train_data <- get(paste0(df_name, "_train"))
  test_data <- get(paste0(df_name, "_test"))
  numvar <- c() 
  lambdavec <- seq(from = 10, to = 100, by = 1)
  #lambdy <- 5
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
  lambda_value <- kneedle(lambdavec, numvar)
  lambda_value <- lambda_value[1]
  
  # STEP 2 Create the final model on the training data using the chosen lambda
  best_model <- glmmLasso(fix = fix_formula,
                          data = train_data,
                          rnd = list(subject_id = ~ 1, range = ~ 1),
                          lambda = lambda_value,  # Use the lambda_value set manually
                          family = gaussian(link = "identity"))
  
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
  
  write.csv(top_features, 
            file = paste0(out_dir, 
                          paste0(df_name, "_gl_long_top_features.csv")), 
            row.names = FALSE)
  
  plot_feat <- top_features %>%
    dplyr::filter(!str_detect(tolower(Feature), "time")) %>%  # Remove rows with "time" in Feature
    slice_head(n = 10)  
  
  # Plot top features
  features <- ggplot(plot_feat, 
                     aes(x = reorder(Feature, Estimate), 
                         y = Estimate)) +
    geom_bar(stat = "identity", 
             fill = "#1C4C98") +
    coord_flip() + 
    theme_bw() +
    ggtitle(paste("Top Features & Coefficients from", df_name)) +
    xlab("Feature") + 
    ylab("Coefficient Estimate") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12))
  
  features 
  
  # Save the plot using ggsave
  ggsave(filename = paste0(out_dir, paste0(df_name, "_top_features.png")),
         plot = features, width = 8, height = 6, units = "in", dpi = 300)
  
  # STEP Predict on test data 
  test_data <- test_data[complete.cases(test_data),]  # Filter complete cases
  pred_risk_scores <- predict(best_model, test_data)
  range_long <- test_data$range
  sub_ids <- test_data$subject_id
  actual_bmi <- test_data$BMI
  pred_df <- data.frame(
    subject_id = as.character(sub_ids),
    time = as.numeric(range_long), 
    actual = as.numeric(actual_bmi),
    predicted = as.numeric(pred_risk_scores))
  
  # Combine predicted and actual values
  #pred_df <- as.data.frame(cbind(test_data$subject_id, test_data$range, test_data$BMI, pred_risk_scores))
  #colnames(pred_df) <- c("subject_id", "time", "actual", "predicted")
  assign(paste0(df_name, "_pred_df"), pred_df) # Dynamically assign the prediction results
  
  # Calculate R-squared value
  actual <- test_data$BMI
  predicted <- pred_risk_scores
  pred_plot <- plot(predicted, 
                    actual, 
                    xlab = paste0(df_name, 
                                  " Predicted BMI"), ylab = "Actual BMI",
       pch = 21, col = "blue", bg = "lightblue") +
  abline(lm(actual ~ predicted), col = "red", lwd = 2)
  
  pred_plot
  
  sse <- sum((actual - predicted)^2)
  sst <- sum((actual - mean(actual))^2)
  testing_rsq <- 1 - sse / sst

  mean_actual <- mean(actual)
  ss_total <- sum((actual - mean_actual)^2)
  ss_residual <- sum((actual - predicted)^2)
  r_squared <- 1 - (ss_residual / ss_total)
  assign(paste0(df_name, "_r_squared"), r_squared)

  # Print R-squared for the current dataset
  print(paste(df_name, "R-squared:", round(r_squared, 3)))
  
  # Save predictions to a CSV file
  file_path <- file.path(out_dir, paste0(df_name, "_predictions.csv"))
  write.csv(pred_df, file_path, row.names = FALSE)
}

# Participants with very odd predictions:
#AHE-055, 	TFA-016

################################################################################

##### Compare models
# Convert the first 2 columns to factors
basic_pred_df[, head(names(basic_pred_df), 2)] <- 
  lapply(basic_pred_df[, head(names(basic_pred_df), 2)], as.factor) 
basic_pred_df <- basic_pred_df %>% rename(y_new_basic = predicted)

meta_pred_df[, head(names(meta_pred_df), 2)] <- 
  lapply(meta_pred_df[, head(names(meta_pred_df), 2)], as.factor) 
meta_pred_df <- meta_pred_df%>% dplyr::rename(y_new_meta_only = predicted)

taxa_pred_df[, head(names(taxa_pred_df), 2)] <- 
  lapply(taxa_pred_df[, head(names(taxa_pred_df), 2)], as.factor) 
taxa_pred_df <- taxa_pred_df%>% dplyr::rename(y_new_tax_only = predicted)

micom_pred_df[, head(names(micom_pred_df), 2)] <- 
  lapply(micom_pred_df[, head(names(micom_pred_df), 2)], as.factor) 
micom_pred_df <- micom_pred_df %>% dplyr::rename(y_new_micom_only = predicted)

pathway_pred_df[, head(names(pathway_pred_df), 2)] <- 
  lapply(pathway_pred_df[, head(names(pathway_pred_df), 2)], as.factor) 
pathway_pred_df <- pathway_pred_df %>% dplyr::rename(y_new_path_only = predicted)

metabo_pred_df[, head(names(metabo_pred_df), 2)] <- 
  lapply(metabo_pred_df[, head(names(metabo_pred_df), 2)], as.factor) 
metabo_pred_df <- metabo_pred_df %>% dplyr::rename(y_new_metabo_only = predicted)

grs_pred_df[, head(names(grs_pred_df), 2)] <- 
  lapply(grs_pred_df[, head(names(grs_pred_df), 2)], as.factor) 
grs_pred_df <- grs_pred_df %>% dplyr::rename(y_new_grs_only = predicted)

all_pred_df[, head(names(all_pred_df), 2)] <- 
  lapply(all_pred_df[, head(names(all_pred_df), 2)], as.factor) 
all_pred_df <- all_pred_df %>% dplyr::rename(y_new_all_only = predicted)

met_basic <- merge(basic_pred_df, meta_pred_df, by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

met_tax <- merge(met_basic, taxa_pred_df, 
                 by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

met_tax_micom <- merge(met_tax, micom_pred_df, 
                       by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

met_tax_micom_path <- merge(met_tax_micom, pathway_pred_df, 
                            by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

met_tax_micom_path_tabo <- merge(met_tax_micom_path, metabo_pred_df, 
                            by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

all_but_all <- merge(met_tax_micom_path_tabo, grs_pred_df, 
                  by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

all_omic <- merge(all_but_all, all_pred_df, 
                  by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x) %>% 
  unique()

# Center and scale the last 5 columns of the dataframe
#all_omic[, 3:8] <- scale(all_omic[, 3:8])
head(all_omic)

########################################################################################
mod_dat = all_omic %>% rename(bmi = actual, Time = time, Cluster = subject_id)
write.csv(mod_dat, file = paste0(out_dir, "june_lasso_long_predictions_df.csv"))

### Single plus omic including time 
lmer_basic <- lmer(bmi ~ y_new_basic + Time + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_meta_b <- lmer(bmi ~ y_new_basic + y_new_meta_only + Time + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_grs_b <- lmer(bmi ~ y_new_basic + y_new_grs_only + Time + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_micom_b <- lmer(bmi ~ y_new_basic + y_new_micom_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_path_b <- lmer(bmi ~ y_new_basic + y_new_path_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_tax_b <- lmer(bmi ~ y_new_basic + y_new_tax_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_metabo_b <- lmer(bmi ~ y_new_basic + y_new_metabo_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_all_b <- lmer(bmi ~ y_new_basic + y_new_all_only + (1|Cluster), data = mod_dat, REML = FALSE)

anova(lmer_basic, lmer_meta_b, test = "LRT")
anova(lmer_basic, lmer_grs_b, test = "LRT")
anova(lmer_basic, lmer_micom_b, test = "LRT")
anova(lmer_basic, lmer_tax_b, test = "LRT")
anova(lmer_basic, lmer_path_b, test = "LRT")
anova(lmer_basic, lmer_metabo_b, test = "LRT")
anova(lmer_basic, lmer_all_b, test = "LRT")

glmmlong_models <- list(
  c("lmer_basic", "lmer_meta_b"),
  c("lmer_basic", "lmer_grs_b"),
  c("lmer_basic", "lmer_micom_b"),
  c("lmer_basic", "lmer_tax_b"),
  c("lmer_basic", "lmer_path_b"),
  c("lmer_basic", "lmer_metabo_b"),
  c("lmer_basic", "lmer_all_b"))

library(lme4)  # Make sure the lme4 package is loaded=
anova_results <- list()
for (model_pair in glmmlong_models) {
  model_1 <- get(model_pair[1])
  model_2 <- get(model_pair[2])
  # Perform ANOVA and check if both models are of class lmerMod
  if (inherits(model_1, "lmerMod") && inherits(model_2, "lmerMod")) {
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
writeLines(html_table, paste0(out_dir, "lasso_long_anova_table_time.html"))

# Repeat without time 
### Single plus omic not including time 
lmer_basic_noT <- lmer(bmi ~ y_new_basic + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_meta_b_noT <- lmer(bmi ~ y_new_basic + y_new_meta_only + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_grs_b_noT <- lmer(bmi ~ y_new_basic + y_new_grs_only + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_micom_b_noT <- lmer(bmi ~ y_new_basic + y_new_micom_only + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_path_b_noT <- lmer(bmi ~ y_new_basic + y_new_path_only + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_tax_b_noT <- lmer(bmi ~ y_new_basic + y_new_tax_only + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_metabo_b_noT <- lmer(bmi ~ y_new_basic + y_new_metabo_only + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_all_b_noT <- lmer(bmi ~ y_new_basic + y_new_all_only + (1|Cluster), data = mod_dat, REML = FALSE)

anova(lmer_basic_noT, lmer_meta_b_noT, test = "LRT")
anova(lmer_basic_noT, lmer_grs_b_noT, test = "LRT")
anova(lmer_basic_noT, lmer_micom_b_noT, test = "LRT")
anova(lmer_basic_noT, lmer_tax_b_noT, test = "LRT")
anova(lmer_basic_noT, lmer_path_b_noT, test = "LRT")
anova(lmer_basic_noT, lmer_metabo_b_noT, test = "LRT")
anova(lmer_basic_noT, lmer_all_b_noT, test = "LRT")

glmmlong_models_noT <- list(
  c("lmer_basic_noT", "lmer_meta_b_noT"),
  c("lmer_basic_noT", "lmer_grs_b_noT"),
  c("lmer_basic_noT", "lmer_micom_b_noT"),
  c("lmer_basic_noT", "lmer_tax_b_noT"),
  c("lmer_basic_noT", "lmer_path_b_noT"),
  c("lmer_basic_noT", "lmer_metabo_b_noT"),
  c("lmer_basic_noT", "lmer_all_b_noT"))
library(lme4)  # Make sure the lme4 package is loaded=
anova_results_noT <- list()
for (model_pair in glmmlong_models_noT) {
  model_1 <- get(model_pair[1])
  model_2 <- get(model_pair[2])
  # Perform ANOVA and check if both models are of class lmerMod
  if (inherits(model_1, "lmerMod") && inherits(model_2, "lmerMod")) {
    anova_result_noT <- anova(model_1, model_2)  # This works with lmerMod objects
    tidied_result <- tidy(anova_result_noT)  # Tidy the ANOVA result using broom::tidy()
    tidied_result$model_comparison <- paste(model_pair[1], "vs", model_pair[2])
    anova_results_noT[[length(anova_results_noT) + 1]] <- tidied_result
  } else {
    message("One of the models in the pair is not a valid lmerMod object.")
  }
}

# Combine all results into one dataframe
anova_table <- do.call(rbind, anova_results_noT)
anova_table_clean <- anova_table %>%
  filter(!is.na(statistic) & !is.na(p.value)) %>%
  mutate(across(where(is.numeric), round, 3)) 

# View the combined table
print(anova_table_clean)
html_table <- kable(anova_table_clean, 
                    format = "html", 
                    table.attr = "class='table table-striped'")

# Save the table as an HTML file
writeLines(html_table, paste0(out_dir, "glmmlasso_long_anova_table_noT.html"))
