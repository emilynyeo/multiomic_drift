# Trying out long lasso
rm(list = ls())
#source("zc_functions.R") 
library(pacman)
p_load(tools, reticulate, viridis, tidyplots, patchwork, jsonlite, maps, ggvenn, 
       caret, caretEnsemble, glmnet, xgboost, ggplot2, glmmLasso, corrplot,
       readr, plyr, dplyr, tidyr, purrr, tibble, stringr, psych, randomForest,  
       reshape2, scales, gridExtra, plotly, sf, tidyverse, naniar, VIM)
'%ni%' <- Negate('%in%')

#out_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/feb20_long/"
out_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_long/"

#long_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/merf_dfs/5.combined/"
long_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/march_20/"

long <- read.csv(file.path(long_dir, 'long_df_imputed_cent_scale.csv')) %>% 
                 dplyr::select(-c("consent", "record_id", "completer", "Peptide_YY",
                                  "Ghrelin", "Leptin", "time.x")) %>%
  dplyr::mutate(time = as.factor(time),
    subject_id = as.factor(subject_id),
    randomized_group = as.numeric(randomized_group),
    sex = as.numeric(sex),
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
       x = "Range", y = "BMI PRS") + 
  theme(legend.position = "none") +  # Remove the legend
  scale_color_manual(values = rainbow(length(unique(long$subject_id)))) 

# Define the column names based on your lists
basic <- c('subject_id','BMI', 'range','age', 'sex')
meta_keep <- c('subject_id','BMI', 'range', 'randomized_group', 'sex', 'race', 
               'age', 'HbA1c', 'HDL', 'homo_ir', 'insulin', 'LDL', 'Glucose')
only_taxa <- c('subject_id','BMI', 'range', 
               grep("^g__", names(long), value = TRUE))
proton_column <- which(names(long) == "proton")
carbon_dioxide_column <- which(names(long) == "X3.methyl.2.oxopentanoate")
only_micom <- c('subject_id','BMI', 'range', 
                names(long)[proton_column:carbon_dioxide_column])
exclude_columns <- unique(c(meta_keep, only_taxa, only_micom))
only_pathway <- c('subject_id', 'BMI', 'range', 
                  setdiff(names(long), exclude_columns))

# Create data frames based on the columns defined
basic <- long[, basic, drop = FALSE] %>% unique()
meta <- long[, meta_keep, drop = FALSE] %>% unique()
taxa <- long[, only_taxa, drop = FALSE] %>% unique()
micom <- long[, only_micom, drop = FALSE] %>% unique()
pathway <- long[, only_pathway, drop = FALSE] %>% unique()

subject_id_count <- meta %>%
  dplyr::filter(range %in% c("BL","6m", "12m")) %>%
  dplyr::group_by(subject_id) %>%
  dplyr::summarize(range_count = n_distinct(range))  # Count the distinct range values
table(subject_id_count$range_count)

# Filter `subject_id`s with fewer than 3 unique range values (0, 6, 12)
missing_subjects  <- subject_id_count %>% dplyr::filter(range_count != 3)
missing_subjects

# Make train and test set
# Test sample names
test_names <- c("ABR-079", "AGA-071", "AHE-055", "ALI-121", "ALO-163", "AMA-031", "ASO-013", "AWI-167", "BMO-164", "CWA-183", "DSC-024", "EBE-130", "EHI-177", "EJO-092", "GFU-188", "HGI-010", "JCA-109", "JGO-100", "KBU-085", "KCE-034", "KHE-170", "LDO-148", "LST-186", "LZD-142", "MAR-119", "MCA-088", "MJA-153", "MWE-112", "NPO-149", "RAE-114", "SBO-020", "SEG-080", "SKA-195", "SLO-178", "SSH-028", "TDU-086","TFA-016", "VCA-041")

# Train sample names
train_names <- c("AAL-144", "ACO-053", "ADA-105", "AKE-009", "AKI-011", "AKO-139", "AMC-155", "AME-128", "AME-157", "ATA-129", "AWA-052", "AWA-083", "BAN-193", "BHO-014", "BIN-201", "BKN-104", "BMI-156", "BSA-174", "CAM-057", "CCO-189", "CED-026", "CEL-073", "CGA-134", "CIS-077", "CKR-078", "CLE-049", "COW-066", "CRO-108", "CWA-161", "EBE-051", "EKA-135", "EKR-045", "ELA-159", "EPO-182", "EVO-184", "FWI-098", "GHA-035", "HDE-154", "IBE-120", "JDI-140", "JER-110", "JFU-027", "JJO-093", "JKN-127", "JPO-022", "JUG-116", "JUT-032", "JVE-126", "KAN-138", "KBR-162", "KEL-185", "KEL-199", "KGI-029", "KHU-196", "KPA-042", "KRI-072", "KVA-038", "KWA-122", "KWA-141", "LBL-047", "LBU-015", "LEL-147", "LFI-003", "LJA-101", "LMC-111", "LPF-198", "LVA-017", "MBA-187", "MCW-065", "MDI-107", "MES-068", "MFB-118", "MGA-076", "MHO-117", "MKE-192", "MMA-036", "MRT-179", "MSH-091", "MST-039", "MWE-143", "MWO-133", "MWY-152", "NAR-099", "NBI-048", "NBI-069", "NCO-171", "NDI-067", "NEL-094", "NKA-090", "NMO-151", "NTA-021", "PBE-123", "QNG-166", "RAF-125", "RAM-050", "RHP-023", "RLA-132", "ROL-006", "SAB-160", "SCA-043", "SCR-061", "SDA-150", "SGA-062", "SKA-087", "SRO-194", "TBU-115", "TFA-172", "TRO-113", "TSH-146", "TSL-056", "WPE-005", "YOR-103", "YSU-097", "ZVU-096")

cat("Length of test names:", length(test_names), "\n")
cat("Length of train names:", length(train_names), "\n")

# List of data frames
data_frames <- c("basic", "meta", "micom", "pathway", "taxa")
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

#### ATTEMPT TO LOOP ALL GLMMLASSO ############################################
#"basic", "meta",
#, "", "taxa", micom, taxa
data_frames <- c("pathway")
for (df_name in data_frames) {
  train_data <- get(paste0(df_name, "_train"))
  test_data <- get(paste0(df_name, "_test"))
  numvar <- c() 
  lambdavec <- seq(from = 10, to = 50, by = 1)
  
  # Loop through each lambda value to perform Lasso regression
  for (lambdy in lambdavec) {
    predictors <- setdiff(names(train_data), c("BMI", "subject_id"))   
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
    lassoFeatures <- lassoFeatures[grep("as.factor", lassoFeatures, invert = TRUE)]  # Exclude factor variables
    lassoFeatures <- unique(c(lassoFeatures))
    numvar <- c(numvar, length(lassoFeatures)) # Store the number of variables included for this lambda
  }
  
  # Plot lambda vs number of variables for each data frame
  ggplot() +
    geom_point(aes(x = lambdavec, y = numvar)) +
    geom_vline(xintercept = 25, color = "blue", linetype = "dashed") +
    theme_bw() +
    ggtitle(paste("Lasso Results for", df_name)) +
    xlab("Penalty Coefficient (Lambda)") +
    ylab("Number of Included Variables")
  
  # Save the plot to a file with the corresponding dataset name
  ggsave(paste0(out_dir, paste0(df_name, "_lambda_elbow.png")), width = 6, height = 4, units = "in", dpi = 320)
  
  # Create the final model on the training data using lambda = 30
  best_model <- glmmLasso(fix = fix_formula,
                          data = train_data,
                          rnd = list(subject_id = ~ 1, range = ~ 1),
                          lambda = 30,
                          family = gaussian(link = "identity"))
  
  # Extract lasso features
  lassoFeatures <- names(best_model$coefficients)
  lassoFeatures <- lassoFeatures[lassoFeatures %ni% c("(Intercept)")]
  lassoFeatures <- lassoFeatures[grep("as.factor", lassoFeatures, invert = TRUE)]  # Exclude factor variables
  lassoFeatures <- unique(c(lassoFeatures))
  lassoFeatures <- lassoFeatures[grep("subject_id|BMI|range", lassoFeatures, invert = TRUE)]  # Exclude specific features
  lassoFeatures <- unique(c(lassoFeatures, "subject_id", "BMI", "range"))
  assign(paste0(df_name, "_lassoFeatures"), lassoFeatures)
  
  # Prepare the data frame for the final model
  bestglm <- as.data.frame(train_data[, lassoFeatures])
  varlist <- names(bestglm)[which(names(bestglm) %ni% c("BMI", "subject_id", "range"))]
  varstring <- paste0(varlist, collapse = " + ", sep = "")
  
  # Fit the final model
  mymod <- lme4::glmer(as.formula(paste0("BMI ~ ", varstring, " + (1|range) + (1|subject_id)")),
                       data = bestglm,
                       family = gaussian(link = "identity"))
  assign(paste0(df_name, "_mymod"), mymod)
  
  # Extract coefficients and plot the top features
  coef_df <- as.data.frame(summary(mymod)$coefficients)
  coef_df$Feature <- rownames(coef_df)
  coef_df <- coef_df %>% arrange(desc(Estimate))  # Sort by coefficients
  assign(paste0(df_name, "_coef_df"), coef_df)
  
  # Filter top 10 features by absolute coefficient magnitude
  top_features <- coef_df %>%
    mutate(abs_estimate = abs(Estimate)) %>%
    dplyr::filter(Feature != "(Intercept)") %>%
    mutate(Feature = str_to_title(Feature),  # Capitalize the first letter of each word
           Feature = str_replace_all(Feature, "[._]", " ")) %>%
    arrange(desc(abs_estimate)) %>%
    head(10)
  
  # Plot top features
  ggplot(top_features, aes(x = reorder(Feature, Estimate), y = Estimate)) +
    geom_bar(stat = "identity", fill = "#1C4C98") +
    coord_flip() + theme_bw() +
    ggtitle(paste("Top Features & Coefficients from", df_name)) +
    xlab("Feature") + ylab("Coefficient Estimate") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12))
  
  # Save the top features plot for the current dataset
  ggsave(paste0(out_dir, paste0(df_name, "_top_features.png")), width = 8, height = 6, units = "in", dpi = 300)
  
  # Predict on test data and calculate R-squared
  test_data <- test_data[complete.cases(test_data),]  # Filter complete cases
  
  # Prepare test data with the variables selected in the model
  test_vars <- test_data %>% dplyr::select(rownames(coef_df)[2:nrow(coef_df)])
  test_vars$Intercept <- 1  # Add intercept column
  test_vars <- test_vars %>% dplyr::select(Intercept, everything())  # Reorder columns
  
  # Predict risk scores
  pred_risk_scores <- as.matrix(test_vars) %*% as.matrix(coef_df$Estimate)
  
  # Combine predicted and actual values
  pred_df <- as.data.frame(cbind(test_data$subject_id, test_data$range, test_data$BMI, scale(pred_risk_scores)))
  colnames(pred_df) <- c("subject_id", "time", "actual", "predicted")
  assign(paste0(df_name, "_pred_df"), pred_df) # Dynamically assign the prediction results
  
  # Calculate R-squared value
  actual <- test_data$BMI
  predicted <- pred_risk_scores
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


# Old way

################################################################################
### BASIC data set
numvar_basic <- c()
lambdavec <- seq(from = 10, to = 50, by = 1)
for (lambdy in lambdavec) {
  predictors <- setdiff(names(basic_train), c("BMI", "subject_id"))  
  predictors_escaped <- paste0("`", predictors, "`", collapse = " + ")
  fix_formula <- as.formula(paste("BMI ~", predictors_escaped))
  lm1_basic <- glmmLasso(fix = fix_formula,
                         data = basic_train,
                         rnd = list(subject_id = ~ 1, range = ~ 1),
                         lambda = lambdy,
                         family = gaussian(link = "identity"))
  summary(lm1_basic)
  lassoFeatures <- names(lm1_basic$coefficients[which(lm1_basic$coefficients != 0)])
  lassoFeatures <- lassoFeatures[lassoFeatures %ni% c("(Intercept)")]
  lassoFeatures <- lassoFeatures[grep("as.factor",lassoFeatures,invert=T)] ####
  lassoFeatures <- unique(c(lassoFeatures))
  numvar_basic <- c(numvar_basic, length(lassoFeatures))
}
plot(x = lambdavec, y = numvar_basic)
ggplot() +
  geom_point(aes(x = lambdavec, y = numvar_basic)) +
  geom_vline(xintercept = 25, color = "blue", linetype = "dashed") +
  theme_bw() +
  ggtitle("Only BAsic Delta") +
  xlab("Penalty Coefficient (Lambda)") +
  ylab("Number of Included Variables")
ggsave(paste0(out_dir,"lambda_elbow_only_basic.png"), width=6, height=4, units="in", dpi=320)

best_basic <- glmmLasso(fix = fix_formula,
                        data = train_meta,
                        rnd = list(subject_id = ~ 1, range = ~ 1),
                        lambda = 30,
                        family = gaussian(link = "identity"))

lassoFeatures <- names(best_basic$coefficients)
lassoFeatures <- lassoFeatures[lassoFeatures %ni% c("(Intercept)")]
lassoFeatures <- lassoFeatures[grep("as.factor",lassoFeatures,invert=T)] ####
lassoFeatures <- unique(c(lassoFeatures))
#lassoFeatures <- gsub("range0.997860960117438", "range", lassoFeatures)
lassoFeatures <- lassoFeatures[grep("subject_id|BMI|range", lassoFeatures, invert = T)]
lassoFeatures <- unique(c(lassoFeatures, "subject_id", "BMI", "range"))

# Check the updated vector
lassoFeatures

bestglm_basic <- as.data.frame(basic_train[,lassoFeatures])
summary(bestglm_basic$BMI)
varlist_basic <- names(bestglm_basic)[which(names(bestglm_basic) %ni% 
                                              c("BMI", "subject_id", "range"))]
varstring_basic <- paste0(varlist_basic, collapse = " + ", sep = "")

# Use selected variables 
mymod_basic <- lme4::glmer(as.formula(paste0("BMI ~ ",varstring_basic, 
                                             " + (1|range) +  (1|subject_id)")), 
                           data = bestglm_basic, 
                           family = gaussian(link = "identity"))
mymodsum <- summary(mymod_basic)
mod_coef_df <- coef(mymodsum) %>% data.frame()

# Extract the fixed effects coefficients from the final model
coef_df <- as.data.frame(summary(mymod_basic)$coefficients)
coef_df$Feature <- rownames(coef_df)
coef_df <- coef_df %>% arrange(desc(Estimate))  # Sort by coefficients

# Filter top 10 features by absolute coefficient magnitude
top_features <- coef_df %>% 
  mutate(abs_estimate = abs(Estimate)) %>% 
  dplyr::filter(Feature != "(Intercept)") %>% 
  mutate(Feature = str_to_title(Feature),  # Capitalize the first letter of each word
         Feature = str_replace_all(Feature, "[._]", " ")) %>% 
  arrange(desc(abs_estimate)) %>% head(10)  

ggplot(top_features, aes(x = reorder(Feature, Estimate), y = Estimate)) +
  geom_bar(stat = "identity", fill = "#1C4C98") +
  coord_flip() + theme_bw() +
  ggtitle("Top Features & Coefficients from the final Basic Model") +
  xlab("Feature") + ylab("Coefficient Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12)) 

# prediction on reserved validation samples
test_basic <- test_basic[complete.cases(test_basic),] # complete cases

# grab only the columns that we have coefs for (starting at index 2 removes the intercept)
test_basic_vars <- test_basic %>% dplyr::select(rownames(mod_coef_df)[2:nrow(mod_coef_df)])
test_basic_vars$Intercept <- 1 # make a column for the intercept
test_basic_vars <- test_basic_vars %>% dplyr::select(Intercept, everything())

#--predict the risk scores without covariates--########
pred_risk_scores_basic <- as.matrix(test_basic_vars) %*% as.matrix(mod_coef_df$Estimate)

# combine the predicted and actual
pred_df_basic <- as.data.frame(cbind(test_basic$subject_id,
                                     test_basic$range,
                                     test_basic$BMI,
                                     scale(pred_risk_scores_basic))); colnames(pred_df_basic) <- c("subject_id",
                                                                                                   "time", 
                                                                                                   "actual",
                                                                                                   "predicted_basic")
# Calculate R-squared value
actual <- test_basic$BMI
predicted <- pred_risk_scores_basic
mean_actual <- mean(actual)
ss_total <- sum((actual - mean_actual)^2)
ss_residual <- sum((actual - predicted)^2)
r_squared <- 1 - (ss_residual / ss_total)

# Print R-squared value
r_squared
# Via function R2
#r2_general(predicted,actual)

# Extract the coefficients (excluding intercept) and calculate absolute values for feature importance
coef_df <- as.data.frame(best_basic$coefficients)
#coef_df <- as.data.frame(coef_df[rownames(coef_df) != "(Intercept)", ]) # Remove intercept
coef_df$Feature <- rownames(coef_df) 
coef_df$Feature <- gsub("range3.99098082279631", "time", coef_df$Feature)
coef_df$Importance <- coef_df$`best_basic$coefficients`

# Filter out rows where Importance is 0
#coef_df <- coef_df[coef_df$Importance != 0, ]

# Plot the feature importances
library(ggplot2)
ggplot(coef_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Basic Feature Importances from GLM Lasso Model") +
  xlab("Features") +
  ylab("Importance (Absolute Coefficients)") +
  # Add extra space around the plot using plot.margin
  #theme(plot.margin = margin(t = 20, r = 20, b = 20, l = 40)) +  # top, right, bottom, left margins
  annotate("text", 
           x = 1, 
           y = max(coef_df$Importance) * 0.95, 
           label = paste("R² =", round(r_squared, 3)), 
           size = 5, 
           color = "red", 
           hjust = 0.9)  # hjust = 0 will align the label to the left

# Save the plot to a file
ggsave(paste0(out_dir, "basic_feature_importances_with_r2.png"), 
       width = 8, height = 6, units = "in", dpi = 300)

file_path <- file.path(out_dir, "long_basic_predictions_feb20.csv")
#write.csv(pred_df_basic, file_path, row.names = FALSE)

#############################################################################################

### Meta data set
train_meta_filtered <- train_meta %>%
  dplyr::filter(complete.cases(.)) %>%  # Keep only rows with complete cases
  dplyr::group_by(subject_id) %>%       # Group by subject_id
  dplyr::filter(all(c("0", "6", "12") %in% range)) %>%  # Ensure all three ranges are present
  dplyr::ungroup()
numvar_meta <- c()
lambdavec <- seq(from = 10, to = 50, by = 1)
for (lambdy in lambdavec) {
  predictors <- setdiff(names(train_meta_filtered), c("BMI", "subject_id"))  
  predictors_escaped <- paste0("`", predictors, "`", collapse = " + ")
  fix_formula <- as.formula(paste("BMI ~", predictors_escaped))
  lm1_meta <- glmmLasso(fix = fix_formula,
                        data = train_meta_filtered,
                        rnd = list(subject_id = ~ 1, range = ~ 1),
                        lambda = lambdy,
                        family = gaussian(link = "identity"))
  summary(lm1_meta)
  lassoFeatures <- names(lm1_meta$coefficients[which(lm1_meta$coefficients != 0)])
  lassoFeatures <- lassoFeatures[lassoFeatures %ni% c("(Intercept)")]
  lassoFeatures <- lassoFeatures[grep("as.factor",lassoFeatures,invert=T)] ####
  lassoFeatures <- unique(c(lassoFeatures))
  numvar_meta <- c(numvar_meta, length(lassoFeatures))
}
plot(x = lambdavec, y = numvar_meta)
ggplot() +
  geom_point(aes(x = lambdavec, y = numvar_meta)) +
  geom_vline(xintercept = 38, color = "blue", linetype = "dashed") +
  theme_bw() +
  ggtitle("Only Meta Delta") +
  xlab("Penalty Coefficient (Lambda)") +
  ylab("Number of Included Variables")
ggsave(paste0(out_dir,"long_lambda_elbow_only_meta_feb20.png"), 
       width=6, height=4, units="in", dpi=320)

best_meta <- glmmLasso(fix = fix_formula,
                       data = train_meta_filtered,
                       rnd = list(subject_id = ~ 1, range = ~ 1),
                       lambda = 38,
                       family = gaussian(link = "identity"))

lassoFeatures <- names(best_meta$coefficients[which(best_meta$coefficients != 0)])
lassoFeatures <- lassoFeatures[lassoFeatures %ni% c("(Intercept)")]
lassoFeatures <- lassoFeatures[grep("as.factor",lassoFeatures,invert=T)] ####
lassoFeatures <- unique(c(lassoFeatures))
lassoFeatures <- gsub("range0.997860960117438", "range", lassoFeatures)
lassoFeatures <- lassoFeatures[grep("subject_id|BMI|range", lassoFeatures, invert = T)]
lassoFeatures <- unique(c(lassoFeatures, "subject_id", "BMI", "range"))

# Check the updated vector
lassoFeatures

bestglm_meta <- as.data.frame(train_meta_filtered[,lassoFeatures])
summary(bestglm_meta$BMI)
varlist_meta <- names(bestglm_meta)[which(names(bestglm_meta) %ni% 
                                            c("BMI", "subject_id", "range"))]
varstring_meta <- paste0(varlist_meta, collapse = " + ", sep = "")

# Use selected variables 
mymod_meta <- lme4::glmer(as.formula(paste0("BMI ~ ",varstring_meta, 
                                            " + (1|range) +  (1|subject_id)")), 
                          data = bestglm_meta, 
                          family = gaussian(link = "identity"))
mymodsum <- summary(mymod_meta)
mod_coef_df <- coef(mymodsum) %>% data.frame()

# Extract the fixed effects coefficients from the final model
coef_df <- as.data.frame(summary(mymod_meta)$coefficients)
coef_df$Feature <- rownames(coef_df)
coef_df <- coef_df %>% arrange(desc(Estimate))  # Sort by coefficients

# Filter top 10 features by absolute coefficient magnitude
top_features <- coef_df %>% 
  mutate(abs_estimate = abs(Estimate)) %>% 
  dplyr::filter(Feature != "(Intercept)") %>% 
  mutate(Feature = str_to_title(Feature),  # Capitalize the first letter of each word
         Feature = str_replace_all(Feature, "[._]", " ")) %>% 
  arrange(desc(abs_estimate)) %>% head(10)  

ggplot(top_features, aes(x = reorder(Feature, Estimate), y = Estimate)) +
  geom_bar(stat = "identity", fill = "#1C7C54") +
  coord_flip() + theme_bw() +
  ggtitle("Top Features & Coefficients from the final Meta Model") +
  xlab("Feature") + ylab("Coefficient Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12)) 

# prediction on reserved validation samples
test_meta <- test_meta[complete.cases(test_meta),] # complete cases
test_meta_filtered <- test_meta %>%
  dplyr::filter(complete.cases(.)) %>%  # Keep only rows with complete cases
  dplyr::group_by(subject_id) %>%       # Group by subject_id
  dplyr::filter(all(c("0", "6", "12") %in% range)) %>%  # Ensure all three ranges are present
  dplyr::ungroup()
# grab only the columns that we have coefs for (starting at index 2 removes the intercept)
test_meta_vars <- test_meta_filtered %>% 
  dplyr::select(rownames(mod_coef_df)[2:nrow(mod_coef_df)])
test_meta_vars$Intercept <- 1 # make a column for the intercept
test_meta_vars <- test_meta_vars %>% dplyr::select(Intercept, everything())

#--predict the risk scores without covariates--########
pred_risk_scores_meta <- as.matrix(test_meta_vars) %*% as.matrix(mod_coef_df$Estimate)

# combine the predicted and actual
pred_df_meta <- as.data.frame(cbind(test_meta_filtered$subject_id,
                                    test_meta_filtered$range,
                                    test_meta_filtered$BMI,
                                    scale(pred_risk_scores_meta))); colnames(pred_df_meta) <- c("subject_id", "time", "actual", "predicted_meta")
file_path <- file.path(out_dir, "long_meta_predictions_feb20.csv")
#rite.csv(pred_df_meta, file_path, row.names = FALSE)

#############################################################################################
### Taxa data set

# function to remove one of the two collinear pairs: 
remove_high_corr_vars <- function(data, threshold = 0.8) {
  data_without_subject_range <- data[, !colnames(data) %in% c("subject_id", "range")]
  cor_matrix <- cor(data_without_subject_range)
  high_corr_pairs <- which(abs(cor_matrix) > threshold, arr.ind = TRUE)
  high_corr_pairs <- high_corr_pairs[high_corr_pairs[,1] < high_corr_pairs[,2], ]
  to_remove <- unique(colnames(data_without_subject_range)[high_corr_pairs[,2]])
  data_cleaned <- data[, !colnames(data) %in% to_remove]
  return(data_cleaned)
}

train_taxa_cleaned <- remove_high_corr_vars(train_taxa)
train_taxa_filtered <- train_taxa_cleaned %>%
  dplyr::filter(complete.cases(.)) %>%  # Keep only rows with complete cases
  dplyr::group_by(subject_id) %>%       # Group by subject_id
  dplyr::filter(all(c("0", "6", "12") %in% range)) %>%  # Ensure all three ranges are present
  dplyr::ungroup()

numvar_taxa <- c()
lambdavec <- seq(from = 32, to = 50, by = 1)
for (lambdy in lambdavec) {
  predictors <- setdiff(names(train_taxa_filtered), c("BMI", "subject_id", "range"))
  predictors_escaped <- paste0("`", predictors, "`", collapse = " + ")
  fix_formula <- as.formula(paste("BMI ~", predictors_escaped))
  lm1_taxa <- glmmLasso(fix = fix_formula,
                        data = train_taxa_filtered,
                        rnd = list(subject_id = ~ 1, range = ~ 1),
                        lambda = lambdy,
                        family = gaussian(link = "identity"))
  summary(lm1_taxa)
  lassoFeatures <- names(lm1_taxa$coefficients[which(lm1_taxa$coefficients != 0)])
  lassoFeatures <- lassoFeatures[lassoFeatures %ni% c("(Intercept)")]
  lassoFeatures <- lassoFeatures[grep("as.factor",lassoFeatures,invert=T)] ####
  lassoFeatures <- unique(c(lassoFeatures))
  numvar_taxa <- c(numvar_taxa, length(lassoFeatures))
}

plot(x = lambdavec, y = numvar_taxa)
ggplot() +
  geom_point(aes(x = lambdavec, y = numvar_taxa)) +
  geom_vline(xintercept = 42, color = "blue", linetype = "dashed") +
  theme_bw() +
  ggtitle("Only Taxa Delta") +
  xlab("Penalty Coefficient (Lambda)") +
  ylab("Number of Included Variables")
ggsave(paste0(out_dir,"long_lambda_elbow_only_taxa_feb20.png"), 
       width=6, height=4, units="in", dpi=320)

best_taxa <- glmmLasso(fix = fix_formula,
                       data = train_taxa_filtered,
                       rnd = list(subject_id = ~ 1, range = ~ 1),
                       lambda = 43,
                       family = gaussian(link = "identity"))

lassoFeatures <- names(best_taxa$coefficients[which(best_taxa$coefficients != 0)])
lassoFeatures <- lassoFeatures[lassoFeatures %ni% c("(Intercept)")]
lassoFeatures <- lassoFeatures[grep("as.factor",lassoFeatures,invert=T)] ####
lassoFeatures <- unique(c(lassoFeatures))
lassoFeatures <- gsub("range0.997860960117438", "range", lassoFeatures)
lassoFeatures <- lassoFeatures[grep("subject_id|BMI|range", lassoFeatures, invert = T)]
lassoFeatures <- unique(c(lassoFeatures, "subject_id", "BMI", "range"))
lassoFeatures <- gsub("`", "", lassoFeatures)
# Check the updated vector
lassoFeatures
bestglm_taxa <- as.data.frame(train_taxa_filtered[,lassoFeatures])
summary(bestglm_taxa$BMI)
varlist_taxa <- names(bestglm_taxa)[which(names(bestglm_taxa) %ni% 
                                            c("BMI", "subject_id", "range"))]
varstring_taxa <- paste0(varlist_taxa, collapse = " + ", sep = "")
varstring_taxa <- paste0("`", varlist_taxa, "`", collapse = " + ")
# Use selected variables 
mymod_taxa <- lme4::glmer(as.formula(paste0("BMI ~ ",varstring_taxa, 
                                            " + (1|range) +  (1|subject_id)")), 
                          data = bestglm_taxa, 
                          family = gaussian(link = "identity"))
mymodsum <- summary(mymod_taxa)
mod_coef_df <- coef(mymodsum) %>% data.frame()

# prediction on reserved validation samples
test_taxa <- test_taxa[complete.cases(test_taxa),] # complete cases
test_taxa_filtered <- test_taxa %>%
  dplyr::filter(complete.cases(.)) %>%  # Keep only rows with complete cases
  dplyr::group_by(subject_id) %>%       # Group by subject_id
  dplyr::filter(all(c("0", "6", "12") %in% range)) %>%  # Ensure all three ranges are present
  dplyr::ungroup()
# grab only the columns that we have coefs for (starting at index 2 removes the intercept)
rownames(mod_coef_df) <- gsub("`", "", rownames(mod_coef_df))
test_taxa_vars <- test_taxa_filtered %>% 
  dplyr::select(rownames(mod_coef_df)[2:nrow(mod_coef_df)])
test_taxa_vars$Intercept <- 1 # make a column for the intercept
test_taxa_vars <- test_taxa_vars %>% dplyr::select(Intercept, everything())

#--predict the risk scores without covariates--########
pred_risk_scores_taxa <- as.matrix(test_taxa_vars) %*% as.matrix(mod_coef_df$Estimate)

# combine the predicted and actual
pred_df_taxa <- as.data.frame(cbind(test_taxa_filtered$subject_id,
                                    test_taxa_filtered$range,
                                    test_taxa_filtered$BMI,
                                    scale(pred_risk_scores_taxa))); colnames(pred_df_taxa) <- c("subject_id",
                                                                                                "time", 
                                                                                                "actual",
                                                                                                "predicted_taxa")
file_path <- file.path(out_dir, "long_taxa_predictions_feb20.csv")
#write.csv(pred_df_taxa, file_path, row.names = FALSE)

##########################################################################################
### Micom data set
train_micom_cleaned <- remove_high_corr_vars(train_micom, threshold = 0.9)
train_micom_filtered <- train_micom_cleaned %>%
  dplyr::filter(complete.cases(.)) %>%  # Keep only rows with complete cases
  dplyr::group_by(subject_id) %>%       # Group by subject_id
  dplyr::filter(all(c("0", "6", "12") %in% range)) %>%  # Ensure all three ranges are present
  dplyr::ungroup()

numvar_micom <- c()
lambdavec <- seq(from = 10, to = 50, by = 1)
for (lambdy in lambdavec) {
  predictors <- setdiff(names(train_micom_filtered), c("BMI", "subject_id", "range"))
  predictors_escaped <- paste0("`", predictors, "`", collapse = " + ")
  fix_formula <- as.formula(paste("BMI ~", predictors_escaped))
  lm1_micom <- glmmLasso(fix = fix_formula,
                         data = train_micom_filtered,
                         rnd = list(subject_id = ~ 1, range = ~ 1),
                         lambda = lambdy,
                         family = gaussian(link = "identity"))
  summary(lm1_micom)
  lassoFeatures <- names(lm1_micom$coefficients[which(lm1_micom$coefficients != 0)])
  lassoFeatures <- lassoFeatures[lassoFeatures %ni% c("(Intercept)")]
  lassoFeatures <- lassoFeatures[grep("as.factor",lassoFeatures,invert=T)] ####
  lassoFeatures <- unique(c(lassoFeatures))
  numvar_micom <- c(numvar_micom, length(lassoFeatures))
}
plot(x = lambdavec, y = numvar_micom)
ggplot() +
  geom_point(aes(x = lambdavec, y = numvar_micom)) +
  geom_vline(xintercept = 35, color = "blue", linetype = "dashed") +
  theme_bw() +
  ggtitle("Only Micom Delta") +
  xlab("Penalty Coefficient (Lambda)") +
  ylab("Number of Included Variables")
ggsave(paste0(out_dir,"long_lambda_elbow_only_micom_feb20.png"), 
       width=6, height=4, units="in", dpi=320)

best_micom <- glmmLasso(fix = fix_formula,
                        data = train_micom_filtered,
                        rnd = list(subject_id = ~ 1, range = ~ 1),
                        lambda = 35,
                        family = gaussian(link = "identity"))

lassoFeatures <- names(best_micom$coefficients[which(best_micom$coefficients != 0)])
lassoFeatures <- lassoFeatures[lassoFeatures %ni% c("(Intercept)")]
lassoFeatures <- lassoFeatures[grep("as.factor",lassoFeatures,invert=T)] ####
lassoFeatures <- unique(c(lassoFeatures))
lassoFeatures <- gsub("range0.997860960117438", "range", lassoFeatures)
lassoFeatures <- lassoFeatures[grep("subject_id|BMI|range", lassoFeatures, invert = T)]
lassoFeatures <- unique(c(lassoFeatures, "subject_id", "BMI", "range"))
lassoFeatures <- gsub("`", "", lassoFeatures)
# Check the updated vector
lassoFeatures
bestglm_micom <- as.data.frame(train_micom_filtered[,lassoFeatures])
summary(bestglm_micom$BMI)
varlist_micom <- names(bestglm_micom)[which(names(bestglm_micom) %ni% 
                                              c("BMI", "subject_id", "range"))]
varstring_micom <- paste0(varlist_micom, collapse = " + ", sep = "")
varstring_micom <- paste0("`", varlist_micom, "`", collapse = " + ")
# Use selected variables 
mymod_micom <- lme4::glmer(as.formula(paste0("BMI ~ ",varstring_micom, 
                                             " + (1|range) +  (1|subject_id)")), 
                           data = bestglm_micom, 
                           family = gaussian(link = "identity"))
mymodsum <- summary(mymod_micom)
mod_coef_df <- coef(mymodsum) %>% data.frame()

# prediction on reserved validation samples
test_micom <- test_micom[complete.cases(test_micom),] # complete cases
test_micom_filtered <- test_micom %>%
  dplyr::filter(complete.cases(.)) %>%  # Keep only rows with complete cases
  dplyr::group_by(subject_id) %>%       # Group by subject_id
  dplyr::filter(all(c("0", "6", "12") %in% range)) %>%  # Ensure all three ranges are present
  dplyr::ungroup()
# grab only the columns that we have coefs for (starting at index 2 removes the intercept)
rownames(mod_coef_df) <- gsub("`", "", rownames(mod_coef_df))
test_micom_vars <- test_micom_filtered %>% 
  dplyr::select(rownames(mod_coef_df)[2:nrow(mod_coef_df)])
test_micom_vars$Intercept <- 1 # make a column for the intercept
test_micom_vars <- test_micom_vars %>% dplyr::select(Intercept, everything())

#--predict the risk scores without covariates--########
pred_risk_scores_micom <- as.matrix(test_micom_vars) %*% as.matrix(mod_coef_df$Estimate)

# combine the predicted and actual
pred_df_micom<- as.data.frame(cbind(test_micom_filtered$subject_id,
                                    test_micom_filtered$range,
                                    test_micom_filtered$BMI,
                                    scale(pred_risk_scores_micom))); colnames(pred_df_micom) <- c("subject_id",
                                                                                                  "time", 
                                                                                                  "actual",
                                                                                                  "predicted_micom")
file_path <- file.path(out_dir, "long_micom_predictions_feb20.csv")
#write.csv(pred_df_micom, file_path, row.names = FALSE)

##########################################################################################
### Pathway data set
train_pathway_cleaned <- remove_high_corr_vars(train_pathway, threshold = 0.9) %>% 
  dplyr::select(-c("cohort_number","time_y", "bmi_prs" ))
train_pathway_filtered <- train_pathway_cleaned %>%
  dplyr::filter(complete.cases(.)) %>%  # Keep only rows with complete cases
  dplyr::group_by(subject_id) %>%       # Group by subject_id
  dplyr::filter(all(c("0", "6", "12") %in% range)) %>%  # Ensure all three ranges are present
  dplyr::ungroup()
numvar_path <- c()
lambdavec <- seq(from = 10, to = 50, by = 1)
for (lambdy in lambdavec) {
  predictors <- setdiff(names(train_pathway_filtered), c("BMI", "subject_id", "range"))
  predictors_escaped <- paste0("`", predictors, "`", collapse = " + ")
  fix_formula <- as.formula(paste("BMI ~", predictors_escaped))
  lm1_pathway <- glmmLasso(fix = fix_formula,
                           data = train_pathway_filtered,
                           rnd = list(subject_id = ~ 1, range = ~ 1),
                           lambda = lambdy,
                           family = gaussian(link = "identity"))
  summary(lm1_pathway)
  lassoFeatures <- names(lm1_pathway$coefficients[which(lm1_pathway$coefficients != 0)])
  lassoFeatures <- lassoFeatures[lassoFeatures %ni% c("(Intercept)")]
  lassoFeatures <- lassoFeatures[grep("as.factor",lassoFeatures,invert=T)] ####
  lassoFeatures <- unique(c(lassoFeatures))
  numvar_path <- c(numvar_path, length(lassoFeatures))
}
plot(x = lambdavec, y = numvar_path)
ggplot() +
  geom_point(aes(x = lambdavec, y = numvar_path)) +
  geom_vline(xintercept = 30, color = "blue", linetype = "dashed") +
  theme_bw() +
  ggtitle("Only Pathway Delta") +
  xlab("Penalty Coefficient (Lambda)") +
  ylab("Number of Included Variables")
ggsave(paste0(out_dir,"long_lambda_elbow_only_pathway_feb20.png"), 
       width=6, height=4, units="in", dpi=320)

best_path <- glmmLasso(fix = fix_formula,
                       data = train_pathway_filtered,
                       rnd = list(subject_id = ~ 1, range = ~ 1),
                       lambda = 32,
                       family = gaussian(link = "identity"))

lassoFeatures <- names(best_path$coefficients[which(best_path$coefficients != 0)])
lassoFeatures <- lassoFeatures[lassoFeatures %ni% c("(Intercept)")]
lassoFeatures <- lassoFeatures[grep("as.factor",lassoFeatures,invert=T)] ####
lassoFeatures <- unique(c(lassoFeatures))
lassoFeatures <- gsub("range0.997860960117438", "range", lassoFeatures)
lassoFeatures <- lassoFeatures[grep("subject_id|BMI|range", lassoFeatures, invert = T)]
lassoFeatures <- unique(c(lassoFeatures, "subject_id", "BMI", "range"))
lassoFeatures <- gsub("`", "", lassoFeatures)

# Check the updated vector
lassoFeatures
bestglm_path <- as.data.frame(train_pathway_filtered[,lassoFeatures])
summary(bestglm_path$BMI)
varlist_path <- names(bestglm_path)[which(names(bestglm_path) %ni% 
                                            c("BMI", "subject_id", "range"))]
varstring_path <- paste0(varlist_path, collapse = " + ", sep = "")
varstring_path <- paste0("`", varlist_path, "`", collapse = " + ")

# Use selected variables 
mymod_path <- lme4::glmer(as.formula(paste0("BMI ~ ",varstring_path, 
                                            " + (1|range) +  (1|subject_id)")), 
                          data = bestglm_path, 
                          family = gaussian(link = "identity"))
mymodsum <- summary(mymod_path)
mod_coef_df <- coef(mymodsum) %>% data.frame()

# Plot pathway features
# Extract the fixed effects coefficients from the final model
coef_df <- as.data.frame(summary(mymod_path)$coefficients)
coef_df$Feature <- rownames(coef_df)
coef_df <- coef_df %>% arrange(desc(Estimate))  # Sort by coefficients

# Filter top 10 features by absolute coefficient magnitude
top_features <- coef_df %>% 
  mutate(abs_estimate = abs(Estimate)) %>% 
  dplyr::filter(Feature != "(Intercept)") %>% 
  mutate(Feature = str_to_title(Feature),  # Capitalize the first letter of each word
         Feature = str_replace_all(Feature, "[._]", " ")) %>% 
  arrange(desc(abs_estimate)) %>% head(10)  

ggplot(top_features, aes(x = reorder(Feature, Estimate), y = Estimate)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() + theme_bw() +
  ggtitle("Top Features & Coefficients from the final Pathway Model") +
  xlab("Feature") + ylab("Coefficient Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12)) 

# prediction on reserved validation samples
test_path <- test_pathway[complete.cases(test_pathway),] # complete cases
test_path_filtered <- test_path %>%
  dplyr::filter(complete.cases(.)) %>%  # Keep only rows with complete cases
  dplyr::group_by(subject_id) %>%       # Group by subject_id
  dplyr::filter(all(c("0", "6", "12") %in% range)) %>%  # Ensure all three ranges are present
  dplyr::ungroup()

# grab only the columns that we have coefs for (starting at index 2 removes the intercept)
rownames(mod_coef_df) <- gsub("`", "", rownames(mod_coef_df))
test_path_vars <- test_path_filtered %>% 
  dplyr::select(rownames(mod_coef_df)[2:nrow(mod_coef_df)])
test_path_vars$Intercept <- 1 # make a column for the intercept
test_path_vars <- test_path_vars %>% dplyr::select(Intercept, everything())

#--predict the risk scores without covariates--########
pred_risk_scores_path <- as.matrix(test_path_vars) %*% as.matrix(mod_coef_df$Estimate)

# combine the predicted and actual
pred_df_path<- as.data.frame(cbind(test_path_filtered$subject_id,
                                   test_path_filtered$range,
                                   test_path_filtered$BMI,
                                   scale(pred_risk_scores_path))); colnames(pred_df_path) <- c("subject_id",
                                                                                               "time", 
                                                                                               "actual",
                                                                                               "predicted_pathway")
file_path <- file.path(out_dir, "long_path_predictions_feb20.csv")
write.csv(pred_df_path, file_path, row.names = FALSE)

##### Compare models
# Convert the first 2 columns to factors
pred_df_basic[, head(names(pred_df_basic), 2)] <- 
  lapply(pred_df_basic[, head(names(pred_df_basic), 2)], as.factor)
pred_df_meta[, head(names(pred_df_meta), 2)] <- 
  lapply(pred_df_meta[, head(names(pred_df_meta), 2)], as.factor)
pred_df_taxa[, head(names(pred_df_taxa), 2)] <- 
  lapply(pred_df_taxa[, head(names(pred_df_taxa), 2)], as.factor)
pred_df_micom[, head(names(pred_df_micom), 2)] <- 
  lapply(pred_df_micom[, head(names(pred_df_micom), 2)], as.factor)
pred_df_path[, head(names(pred_df_path), 2)] <- 
  lapply(pred_df_path[, head(names(pred_df_path), 2)], as.factor)

met_basic <- merge(pred_df_basic, pred_df_meta, by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

met_tax <- merge(met_basic, pred_df_taxa, 
                 by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

met_tax_micom <- merge(met_tax, pred_df_micom, 
                       by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

met_tax_micom_path <- merge(met_tax_micom, pred_df_path, 
                            by = c("subject_id", "time")) %>% 
  dplyr::select(-c(actual.y)) %>% rename(actual = actual.x)

all_omic <- unique(met_tax_micom_path)

# Center and scale the last 5 columns of the dataframe
all_omic[, 3:8] <- scale(all_omic[, 3:8])
head(all_omic)

########################################################################################
mod_dat = all_omic %>% rename(bmi = actual, 
                              Time = time,
                              Cluster = subject_id,
                              y_new_meta_only = predicted_meta,
                              y_new_basic = predicted_basic,
                              y_new_micom_only = predicted_micom,
                              y_new_path_only = predicted_pathway,
                              y_new_tax_only = predicted_taxa)
### Change time to contimuous
lmer_basic <- lmer(bmi ~ y_new_basic + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_meta <- lmer(bmi ~ y_new_meta_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_micom <- lmer(bmi ~ y_new_micom_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_path <- lmer(bmi ~ y_new_path_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_tax <- lmer(bmi ~ y_new_tax_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_time <- lmer(bmi ~ Time+ (1|Cluster), data = mod_dat, REML = FALSE)

### Single plus omic 
lmer_basic <- lmer(bmi ~ y_new_basic + Time + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_meta_b <- lmer(bmi ~ y_new_basic + y_new_meta_only + Time + (1|Cluster), data = mod_dat, REML = FALSE)
lmer_micom_b <- lmer(bmi ~ y_new_basic + y_new_micom_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_path_b <- lmer(bmi ~ y_new_basic + y_new_path_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_tax_b <- lmer(bmi ~ y_new_basic + y_new_tax_only + Time+ (1|Cluster), data = mod_dat, REML = FALSE)
lmer_time_b <- lmer(bmi ~ y_new_basic + Time+ (1|Cluster), data = mod_dat, REML = FALSE)

### Combined PTEV models 
basic <- lmer(bmi ~ y_new_basic + Time+ (1|Cluster), data = mod_dat)
meta_basic <- lmer(bmi ~ y_new_basic + y_new_meta_only + Time + (1|Cluster), data = mod_dat)
meta_basic_tax <- lmer(bmi ~ y_new_basic + y_new_meta_only + 
                         y_new_tax_only + Time + (1|Cluster), data = mod_dat)
meta_basic_tax_path <- lmer(bmi ~ y_new_basic + y_new_meta_only + y_new_tax_only + 
                              y_new_path_only + Time + (1|Cluster), data = mod_dat)
meta_basic_tax_path_micom <- lmer(bmi ~ y_new_basic + y_new_meta_only + 
                                    y_new_tax_only + y_new_path_only + y_new_micom_only + 
                                    Time+ (1|Cluster), data = mod_dat)

# Combined PTEV models sequential 
sjPlot::tab_model(basic, meta_basic, 
                  meta_basic_tax, meta_basic_tax_path, 
                  meta_basic_tax_path_micom,
                  title = "glmlasso long sequential lmer models",
                  string.pred = "Predictors",
                  string.est = "Estimate",
                  string.std = "std. Beta",
                  string.ci = "95% CI",
                  string.se = "std. Error",
                  p.style = c("numeric"), 
                  p.threshold = c(0.05),
                  #dv.labels = model_titles,
                  auto.label = FALSE)

# https://joshuawiley.com/MonashHonoursStatistics/LMM_Comparison.html#nested-models-in-r
# Basic vs Basic + meta
nobs(lmer_basic)
nobs(lmer_meta_b)
logLik(lmer_basic)
logLik(lmer_meta_b)
anova(lmer_basic, lmer_meta_b, test = "LRT")
anova(lmer_basic, lmer_micom_b, test = "LRT")
anova(lmer_basic, lmer_tax_b, test = "LRT")
anova(lmer_basic, lmer_path_b, test = "LRT")

glmmlong_models <- list(
  c("lmer_basic", "lmer_meta_b"),
  c("lmer_basic", "lmer_micom_b"),
  c("lmer_basic", "lmer_tax_b"),
  c("lmer_basic", "lmer_path_b"))
library(lme4)  # Make sure the lme4 package is loaded=
anova_results <- list()
for (model_pair in glmmlong_models) {
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
writeLines(html_table, "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/feb20/lasso_anova_table.html")
########################################################################################

# Make linear models MSE
# Combined models 
lm_meta_tax_time <- lme(actual ~ predicted_meta + predicted_taxa + time, 
                        random = ~1|subject_id, data = all_omic)
summary(lm_meta_tax_time)

### Checks before models 
library(car)
vif(lm(actual ~ predicted_meta + predicted_taxa + predicted_micom + time, data = all_omic))
summary(all_omic)
table(all_omic$subject_id)
cor(all_omic[, c("predicted_meta", "predicted_taxa", "predicted_micom")])
hist(all_omic$actual)
anyDuplicated(all_omic)

# Fit a simpler model first
lm_path <- lme(actual ~ predicted_pathway, random = ~1|subject_id, data = all_omic)
lm_micom <- lme(actual ~ predicted_micom, random = ~1|subject_id, data = all_omic)
lm_tax <- lme(actual ~ predicted_taxa, random = ~1|subject_id, data = all_omic)
lm_meta <- lme(actual ~ predicted_meta, random = ~1|subject_id, data = all_omic)

anova(lm_meta, lm_path)
anova(lm_meta, lm_micom)

model_titles <- c("Model 1: Meta", 
                  "Model 2: Taxa", 
                  "Model 3: Micom", 
                  "Model 4: Pathway")

sjPlot::tab_model(lm_meta, lm_tax, 
                  lm_micom, lm_path, 
                  title = "Comparing single omic glmLASSO models 6m, BL, 12m",
                  string.pred = "Predictors",
                  string.est = "Estimate",
                  string.std = "std. Beta",
                  string.ci = "95% CI",
                  string.se = "std. Error",
                  p.style = c("numeric"), 
                  p.threshold = c(0.05),
                  dv.labels = model_titles,
                  auto.label = FALSE)
### Combined 
lm_meta_tax <- lme(actual ~ predicted_meta + predicted_taxa, 
                   random = ~1|subject_id, data = all_omic)
summary(lm_meta_tax)
lm_meta_tax_micom <- lme(actual ~ predicted_meta + predicted_taxa + 
                           predicted_micom, random = ~1|subject_id, data = all_omic)
summary(lm_meta_tax_micom)
lm_meta_tax_micom_path <- lme(actual ~ predicted_meta + predicted_taxa + 
                                predicted_micom + predicted_pathway, 
                              random = ~1|subject_id, data = all_omic)
lm_no_meta <- lme(actual ~ predicted_taxa + 
                    predicted_micom + predicted_pathway, 
                  random = ~1|subject_id, data = all_omic)
summary(lm_meta_tax_micom_path)
anova(lm_meta, lm_meta_tax)
anova(lm_meta_tax, lm_meta_tax_micom)
anova(lm_meta_tax_micom, lm_meta_tax_micom_path)
anova(lm_meta_tax, lm_meta_tax_micom_path)
anova(lm_meta, lm_meta_tax_micom_path)

model_titles <- c("Model 1: Meta", 
                  "Model 2: Meta + Taxa", 
                  "Model 3: Meta + Taxa + Micom", 
                  "Model 4: Taxa + Micom + Pathway", 
                  "Model 5: Meta + Taxa + Micom + Pathway")

sjPlot::tab_model(lm_meta, lm_meta_tax, 
                  lm_meta_tax_micom, lm_no_meta, lm_meta_tax_micom_path, 
                  title = "Comparing combined omic models 6m, BL, 12m",
                  string.pred = "Predictors",
                  string.est = "Estimate",
                  string.std = "std. Beta",
                  string.ci = "95% CI",
                  string.se = "std. Error",
                  p.style = c("numeric"), 
                  p.threshold = c(0.05),
                  dv.labels = model_titles,
                  auto.label = FALSE)

#### PLOT AIC values
# Extract AIC and BIC for each model
aic_values <- c(AIC(lm_meta), AIC(lm_meta_tax), 
                AIC(lm_meta_tax_micom), AIC(lm_meta_tax_micom_path))
bic_values <- c(BIC(lm_meta), BIC(lm_meta_tax), 
                BIC(lm_meta_tax_micom), BIC(lm_meta_tax_micom_path))

# Create a data frame for plotting
model_names <- c("lm_meta","lm_meta_tax", "lm_meta_tax_micom", "lm_meta_tax_micom_path")
df <- data.frame(Model = model_names, AIC = aic_values, BIC = bic_values)

# Plot AIC and BIC side by side
ggplot(df, aes(x = Model)) +
  geom_bar(aes(y = AIC), stat = "identity", position = "dodge", fill = "blue", alpha = 0.7) +
  geom_bar(aes(y = BIC), stat = "identity", position = "dodge", fill = "red", alpha = 0.7) +
  ylab("Value") +
  ggtitle("Comparison of AIC and BIC for Models") +
  theme_minimal() +
  scale_y_continuous(sec.axis = sec_axis(~., name = "BIC")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))


### Check Correlations
cor_matrix <- cor((all_omic)[3:7], 
                  use = "pairwise.complete.obs", method = "pearson")
melted_cor_matrix <- melt(cor_matrix, na.rm = TRUE)
ggplot(melted_cor_matrix, 
       aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", 
                       high = "red") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 3),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 3)) + 
  labs(title = "Correlations btwn single omics & actual BMI", fill = "Correlation")
