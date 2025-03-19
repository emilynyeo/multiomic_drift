# Trying out long lasso
rm(list = ls())
#source("zc_functions.R") 
library(pacman)
p_load(tools, reticulate, viridis, tidyplots, patchwork, jsonlite, maps, ggvenn, 
       caret, caretEnsemble, glmnet, xgboost, ggplot2, glmmLasso, corrplot,
       readr, plyr, dplyr, tidyr, purrr, tibble, stringr, psych, randomForest,  
       reshape2, scales, gridExtra, plotly, sf, tidyverse, naniar, VIM)
'%ni%' <- Negate('%in%')

out_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/feb20_long/"

long_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/merf_dfs/5.combined/"

test <- read.csv(file.path(long_dir, 'feb20_test_merged_all_omics_raw_meta.csv'))
train <- read.csv(file.path(long_dir, 'feb20_training_merged_all_omics_raw_meta.csv'))

train_slim <- train %>% dplyr::select(-c("Unnamed..0_merged_data", "sample_id", 
                                         "subject_id", "record_id"))
train_slim$time <- as.factor(train_slim$time)
train_slim$all_samples <- as.factor(train_slim$all_samples)
train_slim$randomized_group <- as.numeric(train_slim$randomized_group)
train_slim$sex <- as.numeric(train_slim$sex)
train_slim$race <- as.numeric(train_slim$race)
test_slim <- test %>% 
  dplyr::select(-c("Unnamed..0_merged_data", "sample_id", "subject_id",
                   "record_id"))
test_slim$time <- as.factor(test_slim$time)
test_slim$all_samples <- as.factor(test_slim$all_samples)
test_slim$randomized_group <- as.numeric(test_slim$randomized_group)
test_slim$sex <- as.numeric(test_slim$sex)
test_slim$race <- as.numeric(test_slim$race)
train_slim <- train_slim %>% rename(subject_id = all_samples,
                                    BMI = outcome_BMI_fnl,
                                    range = time,
                                    homo_ir = HOMA_IR,
                                    insulin = Insulin_endo,
                                    LDL = LDL_Calculated,
                                    HDL = HDL_Total_Direct_lipid,
                                    HbA1C = Glucose)
test_slim <- test_slim %>% rename(subject_id = all_samples,
                                  BMI = outcome_BMI_fnl,
                                  range = time,
                                  homo_ir = HOMA_IR,
                                  insulin = Insulin_endo,
                                  LDL = LDL_Calculated,
                                  HDL = HDL_Total_Direct_lipid,
                                  HbA1C = Glucose)

# Define the column names based on your lists
basic <- c('subject_id','BMI', 'range','age', 'sex')
meta_keep <- c('subject_id','BMI', 'range', 'randomized_group', 'sex', 'race', 
               'age', 'HbA1C', 'HDL', 'homo_ir', 'insulin', 'LDL')
only_taxa <- c('subject_id','BMI', 'range', 
               grep("^g__", names(train_slim), value = TRUE))
proton_column <- which(names(train_slim) == "proton")
carbon_dioxide_column <- which(names(train_slim) == "Carbon.dioxide")
only_micom <- c('subject_id','BMI', 'range', 
                names(train_slim)[proton_column:carbon_dioxide_column])
exclude_columns <- unique(c(meta_keep, only_taxa, only_micom,
                            "leptin", "CRP", "ghrelin", "peptide_yy",
                            "cholesterol","tgcyd","Weight"))
only_pathway <- c('subject_id', 'BMI', 'range', 
                  setdiff(names(train_slim), exclude_columns))

# Create data frames based on the columns defined
train_basic <- train_slim[, basic, drop = FALSE]
train_meta <- train_slim[, meta_keep, drop = FALSE]
train_taxa <- train_slim[, only_taxa, drop = FALSE]
train_micom <- train_slim[, only_micom, drop = FALSE]
train_pathway <- train_slim[, only_pathway, drop = FALSE]

test_basic <- test_slim[, basic, drop = FALSE]
test_meta <- test_slim[, meta_keep, drop = FALSE]
test_taxa <- test_slim[, only_taxa, drop = FALSE]
test_micom <- test_slim[, only_micom, drop = FALSE]
test_pathway <- test_slim[, only_pathway, drop = FALSE]

train_basic <- unique(train_basic)
train_meta <- unique(train_meta)
train_micom <- unique(train_micom)
train_pathway <- unique(train_pathway)
train_taxa <- unique(train_taxa)

test_basic <- unique(test_basic)
test_meta <- unique(test_meta)
test_micom <- unique(test_micom)
test_pathway <- unique(test_pathway)
test_taxa <- unique(test_taxa)

subject_id_count <- train_meta %>%
  dplyr::filter(range %in% c(0, 6, 12)) %>%
  dplyr::group_by(subject_id) %>%
  dplyr::summarize(range_count = n_distinct(range))  # Count the distinct range values

# Filter `subject_id`s with fewer than 3 unique range values (0, 6, 12)
missing_subjects <- subject_id_count %>%
  dplyr::filter(range_count < 3) %>%
  pull(subject_id)

# Print the result
missing_subjects
# Remove rows from each dataframe where `subject_id` is in `missing_subjects`
train_basic <- train_basic %>%
  dplyr::filter(!subject_id %in% missing_subjects)

train_meta <- train_meta %>%
  dplyr::filter(!subject_id %in% missing_subjects)

train_micom <- train_micom %>%
  dplyr::filter(!subject_id %in% missing_subjects)

train_pathway <- train_pathway %>%
  dplyr::filter(!subject_id %in% missing_subjects)

train_taxa <- train_taxa %>%
  dplyr::filter(!subject_id %in% missing_subjects)

################################################################################
### BASIC data set
numvar_basic <- c()
lambdavec <- seq(from = 10, to = 50, by = 1)
for (lambdy in lambdavec) {
  predictors <- setdiff(names(train_basic), c("BMI", "subject_id"))  
  predictors_escaped <- paste0("`", predictors, "`", collapse = " + ")
  fix_formula <- as.formula(paste("BMI ~", predictors_escaped))
  lm1_basic <- glmmLasso(fix = fix_formula,
                         data = train_basic,
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

bestglm_basic <- as.data.frame(train_basic[,lassoFeatures])
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
train_pathway_cleaned <- remove_high_corr_vars(train_pathway, threshold = 0.9)
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
#write.csv(pred_df_path, file_path, row.names = FALSE)

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
