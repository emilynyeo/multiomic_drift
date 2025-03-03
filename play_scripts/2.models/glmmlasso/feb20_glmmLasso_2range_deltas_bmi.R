# Trying out long lasso
rm(list = ls())
#source("zc_functions.R") 
#install.packages("rsq")  
library(rsq)
library(pacman)
p_load(tools, reticulate, viridis, tidyplots, patchwork, jsonlite, maps, ggvenn, 
       caret, caretEnsemble, glmnet, xgboost, ggplot2, glmmLasso, corrplot,
       readr, plyr, dplyr, tidyr, purrr, tibble, stringr, psych, randomForest,  
       reshape2, scales, gridExtra, plotly, sf, tidyverse, naniar, VIM)
library(nlme)
library(gridExtra)
library(sjPlot)
library(htmltools)
library(officer)
library(flextable)
library(webshot)
library(apaTables)
library(MuMIn)
'%ni%' <- Negate('%in%')
r2_general <-function(preds,actual){ 
  return(1- sum((preds - actual) ^ 2)/sum((actual - mean(actual))^2))
}

out_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/feb20/"
data_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/zachs_rerun/drift_fs/csv/all_omic_processed_data/deltas/"
train <- read_csv(paste0(data_dir, "feb20_all_delta_train.csv"))
test <- read_csv(paste0(data_dir, "feb20_all_delta_test.csv"))

train_slim <- train %>% dplyr::select(-c("leptin", "peptide_yy","Weight",
                                         "cholesterol","ghrelin","CRP","tgcyd",
                                         "LDL_Calculated", "HDL_Total_Direct_lipid", 
                                         "Triglyceride_lipid", "outcome_wt_fnl","outcome_BMI_fnl",
                                         "Insulin_endo", "Glucose", "HOMA_IR")) %>% rename(age = age.x)
train_slim$range <- as.factor(train_slim$range)

test_slim <- test %>% dplyr::select(-c("leptin", "peptide_yy","Weight",
                                       "cholesterol","ghrelin","CRP","tgcyd",
                                       "LDL_Calculated", "HDL_Total_Direct_lipid", 
                                       "Triglyceride_lipid", "outcome_wt_fnl","outcome_BMI_fnl",
                                       "Insulin_endo", "Glucose", "HOMA_IR")) %>% rename(age = age.x)
test_slim$range <- as.factor(test_slim$range)

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
train_slim$subject_id <- as.factor(train_slim$subject_id)
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
lassoFeatures <- gsub("range0.997860960117438", "range", lassoFeatures)
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
r2_general(predicted,actual)

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
           label = paste("R� =", round(r_squared, 3)), 
           size = 5, 
           color = "red", 
           hjust = 0.9)  # hjust = 0 will align the label to the left

# Save the plot to a file
ggsave(paste0(out_dir, "basic_feature_importances_with_r2.png"), 
       width = 8, height = 6, units = "in", dpi = 300)

file_path <- file.path(out_dir, "basic_predictions_feb20.csv")
# write.csv(pred_df_basic, file_path, row.names = FALSE)

#############################################################################################

### Meta data set
numvar_meta <- c()
lambdavec <- seq(from = 10, to = 50, by = 1)
for (lambdy in lambdavec) {
  predictors <- setdiff(names(train_meta), c("BMI", "subject_id"))  
  predictors_escaped <- paste0("`", predictors, "`", collapse = " + ")
  fix_formula <- as.formula(paste("BMI ~", predictors_escaped))
  lm1_meta <- glmmLasso(fix = fix_formula,
                   data = train_meta,
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
  geom_vline(xintercept = 25, color = "blue", linetype = "dashed") +
  theme_bw() +
  ggtitle("Only Meta Delta") +
  xlab("Penalty Coefficient (Lambda)") +
  ylab("Number of Included Variables")
ggsave(paste0(out_dir,"lambda_elbow_only_meta.png"), width=6, height=4, units="in", dpi=320)

best_meta <- glmmLasso(fix = fix_formula,
                       data = train_meta,
                       rnd = list(subject_id = ~ 1, range = ~ 1),
                       lambda = 30,
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

bestglm_meta <- as.data.frame(train_meta[,lassoFeatures])
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

# grab only the columns that we have coefs for (starting at index 2 removes the intercept)
test_meta_vars <- test_meta %>% dplyr::select(rownames(mod_coef_df)[2:nrow(mod_coef_df)])
test_meta_vars$Intercept <- 1 # make a column for the intercept
test_meta_vars <- test_meta_vars %>% dplyr::select(Intercept, everything())

#--predict the risk scores without covariates--########
pred_risk_scores_meta <- as.matrix(test_meta_vars) %*% as.matrix(mod_coef_df$Estimate)

# combine the predicted and actual
pred_df_meta <- as.data.frame(cbind(test_meta$subject_id,
                               test_meta$range,
                               test_meta$BMI,
                               scale(pred_risk_scores_meta))); colnames(pred_df_meta) <- c("subject_id",
                                                                                 "time", 
                                                                                 "actual",
                                                                                 "predicted_meta")
# Calculate R-squared value
actual <- test_meta$BMI
predicted <- pred_risk_scores_meta
mean_actual <- mean(actual)
ss_total <- sum((actual - mean_actual)^2)
ss_residual <- sum((actual - predicted)^2)
r_squared <- 1 - (ss_residual / ss_total)

# Print R-squared value
r_squared

# Extract the coefficients (excluding intercept) and calculate absolute values for feature importance
coef_df <- as.data.frame(best_meta$coefficients)
#coef_df <- as.data.frame(coef_df[rownames(coef_df) != "(Intercept)", ]) # Remove intercept
coef_df$Feature <- rownames(coef_df) 
coef_df$Feature <- gsub("range3.99098082279631", "time", coef_df$Feature)
coef_df$Importance <- coef_df$`best_meta$coefficients`

# Filter out rows where Importance is 0
coef_df <- coef_df[coef_df$Importance != 0, ]

# Plot the feature importances
library(ggplot2)
ggplot(coef_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Meta Feature Importances from GLM Lasso Model") +
  xlab("Features") +
  ylab("Importance (Absolute Coefficients)") +
  # Add extra space around the plot using plot.margin
  #theme(plot.margin = margin(t = 20, r = 20, b = 20, l = 40)) +  # top, right, bottom, left margins
  annotate("text", 
           x = 1, 
           y = max(coef_df$Importance) * 0.95, 
           label = paste("R� =", round(r_squared, 3)), 
           size = 5, 
           color = "red", 
           hjust = 0.9)  # hjust = 0 will align the label to the left

# Save the plot to a file
ggsave(paste0(out_dir, "meta_feature_importances_with_r2.png"), 
       width = 8, height = 6, units = "in", dpi = 300)

file_path <- file.path(out_dir, "meta_predictions_feb20.csv")
#write.csv(pred_df_meta, file_path, row.names = FALSE)

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
numvar_taxa <- c()
lambdavec <- seq(from = 32, to = 50, by = 1)
for (lambdy in lambdavec) {
  predictors <- setdiff(names(train_taxa_cleaned), c("BMI", "subject_id", "range"))
  predictors_escaped <- paste0("`", predictors, "`", collapse = " + ")
  fix_formula <- as.formula(paste("BMI ~", predictors_escaped))
  lm1_taxa <- glmmLasso(fix = fix_formula,
                      data = train_taxa_cleaned,
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
  geom_vline(xintercept = 43, color = "blue", linetype = "dashed") +
  theme_bw() +
  ggtitle("Only Taxa Delta") +
  xlab("Penalty Coefficient (Lambda)") +
  ylab("Number of Included Variables")
ggsave(paste0(out_dir,"lambda_elbow_only_taxa.png"), width=6, height=4, units="in", dpi=320)

best_taxa <- glmmLasso(fix = fix_formula,
                       data = train_taxa_cleaned,
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
bestglm_taxa <- as.data.frame(train_taxa_cleaned[,lassoFeatures])
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

# grab only the columns that we have coefs for (starting at index 2 removes the intercept)
rownames(mod_coef_df) <- gsub("`", "", rownames(mod_coef_df))
test_taxa_vars <- test_taxa %>% dplyr::select(rownames(mod_coef_df)[2:nrow(mod_coef_df)])
test_taxa_vars$Intercept <- 1 # make a column for the intercept
test_taxa_vars <- test_taxa_vars %>% dplyr::select(Intercept, everything())

#--predict the risk scores without covariates--########
pred_risk_scores_taxa <- as.matrix(test_taxa_vars) %*% as.matrix(mod_coef_df$Estimate)

# combine the predicted and actual
pred_df_taxa <- as.data.frame(cbind(test_taxa$subject_id,
                                    test_taxa$range,
                                    test_taxa$BMI,
                                    scale(pred_risk_scores_taxa))); colnames(pred_df_taxa) <- c("subject_id",
                                                                                      "time", 
                                                                                      "actual",
                                                                                      "predicted_taxa")
# Calculate R-squared value
actual <- test_taxa$BMI
predicted <- pred_risk_scores_taxa
mean_actual <- mean(actual)
ss_regression <- sum((predicted - mean_actual)^2)

ss_total <- sum((actual - mean_actual)^2)
ss_regression/ss_total

ss_residual <- sum((actual - predicted)^2)
r_squared <- 1 - (ss_residual / ss_total)

# OTHER FORMULAS:
rsq <- function (x, y) cor(x, y) ^ 2
rsq(actual, predicted)

# Print R-squared value
r_squared

# Extract the coefficients (excluding intercept) and calculate absolute values for feature importance
coef_df <- as.data.frame(best_taxa$coefficients)
#coef_df <- as.data.frame(coef_df[rownames(coef_df) != "(Intercept)", ]) # Remove intercept
coef_df$Feature <- rownames(coef_df) 
coef_df$Feature <- gsub("range3.99098082279631", "time", coef_df$Feature)
coef_df$Importance <- coef_df$`best_taxa$coefficients`

# Filter out rows where Importance is 0
coef_df <- coef_df[coef_df$Importance != 0, ]

# Plot the feature importances
ggplot(coef_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Taxa Feature Importances from GLM Lasso Model") +
  xlab("Features") +
  ylab("Importance (Absolute Coefficients)") +
  # Add extra space around the plot using plot.margin
  #theme(plot.margin = margin(20, 20, 20, 40)) +  # top, right, bottom, left margins
  annotate("text", 
           x = 1, 
           y = max(coef_df$Importance) * 0.95, 
           label = paste("R� =", round(r_squared, 3)), 
           size = 5, 
           color = "red", 
           hjust = 0.9)  # hjust = 0 will align the label to the left

# Save the plot to a file
ggsave(paste0(out_dir, "taxa_feature_importances_with_r2.png"), 
       width = 8, height = 6, units = "in", dpi = 300)

file_path <- file.path(out_dir, "taxa_predictions_feb20.csv")
#write.csv(pred_df_taxa, file_path, row.names = FALSE)

##########################################################################################
### Micom data set
train_micom_cleaned <- remove_high_corr_vars(train_micom, threshold = 0.9)
numvar_micom <- c()
lambdavec <- seq(from = 10, to = 50, by = 1)
for (lambdy in lambdavec) {
  predictors <- setdiff(names(train_micom_cleaned), c("BMI", "subject_id", "range"))
  predictors_escaped <- paste0("`", predictors, "`", collapse = " + ")
  fix_formula <- as.formula(paste("BMI ~", predictors_escaped))
  lm1_micom <- glmmLasso(fix = fix_formula,
                        data = train_micom_cleaned,
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
ggsave(paste0(out_dir,"lambda_elbow_only_micom.png"), width=6, height=4, units="in", dpi=320)

best_micom <- glmmLasso(fix = fix_formula,
                        data = train_micom_cleaned,
                        rnd = list(subject_id = ~ 1, range = ~ 1),
                        lambda = 30,
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
bestglm_micom <- as.data.frame(train_micom_cleaned[,lassoFeatures])
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

# grab only the columns that we have coefs for (starting at index 2 removes the intercept)
rownames(mod_coef_df) <- gsub("`", "", rownames(mod_coef_df))
test_micom_vars <- test_micom %>% dplyr::select(rownames(mod_coef_df)[2:nrow(mod_coef_df)])
test_micom_vars$Intercept <- 1 # make a column for the intercept
test_micom_vars <- test_micom_vars %>% dplyr::select(Intercept, everything())

#--predict the risk scores without covariates--########
pred_risk_scores_micom <- as.matrix(test_micom_vars) %*% as.matrix(mod_coef_df$Estimate)

# combine the predicted and actual
pred_df_micom<- as.data.frame(cbind(test_micom$subject_id,
                                    test_micom$range,
                                    test_micom$BMI,
                                    scale(pred_risk_scores_micom))); colnames(pred_df_micom) <- c("subject_id",
                                                                                           "time", 
                                                                                           "actual",
                                                                                           "predicted_micom")
# Calculate R-squared value
actual <- test_micom$BMI
predicted <- pred_risk_scores_micom
mean_actual <- mean(actual)
ss_total <- sum((actual - mean_actual)^2)
ss_residual <- sum((actual - predicted)^2)
r_squared <- 1 - (ss_residual / ss_total)

# Print R-squared value
r_squared

# Extract the coefficients (excluding intercept) and calculate absolute values for feature importance
coef_df <- as.data.frame(best_micom$coefficients)
#coef_df <- as.data.frame(coef_df[rownames(coef_df) != "(Intercept)", ]) # Remove intercept
coef_df$Feature <- rownames(coef_df) 
coef_df$Feature <- gsub("range3.99098082279631", "time", coef_df$Feature)
coef_df$Importance <- coef_df$`best_micom$coefficients`

# Filter out rows where Importance is 0
coef_df <- coef_df[coef_df$Importance != 0, ]

# Plot the feature importances
ggplot(coef_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Micom Feature Importances from GLM Lasso Model") +
  xlab("Features") +
  ylab("Importance (Absolute Coefficients)") +
  # Add extra space around the plot using plot.margin
  #theme(plot.margin = margin(20, 20, 20, 40)) +  # top, right, bottom, left margins
  annotate("text", 
           x = 1, 
           y = max(coef_df$Importance) * 0.95, 
           label = paste("R� =", round(r_squared, 3)), 
           size = 5, 
           color = "red", 
           hjust = 0.9)  # hjust = 0 will align the label to the left

# Save the plot to a file
ggsave(paste0(out_dir, "micom_feature_importances_with_r2.png"), 
       width = 8, height = 6, units = "in", dpi = 300)

file_path <- file.path(out_dir, "micom_predictions_feb20.csv")
#write.csv(pred_df_micom, file_path, row.names = FALSE)

##########################################################################################
### Pathway data set
train_pathway_cleaned <- remove_high_corr_vars(train_pathway, threshold = 0.9)
numvar_path <- c()
lambdavec <- seq(from = 10, to = 50, by = 1)
for (lambdy in lambdavec) {
  predictors <- setdiff(names(train_pathway_cleaned), c("BMI", "subject_id", "range"))
  predictors_escaped <- paste0("`", predictors, "`", collapse = " + ")
  fix_formula <- as.formula(paste("BMI ~", predictors_escaped))
  lm1_pathway <- glmmLasso(fix = fix_formula,
                         data = train_pathway_cleaned,
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
ggsave(paste0(out_dir,"lambda_elbow_only_pathway.png"), width=6, height=4, units="in", dpi=320)

best_path <- glmmLasso(fix = fix_formula,
                       data = train_pathway_cleaned,
                       rnd = list(subject_id = ~ 1, range = ~ 1),
                       lambda = 26,
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
bestglm_path <- as.data.frame(train_pathway_cleaned[,lassoFeatures])
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

# grab only the columns that we have coefs for (starting at index 2 removes the intercept)
rownames(mod_coef_df) <- gsub("`", "", rownames(mod_coef_df))
test_path_vars <- test_path %>% dplyr::select(rownames(mod_coef_df)[2:nrow(mod_coef_df)])
test_path_vars$Intercept <- 1 # make a column for the intercept
test_path_vars <- test_path_vars %>% dplyr::select(Intercept, everything())

#--predict the risk scores without covariates--########
pred_risk_scores_path <- as.matrix(test_path_vars) %*% as.matrix(mod_coef_df$Estimate)

# combine the predicted and actual
pred_df_path<- as.data.frame(cbind(test_path$subject_id,
                                    test_path$range,
                                    test_path$BMI,
                                    scale(pred_risk_scores_path))); colnames(pred_df_path) <- c("subject_id",
                                                                                            "time", 
                                                                                            "actual",
                                                                                            "predicted_pathway")


# Calculate R-squared value
actual <- test_path$BMI
predicted <- pred_risk_scores_path
mean_actual <- mean(actual)
ss_total <- sum((actual - mean_actual)^2)
ss_residual <- sum((actual - predicted)^2)
r_squared <- 1 - (ss_residual / ss_total)

# Print R-squared value
r_squared

# Extract the coefficients (excluding intercept) and calculate absolute values for feature importance
coef_df <- as.data.frame(best_path$coefficients)
#coef_df <- as.data.frame(coef_df[rownames(coef_df) != "(Intercept)", ]) # Remove intercept
coef_df$Feature <- rownames(coef_df) 
coef_df$Feature <- gsub("range3.99098082279631", "time", coef_df$Feature)
coef_df$Importance <- coef_df$`best_path$coefficients`

# Filter out rows where Importance is 0
coef_df <- coef_df[coef_df$Importance != 0, ]

# Plot the feature importances
ggplot(coef_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Micom Feature Importances from GLM Lasso Model") +
  xlab("Features") +
  ylab("Importance (Absolute Coefficients)") +
  # Add extra space around the plot using plot.margin
  #theme(plot.margin = margin(20, 20, 20, 40)) +  # top, right, bottom, left margins
  annotate("text", 
           x = 1, 
           y = max(coef_df$Importance) * 0.95, 
           label = paste("R� =", round(r_squared, 3)), 
           size = 5, 
           color = "red", 
           hjust = 0.9)  # hjust = 0 will align the label to the left
# Save the plot to a file
ggsave(paste0(out_dir, "path_feature_importances_with_r2.png"), 
       width = 8, height = 6, units = "in", dpi = 300)

file_path <- file.path(out_dir, "path_predictions_feb20.csv")
#write.csv(pred_df_path, file_path, row.names = FALSE)

### Basic model ###
numvar_meta <- c()
lambdavec <- seq(from = 10, to = 50, by = 1)
for (lambdy in lambdavec) {
  predictors <- c("age", "sex")
  predictors_escaped <- paste0("`", predictors, "`", collapse = " + ")
  fix_formula <- as.formula(paste("BMI ~", predictors_escaped))
  lm1_meta <- glmmLasso(fix = fix_formula,
                        data = train_meta,
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
  geom_vline(xintercept = 25, color = "blue", linetype = "dashed") +
  theme_bw() +
  ggtitle("Only Meta Delta") +
  xlab("Penalty Coefficient (Lambda)") +
  ylab("Number of Included Variables")
ggsave(paste0(out_dir,"lambda_elbow_only_basic.png"), width=6, height=4, units="in", dpi=320)

best_meta <- glmmLasso(fix = fix_formula,
                       data = train_meta,
                       rnd = list(subject_id = ~ 1, range = ~ 1),
                       lambda = 30,
                       family = gaussian(link = "identity"))

lassoFeatures <- names(best_meta$coefficients)
lassoFeatures <- lassoFeatures[lassoFeatures %ni% c("(Intercept)")]
lassoFeatures <- lassoFeatures[grep("as.factor",lassoFeatures,invert=T)] ####
lassoFeatures <- unique(c(lassoFeatures))
lassoFeatures <- gsub("range0.997860960117438", "range", lassoFeatures)
lassoFeatures <- lassoFeatures[grep("subject_id|BMI|range", lassoFeatures, invert = T)]
lassoFeatures <- unique(c(lassoFeatures, "subject_id", "BMI", "range"))

# Check the updated vector
lassoFeatures

bestglm_meta <- as.data.frame(train_meta[,lassoFeatures])
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
test_basic <- test_meta[complete.cases(test_meta),] # complete cases

# grab only the columns that we have coefs for (starting at index 2 removes the intercept)
test_meta_vars <- test_basic %>% dplyr::select(rownames(mod_coef_df)[2:nrow(mod_coef_df)])
test_meta_vars$Intercept <- 1 # make a column for the intercept
test_meta_vars <- test_meta_vars %>% dplyr::select(Intercept, everything())

#--predict the risk scores without covariates--########
pred_risk_scores_basic <- as.matrix(test_meta_vars) %*% as.matrix(mod_coef_df$Estimate)

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
predicted <- pred_risk_scores_path
mean_actual <- mean(actual)
ss_total <- sum((actual - mean_actual)^2)
ss_residual <- sum((actual - predicted)^2)
r_squared <- 1 - (ss_residual / ss_total)

# Print R-squared value
r_squared

# Extract the coefficients (excluding intercept) and calculate absolute values for feature importance
coef_df <- as.data.frame(best_basic$coefficients)
#coef_df <- as.data.frame(coef_df[rownames(coef_df) != "(Intercept)", ]) # Remove intercept
coef_df$Feature <- rownames(coef_df) 
coef_df$Feature <- gsub("range3.99098082279631", "time", coef_df$Feature)
coef_df$Importance <- coef_df$`best_basic$coefficients`

# Plot the feature importances
ggplot(coef_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Basic Feature Importances from GLM Lasso Model") +
  xlab("Features") +
  ylab("Importance (Absolute Coefficients)") +
  # Add extra space around the plot using plot.margin
  #theme(plot.margin = margin(20, 20, 20, 40)) +  # top, right, bottom, left margins
  annotate("text", 
           x = 1, 
           y = max(coef_df$Importance) * 0.95, 
           label = paste("R� =", round(r_squared, 3)), 
           size = 5, 
           color = "red", 
           hjust = 0.9)  # hjust = 0 will align the label to the left

# Save the plot to a file
ggsave(paste0(out_dir, "basic_feature_importances_with_r2.png"), 
       width = 8, height = 6, units = "in", dpi = 300)

file_path <- file.path(out_dir, "basic_predictions_feb20.csv")
#write.csv(pred_df_basic, file_path, row.names = FALSE)

################################################################################

##### Compare models

# Convert the last 3 columns to numeric
pred_df_basic[, tail(names(pred_df_basic), 3)] <- 
  lapply(pred_df_basic[, tail(names(pred_df_basic), 3)], as.numeric)
pred_df_meta[, tail(names(pred_df_meta), 3)] <- 
  lapply(pred_df_meta[, tail(names(pred_df_meta), 3)], as.numeric)
pred_df_taxa[, tail(names(pred_df_taxa), 3)] <- 
  lapply(pred_df_taxa[, tail(names(pred_df_taxa), 3)], as.numeric)
pred_df_micom[, tail(names(pred_df_micom), 3)] <- 
  lapply(pred_df_micom[, tail(names(pred_df_micom), 3)], as.numeric)
pred_df_path[, tail(names(pred_df_path), 3)] <- 
  lapply(pred_df_path[, tail(names(pred_df_path), 3)], as.numeric)

met_basic <- merge(pred_df_meta, pred_df_basic, 
                 by = c("subject_id", "time")) %>% 
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
all_omic$predicted_taxa <- as.numeric(all_omic$predicted_taxa)
all_omic$predicted_micom <- as.numeric(all_omic$predicted_micom)
all_omic$predicted_pathway <- as.numeric(all_omic$predicted_pathway)
all_omic[, 3:8] <- scale(all_omic[, 3:8])

########################################################################################
mod_dat = all_omic %>% rename(bmi = actual, 
                              Time = time,
                              Cluster = subject_id,
                              y_new_meta_only = predicted_meta,
                              y_new_basic = predicted_basic,
                              y_new_micom_only = predicted_micom,
                              y_new_path_only = predicted_pathway,
                              y_new_tax_only = predicted_taxa)

mod_dat$Time <- as.numeric(mod_dat$Time)
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
                  title = "glmlasso delta sequential lmer models",
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
########################################################################################
# Make linear models 
basic_model <- lme(actual ~ predicted_basic + time, 
                   random = ~1|subject_id, data = all_omic)


# Combined models 
lm_meta_tax_time <- lme(actual ~ predicted_meta + predicted_taxa + time, 
                   random = ~1|subject_id, data = all_omic)
summary(lm_meta_tax_time)

### Checks before models 
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
lm_basic <- lme(actual ~ predicted_basic, random = ~1|subject_id, data = all_omic)

model_titles <- c("1: Basic",
                  "2: Meta", 
                  "3: Taxa", 
                  "4: Micom", 
                  "5: Pathway")

# Create the table as a grid object
sjPlot::tab_model(lm_basic, lm_meta, lm_tax, 
                  lm_micom, lm_path, 
                  title = "Comparing single omic glmLASSO models deltas",
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
lm_meta_basic <- lme(actual ~ predicted_basic + predicted_meta, 
                   random = ~1|subject_id, data = all_omic)
lm_meta_tax <- lme(actual ~ predicted_meta + predicted_taxa, 
                                random = ~1|subject_id, data = all_omic)
lm_basic_tax <- lme(actual ~ predicted_basic + predicted_taxa, 
                   random = ~1|subject_id, data = all_omic)
lm_basic_micom <- lme(actual ~ predicted_basic + predicted_micom, 
                    random = ~1|subject_id, data = all_omic)
lm_basic_path <- lme(actual ~ predicted_basic + predicted_pathway, 
                      random = ~1|subject_id, data = all_omic)

lm_meta_tax <- lme(actual ~ predicted_meta + predicted_taxa, 
                   random = ~1|subject_id, data = all_omic)
lm_basic_tax <- lme(actual ~ predicted_basic + predicted_taxa, 
                    random = ~1|subject_id, data = all_omic)
lm_basic_micom <- lme(actual ~ predicted_basic + predicted_micom, 
                      random = ~1|subject_id, data = all_omic)
lm_basic_path <- lme(actual ~ predicted_basic + predicted_pathway, 
                     random = ~1|subject_id, data = all_omic)

lm_basic_all_omic <- lme(actual ~ predicted_basic + predicted_meta + 
                         predicted_taxa + predicted_pathway + predicted_micom, 
                         random = ~1|subject_id, data = all_omic)
# No basic 
lm_meta_tax_micom <- lme(actual ~ predicted_meta + predicted_taxa + 
                           predicted_micom, random = ~1|subject_id, data = all_omic)
lm_meta_tax_micom_path <- lme(actual ~ predicted_meta + predicted_taxa + 
                                predicted_micom + predicted_pathway, 
                              random = ~1|subject_id, data = all_omic)
lm_no_meta <- lme(actual ~ predicted_taxa + 
                    predicted_micom + predicted_pathway, 
                  random = ~1|subject_id, data = all_omic)
summary(lm_meta_tax_micom_path)



# Fit both models using Maximum Likelihood (ML)
lm_basic_ml <- lme(fixed = actual ~ predicted_basic,  # Fixed effects
                   random = ~ 1 | subject_id,  # Random effects
                   method = "ML",  # Use Maximum Likelihood
                   data = all_omic)

lm_basic_path_ml <- lme(fixed = actual ~ predicted_basic + predicted_taxa,  # Different fixed effects
                        random = ~ 1 | subject_id,  # Random effects
                        method = "ML",  # Use Maximum Likelihood
                        data = all_omic)

# Now perform ANOVA (ML comparison)
anova(lm_basic_ml, lm_basic_path_ml)


# Define model titles for clarity
model_titles <- c("Model 0: Basic" ,
                  "Model 1: Basic + Meta", 
                  "Model 2: Meta + Taxa", 
                  "Model 3: Meta + Taxa + Micom", 
                  "Model 4: Taxa + Micom + Pathway", 
                  "Model 5: Meta + Taxa + Micom + Pathway")

sjPlot::tab_model(lm_basic, lm_meta_basic, lm_meta_tax, 
                  lm_meta_tax_micom, lm_no_meta, lm_meta_tax_micom_path, 
                  title = "Comparing combined omic delta models",
                  string.pred = "Predictors",
                  string.est = "Estimates",
                  string.std = "std. Beta",
                  string.ci = "95% CI",
                  string.se = "std. Error",
                  p.style = c("numeric"), 
                  p.threshold = c(0.05),
                  dv.labels = model_titles,
                  auto.label = FALSE)

# Save the plot
#output_file <- file.path(out_dir, "delta_combined_omics_models_feb20.png")
#ggsave(output_file, plot_model, width = 10, height = 6)

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
  labs(title = "Correlations btwn HIGH corr pathway2", fill = "Correlation")

# Run the ANOVA tests
anova_basic_meta <- anova(lm_basic, lm_meta_basic)
anova_basic_tax <- anova(lm_basic, lm_basic_tax)
anova_basic_micom <- anova(lm_basic, lm_basic_micom)
anova_basic_path <- anova(lm_basic, lm_basic_path)
anova_basic_all <- anova(lm_basic, lm_basic_all_omic)

meta_r <- rsq.lmm(lm_meta,adj=TRUE)
basic_r <- rsq.lmm(lm_basic,adj=TRUE)
meta_basic_r <- rsq.lmm(lm_meta_basic,adj=TRUE)
tax_basic_r <- rsq.lmm(lm_basic_tax,adj=TRUE)
micom_basic_r <- rsq.lmm(lm_basic_micom,adj=TRUE)
path_basic_r <- rsq.lmm(lm_basic_path,adj=TRUE)
all_r <- rsq.lmm(lm_basic_all_omic,adj=TRUE)
  
### Plot r-squared 
# Create a data frame with model names and their corresponding R-squared values
r_squared_data <- data.frame(
  Model = rep(c("lm_basic", "lm_basic_meta", "lm_basic_tax", 
                "lm_basic_micom", "lm_basic_path", "basic_all_omic"), each = 3),
  Type = rep(c("Model_Total", "Fixed_Effects", "Random_Effects"), times = 6),
  R_squared = c(path_basic_r$model, path_basic_r$fixed, path_basic_r$random,
                meta_basic_r$model, meta_basic_r$fixed, meta_basic_r$random,
                tax_basic_r$model, tax_basic_r$fixed, tax_basic_r$random,
                micom_basic_r$model, micom_basic_r$fixed, micom_basic_r$random,
                basic_r$model, basic_r$fixed, basic_r$random,
                all_r$model, all_r$fixed, all_r$random))

# Plot the data using ggplot2
r2_plot <- ggplot(r_squared_data, aes(x = Model, y = R_squared, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(x = "Model", y = "R Squared", 
       title = "R Squared for Different Models (glmlasso delta)") +
  scale_fill_manual(values = c("lightblue", "lightgreen", "lightcoral")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 0.6)

output_file <- file.path(out_dir, "glm_delta_r2_models_feb20.png")
ggsave(output_file, r2_plot, width = 10, height = 6)

### PLOT ANOVA RESULTS
anova_results <- data.frame(
  Model_Comparison = c("basic vs meta + basic", "basic vs basic + tax", 
                       "basic vs basic + micom", "basic vs basic + path",
                       "basic vs basic + all omic"),
  L_Ratio = c(anova_basic_meta$L.Ratio[2], anova_basic_tax$L.Ratio[2], 
              anova_basic_micom$L.Ratio[2], anova_basic_path$L.Ratio[2],
              anova_basic_all$L.Ratio[2]),
  p_value = c(anova_basic_meta$`p-value`[2], anova_basic_tax$`p-value`[2], 
              anova_basic_micom$`p-value`[2], anova_basic_path$`p-value`[2],
              anova_basic_all$`p-value`[2]))

# Plot ANOVA
anova_plot <- ggplot(anova_results, aes(x = Model_Comparison, y = L_Ratio)) +
  geom_bar(stat = "identity", fill = "skyblue") +  # Create bars
  geom_text(aes(label = round(p_value, 3)), vjust = -0.5) +  # Add p_value above the bars
  theme_minimal() +  # Clean theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  labs(x = "Model Comparison", y = "L Ratio", 
       title = "L Ratio for Each Model Comparison (glmlasso delta models)") +
  ylim(0, 7)

output_file <- file.path(out_dir, "glm_delta_anova_models_feb20.png")
ggsave(output_file, anova_plot, width = 10, height = 6)
