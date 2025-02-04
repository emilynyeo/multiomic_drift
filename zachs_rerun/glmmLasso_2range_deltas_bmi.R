# Trying out long lasso
rm(list = ls())
source("zc_functions.R") 
library(pacman)
p_load(tools, reticulate, viridis, tidyplots, patchwork, jsonlite, maps, ggvenn, 
       caret, caretEnsemble, glmnet, xgboost, ggplot2, glmmLasso, corrplot,
       readr, plyr, dplyr, tidyr, purrr, tibble, stringr, psych, randomForest,  
       reshape2, scales, gridExtra, plotly, sf, tidyverse, naniar, VIM)
'%ni%' <- Negate('%in%')

data_dir <- "drift_fs/csv/all_omic_processed_data/deltas/"
train <- read_csv(paste0(data_dir, "jan30_all_delta_train_imp_varcheck.csv"))
test <- read_csv(paste0(data_dir, "jan30_all_delta_test_imp_varcheck.csv"))

### LASSO Lambda Search ###
train_bmi <- train %>% select(-c("Weight"))
names(train_bmi) <- sapply(names(train_bmi), 
                           function(x) gsub("[[:space:].-]", "_", x))
names(train_bmi)
# options(expressions = 50000) 
train_bmi$subject_id <- as.factor(train_bmi$subject_id)
train_bmi$range <- as.factor(train_bmi$range)
#train_bmi$sex <- as.factor(train_bmi$sex)
#train_bmi$randomized_group <- as.factor(train_bmi$randomized_group)
#train_bmi$race <- as.factor(train_bmi$race)

### Check Correlations
cor_matrix <- cor((train_bmi)[2:918])
highly_correlated <- findCorrelation(cor_matrix, cutoff = 0.9)
train_bmi_clean <- train_bmi[, -highly_correlated]
train_bmi_correlated <- train_bmi[, highly_correlated]

cor_matrix_meta <- cor((train_bmi)[2:17], 
                  use = "pairwise.complete.obs", method = "pearson")
cor_matrix_path <- cor((train_bmi)[18:374], 
                       use = "pairwise.complete.obs", method = "pearson")
cor_matrix_path1 <- cor((train_bmi)[18:196], 
                       use = "pairwise.complete.obs", method = "pearson")
cor_matrix_path2 <- cor((train_bmi)[197:374], 
                        use = "pairwise.complete.obs", method = "pearson")
cor_matrix_tax <- cor((train_bmi)[375:827], 
                       use = "pairwise.complete.obs", method = "pearson")
cor_matrix_micom <- cor((train_bmi)[828:918], 
                      use = "pairwise.complete.obs", method = "pearson")

# Melt the correlation matrix into a long format (use only upper triangle)
melted_cor_matrix <- melt(cor_matrix_path2, na.rm = TRUE)
melted_cor_matrix$Var1 <- str_sub(melted_cor_matrix$Var1, 1, 12)
melted_cor_matrix$Var2 <- str_sub(melted_cor_matrix$Var2, 1, 12)

filtered_melted_cor_matrix <- melted_cor_matrix %>%
                              filter(value >= 0.95 & value < 0.99999999999999)

# Create the plot with trimmed labels and smaller axis text
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

### Actual ###
varlist_nocolin <- c(colnames(train_bmi[2:15]))
varlist_nocolin <- setdiff(colnames(train_bmi)[1:15], 
                           c("subject_id", "BMI", "range"))
varstring_nocolin <- paste0(varlist_nocolin, collapse = " + ", sep = "")

# Create the formula
nchar(varstring_nocolin)
options(expressions = 15000)
formula <- as.formula(paste0("BMI ~ ", varstring_nocolin))
numvariables <- c()
lambdavec <- seq(from = 10, to = 50, by = 1)
for(lambdy in lambdavec){
  lm1 <- glmmLasso(fix = as.formula(paste0("BMI ~ ", varstring_nocolin)),
                   data = train_bmi_clean,
                   rnd = list(subject_id = ~ 1), 
                              # range = ~ 1),  # time point variable
                   lambda=lambdy,
                   family = gaussian(link = "identity"))
}
  summary(lm1)
  lassoFeatures <- names(lm1$coefficients[which(lm1$coefficients != 0)])
  lassoFeatures <- lassoFeatures[lassoFeatures %ni% c("(Intercept)")]
  lassoFeatures <- lassoFeatures[grep("as.factor",lassoFeatures,invert=T)] ####
  lassoFeatures <- unique(c(lassoFeatures, "Participant_ID", "site_name", 
                            "diagnosis", "consent_age", "sex", "race", "Antibiotics"))
  numvariables <- c(numvariables, length(lassoFeatures))

plot(x = lambdavec, y = numvariables)

numvariables <- numvariables - 7; #numvariables[which(numvariables < 0)] <- 0 ####

ggplot() +
  geom_point(aes(x = lambdavec, y = numvariables)) +
  geom_vline(xintercept = 25, color = "blue", linetype = "dashed") +
  theme_bw() +
  xlab("Penalty Coefficient (Lambda)") +
  ylab("Number of Included Variables")
ggsave("lambda_elbow.png", width=6, height=4, units="in", dpi=320)

lassoFeatures