# BL predicting BMI change 
source("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/zachs_rerun/zc_functions.R") 
library(dplyr)
library(tidyr)
library(tools)
library(broom)
library(stargazer)
library(sjPlot)
library(lme4)
library(broom)
library(knitr)

train_and_save_models <- function(data, target_var, 
                                  train_control, result_prefix, 
                                  test_size = 0.3) {
  # Split data into training and testing sets
  set.seed(123) # Ensure reproducibility
  train_indices <- sample(seq_len(nrow(data)), size = (1 - test_size) * nrow(data))
  train_data <- data[train_indices, ]
  test_data <- data[-train_indices, ]
  
  # Train models
  results <- train_all_models(train_data, target_var, train_control)
  
  # Combine importances
  feature_importance <- combine_importances(
    list(results$rf_model, results$lasso_model, 
         results$ridge_model, results$elastic_net_model, results$xgb_model),
    c("RF_Importance", "Lasso_Importance", 
      "Ridge_Importance", "Enet_Importance", "XGBoost_Importance")
  )
  # ensure output directory exists
  if (!dir.exists("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/april_aim1/")) {
    dir.create("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/april_aim1/", recursive = TRUE)
  }
  
  write.csv(feature_importance, paste0("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/april_aim1/", 
                                       result_prefix, "_feature_importance.csv"), 
            row.names = FALSE)
  
  # Extract best betas
  beta_coefficients <- extract_best_betas(
    list(results$lasso_model, results$ridge_model, results$elastic_net_model),
    c("Lasso_Beta", "Ridge_Beta", "Enet_Beta")
  )
  write.csv(beta_coefficients, paste0("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/april_aim1/", result_prefix, "_beta.csv"), row.names = FALSE)
  
  # Initialize an empty DataFrame to store performance metrics
  metrics_df <- data.frame(Model = character(), 
                           DataType = character(), 
                           R2 = numeric(), MAE = numeric(), 
                           RMSE = numeric(), stringsAsFactors = FALSE)
  
  # Calculate and store metrics for each model on both training and testing data
  for (model_name in names(results[-length(results)])) {
    model <- results[[model_name]]
    
    # Training metrics
    train_metrics <- calculate_metrics(model, train_data, target_var, model_name, "Train")
    metrics_df <- rbind(metrics_df, train_metrics)
    
    # Testing metrics
    test_metrics <- calculate_metrics(model, test_data, target_var, model_name, "Test")
    metrics_df <- rbind(metrics_df, test_metrics)
  }
  
  # Save the metrics DataFrame as CSV
  write.csv(metrics_df, paste0("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/april_aim1/", result_prefix, "_metrics.csv"), row.names = FALSE)
  
  # ensure the model directory exists
  if (!dir.exists("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/april_aim1/models/")) {
    dir.create("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/april_aim1/models/", recursive = TRUE)
  }
  # Save model results
  saveRDS(results, paste0("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/april_aim1/models/", result_prefix, "_results.rds"))
  
  return(results)
}

create_plots <- function(data_list, max_r2, titles) {
  custom_colors <- c("#960018", "#8B3A3A", "#996759", "#D8A39D", "#C74375", "#582F2B")
  
  plots <- lapply(seq_along(data_list), function(i) {
    ggplot(data_list[[i]], aes(x = R2, y = Model, fill = Model)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(
        title = titles[i],
        x = "R²",  # fixed the character
        y = "Model",
        fill = "Model"
      ) +
      coord_cartesian(xlim = c(0, max_r2)) +
      theme_minimal() +
      scale_fill_manual(values = custom_colors) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15, angle = 45),
        legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 20)
      )
  })
  
  # Combine the plots vertically (you might want to use patchwork or cowplot for this)
  Reduce(`/`, plots)
}
# Set the input directory
long_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/april_processing/"

# Read the data from the CSV file
long <- read.csv(file.path(long_dir, "long.csv")) %>% dplyr::filter(!time == 12) 

bmi <- long %>% dplyr::select(c(subject_id, time, outcome_BMI_fnl)) %>% 
  tidyr::pivot_wider(names_from = time,
    values_from = outcome_BMI_fnl,names_prefix = "BMI_time_") %>% 
  mutate(bmi_0_6 = BMI_time_6 - BMI_time_0) %>% 
  dplyr::select(-c(BMI_time_0, BMI_time_6))

long_bl <- long %>% dplyr::filter(time == 0)
BL <- merge(bmi, long_bl, by = "subject_id") %>%
      dplyr::filter(!is.na(bmi_0_6)) # 26 removed 

# Make train and test sets 
# test sample names
test_names <- c("ABR-079", "AGA-071", "AHE-055", "ALI-121", "ALO-163", "AMA-031", "ASO-013", "AWI-167", "BMO-164", "CWA-183", "DSC-024", "EBE-130", "EHI-177", "EJO-092", "GFU-188", "HGI-010", "JCA-109", "JGO-100","KBU-085", "KCE-034", "KHE-170", "LDO-148", "LST-186", "LZD-142", "MAR-119", "MCA-088", "MJA-153", "MWE-112", "NPO-149", "RAE-114", "SBO-020", "SEG-080", "SKA-195", "SLO-178", "SSH-028", "TDU-086", "TFA-016", "VCA-041")

# train sample names
train_names <- c("AAL-144", "ACO-053", "ADA-105", "AKE-009", "AKI-011", "AKO-139", "AMC-155", "AME-128",  "AME-157", "ATA-129", "AWA-052", "AWA-083", "BAN-193", "BHO-014", "BIN-201", "BKN-104", "BMI-156", "BSA-174", "CAM-057", "CCO-189", "CED-026", "CEL-073", "CGA-134", "CIS-077", "CKR-078", "CLE-049", "COW-066", "CRO-108", "CWA-161", "EBE-051", "EKA-135", "EKR-045", "ELA-159", "EPO-182", "EVO-184", "FWI-098", "GHA-035", "HDE-154", "IBE-120", "JDI-140", "JER-110", "JFU-027", "JJO-093", "JKN-127", "JPO-022", "JUG-116", "JUT-032", "JVE-126", "KAN-138", "KBR-162", "KEL-185", "KEL-199", "KGI-029", "KHU-196", "KPA-042", "KRI-072", "KVA-038", "KWA-122", "KWA-141", "LBL-047", "LBU-015", "LEL-147", "LFI-003", "LJA-101",  "LMC-111", "LPF-198", "LVA-017", "MBA-187", "MCW-065", "MDI-107", "MES-068", "MFB-118", "MGA-076", "MHO-117", "MKE-192", "MMA-036", "MRT-179", "MSH-091", "MST-039", "MWE-143",  "MWO-133", "MWY-152", "NAR-099", "NBI-048", "NBI-069", "NCO-171", "NDI-067", "NEL-094", "NKA-090", "NMO-151", "NTA-021", "PBE-123", "QNG-166", "RAF-125", "RAM-050", "RHP-023",  "RLA-132", "ROL-006", "SAB-160", "SCA-043", "SCR-061", "SDA-150", "SGA-062", "SKA-087", "SRO-194", "TBU-115", "TFA-172", "TRO-113", "TSH-146", "TSL-056", "WPE-005", "YOR-103",  "YSU-097", "ZVU-096")

cat("Length of test names:", length(test_names), "\n")
cat("Length of train names:", length(train_names), "\n")

# Create train and test sets for the current subset
#train_set <- BL[BL[["subject_id"]] %in% train_names, ]  # Select rows where ID_VAR is in train_names
#test_set <- BL[BL[["subject_id"]] %in% test_names, ]    # Select rows where ID_VAR is in test_names

# Define the column names based on your lists
basic <- c('subject_id', 'bmi_0_6','outcome_BMI_fnl', 'time','age', 'sex')
meta_keep <- c('subject_id', 'bmi_0_6','outcome_BMI_fnl', 'time', 'randomized_group', 
               'cohort_number', 'sex', 'race', 'age', 'Hemoglobin_A1C', 'HDL_Total_Direct_lipid', 
               'HOMA_IR', 'Insulin_endo', 'LDL_Calculated', 'Glucose.x')
only_grs <- c('subject_id', 'bmi_0_6','outcome_BMI_fnl', 'time','bmi_prs')
only_taxa <- c('subject_id', 'bmi_0_6', 'outcome_BMI_fnl', 'time', grep("^g__", names(long), value = TRUE))

micom_start <- which(names(long) == "Diacetyl")
micom_end <- which(names(long) == "aldehydo.D.xylose")
only_micom <- c('subject_id', 'bmi_0_6','outcome_BMI_fnl', 'time', names(long)[micom_start:micom_end])

path_start <- which(names(long) == "arginine..ornithine.and.proline.interconversion")
path_end <- which(names(long) == "UDP.N.acetyl.D.glucosamine.biosynthesis.I")
only_pathway <- c('subject_id', 'bmi_0_6','outcome_BMI_fnl', 'time', names(long)[path_start:path_end])

tabo_start <- which(names(long) == "non_HDL_C")
tabo_end <- which(names(long) == "IDL_TG_pct")
only_tabo <- c('subject_id', 'bmi_0_6', 'outcome_BMI_fnl', 'time', names(long)[tabo_start:tabo_end])

# Create data frames based on the columns defined
basic <- BL[, basic, drop = FALSE] %>% unique()
meta <- BL[, meta_keep, drop = FALSE] %>% unique()
grs <- BL[, only_grs, drop = FALSE] %>% unique()
taxa <- BL[, only_taxa, drop = FALSE] %>% unique()
micom <- BL[, only_micom, drop = FALSE] %>% unique()
pathway <- BL[, only_pathway, drop = FALSE] %>% unique()
metabo <- BL[, only_tabo, drop = FALSE] %>% unique()

# Make test and tain sets for each omic 
data_frames <- c("basic", "meta", "grs", "micom", "pathway", "taxa", "metabo")
models <- list()
for (df in data_frames) {
  df_data <- get(df)  # Get the data frame by name
  # Split into training and test sets
  train_set <- df_data[df_data$subject_id %in% train_names, ]
  test_set <- df_data[df_data$subject_id %in% test_names, ]
  assign(paste0(df, "_train"), train_set)
  assign(paste0(df, "_test"), test_set)
  
  set.seed(123)
  train_control <- trainControl(method = "cv", number = 5, search = "grid")
  
  m0_6_results <- train_and_save_models(
    train_set, "bmi_0_6", train_control,
    paste0("BL_pred_0_6m_bmi_", df))
  
  models[[df]] <- m0_6_results # store or print
}

# Read in data 
base_path <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/april_aim1"
df_names = list(basic = basic, meta = meta, grs = grs, taxa = taxa, 
                micom = micom, pathway = pathway, metabo = metabo)

made_dfs <- list()
for (df_name in names(df_names)) {
  df <- df_names[[df_name]] 
  
  file_paths <- list(
    beta = paste0("BL_pred_0_6m_bmi_", df_name, "_beta.csv"),
    feature_importance = paste0("BL_pred_0_6m_bmi_", df_name,"_feature_importance.csv"),
    metrics = paste0("BL_pred_0_6m_bmi_", df_name, "_metrics.csv"))
  
  data_list <- lapply(file_paths, function(path) read.csv(file.path(base_path, path)))
  made_dfs[[df_name]] <- data_list # store or print
  
  top_20_features <- get_top_n_features_all_models(data_list$feature_importance, 20)
  print(top_20_features)
  
  # Extract metrics and max RÃƒâ€šÃ‚Â² for all datasets
  results_all <- extract_metrics(made_dfs[[df_name]])
  
  # Calculate overall max RÃƒâ€šÃ‚Â² for each main category
  max_r2 <- results_all$max_r2
  max_r2 <- 0.5
  
  # Prepare data and titles for genus
  all_data_list <- list(results_all$metrics)
  
  # rename the models
  all_data_list[[1]]$Model <- c("Lasso", "Ridge", "Elastic Net", "Random Forest", "XGBoost")
  titles <- c(paste0(df_name, " BL omic predicting BMI 0-6m - Model Testing R2"))
  
  # Generate combined genus plot
  omic_plot_BL_0_6bmi <- create_plots(all_data_list, max_r2, titles)
  pdf(file = file.path(base_path, 
                       paste0("figures/BL_pred_0_6m_bmi_", df_name, ".pdf")), 
      width = 7, height = 7)
  print(omic_plot_BL_0_6bmi)
  dev.off()
}


### Make predictions on test set 

library(purrr)
# Create a list to store predictions
all_predictions <- list()
for (omic_name in names(models)) {
  omic_models <- models[[omic_name]]
  test_data_name <- paste0(omic_name, "_test")
  if (!exists(test_data_name)) {
    warning(paste("Test data", test_data_name, "not found"))
    next
  }
  test_data <- get(test_data_name)
  omic_predictions <- list()
  for (model_name in names(omic_models)) {
    model <- omic_models[[model_name]]
    model_vars <- model$finalModel$xNames
    cat("\n[INFO] Omic:", omic_name, "| Model:", model_name, "\n")
    cat("Expected predictors:", paste(model_vars, collapse = ", "), "\n")
    cat("Test data columns:", paste(colnames(test_data), collapse = ", "), "\n")
    if (!all(model_vars %in% colnames(test_data))) {
      warning(paste("Missing predictors for", model_name, "in", omic_name))
      next
    }
    complete_test_data <- test_data[complete.cases(test_data[, model_vars]), ]
    if (nrow(complete_test_data) > 0) {
      tryCatch({
        preds <- predict(model, newdata = complete_test_data)
        omic_predictions[[model_name]] <- preds
      }, error = function(e) {
        warning(paste("Prediction failed for", model_name, "in", omic_name, ":", e$message))
        omic_predictions[[model_name]] <- NA
      })
    } else {
      warning(paste("No complete test cases for", model_name, "in", omic_name))
      omic_predictions[[model_name]] <- NA
    }
  }
  
  all_predictions[[omic_name]] <- omic_predictions
}

########### Merge ##########
final_list <- list()
for (omic_name in names(all_predictions)) {
  preds <- all_predictions[[omic_name]]
  if (!all(c("elastic_net_model", "rf_model", "ridge_model") %in% names(preds))) next
  lasso <- preds$elastic_net_model
  rf <- preds$rf_model
  ridge <- preds$ridge_model
  ids <- intersect(names(rf), names(ridge))
  if (length(ids) == 0) next
  test_data_name <- paste0(omic_name, "_test")
  if (!exists(test_data_name)) next
  test_data <- get(test_data_name)
  
  # Convert rownames to a column so we can match against prediction IDs
  test_data$pred_id <- rownames(test_data)
  test_subset <- test_data[test_data$pred_id %in% ids, c("pred_id", "bmi_0_6", "subject_id")]
  
  # Create prediction dataframe
  merged_df <- data.frame(
    pred_id = ids,
    lasso_pred = as.numeric(lasso[ids]),
    rf_pred = as.numeric(rf[ids]),
    ridge_pred = as.numeric(ridge[ids]),
    omic = omic_name)
  # Merge on pred_id
  final_list[[omic_name]] <- merge(merged_df, test_subset, by = "pred_id")
}

# Combine and clean
final_df <- do.call(rbind, final_list)
rownames(final_df) <- NULL
head(final_df)


### Spliot bby model RF
rf <- final_df %>% dplyr::select(-c(lasso_pred, ridge_pred)) %>% 
  pivot_wider(names_from = omic,
              values_from = rf_pred,
              names_prefix = "prediction_")

#### ANOVA for RF 
basic_rf <- lm(bmi_0_6 ~ prediction_basic, data = rf)
meta_rf <- lm(bmi_0_6 ~ prediction_basic + prediction_meta, data = rf)
grs_rf <- lm(bmi_0_6 ~ prediction_basic + prediction_grs, data = rf)
taxa_rf <- lm(bmi_0_6 ~ prediction_basic + prediction_taxa, data = rf)
pathway_rf <- lm(bmi_0_6 ~ prediction_basic + prediction_pathway, data = rf)
micom_rf <- lm(bmi_0_6 ~ prediction_basic + prediction_micom, data = rf)
metabo_rf <- lm(bmi_0_6 ~ prediction_basic + prediction_metabo, data = rf)

anova(basic_rf, meta_rf)
anova(basic_rf, grs_rf)
anova(basic_rf, taxa_rf)
anova(basic_rf, pathway_rf)
anova(basic_rf, micom_rf)
anova(basic_rf, metabo_rf)

summary(basic_rf)
summary(meta_rf)

rf_models <- list(
  c("basic_rf", "meta_rf"),
  c("basic_rf", "grs_rf"),
  c("basic_rf", "taxa_rf"),
  c("basic_rf", "pathway_rf"),
  c("basic_rf", "micom_rf"),
  c("basic_rf", "metabo_rf"))

anova_results_rf <- list()
for (model_pair in rf_models) {
  model_1 <- get(model_pair[1])
  model_2 <- get(model_pair[2])
  # Perform ANOVA without checking the model crf
  anova_result_rf <- anova(model_1, model_2)  # This will work with any model objects that anova can handle
  tidied_result <- tidy(anova_result_rf)  # Tidy the ANOVA result using broom::tidy()
  tidied_result$model_comparison <- paste(model_pair[1], "vs", model_pair[2])
  anova_results_rf[[length(anova_results_rf) + 1]] <- tidied_result
}

# Combine all results into one dataframe
anova_table_rf <- do.call(rbind, anova_results_rf)
anova_table_rf_clean <- anova_table_rf %>%
  filter(!is.na(statistic) & !is.na(p.value)) %>%
  mutate(across(where(is.numeric), round, 3)) 

print(anova_table_rf_clean)

## Anova table with R2:

# Initialize lists to store adjusted R² values
adj_r_squared_list <- list()

# Loop through each model comparison
for (model_pair in rf_models) {
  model_1_name <- model_pair[1]
  model_2_name <- model_pair[2]
  model_1 <- get(model_1_name)
  model_2 <- get(model_2_name)
  
  # Perform ANOVA and tidy the result
  anova_result_rf <- anova(model_1, model_2)
  tidied_result <- broom::tidy(anova_result_rf)
  tidied_result$model_comparison <- paste(model_1_name, "vs", model_2_name)
  
  # Store the adjusted R² values
  adj_r_squared_list[[length(adj_r_squared_list) + 1]] <- tibble(
    model_comparison = paste(model_1_name, "vs", model_2_name),
    adj_r2_model_1 = summary(model_1)$adj.r.squared,
    adj_r2_model_2 = summary(model_2)$adj.r.squared)
  
  anova_results_rf[[length(anova_results_rf) + 1]] <- tidied_result # Store 
}

# Combine all results
anova_table_rf <- do.call(rbind, anova_results_rf)
adj_r_squared_df <- bind_rows(adj_r_squared_list)

# Clean and round ANOVA table
anova_table_rf_clean <- anova_table_rf %>%
  filter(!is.na(statistic) & !is.na(p.value)) %>%
  mutate(across(where(is.numeric), round, 3))

# Merge adjusted R² into ANOVA table
anova_table_rf_clean <- left_join(anova_table_rf_clean, adj_r_squared_df, by = "model_comparison")
print(anova_table_rf_clean)

#### TRY PLOTTING IT ####

# Prepare data: One bar for Basic model, one for each Extended
plot_data_single_basic <- anova_table_rf_clean %>%
  distinct(model_comparison, adj_r2_model_1, adj_r2_model_2, p.value) %>%
  dplyr::mutate(extended_model = gsub("basic_rf vs ", "", model_comparison))

# Define custom colors
custom_colors <- c(
  "Basic"            = "#582F2B",
  "Basic + Meta"     = "#960018",
  "Basic + Taxa"     = "#8B3A5A",
  "Basic + Micom"    = "#996759",
  "Basic + GRS"      = "#F6D1D1",
  "Basic + Pathway"  = "#FF9E80",
  "Basic + Metabo"   = "#C74375")

show_col(c("#FF9E80", "#E4717A", "#F6D1D1", "#EFDECD"))

# Prepare data with friendly labels
plot_data_long <- bind_rows(tibble(model = "Basic", 
                                   adj_r2 = plot_data_single_basic$adj_r2_model_1[1], 
                                   p.value = NA),
  tibble(model = paste("Basic +", plot_data_single_basic$extended_model),
         adj_r2 = plot_data_single_basic$adj_r2_model_2,
         p.value = plot_data_single_basic$p.value)) %>% 
  dplyr::mutate(model = dplyr::recode(model,
                               "Basic + meta_rf"    = "Basic + Meta",
                               "Basic + grs_rf"     = "Basic + GRS",
                               "Basic + taxa_rf"    = "Basic + Taxa",
                               "Basic + micom_rf"   = "Basic + Micom",
                               "Basic + pathway_rf" = "Basic + Pathway",
                               "Basic + metabo_rf"  = "Basic + Metabo",
                               "Basic"              = "Basic")) %>%
  dplyr::mutate(model = factor(model, levels = c("Basic",
                     plot_data_long %>%
                       dplyr::filter(model != "Basic") %>%
                       arrange(desc(adj_r2)) %>%
                       pull(model))))

# Plot
ggplot(plot_data_long, aes(x = model, y = adj_r2, fill = model)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(data = filter(plot_data_long, !is.na(p.value)),
    aes(label = paste0("p = ", round(p.value, 3))),
    vjust = -1, hjust = -0.1, size = 3.5, angle = 45) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(limits = c(0, 0.21)) +
  labs(title = "Adjusted R²: Basic vs Extended Models",
    x = "Model", y = "Adjusted R²",
    fill = "Model") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#### RIDGE #################################################################### ####
