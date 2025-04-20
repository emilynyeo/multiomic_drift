library(dplyr)
library(tidyr)
library(tools)
library(broom)
library(stargazer)
library(sjPlot)
library(lme4)
library(broom)
library(knitr)

# Set the input directory
delta_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/april_processing/"

# Read the data from the CSV file
delta <- read.csv(file.path(delta_dir, "all_delta.csv"))

# Make train and test sets 
# test sample names
test_names <- c("ABR-079", "AGA-071", "AHE-055", "ALI-121", "ALO-163", "AMA-031", "ASO-013", "AWI-167", "BMO-164", "CWA-183", "DSC-024", "EBE-130", "EHI-177", "EJO-092", "GFU-188", "HGI-010", "JCA-109", "JGO-100","KBU-085", "KCE-034", "KHE-170", "LDO-148", "LST-186", "LZD-142", "MAR-119", "MCA-088", "MJA-153", "MWE-112", "NPO-149", "RAE-114", "SBO-020", "SEG-080", "SKA-195", "SLO-178", "SSH-028", "TDU-086", "TFA-016", "VCA-041")

# train sample names
train_names <- c("AAL-144", "ACO-053", "ADA-105", "AKE-009", "AKI-011", "AKO-139", "AMC-155", "AME-128",  "AME-157", "ATA-129", "AWA-052", "AWA-083", "BAN-193", "BHO-014", "BIN-201", "BKN-104", "BMI-156", "BSA-174", "CAM-057", "CCO-189", "CED-026", "CEL-073", "CGA-134", "CIS-077", "CKR-078", "CLE-049", "COW-066", "CRO-108", "CWA-161", "EBE-051", "EKA-135", "EKR-045", "ELA-159", "EPO-182", "EVO-184", "FWI-098", "GHA-035", "HDE-154", "IBE-120", "JDI-140", "JER-110", "JFU-027", "JJO-093", "JKN-127", "JPO-022", "JUG-116", "JUT-032", "JVE-126", "KAN-138", "KBR-162", "KEL-185", "KEL-199", "KGI-029", "KHU-196", "KPA-042", "KRI-072", "KVA-038", "KWA-122", "KWA-141", "LBL-047", "LBU-015", "LEL-147", "LFI-003", "LJA-101",  "LMC-111", "LPF-198", "LVA-017", "MBA-187", "MCW-065", "MDI-107", "MES-068", "MFB-118", "MGA-076", "MHO-117", "MKE-192", "MMA-036", "MRT-179", "MSH-091", "MST-039", "MWE-143",  "MWO-133", "MWY-152", "NAR-099", "NBI-048", "NBI-069", "NCO-171", "NDI-067", "NEL-094", "NKA-090", "NMO-151", "NTA-021", "PBE-123", "QNG-166", "RAF-125", "RAM-050", "RHP-023",  "RLA-132", "ROL-006", "SAB-160", "SCA-043", "SCR-061", "SDA-150", "SGA-062", "SKA-087", "SRO-194", "TBU-115", "TFA-172", "TRO-113", "TSH-146", "TSL-056", "WPE-005", "YOR-103",  "YSU-097", "ZVU-096")

cat("Length of test names:", length(test_names), "\n")
cat("Length of train names:", length(train_names), "\n")

# Create train and test sets for the current subset
train_set <- delta[delta[["subject_id"]] %in% train_names, ]  # Select rows where ID_VAR is in train_names
test_set <- delta[delta[["subject_id"]] %in% test_names, ]    # Select rows where ID_VAR is in test_names

test <- test_set %>% dplyr::select(c(subject_id, outcome_BMI_fnl, time))

predict_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/april/final_merf_dfs"

basic <- read.csv(file.path(predict_dir, "basic_delta.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% 
  dplyr::rename(y_new_basic = y_hat_new)

meta <- read.csv(file.path(predict_dir, "meta_keep_delta.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% 
  dplyr::rename(y_new_meta = y_hat_new)

grs <- read.csv(file.path(predict_dir, "only_grs_delta.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% 
  dplyr::rename(y_new_grs = y_hat_new)

taxa <- read.csv(file.path(predict_dir, "only_taxa_delta.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% 
  dplyr::rename(y_new_taxa = y_hat_new)

pathway <- read.csv(file.path(predict_dir, "only_pathway_delta.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% 
  dplyr::rename(y_new_pathway = y_hat_new)

micom <- read.csv(file.path(predict_dir, "only_micom_delta.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% 
  dplyr::rename(y_new_micom = y_hat_new)

metabo <- read.csv(file.path(predict_dir, "only_metabo_delta.csv")) %>% 
  dplyr::select(-starts_with("Top_15_"), -starts_with("R_squared")) %>% 
  distinct() %>% 
  dplyr::rename(y_new_metabo = y_hat_new)

# Merge the first two data frames (grs and meta)
merged_basic_meta <- merge(basic, meta, by = c("Model", "Cluster", "Time"))
merged_grs_meta <- merge(grs, merged_basic_meta, by = c("Model", "Cluster", "Time"))
merged_grs_meta_micom <- merge(merged_grs_meta, micom, 
                               by = c("Model", "Cluster", "Time"))
merged_grs_meta_micom_pathway <- merge(merged_grs_meta_micom, pathway, 
                                       by = c("Model", "Cluster", "Time"))
merged_grs_meta_micom_pathway_metab <- merge(merged_grs_meta_micom_pathway, metabo, 
                                             by = c("Model", "Cluster", "Time"))
merged_all <- merge(merged_grs_meta_micom_pathway_metab, taxa, 
                    by = c("Model", "Cluster", "Time"))

# Merge merged_all with test by matching 'Cluster' with 'all_samples' and 'Time' with 'time'
merged <- merge(merged_all, test, 
                by.x = c("Cluster", "Time"), 
                by.y = c("subject_id", "time")) %>% 
  dplyr::rename(bmi = outcome_BMI_fnl)

# For each prediction column, compute R-squared and store as new column
prediction_cols <- grep("^y_new", names(merged), value = TRUE)
for (col in prediction_cols) {
  r_squared <- cor(merged[[col]], merged$bmi, use = "complete.obs")^2
  new_col_name <- paste0("R2_", col)
  merged[[new_col_name]] <- r_squared
}

# Select just the R2 columns and Model
rm(merged_df, merged_df_small, test_all, merged_grs_meta,
   merged_grs_meta_micom, merged_grs_meta_micom_pathway)

### Filter by model
df_mse <- merged %>%
  dplyr::filter(Model == "MSE Model")
df_prev <- merged %>%
  dplyr::filter(Model == "Prev Model")
df_ptev <- merged %>%
  dplyr::filter(Model == "PTEV Model")
df_oob <- merged %>%
  dplyr::filter(Model == "OOB Model")

# Remove duplicates 
df_mse <- distinct(df_mse)
df_prev <- distinct(df_prev)
df_ptev <- distinct(df_ptev)
df_oob <- distinct(df_oob)

# Create a named list of dataframes
dataframes <- list(
  df_mse = df_mse,
  df_prev = df_prev,
  df_ptev = df_ptev,
  df_oob = df_oob
)

# Loop over each dataframe
for (df_name in names(dataframes)) {
  mod_dat <- dataframes[[df_name]]
  # Fit models
  lmer_basic <- lmer(bmi ~ y_new_basic + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_meta_b <- lmer(bmi ~ y_new_basic + y_new_meta + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_micom_b <- lmer(bmi ~ y_new_basic + y_new_micom + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_path_b <- lmer(bmi ~ y_new_basic + y_new_pathway + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_tax_b <- lmer(bmi ~ y_new_basic + y_new_taxa + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_metabo_b <- lmer(bmi ~ y_new_basic + y_new_metabo + (1|Cluster), data = mod_dat, REML = FALSE)
  
  # Optional: ANOVA tests
  anova(lmer_basic, lmer_meta_b)
  anova(lmer_basic, lmer_micom_b)
  anova(lmer_basic, lmer_path_b)
  anova(lmer_basic, lmer_tax_b)
  anova(lmer_basic, lmer_metabo_b)
  
  # Table for model comparisons (optional, adjust to your models)
  sjPlot::tab_model(
    lmer_basic, lmer_meta_b, lmer_tax_b, lmer_path_b, lmer_micom_b,
    title = paste("glmlasso delta sequential lmer models -", df_name),
    string.pred = "Predictors",
    string.est = "Estimate",
    string.std = "std. Beta",
    string.ci = "95% CI",
    string.se = "std. Error",
    p.style = c("numeric"),
    p.threshold = c(0.05),
    auto.label = FALSE
  )
  
  # Define model pairs
  glmmlass_lmer_models <- list(
    c("lmer_basic", "lmer_meta_b"),
    c("lmer_basic", "lmer_tax_b"),
    c("lmer_basic", "lmer_micom_b"),
    c("lmer_basic", "lmer_path_b"),
    c("lmer_basic", "lmer_metabo_b")
  )
  
  anova_results <- list()
  for (model_pair in glmmlass_lmer_models) {
    model_1 <- get(model_pair[1])
    model_2 <- get(model_pair[2])
    
    if (inherits(model_1, "lmerMod") && inherits(model_2, "lmerMod")) {
      anova_result <- anova(model_1, model_2)
      tidied_result <- broom::tidy(anova_result)
      tidied_result$model_comparison <- paste(model_pair[1], "vs", model_pair[2])
      anova_results[[length(anova_results) + 1]] <- tidied_result
    } else {
      message("One of the models in the pair is not a valid lmerMod object.")
    }
  }
  
  # Combine and clean ANOVA results
  anova_table <- do.call(rbind, anova_results)
  anova_table_clean <- anova_table %>%
    filter(!is.na(statistic) & !is.na(p.value)) %>%
    mutate(across(where(is.numeric), round, 3))
  
  print(anova_table_clean)
  html_table <- kable(anova_table_clean, format = "html", table.attr = "class='table table-striped'")
  
  # Save the table
  output_path <- paste0("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/april/anova_results/merf_delta_anova_table_", df_name, ".html")
  #writeLines(html_table, output_path)
}
