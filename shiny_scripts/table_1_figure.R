# MAKE TABLE 1 & Supplementary :
rm(list = ls())
# Combined all input processing script
omic_out_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/april/anova_results/plots/"
#sink(paste0(omic_out_dir,"processing_output.txt"), split = TRUE)
#sink()
pacman::p_load(knitr, data.table, dplyr, tidyr, tableone, kableExtra, readxl,
               readr, car, RColorBrewer, gridExtra, mlbench, earth, ggplot2, missForest,
               AppliedPredictiveModeling, caret, reshape2, corrplot, stringr,
               summarytools, grid, mice, plyr, mlmRev, cowplot, ape, e1071,
               jtools, broom, patchwork, phyloseq, microbiome, glmnet, ISLR,
               MicrobiomeStat, ANCOMBC, ape, vegan, zCompositions, janitor, MASS,
               RColorBrewer, DT, ggpubr, microbiomeutilities, compositions, VIM)
library(naniar)
library(tibble)
library(table1)


meta <- read_csv("~/projects/research/Stanislawski/BMI_risk_scores/data/correct_meta_files/ashleys_meta/DRIFT_working_dataset_meta_deltas_filtered_05.21.2024.csv")
# Replace spaces with underscores in column names
colnames(meta) <- gsub(" ", "_", colnames(meta))
length(unique(meta$subject_id))
meta <- meta[meta$consent != "no", ]
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

vis_miss(a2_extra)
age <- a2_extra %>% dplyr::select(c(subject_id, age))
                                  
# long 
long_imputed <- read_csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/april_processing/long_april29.csv")

# deltas
all_delta <- read_csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/april_processing/all_delta_april29.csv")

### Supplementary figure of BMI change:

# Create the spaghetti plot
ggplot(long_imputed, aes(x = as.factor(time), y = outcome_BMI_fnl, 
                         group = subject_id, color = as.factor(randomized_group))) +
  facet_wrap(vars(randomized_group)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 2) +
  geom_smooth(aes(group = randomized_group), method = "loess", se = TRUE, 
              size = 1.2, linetype = "solid") +
  scale_color_manual(
    values = c("0" = "#218a93", "1" = "#a44d77"),  # custom colors
    labels = c("DCR", "IMF")) +
  labs(title = "BMI Trajectories Across Study Timepoints",
       x = "Time Point",
       y = expression("BMI (kg/m"^2*")"),
       color = "Diet Group") +
  theme_minimal() +
  theme(
    strip.text = element_blank(),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14))

### for table 1 ####

# Exploring other splits:
unique_subjects <- unique(long_imputed$subject_id)
# 2. 
set.seed(50) #13, 25
train_subjects <- sample(unique_subjects, size = 0.75 * length(unique_subjects))
test_subjects <- setdiff(unique_subjects, train_subjects)

# 4. Assign split
# on long 
long_new_split <- long_imputed %>%
  mutate(split_group = case_when(
    subject_id %in% train_subjects ~ "train",
    subject_id %in% test_subjects ~ "test",
    TRUE ~ NA_character_)) %>% 
  dplyr::select(-c("...1", "consent", "completer", 
                   "Peptide_YY", "Ghrelin", "Leptin")) %>%
  dplyr::mutate(time = as.factor(time),
                subject_id = as.factor(subject_id),
                randomized_group = as.factor(randomized_group),
                sex = as.numeric(sex),
                race = as.numeric(race)) %>% 
  dplyr::rename(BMI = outcome_BMI_fnl,
                range = time,
                homo_ir = HOMA_IR,
                insulin = Insulin_endo,
                LDL = LDL_Calculated,
                HDL = HDL_Total_Direct_lipid,
                HbA1c = Hemoglobin_A1C)

long_new_split_age <- merge(long_new_split, age, by = "subject_id", all.x = TRUE)

# Define the column names based on your lists
basic <- c('subject_id','BMI', 'range', 'sex', 'age') #
meta_keep <- c('subject_id','BMI', 'range', 'randomized_group', 'sex', 'race', 
               'age', 'HbA1c', 'HDL', 'homo_ir', 'insulin', 'LDL', 'Glucose.x') # 
only_taxa <- c('subject_id','BMI', 'range', grep("^g__", names(long_new_split_age), value = TRUE))

micom_start <- which(names(long_new_split_age) == "Diacetyl")
micom_end <- which(names(long_new_split_age) == "aldehydo.D.xylose")
only_micom <- c('subject_id','BMI', 'range', names(long_new_split_age)[micom_start:micom_end])

path_start <- which(names(long_new_split_age) == "arginine..ornithine.and.proline.interconversion")
path_end <- which(names(long_new_split_age) == "UDP.N.acetyl.D.glucosamine.biosynthesis.I")
only_pathway <- c('subject_id','BMI', 'range', names(long_new_split_age)[path_start:path_end])

tabo_start <- which(names(long_new_split_age) == "non_HDL_C")
tabo_end <- which(names(long_new_split_age) == "IDL_TG_pct")
only_tabo <- c('subject_id','BMI', 'range', names(long_new_split_age)[tabo_start:tabo_end])

all_col <- c('subject_id','BMI', 'range',
             'randomized_group', 'sex', 'race', 
             'age', 'HbA1c', 'HDL', 'homo_ir', 'insulin', 'LDL', 'Glucose.x',
             grep("^g__", names(long_new_split), value = TRUE),
             names(long_new_split)[micom_start:micom_end],
             names(long_new_split)[path_start:path_end],
             names(long_new_split)[tabo_start:tabo_end])

long_feature <- data.frame(Group = c("Basic", "Meta", "Taxa", 
                                     "Micom", "Pathway", "Metabo", "All"),
                           NumColumns = c((length(basic)-3), (length(meta_keep)-3), (length(only_taxa)-3),
                                          (length(only_micom)-3), (length(only_pathway)-3),
                                          (length(only_tabo)-3), (length(all_col)-3)))
colnames(long_feature) <- c("Data Type", "Number of Features")
ft <- flextable(long_feature)
ft <- autofit(ft)
output_path <- file.path(omic_out_dir, "paper_render/basic_plus_diet/Long_Feature_Count_Table.docx")


# save to word
doc <- read_docx()
doc <- body_add_flextable(doc, ft)
print(doc, target = output_path)

long_new_split_table <- long_new_split_age %>% 
  dplyr::rename(`Diet Group IMF` = randomized_group,
                `Cohort Number` = cohort_number,
                Sex = sex,
                Race = race,
                Age = age.y, 
                Glucose = Glucose.x,
                `Train / Test Subset` = split_group)

# Rename factor levels first
long_new_split_table$Timepoint <- factor(long_new_split_table$range,
                                         levels = c(0, 6, 12),
                                         labels = c("Baseline", "6 Months", "12 Months"))
long_new_split_table$Race <- factor(long_new_split_table$Race,
                                    levels = c(2, 3, 5, 7),
                                    labels = c("Asian", 
                                               "Black or African American",
                                               "White",
                                               "Other"))
long_new_split_table$Sex <- factor(long_new_split_table$Sex,
                                   levels = c(0,1),
                                   labels = c("Female", "Male"))

long_new_split_table$`Diet Group IMF` <- factor(long_new_split_table$`Diet Group IMF`,
                                                levels = c(0,1),
                                                labels = c("DCR", "IMF"))

long_new_split_table0 <- long_new_split_table %>% dplyr::filter(Timepoint == "Baseline")
long_new_split_table6 <- long_new_split_table %>% dplyr::filter(Timepoint == "6 Months")
long_new_split_table12 <- long_new_split_table %>% dplyr::filter(Timepoint == "12 Months")

### FLEXTABLES ###

# Stats for the tables below:
compare_groups <- function(data, group_var, variables, factor_vars) {
  results <- list()
  for (var in variables) {
    if (var %in% factor_vars) {
      # Chi-squared test for categorical variables
      tbl <- table(data[[group_var]], data[[var]])
      test <- chisq.test(tbl)
      results[[var]] <- data.frame(
        Variable = var,
        Test = "Chi-squared",
        P_value = test$p.value)
    } else {
      # t-test for continuous variables
      test <- t.test(data[[var]] ~ data[[group_var]])
      results[[var]] <- data.frame(
        Variable = var,
        Test = "t-test",
        P_value = test$p.value)
    }
  }
  do.call(rbind, results)
}

# List of timepoints
timepoints <- unique(long_new_split_table$Timepoint)

# Variables to compare
vars_to_check <- c("Sex", "Age", "Race", "Cohort Number", "BMI")
factor_vars <- c("Sex", "Race", "Cohort Number")

# Store all results
all_results <- list()
for (tp in timepoints) {
  subset_data <- subset(long_new_split_table, Timepoint == tp)
  
  test_results <- compare_groups(
    data = subset_data,
    group_var = "Diet Group IMF", 
    variables = vars_to_check,
    factor_vars = factor_vars)
  test_results$Timepoint <- tp
  all_results[[tp]] <- test_results
}

final_results <- do.call(rbind, all_results)
write.csv(final_results, paste0(omic_out_dir, "paper_render/basic_plus_diet/group_comparison_by_timepoint.csv"), row.names = FALSE)

# all long
table_long1 <- CreateTableOne(
  data = long_new_split_table[c("Timepoint", "Sex" , "Age", "Diet Group IMF", 
                                "Race", "Cohort Number", "BMI")],
  strata = "Timepoint",
  factorVars = c("Sex", "Race", "Cohort Number", "Diet Group IMF"),
  addOverall = FALSE,
  smd = TRUE) %>%
  print(showAllLevels = TRUE) %>%
  as.data.frame() %>% 
  dplyr::rename(`6-Months` = `6 Months`,
                `12-Months` = `12 Months`,
                `P-value` = p,
                `Level` = level) #

table_long1 <- table_long1 %>%
  mutate(across(where(is.character), ~ str_replace_all(., "\\(\\s+", "("))) %>% 
  mutate(across(where(is.numeric), ~ round(., 1))) %>% 
  mutate(RowName = rownames(.)) %>%
  filter(!RowName %in% c("Timepoint....", "X", "X.1")) %>%
  dplyr::select(-test) %>%
  relocate(RowName) %>%
  mutate(RowName = RowName |> 
           str_remove("^X\\.\\.\\.") |> 
           str_replace_all("\\.\\.\\.\\.", " (%)") |>
           str_replace_all("\\.", " ") |>
           str_replace_all("..mean..SD..", " (mean (SD))"), 
         RowName = if_else(str_starts(RowName, "X "), "", RowName))

set_flextable_defaults(
  border.color = "grey30", 
  font.family = "Arial",
  font.size = 10, 
  font.color = "black", 
  padding = 1,
  line_spacing = 1.5)

ft <- table_long1 %>%
  flextable() %>%
  autofit() %>%
  align(j = c(2:6), align = c("center"), part = "body") %>% 
  set_header_labels(RowName = "", 
                    Level = "Group") %>% 
  hline(i = c(1, 3, 4, 6, 10, 15)) %>% 
  vline(j = c('Level', '12-Months'), part = "body") %>% 
  set_caption(caption = as_paragraph(as_b("Characteristics of Included Participants per Timepoint"))) %>% 
  add_header_row(values = c("Characteristic", "Timepoint"), 
                 colwidths = c(3, 3)) %>% 
  bold(bold = TRUE, i = 2, part = "header") %>% 
  bg(i = 2, part = "header", bg = "grey80") 

ft

# Export to Word
doc <- read_docx() %>%
  body_add_flextable(value = ft) %>%
  body_add_par("", style = "Normal")  # Optional: adds a blank line after the table

print(doc, target = paste0(omic_out_dir, "paper_render/basic_plus_diet/table1_flextable_output.docx"))

# long 0m
table0 <- CreateTableOne(data = long_new_split_table0[c(4:8, 397:398)], 
                         strata = c("Train / Test Subset"),
                         factorVars = c("Sex", "Race", "Cohort Number"),
                         addOverall = TRUE) %>%
  print(showAllLevels = TRUE) %>% 
  as.data.frame()
colnames(table0) <- make.names(colnames(table0), unique = TRUE)
table0 <- subset(table0, select = -test.1)  %>% 
  dplyr::rename(`Overall Baeline` = Overall,
                `Test Set Baseline` = test,
                `Train Set Baeline` = train,
                `Baseline P-value` = p)

# long 6m
table6 <- CreateTableOne(data = long_new_split_table6[c(4:8, 397:398)], 
                         strata = c("Train / Test Subset"),
                         factorVars = c("Sex", "Race", "Cohort Number"),
                         addOverall = FALSE) %>% 
  print(showAllLevels = TRUE) %>% 
  as.data.frame()
colnames(table6) <- make.names(colnames(table6), unique = TRUE)
table6 <- subset(table6, select = -test.1) %>% 
  dplyr::rename(`Test Set 6m` = test,
                `Train Set 6m` = train,
                `P-value 6m` = p)

# long 12m
table12 <- CreateTableOne(data = long_new_split_table12[c(4:8, 397:398)], 
                          strata = c("Train / Test Subset"),
                          factorVars = c("Sex", "Race", "Cohort Number"),
                          addOverall = FALSE) %>% 
  print(showAllLevels = TRUE) %>% 
  as.data.frame() 
colnames(table12) <- make.names(colnames(table12), unique = TRUE)
table12 <- subset(table12, select = -test.1) %>% 
  dplyr::rename(`Test Set 12m` = test,
                `Train Set 12m` = train,
                `P-value 12m` = p)

table0$RowName <- rownames(table0)
table6$RowName <- rownames(table6)
table12$RowName <- rownames(table12)
t0_6 <- merge(table0, table6, by = "RowName", all = TRUE, sort = FALSE)
t0_6_12 <- merge(t0_6, table12, by = "RowName", all = TRUE, sort = FALSE)

t0_6_12$RowName <- t0_6_12$RowName |>
  str_remove("^X\\.\\.\\.") |>           # remove prefix X...
  str_replace_all("\\.\\.\\.", "") |>   # remove triple dots
  str_replace_all("\\.", " ")           # replace remaining dots with spaces

t0_6_12 <- t0_6_12 %>% 
  dplyr::filter(RowName != "TrainTest Subset ") %>% 
  dplyr::filter(RowName != "X 9") %>% 
  dplyr::select(-c("level", "level.y")) %>%
  mutate(RowName = case_when(
    RowName == "Diet Group IMF " ~ "Diet Group (%)",
    RowName == "BMI  mean  SD  " ~ "BMI (mean ± SD)",
    RowName == "Age  mean  SD  " ~ "Age (mean ± SD)",
    RowName == "Sex " ~ "Sex (%)",
    RowName == "Cohort Number " ~ "Cohort Number (%)",
    RowName == "Race " ~ "Race (%)",
    TRUE ~ RowName))

t0_6_12_fixed <- t0_6_12

# Define which rows are factor levels (based on how you've structured the table)
# Here, we blank out the RowName for rows that are levels
t0_6_12_fixed$RowName <- ifelse(
  grepl("^X|^\\s*$", t0_6_12_fixed$RowName), 
  "", t0_6_12_fixed$RowName)

library(flextable)
styled_table <- t0_6_12_fixed %>%
  flextable() %>%
  set_header_labels(
    RowName = "Measure",
    level.x = "Level",
    `Overall Baeline` = "Overall N",
    `Test Set Baseline` = "Baseline Test",
    `Train Set Baeline` = "Baseline Train",
    `Baseline P-value` = "P (Base)",
    `Test Set 6m` = "6m Test",
    `Train Set 6m` = "6m Train",
    `P-value 6m` = "P (6m)",
    `Test Set 12m` = "12m Test",
    `Train Set 12m` = "12m Train",
    `P-value 12m` = "P (12m)") %>%
  autofit() %>%
  theme_box() %>%
  fontsize(size = 8) %>% 
  padding(padding = 1) %>% 
  padding(i = ~ RowName == "", j = "level.x", padding.left = 2) %>%
  #color(i = ~ RowName %in% c("Sex", "Race", "Diet Group", "Cohort Number"), color = "Black") %>%
  bg(i = ~ RowName %in% c("Sex (%)", "Race (%)", "Diet Group (%)", "Cohort Number (%)"), bg = "white") %>%
  bg(j = c("Baseline P-value", "P-value 6m", "P-value 12m"), bg = "#cfe2f3") %>% 
  width(width = 0.9) %>%
  width(j = c("Baseline P-value", "P-value 6m", "P-value 12m"), width = 0.7)


long_t1 <- paste0(omic_out_dir, "paper_render/basic_plus_diet/long_table1.docx")
props_long <- officer::prop_section(
  page_size = officer::page_size(width = 21 / 2.54, height = 29.7 / 2.54, 
                                 orient = "landscape", unit = "in"),
  page_margins = officer::page_mar(top = 0.5, bottom = 0.2, left = 0.2, right = 0.2))

save_as_docx(`long_table_one` = styled_table, 
             pr_section = props_long,
             path = long_t1)

# You are comparing 2 × 3 = 6 groups
# Each p-value is testing for differences across those 6 groups
table1 <- CreateTableOne(data = long_new_split_table[c(4:8, 397:399)], 
                         strata = c("Train / Test Subset", "Timepoint"),
                         factorVars = c("Sex", "Race", "Cohort Number"),
                         addOverall = TRUE)

# Convert to data frame
table1_df <- table1 %>%
  print() %>% 
  as.data.frame()

rows_to_remove <- c("Train / Test Subset = train (%)", "Timepoint (%)", 
                    "   Baseline", "   6 Months", "   12 Months")
table1_df <- table1_df[!rownames(table1_df) %in% rows_to_remove, ]
rownames(table1_df) <- gsub("^X\\.\\.\\.|^X\\.\\.\\.", "", rownames(table1_df))
rownames(table1_df) <- gsub("\\.\\.\\.", "", rownames(table1_df))
rownames(table1_df) <- gsub("\\.", " ", rownames(table1_df))
rows_to_remove <- c("TrainTest Subsettrain ", "Timepoint ", 
                    "Baseline", "6 Months", "12 Months")
table1_df_edit <- table1_df[!rownames(table1_df) %in% rows_to_remove, ] 

# Create Table 2: Number of features per group
feature_table <- data.frame(
  Variable = c("basic", "meta", "taxa", "micom", "pathway", "metabo", "all"),
  Test = rep("", 7),  # To match column structure
  Train = c(length(basic), length(meta_keep), length(only_taxa),
            length(only_micom), length(only_pathway), 
            length(only_tabo), length(all_col)),
  p.value = rep("", 7))

#################### TRYING ASHLEYS WYA OF COUTING FOR LONG 

subject_counts <- long_new_split_table %>% 
  dplyr::select(as.character("subject_id"), range) %>% 
  unique() %>% 
  dplyr::select(range) %>% 
  table() %>% 
  as.data.frame() %>% 
  dplyr::rename("All Subjects" = "Freq") %>% t()

train_counts <- long_new_split_table %>% 
  dplyr::filter(`Train / Test Subset` == "train") %>% 
  dplyr::select(range) %>% 
  table() %>% 
  as.data.frame() %>%
  dplyr::rename("Train Set" = "Freq") %>% t()

test_counts <- long_new_split_table %>% 
  dplyr::filter(`Train / Test Subset` == "test") %>% 
  dplyr::select(range) %>% 
  table() %>% 
  as.data.frame() %>%
  dplyr::rename("Test Set" = "Freq") %>% t()


subject_counts <- as.data.frame(subject_counts)
names(subject_counts) <- subject_counts[1,]
subject_counts <- subject_counts[-1,]

train_counts <- as.data.frame(train_counts)
names(train_counts) <- train_counts[1,]
train_counts <- train_counts[-1,]

test_counts <- as.data.frame(test_counts)
names(test_counts) <- test_counts[1,]
test_counts <- test_counts[-1,]
total_long <- rbind(subject_counts, train_counts, test_counts) 
knitr::kable(total_long)

# try both time and split

### Summary for just train and test
summary_train_test <- long_new_split_table %>%
  dplyr::rename(Subset = `Train / Test Subset`) %>%
  dplyr::group_by(Subset, range) %>%
  dplyr::summarise(
    Mean_Age = mean(Age, na.rm = TRUE),
    `Proportion Male` = mean(Sex == 1, na.rm = TRUE),
    `Proportion IMF` = mean(`Diet Group IMF` == 1, na.rm = TRUE),
    .groups = "drop")

### Summary for all
summary_all <- long_new_split_table %>%
  dplyr::group_by(range) %>%
  dplyr::summarise(
    Subset = "all",
    Mean_Age = mean(Age, na.rm = TRUE),
    `Proportion Male` = mean(Sex == 1, na.rm = TRUE),
    `Proportion IMF` = mean(`Diet Group IMF` == 1, na.rm = TRUE),
    .groups = "drop") %>%
  dplyr::select(Subset, range, everything())

### Combine both
summary_combined <- bind_rows(summary_train_test, summary_all)

summary_table <- summary_combined %>%
  tidyr::pivot_longer(
    cols = -c(Subset, range),
    names_to = "Variable",
    values_to = "Value") %>%
  tidyr::pivot_wider(
    names_from = range,
    values_from = Value)

knitr::kable(summary_table)

library(skimr)
group_by(long_new_split_table[1:10], range) %>% skim()

# now on delta ###################################################################

delta_new_split <- all_delta %>%
  mutate(split_group = case_when(
    subject_id %in% train_subjects ~ "train",
    subject_id %in% test_subjects ~ "test",
    TRUE ~ NA_character_)) %>% 
  dplyr::select(-c("...1", "consent", "completer", 
                   "Peptide_YY", "Ghrelin", "Leptin")) %>%
  dplyr::mutate(time = as.factor(time),
                subject_id = as.factor(subject_id),
                randomized_group = as.factor(randomized_group),
                sex = as.numeric(sex),
                race = as.numeric(race)) %>% 
  dplyr::rename(BMI = outcome_BMI_fnl,
                range = time,
                homo_ir = HOMA_IR,
                insulin = Insulin_endo,
                LDL = LDL_Calculated,
                HDL = HDL_Total_Direct_lipid,
                HbA1c = Hemoglobin_A1C)

delta_new_split_table <- delta_new_split %>% 
  dplyr::rename(`Diet Group IMF` = randomized_group,
                `Cohort Number` = cohort_number,
                Sex = sex,
                Race = race,
                Age = age, 
                Glucose = Glucose.x,
                `Train / Test Subset` = split_group)

# Rename factor levels first
delta_new_split_table$Timepoint <- factor(delta_new_split_table$range,
                                          levels = c(0, 6, 12),
                                          labels = c("Baseline", "6 Months", "12 Months"))

# T1 with omic feature count 
basic <- c('subject_id','BMI', 'range', 'sex', 'age') #
meta_keep <- c('subject_id','BMI', 'range', 'randomized_group', 'sex', 'race', 
               'age', 'HbA1c', 'HDL', 'homo_ir', 'insulin', 'LDL', 'Glucose.x') # 
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

# Create data frames based on the columns defined
basic <- all_deltas[, basic, drop = FALSE] %>% unique()
meta <- all_deltas[, meta_keep, drop = FALSE] %>% unique()
taxa <- all_deltas[, only_taxa, drop = FALSE] %>% unique()
micom <- all_deltas[, only_micom, drop = FALSE] %>% unique()
pathway <- all_deltas[, only_pathway, drop = FALSE] %>% unique()
metabo <- all_deltas[, only_tabo, drop = FALSE] %>% unique()
all <- all_deltas[, all_col, drop = FALSE] %>% unique() %>% 
  dplyr::mutate(randomized_group = as.numeric(randomized_group))


#################### TRYING ASHLEYS WYA OF COUTING FOR DELTA

subject_counts_delta <- delta_new_split_table %>% 
  dplyr::select(as.character("subject_id"), range) %>% 
  unique() %>% 
  dplyr::select(range) %>% 
  table() %>% 
  as.data.frame() %>% 
  dplyr::rename("All Subjects" = "Freq") %>% t()

train_counts_delta <- delta_new_split_table %>% 
  dplyr::filter(`Train / Test Subset` == "train") %>% 
  dplyr::select(range) %>% 
  table() %>% 
  as.data.frame() %>%
  dplyr::rename("Train Set" = "Freq") %>% t()

test_counts_delta <- delta_new_split_table %>% 
  dplyr::filter(`Train / Test Subset` == "test") %>% 
  dplyr::select(range) %>% 
  table() %>% 
  as.data.frame() %>%
  dplyr::rename("Test Set" = "Freq") %>% t()


subject_counts_delta <- as.data.frame(subject_counts_delta)
names(subject_counts_delta) <- subject_counts_delta[1,]
subject_counts_delta <- subject_counts_delta[-1,]

train_counts_delta <- as.data.frame(train_counts_delta)
names(train_counts_delta) <- train_counts_delta[1,]
train_counts_delta <- train_counts_delta[-1,]

test_counts_delta <- as.data.frame(test_counts_delta)
names(test_counts_delta) <- test_counts_delta[1,]
test_counts_delta <- test_counts_delta[-1,]
total_delta <- rbind(subject_counts_delta, train_counts_delta, test_counts_delta) 
knitr::kable(total_delta)

###### GRS PLOT
library(geomtextpath)
library(hrbrthemes)

ggplot(long_imputed, 
       aes(x = bmi_prs, 
           y = outcome_BMI_fnl)) +
  geom_point(aes(color = sex)) +
  geom_labelsmooth(aes(label = sex, 
                       color = sex, 
                       group = sex),
                   method = "lm", formula = y ~ x,
                   fill = "white", size = 6, 
                   linewidth = 1.5, 
                   boxlinewidth = 0.6) +
  theme_bw() +
  guides(color = 'none')


ggplot(all_delta, 
       aes(x = bmi_prs, 
           y = outcome_BMI_fnl)) +
  geom_point(aes(color = sex)) +
  geom_labelsmooth(aes(label = sex, 
                       color = sex, 
                       group = sex),
                   method = "lm", formula = y ~ x,
                   fill = "white", size = 4, 
                   linewidth = 1.5, 
                   boxlinewidth = 0.6) +
  theme_bw() +
  guides(color = 'none')
