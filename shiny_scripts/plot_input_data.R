# Load and rename Script 

library(pacman)
p_load(tools, reticulate, viridis, tidyplots, patchwork, jsonlite, maps, ggvenn, 
       caret, caretEnsemble, glmnet, xgboost, ggplot2, glmmLasso, corrplot,
       readr, plyr, dplyr, tidyr, purrr, tibble, stringr, psych, randomForest,  
       reshape2, scales, gridExtra, plotly, sf, tidyverse, naniar, VIM, gridExtra,
       sjPlot, htmltools, officer, flextable, webshot, apaTables, MuMIn, lme4, 
       glue, grid, rsq, pheatmap, GGally, VennDiagram, glmmTMB, broom.mixed, gt,
       patchwork, tidyverse, ggbeeswarm, scales, viridis)

source("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/shiny_scripts/plotting_utils.R")

# Load Data from Models ##############################################################################################################

# MERF LONG 
#predict_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/april/final_merf_dfs"
#pred_dir_long <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/april/new_split"
pred_dir_long <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/may_basic_plus/long"
#merf_long <- read.csv('/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/april/anova_results/new_split/april_long_predictions_df_april29.csv') %>% 
merf_long <- read.csv('/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/may_basic_plus/long/merf_long_predictions_df.csv') %>%  
  dplyr::filter(Model == "MSE Model") %>% 
  rename(y_new_meta_noas = y_new_meta_no_age_sex,
         y_new_all_noas = y_new_all_no_age_sex,
         y_new_all_nclin = y_new_all_no_clin)

basic <- read.csv(file.path(pred_dir_long, "basic_long.csv")) %>% dplyr::rename(y_new_basic = y_hat_new)
meta <- read.csv(file.path(pred_dir_long, "meta_keep_long.csv")) %>% dplyr::rename(y_new_meta = y_hat_new)
meta_no_age_sex <- read.csv(file.path(pred_dir_long, "meta_keep_no_age_sex_long.csv")) %>% dplyr::rename(y_new_meta_no_age_sex = y_hat_new)
grs <- read.csv(file.path(pred_dir_long, "only_grs_long.csv")) %>% dplyr::rename(y_new_grs = y_hat_new)
taxa <- read.csv(file.path(pred_dir_long, "only_taxa_long.csv")) %>% dplyr::rename(y_new_taxa = y_hat_new)
pathway <- read.csv(file.path(pred_dir_long, "only_pathway_long.csv")) %>% dplyr::rename(y_new_pathway = y_hat_new)
micom <- read.csv(file.path(pred_dir_long, "only_micom_long.csv")) %>% dplyr::rename(y_new_micom = y_hat_new)
metabo <- read.csv(file.path(pred_dir_long, "only_metabo_long.csv")) %>% dplyr::rename(y_new_metabo = y_hat_new)
all <- read.csv(file.path(pred_dir_long, "only_all_long.csv")) %>% dplyr::rename(y_new_all = y_hat_new)
all_no_age_sex <- read.csv(file.path(pred_dir_long, "only_all_no_age_sex_long.csv")) %>% dplyr::rename(y_new_all_no_age_sex = y_hat_new)
all_no_clin <- read.csv(file.path(pred_dir_long, "only_all_no_clin_long.csv")) %>% dplyr::rename(y_new_all_no_clin = y_hat_new)

# MERF DELTA
merf_delta <- read.csv('/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/may_basic_plus/delta/merf_delta_predictions_df.csv') %>% 
  dplyr::filter(Model == "MSE Model") %>% 
  rename(y_new_meta_noas = y_new_meta_no_age_sex, 
         y_new_all_noas = y_new_all_no_age_sex,
         y_new_all_nclin = y_new_all_no_clin)
predict_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/may_basic_plus/"
basic_md <- read.csv(file.path(predict_dir, "basic_delta_april29.csv")) %>% dplyr::rename(y_new_basic = y_hat_new)
meta_md <- read.csv(file.path(predict_dir, "meta_keep_delta_april29.csv")) %>% dplyr::rename(y_new_meta = y_hat_new)
meta_no_age_sex_md <- read.csv(file.path(predict_dir, "meta_keep_no_age_sex_delta_april29.csv")) %>% dplyr::rename(y_new_meta_no_age_sex = y_hat_new)
grs_md <- read.csv(file.path(predict_dir, "only_grs_delta_april29.csv")) %>% dplyr::rename(y_new_grs = y_hat_new)
taxa_md <- read.csv(file.path(predict_dir, "only_taxa_delta_april29.csv")) %>% dplyr::rename(y_new_taxa = y_hat_new)
pathway_md <- read.csv(file.path(predict_dir, "only_pathway_delta_april29.csv")) %>% dplyr::rename(y_new_pathway = y_hat_new)
micom_md <- read.csv(file.path(predict_dir, "only_micom_delta_april29.csv")) %>% dplyr::rename(y_new_micom = y_hat_new)
metabo_md <- read.csv(file.path(predict_dir, "only_metabo_delta_april29.csv")) %>% dplyr::rename(y_new_metabo = y_hat_new)
all_md <- read.csv(file.path(predict_dir, "only_all_delta_april29.csv")) %>% dplyr::rename(y_new_all = y_hat_new)
all_no_age_sex_md <- read.csv(file.path(predict_dir, "only_all_no_age_sex_delta_april29.csv")) %>% dplyr::rename(y_new_all_no_age_sex = y_hat_new)
all_no_clin_md <- read.csv(file.path(predict_dir, "only_all_no_clin_delta_april29.csv")) %>% dplyr::rename(y_new_all_no_clin = y_hat_new)

# GLMMLASSO LONG
#gl_lasso_long <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_long/new_split_may/april_long_predictions_df_april29.csv") %>% 
gl_lasso_long <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/long/june_lasso_long_predictions_df.csv") %>%
  dplyr::rename(y_new_meta = y_new_meta_only,
                y_new_meta_noas = y_new_meta_nas_only,
                y_new_grs = y_new_grs_only,
                y_new_taxa = y_new_tax_only,
                y_new_micom = y_new_micom_only,
                y_new_pathway = y_new_path_only,
                y_new_metabo = y_new_metabo_only,
                y_new_all = y_new_all_only,
                y_new_all_noas = y_new_all_nas_only,
                y_new_all_nclin = y_new_all_nclin_only)

gl_ftimp_long_basic <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/long/may_basic_plusbasic_gl_long_top_features.csv")
gl_ftimp_long_meta <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/long/may_basic_plusmeta_gl_long_top_features.csv")
gl_ftimp_long_meta_no_age_sex <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/long/meta_no_age_sex_gl_long_top_features.csv")
gl_ftimp_long_grs <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/long/grs_gl_long_top_features.csv")
gl_ftimp_long_taxa <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/long/may_basic_plustaxa_gl_long_top_features.csv")
gl_ftimp_long_micom <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/long/may_basic_plusmicom_gl_long_top_features.csv")
gl_ftimp_long_pathway <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/long/may_basic_pluspathway_gl_long_top_features.csv")
gl_ftimp_long_metabo <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/long/may_basic_plusmetabo_gl_long_top_features.csv")
gl_ftimp_long_all <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/long/may_basic_plusall_gl_long_top_features.csv")
gl_ftimp_long_all_no_age_sex <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/long/all_no_age_sex_gl_long_top_features.csv")
gl_ftimp_long_all_no_clin <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/long/all_no_clin_gl_long_top_features.csv")

# GLMMLASSO DELTA 
gl_lasso_delta <- read.csv('/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/delta/gl_delta_predictions_df.csv') %>% 
  dplyr::rename(y_new_basic = y_new_basic_only,
                y_new_meta = y_new_meta_only,
                y_new_meta_noas = y_new_meta_noas_only,
                y_new_grs = y_new_grs_only,
                y_new_taxa = y_new_tax_only,
                y_new_micom = y_new_micom_only,
                y_new_pathway = y_new_path_only,
                y_new_metabo = y_new_metab_only,
                y_new_all = y_new_all_only, 
                y_new_all_noas = y_new_all_noas_only,
                y_new_all_nclin = y_new_all_nclin_only)

gl_ftimp_delta_basic <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/delta/basic_gl_delta_top_features.csv")
gl_ftimp_delta_meta <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/delta/meta_gl_delta_top_features.csv")
gl_ftimp_delta_meta_no_age_sex <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/delta/meta_no_age_sex_gl_delta_top_features.csv")
gl_ftimp_delta_grs <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/delta/grs_gl_delta_top_features.csv")
gl_ftimp_delta_taxa <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/delta/taxa_gl_delta_top_features.csv")
gl_ftimp_delta_micom <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/delta/micom_gl_delta_top_features.csv")
gl_ftimp_delta_pathway <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/delta/pathway_gl_delta_top_features.csv")
gl_ftimp_delta_metabo <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/delta/metabo_gl_delta_top_features.csv")
gl_ftimp_delta_all <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/delta/all_gl_delta_top_features.csv")
gl_ftimp_delta_all_no_age_sex <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/delta/all_no_age_sex_gl_delta_top_features.csv")
gl_ftimp_delta_all_no_clin <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/may_basic_plus/delta/all_no_clin_gl_delta_top_features.csv")

### Read in MERF SHAP data ####################################################################################

long_merf_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/may_basic_plus/long/"

shap_dfs <- list(
  shap_long_basic  = read_csv(file.path(long_merf_dir, "shap_long_basic_BEST_MSE_Model.csv")),
  shap_long_meta   = read_csv(file.path(long_merf_dir, "shap_long_meta_keep_MSE_Model.csv")),
  shap_long_meta_no_age_sex   = read_csv(file.path(long_merf_dir, "shap_long_meta_keep_no_age_sex_MSE_Model.csv")),
  shap_long_grs    = read_csv(file.path(long_merf_dir, "shap_long_only_grs_MSE_Model.csv")),
  shap_long_taxa   = read_csv(file.path(long_merf_dir, "shap_long_only_taxa_MSE_Model.csv")),
  shap_long_pathway= read_csv(file.path(long_merf_dir, "shap_long_only_pathway_MSE_Model.csv")),
  shap_long_micom  = read_csv(file.path(long_merf_dir, "shap_long_only_micom_MSE_Model.csv")),
  shap_long_metabo = read_csv(file.path(long_merf_dir, "shap_long_only_metabo_MSE_Model.csv")),
  shap_long_all    = read_csv(file.path(long_merf_dir, "shap_long_only_all_MSE_Model.csv")),
  shap_long_all_no_age_sex    = read_csv(file.path(long_merf_dir, "shap_long_only_all_no_age_sex_MSE_Model.csv")),
  shap_long_all_no_clin    = read_csv(file.path(long_merf_dir, "shap_long_only_all_no_clin_MSE_Model.csv")))

delta_merf_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/may_basic_plus/delta"

shap_delta_dfs <- list(
  shap_delta_basic  = read_csv(file.path(delta_merf_dir, "shap_long_basic_MSE_Model.csv")),
  shap_delta_meta   = read_csv(file.path(delta_merf_dir, "shap_long_meta_keep_MSE_Model.csv")),
  shap_delta_meta_no_age_sex   = read_csv(file.path(delta_merf_dir, "shap_long_meta_keep_no_age_sex_MSE_Model.csv")),
  shap_delta_grs    = read_csv(file.path(delta_merf_dir, "shap_long_only_grs_MSE_Model.csv")),
  shap_delta_taxa   = read_csv(file.path(delta_merf_dir, "shap_long_only_taxa_MSE_Model.csv")),
  shap_delta_pathway= read_csv(file.path(delta_merf_dir, "shap_long_only_pathway_MSE_Model.csv")),
  shap_delta_micom  = read_csv(file.path(delta_merf_dir, "shap_long_only_micom_MSE_Model.csv")),
  shap_delta_metabo = read_csv(file.path(delta_merf_dir, "shap_long_only_metabo_MSE_Model.csv")),
  shap_delta_all    = read_csv(file.path(delta_merf_dir, "shap_long_only_all_MSE_Model.csv")),
  shap_delta_all_no_age_sex    = read_csv(file.path(delta_merf_dir, "shap_long_only_all_no_age_sex_MSE_Model.csv")),
  shap_delta_all_no_clin    = read_csv(file.path(delta_merf_dir, "shap_long_only_all_no_clin_MSE_Model.csv")))


###### RENAME VARIABLES ####################################################################################

# For GlmLasso
gl_ftimp_delta_meta <- gl_ftimp_delta_meta %>%
  dplyr::mutate(
    Feature = dplyr::recode(
      Feature,
      "Ldl" = "LDL",
      "Randomized group1" = "IMF Diet Grp.",
      "Sex" = "Sex (Male)",
      "Hdl" = "HDL",
      "Glucose x" = "Glucose",
      "Homo ir" = "Homo-IR"))

gl_ftimp_delta_grs <- gl_ftimp_delta_grs %>%
  dplyr::mutate(
    Feature = dplyr::recode(
      Feature,
      "Bmi prs" = "GRS"))

gl_ftimp_long_grs <- gl_ftimp_long_grs %>%
  dplyr::mutate(
    Feature = dplyr::recode(
      Feature,
      "Bmi prs" = "GRS"))

gl_ftimp_delta_pathway <- gl_ftimp_delta_pathway %>%
  dplyr::mutate(
    Feature = dplyr::recode(
      Feature,
      "Superpathway of sulfolactate degradation" = "Sulfolactate degr.",
      "Superpathway of l phenylalanine biosynthesis" = "L-Phe biosynthesis",
      "Urate biosynthesis inosine 5  Phosphate degradation" = "Urate-IMP pathway",
      "Isoprene biosynthesis ii Engineered" = "Isoprene biosynth. II",
      "Mixed acid fermentation" = "Mixed acid ferm.",
      "Tca cycle vii Acetate producers" = "TCA VII (acetate)",
      "Hexitol fermentation to lactate Formate Ethanol and acetate" = "Hexitol -> Lac/For/Eth/Ace",
      "Superpathway of R r Butanediol biosynthesis" = "(R,R)-Butanediol biosynth.",
      "Heme biosynthesis ii Anaerobic" = "Heme biosynth. II (anaerobic)"))

gl_ftimp_delta_metabo <- gl_ftimp_delta_metabo %>%
  dplyr::mutate(
    Feature = dplyr::recode(
      Feature,
      "Ldl" = "LDL",
      "Glyca" = "Glycoprotein acetyls",
      "Vldl size" = "VLDL-size",
      "Hdl size" = "HDL-size",
      "Gly" = "Glycine",
      "Tyr" = "Tyrosine",
      "Phe" = "Phenylalanine",
      "Vldl tg" = "TGs in VLDL",
      "Tg by pg" = "TGs/Phosphoglycerides",
      "Dha pct" = "DHA (% of FAs)",
      "Vldl l" = "Lipids in VLDL"))

gl_ftimp_delta_all <- gl_ftimp_delta_all %>%
  dplyr::mutate(
    Feature = dplyr::recode(
      Feature,
      "Superpathway of sulfolactate degradation" = "Sulfolactate degr.",
      "Superpathway of l phenylalanine biosynthesis" = "L-Phe biosynthesis",
      "Ldl" = "LDL",
      "Glyca" = "Glycoprotein acetyls",
      "Vldl size" = "VLDL-size",
      "Hdl size" = "HDL-size",
      "Gly" = "Glycine",
      "Vldl tg" = "TGs in VLDL",
      "Tg by pg" = "TGs/Phosphoglycerides",
      "Dha pct" = "DHA (% of FAs)",
      "Vldl l" = "Lipids in VLDL"))

gl_ftimp_long_meta <- gl_ftimp_long_meta %>%
  dplyr::mutate(
    Feature = dplyr::recode(
      Feature,
      "Ldl" = "LDL",
      "Randomized group1" = "IMF Diet Grp.",
      "Sex" = "Sex (Male)",
      "Hdl" = "HDL",
      "Glucose x" = "Glucose",
      "Homo ir" = "Homo-IR"))

gl_ftimp_long_metabo <- gl_ftimp_long_metabo %>%
  dplyr::mutate(
    Feature = dplyr::recode(
      Feature,
      "Ldl" = "LDL",
      "Glyca" = "Glycoprotein acetyls",
      "Bohbutyrate" = "3-Hydroxybutyrate",
      "Ldl size" = "LDL-size",
      "Hdl size" = "HDL-size",
      "Gly" = "Glycine",
      "Gln" = "Glutamine",
      "Vldl tg" = "TGs in VLDL",
      "Tyr" = "Tyrosine",
      "Phe" = "Phenylalanine",
      "Tg by pg" = "TGs/Phosphoglycerides",
      "Dha pct" = "DHA (% of FAs)",
      "Vldl l" = "Lipids in VLDL"))


## MERF RENAMING ####################################################################################

# Define a renaming dictionary
feature_rename_basic <- c(
  "age" = "Age",
  "randomized_group" = "IMF Diet Grp.",
  "sex" = "Sex (Male)")

feature_rename_grs <- c(
  "bmi_prs" = "GRS")

feature_rename_meta <- c(
  "Triglyceride_lipid" = "Triglycerides",
  "HDL_Total_Direct_lipid" = "HDL",
  "LDL_Calculated" = "LDL",
  "Glucose.x" = "Glucose",
  "Insulin_endo" = "Insulin",
  "HOMA_IR" = "Homo-IR",
  "race" = "Race",
  "age" = "Age",
  "sex" = "Sex (Male)",
  "randomized_group" = "IMF Diet Grp.")

feature_rename_metabo <- c(
  "PUFA_by_MUFA" = "PUFA:MUFA ratio",
  "GlycA" = "Glycoprotein Acetyls",
  "DHA_pct" = "DHA (% of FAs)",
  "Gly" = "Glycine",
  "Ala" = "Alanine",
  "Glucose.y" = "Glucose",
  "Gln" = "Glutamine-IR",
  "HDL_CE" = "Chol. esters in HDL",
  "DHA_pct" = "DHA (% of FAs)",
  "IDL_CE_pct" = "Cholesteryl Esters (% of IDL Lipids)",
  "PUFA_pct" = "PUFA (%)",
  "HDL_size" = "HDL-size",
  "VLDL_size" = "VLDL-size",
  "VLDL_TG" = "TGs in VLDL",
  "Phe" = "Phenylalanine", 
  "His" = "Histidine", 
  "LA" = "Linoleic Acid (% of FAs)")

feature_rename_taxa <- c(
  "g__Anaerotruncus" = "G Anaerotruncus",
  "g__Anaerobutyricum" = "G Anaerobutyricum",
  "g__Anaerotignum_189125" = "G Anaerotignum 189125",
  "g__CAG-127" = "G CAG-127",
  "g__Clostridium_Q_135822" = "G Clostridium Q 135822",
  "g__Collinsella" = "G Collinsella",
  "g__Coprenecus" = "G Coprenecus",
  "g__Dorea_A" = "G Dorea A",
  "g__Enterocloster" = "G Enterocloster",
  "g__Eubacterium_I" = "G Eubacterium I",
  "g__Faecalibacillus" = "G Faecalibacillus",
  "g__Faecousia" = "G Faecousia",
  "g__Frisingicoccus" = "G Frisingicoccus",
  "g__Haemophilus_D_735815" = "G Haemophilus D 735815",
  "g__Holdemanella" = "G Holdemanella",
  "g__Lawsonibacter" = "G Lawsonibacter",
  "g__Acetatifactor" = "G Acetatifactor",
  "g__Agathobacter_164117" = "G Agathobacter 164117",
  "g__Agathobaculum" = "G Agathobaculum",
  "g__Akkermansia" = "G Akkermansia",
  "g__Alistipes_A_871400" = "G Alistipes A 871400",
  "g__BX12" = "G BX12",
  "g__Bifidobacterium_388775" = "G Bifidobacterium 388775",
  "g__CAG-127" = "G CAG-127",
  "g__CAG-274" = "G CAG-274",
  "g__Clostridium_Q_135822" = "G Clostridium Q 135822",
  "g__Collinsella" = "G Collinsella",
  "g__Coprococcus_A_187866" = "G Coprococcus A 187866",
  "g__Copromonas" = "G Copromonas",
  "g__Eggerthella" = "G Eggerthella",
  "g__Eisenbergiella" = "G Eisenbergiella",
  "g__Streptococcus" = "G Streptococcus",
  "g__UBA1417" = "g__UBA1417",
  "g__Intestinibacter" = "G Intestinibacter",
  "g__Ruminiclostridium_E" = "G Ruminiclostridium E",
  "g__Ruminococcus_D" = "G Ruminococcus D", 
  "g___1" = "G 1",
  "g__Porcincola" = "G Porcincola",
  "g__Ruminococcus_C_58660" = "G Ruminococcus C58660",
  "g__Longicatena" = "G Longicatena",
  "g__Ruminococcus_D" = "G Ruminococcus D",
  "g__Marvinbryantia" = "G Marvinbryantia"
)

feature_rename_pathway <- c(
  "sulfate.reduction.I..assimilatory." = "Sulfate reduction assimilation",
  "superpathway.of.phylloquinol.biosynthesis" = "Phylloquinol biosyn. pathways",
  "D.glucarate.degradation.I" = "D-glucarate degr. ",
  "L.lysine.biosynthesis.II" = "L-lysine biosyn. II",
  "L.lysine.fermentation.to.acetate.and.butanoate" = "Lys ferm. to acetate/butanoate",
  "L.rhamnose.degradation.I" = "L-rhamnose degr. I",
  "TCA.cycle.VI..obligate.autotrophs. " = "TCA VI (autotrophs)",
  "aromatic.biogenic.amine.degradation..bacteria." = "Aromatic amine degr. (bacteria)",
  "glucose.and.glucose.1.phosphate.degradation" = "Glucose/G1P degr.",
  "incomplete.reductive.TCA.cycle" = "Incomplete reductive TCA",
  "isoprene.biosynthesis.II..engineered." = "Isoprene biosynth. II (eng.)",
  "isopropanol.biosynthesis" = "Isopropanol biosyn.",
  "Bifidobacterium.shunt" = "Bifidobacterium shunt",
  "GDP.mannose.biosynthesis" = "GDP mannose biosyn.",
  "L.histidine.degradation.I" = "L-histadine. degr.I ",
  "UDP.2.3.diacetamido.2.3.dideoxy..alpha..D.mannuronate.biosynthesis" = "UDP-DDDA-ManA biosyn.",
  "X1.4.dihydroxy.6.naphthoate.biosynthesis.I " = "1,4 dihydroxy naphthoate biosyn.I",
  "X1.4.dihydroxy.6.naphthoate.biosynthesis.II" = "1,4 dihydroxy naphthoate biosyn.II",
  "allantoin.degradation.to.glyoxylate.III" = "Allantoin -> Glyoxylate III",
  "arginine..ornithine.and.proline.interconversion" = "Arg Orn Pro interconversion",
  "glycolysis.V..Pyrococcus." = "Glycolysis V Pyro.",
  "heme.biosynthesis.II..anaerobic." = "Heme biosyn.II anaerobic",
  "lactose.and.galactose.degradation.I" = "lactose/galactose degr.I",
  "pyrimidine.deoxyribonucleotides.de.novo.biosynthesis.II" = "Pyr dNTP biosyn.II",
  "hexitol.fermentation.to.lactate..formate..ethanol.and.acetate" = "Hexitol Fermentation",
  "X4.deoxy.L.threo.hex.4.enopyranuronate.degradation" = "Deoxy enopyranuronate degr.",
  "TCA.cycle.VI..obligate.autotrophs." = "TCA VI obligate autotrophs",
  "reductive.acetyl.coenzyme.A.pathway" = "Red. Acetyl-CoA Pathway",
  "superpathway.of.sulfur.oxidation..Acidianus.ambivalens." = "Acidianus Ambivalens (Sulfur Oxi.)",
  "peptidoglycan.biosynthesis.V...beta..lactam.resistance." = "Peptidoglycan Biosyn.V",
  "superpathway.of.2.3.butanediol.biosynthesis" = "2,3 Butanediol Biosynthesis Pathways",
  "superpathway.of..R.R..butanediol.biosynthesis" = "R,R Butanediol Biosynthesis Pathway",
  "NAD.salvage.pathway.II" = "NAD Salvage Pathway II",
  "inosine.5..phosphate.biosynthesis.III" = "Inosine-5-phosphate Biosyn.III",
  "fatty.acid.salvage" = "Fatty Acid Salvage",
  "ketogluconate.metabolism" = "Ketogluconate Metabolism",
  "superpathway.of.sulfolactate.degradation" = "Sulfolactate Degradation Pathways", 
  "superpathway.of.UDP.N.acetylglucosamine.derived.O.antigen.building.blocks.biosynthesis" = "UDPN Acetylglucosamine-O-antigen Biosyn.",
  "mono.trans..poly.cis.decaprenyl.phosphate.biosynthesis" = "Mono-trans-poly-cis-decaprenyl-phosphate Biosyn.",
  "thiazole.biosynthesis.II..Bacillus." = "Thiazole Biosyn.II Bacillus.",
  "L.histidine.degradation.II" = "L-histidine Degradation II",
  "superpathway.of.L.phenylalanine.biosynthesis" = "L-phenylalanine Biosyn. Pathways",
  "superpathway.of.hexitol.degradation..bacteria." = "Bacterial Hexitol Degradation Pathways",
  "UDP.2.3.diacetamido.2.3.dideoxy..alpha..D.mannuronate.biosynthesis" = "UDP-DDDA-ManA biosyn.",
  "X1.4.dihydroxy.6.naphthoate.biosynthesis.I " = "1,4 dihydroxy naphthoate biosyn.I",
  "X1.4.dihydroxy.6.naphthoate.biosynthesis.II" = "1,4 dihydroxy naphthoate biosyn.II",
  "allantoin.degradation.to.glyoxylate.III" = "Allantoin -> Glyoxylate III",
  "arginine..ornithine.and.proline.interconversion" = "Arg Orn Pro interconversion",
  "glycolysis.V..Pyrococcus." = "Glycolysis V Pyro.",
  "heme.biosynthesis.II..anaerobic." = "heme biosyn.II anaerobic",
  "lactose.and.galactose.degradation.I" = "lactose/galactose degr.I",
  "pyrimidine.deoxyribonucleotides.de.novo.biosynthesis.II" = "Pyr dNTP biosyn.II",
  "reductive.acetyl.coenzyme.A.pathway" = "Red. Acetyl-CoA Pathway",
  "superpathway.of.histidine..purine..and.pyrimidine.biosynthesis" = "His. Purine & Pyrimidine Biosyn.",
  "superpathway.of.pyridoxal.5..phosphate.biosynthesis.and.salvage" = "Pyridoxal-5-phosphate Biosyn. & Salvage",
  "purine.nucleotides.degradation.II..aerobic." = "Aerobic Purine Degradation II aerobic.",
  "superpathway.of.L.methionine.biosynthesis..by.sulfhydrylation." = "Sulfhydrylation -> L-methionine Biosyn.",
  "myo.inositol.degradation.I" = "Myo-inositol Degradation I"
)

feature_rename_all <- c(
  "reductive.acetyl.coenzyme.A.pathway" = "Red. Acetyl-CoA Pathway",
  "IDL_CE_pct" = "Cholesteryl Esters (% of IDL Lipids)",
  "aromatic.biogenic.amine.degradation..bacteria." = "Aromatic amine degr. (bacteria)",
  "Dephospho.CoA" = "Dephospho-CoA",
  "GlycA" = "Glycoprotein Acetyls",
  "DHA_pct" = "DHA (% of FAs)",
  "Gly" = "Glycine",
  "Ala" = "Alanine",
  "Glucose.y" = "Glucose",
  "Gln" = "Glutamine-IR",
  "X3.methoxy.acetaminophen" = "3-Methoxy Acetaminophen",
  "DHA_pct" = "DHA (% of FAs)",
  "aldehydo.D.xylose" = "Ald. D-xylose",
  "PUFA_pct" = "PUFA (%)",
  "VLDL_size" = "VLDL-size",
  "L.fucose" = "L-fucose",
  "Thiamin.monophosphate" = "Thiamin monophos.",
  "Insulin_endo" = "Insulin",
  "HOMA_IR" = "Homo-IR",
  "MUFA_pct" = "MUFA (% of FAs)",
  "LA" = "Linoleic acid", 
  "Omega_6" = "Omega6 FA",
  "TCA.cycle.VI..obligate.autotrophs." = "TCA VI (autotrophs)",
  "superpathway.of.sulfur.oxidation..Acidianus.ambivalens." = "Acidianus Ambivalens (Sulfur Oxi.)",
  "NAD.salvage.pathway.II" = "NAD Salvage Pathway II", 
  "g__Acetatifactor" = "G Acetatifactor",
  "indole.3.acetate" = "Indole-3-acetate")

feature_rename_grs_lg <- feature_rename_grs
feature_rename_basic_lg <- feature_rename_basic
feature_rename_meta_lg <- feature_rename_meta
feature_rename_metabo_lg <- feature_rename_metabo
feature_rename_pathway_lg <- feature_rename_pathway
feature_rename_taxa_lg <- feature_rename_taxa
feature_rename_all_lg <- feature_rename_all

feature_rename_grs_delta <- feature_rename_grs
feature_rename_basic_delta <- feature_rename_basic
feature_rename_meta_delta <- feature_rename_meta
feature_rename_metabo_delta <- feature_rename_metabo
feature_rename_pathway_delta <- feature_rename_pathway
feature_rename_taxa_delta <- feature_rename_taxa
feature_rename_all_delta <- feature_rename_all

######### RENAME SHAPS ####################################################################################

rename_maps <- list(
  shap_long_basic = feature_rename_basic,
  shap_long_grs = feature_rename_grs,
  shap_long_meta     = feature_rename_meta,
  shap_long_metabo   = feature_rename_metabo,
  shap_long_taxa     = feature_rename_taxa,
  shap_long_pathway  = feature_rename_pathway,
  shap_long_all = feature_rename_all)

# shap_dfs <- lapply(names(shap_dfs), function(nm) {
#   if (nm %in% names(rename_maps)) {
#     shap_dfs[[nm]] %>%
#       dplyr::mutate(Feature = recode(Feature, !!!rename_maps[[nm]]))
#   } else {
#     shap_dfs[[nm]]
#   }
# }) %>% setNames(names(shap_dfs))

# shap_dfs <- lapply(names(shap_dfs), function(nm) {
#   if (nm %in% names(rename_maps)) {
#     shap_dfs[[nm]] %>%
#       dplyr::mutate(Feature = dplyr::recode(Feature, .replace = rename_maps[[nm]]))
#   } else {
#     shap_dfs[[nm]]
#   }
# }) %>% setNames(names(shap_dfs))

shap_dfs <- lapply(names(shap_dfs), function(nm) {
  df <- shap_dfs[[nm]]
  map <- rename_maps[[nm]]
  
  if (!is.null(map) && "Feature" %in% names(df)) {
    current_features <- unique(df$Feature)
    matched_keys <- intersect(current_features, names(map))
    
    # Only apply recode if there are matches
    if (length(matched_keys) > 0) {
      valid_map <- map[matched_keys]
      df <- df %>%
        dplyr::mutate(Feature = dplyr::recode(Feature, !!!as.list(valid_map)))
    } else {
      message(glue::glue("No matching features to rename in {nm}"))
    }
  }
  return(df)
}) %>% setNames(names(shap_dfs))


# and for delta shaps 
rename_maps_delta <- list(
  shap_delta_basic = feature_rename_basic,
  shap_delta_grs = feature_rename_grs,
  shap_delta_meta     = feature_rename_meta,
  shap_delta_metabo   = feature_rename_metabo,
  shap_delta_taxa     = feature_rename_taxa,
  shap_delta_pathway  = feature_rename_pathway,
  shap_delta_all = feature_rename_all)
# 
# shap_delta_dfs <- lapply(names(shap_delta_dfs), function(nm) {
#   if (nm %in% names(rename_maps_delta)) {
#     shap_delta_dfs[[nm]] %>%
#       dplyr::mutate(Feature = recode(Feature, !!!rename_maps_delta[[nm]]))
#   } else {
#     shap_delta_dfs[[nm]]
#   }
# }) %>% setNames(names(shap_delta_dfs))


shap_delta_dfs <- lapply(names(shap_delta_dfs), function(nm) {
  df <- shap_delta_dfs[[nm]]
  map <- rename_maps_delta[[nm]]
  
  if (!is.null(map) && "Feature" %in% names(df)) {
    current_features <- unique(df$Feature)
    matched_keys <- intersect(current_features, names(map))
    
    if (length(matched_keys) > 0) {
      valid_map <- map[matched_keys]
      df <- df %>%
        dplyr::mutate(Feature = dplyr::recode(Feature, !!!as.list(valid_map)))
    } else {
      message(glue::glue("No matching features to rename in {nm}"))
    }
  }
  
  return(df)
}) %>% setNames(names(shap_delta_dfs))


####################################################################################

# Apply renaming function to the column MERF Long
meta <- meta %>%
  mutate(Top_15_Feature_Importances = map_chr(Top_15_Feature_Importances, 
                                              ~rename_features_in_json(.x, feature_rename_meta)))
#grs <- grs %>%
#  mutate(Top_15_Feature_Importances = map_chr(Top_15_Feature_Importances, 
#                                              ~rename_features_in_json(.x, feature_rename_grs)))

metabo <- metabo %>%
  mutate(Top_15_Feature_Importances = map_chr(Top_15_Feature_Importances, 
                                              ~rename_features_in_json(.x, feature_rename_metabo_lg)))
taxa <- taxa %>%
  mutate(Top_15_Feature_Importances = map_chr(Top_15_Feature_Importances, 
                                              ~rename_features_in_json(.x, feature_rename_taxa_lg)))
pathway <- pathway %>%
  mutate(Top_15_Feature_Importances = map_chr(Top_15_Feature_Importances, 
                                              ~rename_features_in_json(.x, feature_rename_pathway_lg)))
all <- all %>%
  mutate(Top_15_Feature_Importances = map_chr(Top_15_Feature_Importances, 
                                              ~rename_features_in_json(.x, feature_rename_all)))
# Apply to the column MERF Delta
meta_md <- meta_md %>%
  mutate(Top_15_Feature_Importances = map_chr(Top_15_Feature_Importances, 
                                              ~rename_features_in_json(.x, feature_rename_meta)))
#grs_md <- grs_md %>%
# mutate(Top_15_Feature_Importances = map_chr(Top_15_Feature_Importances, 
#                                              ~rename_features_in_json(.x, feature_rename_grs)))
metabo_md <- metabo_md %>%
  mutate(Top_15_Feature_Importances = map_chr(Top_15_Feature_Importances, 
                                              ~rename_features_in_json(.x, feature_rename_metabo)))
taxa_md <- taxa_md %>%
  mutate(Top_15_Feature_Importances = map_chr(Top_15_Feature_Importances, 
                                              ~rename_features_in_json(.x, feature_rename_taxa)))
pathway_md <- pathway_md %>%
  mutate(Top_15_Feature_Importances = map_chr(Top_15_Feature_Importances, 
                                              ~rename_features_in_json(.x, feature_rename_pathway)))
all_md <- all_md %>%
  mutate(Top_15_Feature_Importances = map_chr(Top_15_Feature_Importances, 
                                              ~rename_features_in_json(.x, feature_rename_all)))

####################################################################################



