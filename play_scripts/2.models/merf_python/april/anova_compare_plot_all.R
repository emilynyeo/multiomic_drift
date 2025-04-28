rm(list = ls())
#source("zc_functions.R") 
#install.packages("rsq")  
library(rsq)
library(pacman)
p_load(tools, reticulate, viridis, tidyplots, patchwork, jsonlite, maps, ggvenn, 
       caret, caretEnsemble, glmnet, xgboost, ggplot2, glmmLasso, corrplot,
       readr, plyr, dplyr, tidyr, purrr, tibble, stringr, psych, randomForest,  
       reshape2, scales, gridExtra, plotly, sf, tidyverse, naniar, VIM)
#library(nlme)
library(gridExtra)
library(sjPlot)
library(htmltools)
library(officer)
library(flextable)
library(webshot)
library(apaTables)
library(MuMIn)
library(lme4)  #
library(ggplot2)
library(glue)  # if using glue()
library(grid)
'%ni%' <- Negate('%in%')
r2_general <-function(preds,actual){ 
  return(1- sum((preds - actual) ^ 2)/sum((actual - mean(actual))^2))
}
q_squared <- function(actual, predicted) {
  ss_res <- sum((actual - predicted)^2)
  ss_tot <- sum((actual - mean(actual))^2)
  1 - (ss_res / ss_tot)
}

# ANOVA & Feat Imp Comparisons 
omic_out_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/april/anova_results/plots/"
# MERF LONG 
predict_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/april/final_merf_dfs"
merf_long <- read.csv('/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/april/anova_results/april_long_predictions_df.csv') %>% 
             dplyr::filter(Model == "MSE Model")

basic <- read.csv(file.path(predict_dir, "basic_long.csv")) %>% dplyr::rename(y_new_basic = y_hat_new)
meta <- read.csv(file.path(predict_dir, "meta_keep_long.csv")) %>% dplyr::rename(y_new_meta = y_hat_new)
grs <- read.csv(file.path(predict_dir, "only_grs_long.csv")) %>% dplyr::rename(y_new_grs = y_hat_new)
taxa <- read.csv(file.path(predict_dir, "only_taxa_long.csv")) %>% dplyr::rename(y_new_taxa = y_hat_new)
pathway <- read.csv(file.path(predict_dir, "only_pathway_long.csv")) %>% dplyr::rename(y_new_pathway = y_hat_new)
micom <- read.csv(file.path(predict_dir, "only_micom_long.csv")) %>% dplyr::rename(y_new_micom = y_hat_new)
metabo <- read.csv(file.path(predict_dir, "only_metabo_long.csv")) %>% dplyr::rename(y_new_metabo = y_hat_new)

# MERF DELTA
merf_delta <- read.csv('/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/april/anova_results/april_delta_predictions_df.csv') %>% 
              dplyr::filter(Model == "MSE Model")

basic_md <- read.csv(file.path(predict_dir, "basic_delta.csv")) %>% dplyr::rename(y_new_basic = y_hat_new)
meta_md <- read.csv(file.path(predict_dir, "meta_keep_delta.csv")) %>% dplyr::rename(y_new_meta = y_hat_new)
grs_md <- read.csv(file.path(predict_dir, "only_grs_delta.csv")) %>% dplyr::rename(y_new_grs = y_hat_new)
taxa_md <- read.csv(file.path(predict_dir, "only_taxa_delta.csv")) %>% dplyr::rename(y_new_taxa = y_hat_new)
pathway_md <- read.csv(file.path(predict_dir, "only_pathway_delta.csv")) %>% dplyr::rename(y_new_pathway = y_hat_new)
micom_md <- read.csv(file.path(predict_dir, "only_micom_delta.csv")) %>% dplyr::rename(y_new_micom = y_hat_new)
metabo_md <- read.csv(file.path(predict_dir, "only_metabo_delta.csv")) %>% dplyr::rename(y_new_metabo = y_hat_new)

# GLMMLASSO LONG
gl_lasso_long <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_long/april_long_predictions_df.csv") %>% 
  dplyr::rename(y_new_meta = y_new_meta_only,
                y_new_taxa = y_new_tax_only,
                y_new_micom = y_new_micom_only,
                y_new_pathway = y_new_path_only,
                y_new_metabo = y_new_metabo_only)

gl_ftimp_long_basic <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_long/basic_gl_long_top_features.csv")
gl_ftimp_long_meta <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_long/meta_gl_long_top_features.csv")
gl_ftimp_long_taxa <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_long/taxa_gl_long_top_features.csv")
gl_ftimp_long_micom <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_long/micom_gl_long_top_features.csv")
gl_ftimp_long_pathway <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_long/pathway_gl_long_top_features.csv")
gl_ftimp_long_metabo <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_long/metabo_gl_long_top_features.csv")

# GLMMLASSO DELTA 
gl_lasso_delta <- read.csv('/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_delta/april_delta_predictions_df.csv') %>% 
  dplyr::rename(y_new_basic = y_new_basic_only,
                y_new_meta = y_new_meta_only,
                y_new_taxa = y_new_tax_only,
                y_new_micom = y_new_micom_only,
                y_new_pathway = y_new_path_only,
                y_new_metabo = y_new_metab_only)

gl_ftimp_delta_basic <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_delta/basic_gl_delta_top_features.csv")
gl_ftimp_delta_meta <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_delta/meta_gl_delta_top_features.csv")
gl_ftimp_delta_taxa <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_delta/taxa_gl_delta_top_features.csv")
gl_ftimp_delta_micom <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_delta/micom_gl_delta_top_features.csv")
gl_ftimp_delta_pathway <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_delta/pathway_gl_delta_top_features.csv")
gl_ftimp_delta_metabo <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_delta/metabo_gl_delta_top_features.csv")

### PREDICTED vs ACTUAL BMI 
omic_dfs <- list(
  merf_long = merf_long,
  gl_lasso_long = gl_lasso_long,
  merf_delta = merf_delta, 
  gl_lasso_delta = gl_lasso_delta)

for (omic_name in names(omic_dfs)) {
  omic <- omic_dfs[[omic_name]]
  basic_plot <- ggplot(omic, aes(x = bmi, y = y_new_meta, color = Cluster, shape = as.factor(Time))) +
    geom_point(alpha = 3.3, size = 3) +
    geom_smooth(method = "lm", se = FALSE) +  # optional trend lines
    ylim(23, 52) +
    labs(title = paste0(omic_name, " Predicted BMI Values vs. Actual BMI"),
         x = "Actual BMI",
         y = "Predicted BMI",
         color = "Predictor") +
    theme_minimal()
  basic_plot
  assign(paste0(omic_name, "_basic_plot"), basic_plot)
  
  # Pivot to long format for all omic predictions on test set 
  pred_data <- omic 
  gl_long <- pred_data %>% pivot_longer(cols = starts_with("y_new"),
                                        names_to = "predictor_type", 
                                        values_to = "prediction") %>% 
    mutate(Time = as.factor(Time)) %>% as.data.frame()
  
  long_plot <- ggplot(gl_long, aes(x = bmi, y = prediction, color = predictor_type)) +
    geom_point(alpha = 0.6, size = 3) +
    geom_smooth(method = "lm", se = FALSE) +  # optional trend lines
    #ylim(23, 52) +
    labs(title = paste0(omic_name, " Predicted BMI Values vs. Actual BMI"),
         x = "Actual Test Set BMI",
         y = "Omic Risk Score Predicted BMI",
         color = "Predictor") +
    theme_minimal()
  long_plot
  assign(paste0(omic_name, "_each_omic_plot"), long_plot)
  # Save the long plot
  ggsave(file.path(omic_out_dir, paste0(omic_name, "_each_omic_plot.png")), plot = long_plot, width = 8, height = 6)
  
  ### ANOVA 
  # Single basic model plus each omic 
  mod_dat <- omic
  lmer_basic <- lmer(bmi ~ y_new_basic + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_meta_b <- lmer(bmi ~ y_new_basic + y_new_meta + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_micom_b <- lmer(bmi ~ y_new_basic + y_new_micom  + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_path_b <- lmer(bmi ~ y_new_basic + y_new_pathway + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_tax_b <- lmer(bmi ~ y_new_basic + y_new_taxa + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_metabo_b <- lmer(bmi ~ y_new_basic + y_new_metabo + (1|Cluster), data = mod_dat, REML = FALSE)
  
  # Put your models into a named list
  models <- list(
    Basic   = lmer_basic,
    `Basic + Meta`  = lmer_meta_b,
    `Basic + Micom` = lmer_micom_b,
    `Basic + Pathway` = lmer_path_b,
    `Basic + Taxa`    = lmer_tax_b,
    `Basic + Metabo`  = lmer_metabo_b)
  
  # Extract R2 values for each model
  r2_df <- do.call(rbind, lapply(names(models), function(name) {
    r2_vals <- r.squaredGLMM(models[[name]])
    data.frame(Model = name,
      R2m = r2_vals[1, "R2m"],
      R2c = r2_vals[1, "R2c"])
  }))
  print(r2_df)
  assign(paste0(omic_name, "_r2_df"), r2_df)
  
  glmmlass_lmer_models <- list(
    c("lmer_basic", "lmer_meta_b"),
    c("lmer_basic", "lmer_tax_b"),
    c("lmer_basic", "lmer_micom_b"),
    c("lmer_basic", "lmer_path_b"),
    c("lmer_basic", "lmer_metabo_b"))
  
  anova_results <- list()  # empty list to store ANOVA results
  all_anova_tables <- list()
  for (model_pair in glmmlass_lmer_models) {
    model_1 <- get(model_pair[1])
    model_2 <- get(model_pair[2])
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
  
  print(anova_table_clean)
  all_anova_tables[[omic_name]] <- anova_table_clean
  
  anova_labels <- anova_table_clean %>%
    dplyr::mutate(
      Compared_Model = dplyr::recode(
        model_comparison,
        "lmer_basic vs lmer_meta_b" = "Basic + Meta",
        "lmer_basic vs lmer_tax_b" = "Basic + Taxa",
        "lmer_basic vs lmer_micom_b" = "Basic + Micom",
        "lmer_basic vs lmer_path_b" = "Basic + Pathway",
        "lmer_basic vs lmer_metabo_b" = "Basic + Metabo"
      ),
      Significance = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ ""
      )
    ) %>%
    dplyr::select(Compared_Model, p.value, Significance)
  
  r2_annotated <- r2_df %>%
    left_join(anova_labels, by = c("Model" = "Compared_Model"))
  
  #r2_annotated <- r2_annotated %>%
  #      dplyr::mutate(Model = factor(Model, 
  #                    levels = c("Basic", r2_annotated %>% 
  #                               dplyr::filter(Model != "Basic") %>%
  #      arrange(R2m) %>% 
  #      pull(Model)))) %>%
  #      unique()
  
  model_order <- r2_annotated %>%
    dplyr::filter(Model != "Basic") %>%
    arrange(R2m) %>%
    pull(Model) %>%
    unique() #
  
  r2_annotated <- r2_annotated %>%
    dplyr::mutate(Model = factor(Model, levels = c("Basic", model_order)))
  
  sig_plot <- ggplot(r2_annotated, aes(x = Model, y = R2m, fill = Model)) +
    geom_col(width = 0.7) +
    geom_text(aes(y = R2m + 0.03, label = Significance), size = 18) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = c(
      "Basic + Meta"     = "#960018",
      "Basic + Taxa"     = "#8B3A3A",
      "Basic + Micom"    = "#996759",
      "Basic + Pathway"  = "#D8A39D",
      "Basic + Metabo"   = "#C74375",
      "Basic"            = "#582F2B")) +
    labs(title = glue("Marginal R² by Model – {omic_name}"),
      x = "Model Comparisons",
      y = "R2(m) = variance explained by the fixed effect") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none",
          plot.margin = ggplot2::margin(t = 2, r = 2, b = 2, l = 2, unit = "cm"),
      axis.text.y = element_text(face = "bold", hjust = 1, size = 14),
      axis.text.x = element_text(face = "bold", size = 14, hjust = 1, angle = 45),
      axis.title.x = element_text(margin = ggplot2::margin(1.5,0,0,0, unit = "cm"), size = 16),
      axis.title.y = element_text(margin = ggplot2::margin(0,1.5,0,0, unit = "cm"), size = 16))
  sig_plot
  assign(paste0(omic_name, "_sig_plot"), sig_plot)
  # Save the sig plot
  ggsave(file.path(omic_out_dir, paste0(omic_name, "_sig_plot.png")), plot = sig_plot, 
         width = 10, height = 8)
  
}

library(jsonlite)
library(purrr)

# Assuming `basic` is your dataframe
# Step 1: Parse the JSON strings
ft_dfs <- list(basic_lg = basic, meta_lg = meta, taxa_lg = taxa, 
          micom_lg = micom, pathway_lg = pathway, metabo_lg = metabo)

ft_dfs_delta <- list(basic_delta = basic_md, meta_delta = meta_md, taxa_delta = taxa_md, 
               micom_delta = micom_md, pathway_delta = pathway_md, metabo_delta = metabo_md)

ft_colors <- c(meta_lg = "#960018", taxa_lg = "#8B3A3A", micom_lg = "#996759",     
  pathway_lg = "#D8A39D", metabo_lg  = "#C74375", basic_lg = "#582F2B")

#ft_colors <- c(meta_delta = "#960018", taxa_delta = "#8B3A3A", micom_delta = "#996759",     
#               pathway_delta = "#D8A39D", metabo_delta  = "#C74375", basic_delta = "#582F2B")

# store plots if you want to print/save them later
plots <- list()
for (ft_name in names(ft_dfs)) {
  omic <- ft_dfs[[ft_name]]
  
  omic_parsed <- omic %>%
    dplyr::mutate(Top_15_Feature_Importances = gsub("'", '"', Top_15_Feature_Importances)) %>%
    dplyr::mutate(parsed = map(Top_15_Feature_Importances, fromJSON))
  
  omic_long <- omic_parsed %>%
    unnest(parsed) %>%
    dplyr::select(-Top_15_Feature_Importances) %>% filter(Feature != "time") 
  
  avg_importance <- omic_long %>%
    dplyr::group_by(Feature) %>%
    summarise(avg_importance = mean(Importance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Feature = reorder(Feature, avg_importance)) %>%
    head(10)

  merf_long_ft_imp <- ggplot(avg_importance, aes(x = Feature, y = avg_importance)) +
    geom_col(fill = ft_colors[[ft_name]], width = 0.7) +
    labs(title = paste("Avg. MERF Ft. Importance -", toupper(ft_name)),
      x = "Feature",
      y = "Average Importance") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none",
          plot.margin = ggplot2::margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
          plot.title = element_text(hjust = 1.4, face = "bold", size = 16),
          axis.text.y = element_text(face = "bold", hjust = 1, size = 14),
          axis.text.x = element_text(face = "bold", size = 14),
          axis.title.x = element_text(margin = ggplot2::margin(1,0,0,0, unit = "cm"), size = 16),
          axis.title.y = element_text(margin = ggplot2::margin(0,1,0,0, unit = "cm"), size = 16)) +
    coord_flip()

  plots[[ft_name]] <- merf_long_ft_imp # store or print
  
  ggsave(filename = file.path(omic_out_dir, paste0("ft_imp_merf_long/merf_avg_importance_", ft_name, ".png")),
   #ggsave(filename = file.path(omic_out_dir, paste0("ft_imp_merf_delta/merf_avg_importance_", ft_name, ".png")),     
          plot = merf_long_ft_imp, width = 9, height = 5, dpi = 300)
}


### GLMMLASSO Top features
lg_lass_ft_imp <- list(basic_lg = gl_ftimp_long_basic, meta_lg = gl_ftimp_long_meta, 
                       taxa_lg = gl_ftimp_long_taxa, micom_lg = gl_ftimp_long_micom, 
                       pathway_lg = gl_ftimp_long_pathway, metabo_lg = gl_ftimp_long_metabo)

#delta_lass_ft_imp <- list(basic_delta = gl_ftimp_delta_basic, meta_delta = gl_ftimp_delta_meta, 
#                       taxa_delta = gl_ftimp_delta_taxa, micom_delta = gl_ftimp_delta_micom, 
#                       pathway_delta = gl_ftimp_delta_pathway, metabo_delta = gl_ftimp_delta_metabo)

ft_colors <- c(meta_lg = "#960018", taxa_lg = "#8B3A3A", micom_lg = "#996759",     
               pathway_lg = "#D8A39D", metabo_lg  = "#C74375", basic_lg = "#582F2B")

#ft_colors <- c(meta_delta = "#960018", taxa_delta = "#8B3A3A", micom_delta = "#996759",     
#               pathway_delta = "#D8A39D", metabo_delta  = "#C74375", basic_delta = "#582F2B")

lass_plots <- list()
for (ft_name in names(lg_lass_ft_imp)) {
  top_features <- lg_lass_ft_imp[[ft_name]]
  
  plot_feat <- top_features %>%
    dplyr::filter(!str_detect(tolower(Feature), "time")) %>%  # Remove rows with "time" in Feature
    slice_head(n = 10)
  
  features <- ggplot(plot_feat, aes(x = reorder(Feature, Estimate), y = Estimate)) +
    geom_col(fill = ft_colors[[ft_name]], width = 0.7) +
    coord_flip() + theme_bw() +
    ggtitle(paste("Top Features & Coefficients from", ft_name)) +
    xlab("Feature") + ylab("Coefficient Estimate") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12)) +
    theme(legend.position = "none",
          plot.margin = ggplot2::margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
          plot.title = element_text(hjust = 1, face = "bold", size = 16),
          axis.text.y = element_text(face = "bold", hjust = 1, size = 14),
          axis.text.x = element_text(face = "bold", size = 14),
          axis.title.x = element_text(margin = ggplot2::margin(1,0,0,0, unit = "cm"), size = 16),
          axis.title.y = element_text(margin = ggplot2::margin(0,1,0,0, unit = "cm"), size = 16))
  
  lass_plots[[ft_name]] <- features # store or print
  
  ggsave(filename = file.path(omic_out_dir, paste0("ft_imp_glmmlasso/glm_avg_importance_", ft_name, ".png")),
         plot = features, width = 9, height = 5, dpi = 300)
}


