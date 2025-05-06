rm(list = ls())
#source("zc_functions.R") 
#install.packages("rsq")  
library(pacman)
p_load(tools, reticulate, viridis, tidyplots, patchwork, jsonlite, maps, ggvenn, 
       caret, caretEnsemble, glmnet, xgboost, ggplot2, glmmLasso, corrplot,
       readr, plyr, dplyr, tidyr, purrr, tibble, stringr, psych, randomForest,  
       reshape2, scales, gridExtra, plotly, sf, tidyverse, naniar, VIM, gridExtra,
       sjPlot, htmltools, officer, flextable, webshot, apaTables, MuMIn, lme4, 
       glue, grid, rsq, pheatmap, GGally, VennDiagram, glmmTMB, broom.mixed, gt)
#library(nlme)

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
venn_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/april/anova_results/plots/venn/"
# MERF LONG 
predict_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/april/final_merf_dfs"
merf_long <- read.csv('/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/april/anova_results/april_long_predictions_df_april29.csv') %>% 
             dplyr::filter(Model == "MSE Model")

basic <- read.csv(file.path(predict_dir, "basic_long_april29.csv")) %>% dplyr::rename(y_new_basic = y_hat_new)
meta <- read.csv(file.path(predict_dir, "meta_keep_long_april29.csv")) %>% dplyr::rename(y_new_meta = y_hat_new)
grs <- read.csv(file.path(predict_dir, "only_grs_long_april29.csv")) %>% dplyr::rename(y_new_grs = y_hat_new)
taxa <- read.csv(file.path(predict_dir, "only_taxa_long_april29.csv")) %>% dplyr::rename(y_new_taxa = y_hat_new)
pathway <- read.csv(file.path(predict_dir, "only_pathway_long_april29.csv")) %>% dplyr::rename(y_new_pathway = y_hat_new)
micom <- read.csv(file.path(predict_dir, "only_micom_long_april29.csv")) %>% dplyr::rename(y_new_micom = y_hat_new)
metabo <- read.csv(file.path(predict_dir, "only_metabo_long_april29.csv")) %>% dplyr::rename(y_new_metabo = y_hat_new)
all <- read.csv(file.path(predict_dir, "only_all_long_april29.csv")) %>% dplyr::rename(y_new_all = y_hat_new)

# MERF DELTA
merf_delta <- read.csv('/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/merf_python/april/anova_results/april_delta_predictions_df_april29.csv') %>% 
              dplyr::filter(Model == "MSE Model")

basic_md <- read.csv(file.path(predict_dir, "basic_delta_april29.csv")) %>% dplyr::rename(y_new_basic = y_hat_new)
meta_md <- read.csv(file.path(predict_dir, "meta_keep_delta_april29.csv")) %>% dplyr::rename(y_new_meta = y_hat_new)
grs_md <- read.csv(file.path(predict_dir, "only_grs_delta_april29.csv")) %>% dplyr::rename(y_new_grs = y_hat_new)
taxa_md <- read.csv(file.path(predict_dir, "only_taxa_delta_april29.csv")) %>% dplyr::rename(y_new_taxa = y_hat_new)
pathway_md <- read.csv(file.path(predict_dir, "only_pathway_delta_april29.csv")) %>% dplyr::rename(y_new_pathway = y_hat_new)
micom_md <- read.csv(file.path(predict_dir, "only_micom_delta_april29.csv")) %>% dplyr::rename(y_new_micom = y_hat_new)
metabo_md <- read.csv(file.path(predict_dir, "only_metabo_delta_april29.csv")) %>% dplyr::rename(y_new_metabo = y_hat_new)
all_md <- read.csv(file.path(predict_dir, "only_all_delta_april29.csv")) %>% dplyr::rename(y_new_all = y_hat_new)

# GLMMLASSO LONG
gl_lasso_long <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_long/april_long_predictions_df_april29.csv") %>% 
  dplyr::rename(y_new_meta = y_new_meta_only,
                y_new_taxa = y_new_tax_only,
                y_new_micom = y_new_micom_only,
                y_new_pathway = y_new_path_only,
                y_new_metabo = y_new_metabo_only,
                y_new_all = y_new_all_only)

gl_ftimp_long_basic <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_long/basic_gl_long_top_features_april29.csv")
gl_ftimp_long_meta <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_long/meta_gl_long_top_features_april29.csv")
gl_ftimp_long_taxa <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_long/taxa_gl_long_top_features_april29.csv")
gl_ftimp_long_micom <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_long/micom_gl_long_top_features_april29.csv")
gl_ftimp_long_pathway <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_long/pathway_gl_long_top_features_april29.csv")
gl_ftimp_long_metabo <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_long/metabo_gl_long_top_features_april29.csv")
gl_ftimp_long_all <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_long/all_gl_long_top_features_april29.csv")

# GLMMLASSO DELTA 
gl_lasso_delta <- read.csv('/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_delta/april_delta_predictions_df_april29.csv') %>% 
  dplyr::rename(y_new_basic = y_new_basic_only,
                y_new_meta = y_new_meta_only,
                y_new_taxa = y_new_tax_only,
                y_new_micom = y_new_micom_only,
                y_new_pathway = y_new_path_only,
                y_new_metabo = y_new_metab_only,
                y_new_all = y_new_all_only)

gl_ftimp_delta_basic <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_delta/basic_gl_delta_top_features_april29.csv")
gl_ftimp_delta_meta <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_delta/meta_gl_delta_top_features_april29.csv")
gl_ftimp_delta_taxa <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_delta/taxa_gl_delta_top_features_april29.csv")
gl_ftimp_delta_micom <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_delta/micom_gl_delta_top_features_april29.csv")
gl_ftimp_delta_pathway <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_delta/pathway_gl_delta_top_features_april29.csv")
gl_ftimp_delta_metabo <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_delta/metabo_gl_delta_top_features_april29.csv")
gl_ftimp_delta_all <- read.csv("/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/play_scripts/2.models/glmmlasso/march30_delta/all_gl_delta_top_features_april29.csv")

### PREDICTED vs ACTUAL BMI 
omic_dfs <- list(
  merf_long = merf_long,
  gl_lasso_long = gl_lasso_long,
  merf_delta = merf_delta, 
  gl_lasso_delta = gl_lasso_delta)

omic_name_labels <- c(
  "merf_long"      = "Merf Long",
  "gl_lasso_long"  = "GlmmLasso Long",
  "merf_delta"     = "Merf Delta",
  "gl_lasso_delta" = "GlmmLasso Delta")

for (omic_name in names(omic_dfs)) {
  omic <- omic_dfs[[omic_name]]
  
  # plot linear prediction trend
  cluster_palette <- c("#860018", "#E9A39D", "#C96374", "#D24934", "#582F2B",
    "#A0522D", "#D2691E", "#BD7C5C", "#F08080", "#A12222", "#8B0000", "#FFA07A")
  pred_data <- omic 
  gl_long <- pred_data %>% pivot_longer(cols = starts_with("y_new"),
                                        names_to = "predictor_type", 
                                        values_to = "prediction") %>% 
    mutate(Time = as.factor(Time)) %>% as.data.frame()
  
  long_plot <- ggplot(gl_long, aes(x = bmi, y = prediction, 
                                   color = predictor_type)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_smooth(method = "lm", se = FALSE, size = 1.5) +
    scale_color_manual(values = cluster_palette) +
    labs(title = glue("{omic_name_labels[[omic_name]]} Predicted BMI Values vs. Actual BMI"),
         x = "Actual Test Set BMI",
         y = "Omic Risk Score Predicted BMI",
         color = "Predictor") +
    theme_minimal()
  long_plot
  assign(paste0(omic_name, "_each_omic_plot_april29"), long_plot)
  #ggsave(file.path(omic_out_dir, paste0(omic_name, "_each_omic_plot_april29.png")), plot = long_plot, width = 8, height = 6)
  
  ### ANOVA 
  mod_dat <- omic
  lmer_basic <- lmer(bmi ~ y_new_basic + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_meta_b <- lmer(bmi ~ y_new_basic + y_new_meta + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_micom_b <- lmer(bmi ~ y_new_basic + y_new_micom  + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_path_b <- lmer(bmi ~ y_new_basic + y_new_pathway + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_tax_b <- lmer(bmi ~ y_new_basic + y_new_taxa + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_metabo_b <- lmer(bmi ~ y_new_basic + y_new_metabo + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_all_b <- lmer(bmi ~ y_new_basic + y_new_all + (1|Cluster), data = mod_dat, REML = FALSE)
  
  # Put your models into a named list
  models <- list(
    Basic   = lmer_basic,
    `Basic + Meta`  = lmer_meta_b,
    `Basic + Micom` = lmer_micom_b,
    `Basic + Pathway` = lmer_path_b,
    `Basic + Taxa`    = lmer_tax_b,
    `Basic + Metabo`  = lmer_metabo_b,
    `Basic + All` = lmer_all_b)
  
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
    c("lmer_basic", "lmer_metabo_b"),
    c("lmer_basic", "lmer_all_b"))
  
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
    mutate(across(where(is.numeric), round, 3)) %>% 
    dplyr::select(-c("term", "npar" ))

  # Extract model names for reference
  anova_table_clean$model_1 <- gsub(" vs .*", "", anova_table_clean$model_comparison)
  anova_table_clean$model_2 <- gsub(".* vs ", "", anova_table_clean$model_comparison)
  
  # Add R2c values
  anova_table_with_r2m <- anova_table_clean %>%
    left_join(r2_df %>% rename(R2m_model_1 = R2m), by = c("model_1" = "Model")) %>%
    left_join(r2_df %>% rename(R2m_model_2 = R2m), by = c("model_2" = "Model")) %>%
    dplyr::select(model_comparison, statistic, p.value, R2m_model_1, R2m_model_2)
  
  # Store the enriched table
  all_anova_tables[[omic_name]] <- anova_table_with_r2m
  print(anova_table_with_r2m)
  
  # Annotate significance and select relevant columns only
    anova_labels <- anova_table_clean %>%
    dplyr::mutate(
      Compared_Model = dplyr::recode(
        model_comparison,
        "lmer_basic vs lmer_meta_b" = "Basic + Meta",
        "lmer_basic vs lmer_tax_b" = "Basic + Taxa",
        "lmer_basic vs lmer_micom_b" = "Basic + Micom",
        "lmer_basic vs lmer_path_b" = "Basic + Pathway",
        "lmer_basic vs lmer_metabo_b" = "Basic + Metabo",
        "lmer_basic vs lmer_all_b" = "Basic + All"),
      Significance = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ "")) %>%
    dplyr::select(Compared_Model, p.value, Significance)
  
  r2_annotated <- r2_df %>%
    left_join(anova_labels, by = c("Model" = "Compared_Model"))
  
  model_order <- r2_annotated %>%
    dplyr::filter(Model != "Basic") %>%
    arrange(R2m) %>%
    pull(Model) %>%
    unique() #
  
  r2_annotated <- r2_annotated %>%
    dplyr::mutate(Model = factor(Model, levels = c("Basic", model_order)))
  
  sig_plot <- ggplot(r2_annotated, aes(x = Model, y = R2m, fill = Model)) +
    geom_col(width = 0.7) +
    geom_text(aes(y = R2m + 0.03, label = Significance), size = 12) +
    scale_y_continuous(limits = c(0, 0.55),
                       expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = c(
      "Basic + Meta"     = "#960018",
      "Basic + Taxa"     = "#8B3A3A",
      "Basic + Micom"    = "#996759",
      "Basic + Pathway"  = "#D8A39D",
      "Basic + Metabo"   = "#C76374",
      "Basic + All"      = "#C13427",
      "Basic"            = "#582F2B")) +
    labs(title = glue("Marginal R² by Model – {omic_name_labels[[omic_name]]}"),
      x = "Model Comparisons",
      y = expression(R^2*" = fixed effect variance")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none",
          plot.margin = ggplot2::margin(t = 2, r = 2, b = 2, l = 2, unit = "cm"),
      axis.text.y = element_text(face = "bold", hjust = 1, size = 14),
      axis.text.x = element_text(face = "bold", size = 14, hjust = 1, angle = 45),
      axis.title.x = element_text(margin = ggplot2::margin(0.5,0,0,0, unit = "cm"), size = 14),
      axis.title.y = element_text(margin = ggplot2::margin(0,0.8,0,0, unit = "cm"), size = 14))
  sig_plot
  assign(paste0(omic_name, "_sig_plot"), sig_plot)
  # Save the sig plot
  #ggsave(file.path(omic_out_dir, paste0(omic_name, "_sig_plot_april29.png")), plot = sig_plot,  width = 10, height = 8)
}

# Combine all plots into 1
plot_names <- glue("{names(omic_dfs)}_sig_plot")
sig_plots <- mget(plot_names)
combined_plot <- wrap_plots(sig_plots, ncol = 2) + 
  plot_annotation(title = "Comparison of Marginal R² Across BMI Modeling Scores", 
                  theme = theme(plot.title = element_text(size = 18, face = "bold"))) &
  theme(plot.margin = ggplot2::margin(20, 20, 20, 20, unit = "pt"))
print(combined_plot)


plots_linear <- glue("{names(omic_dfs)}_each_omic_plot_april29")
sig_linear_plot <- mget(plots_linear)
combined_linear_plot <- wrap_plots(sig_linear_plot, ncol = 2) + 
  plot_annotation(title = "BMI prediction accuracy across models", 
                  theme = theme(plot.title = element_text(size = 18, face = "bold"))) &
  theme(plot.margin = ggplot2::margin(20, 20, 20, 20, unit = "pt"))
print(combined_linear_plot)


######## PLOT feature importances for MERF ######################################
ft_dfs <- list(basic_lg = basic, meta_lg = meta, 
               taxa_lg = taxa, micom_lg = micom, 
               pathway_lg = pathway, metabo_lg = metabo, 
               all_lg = all)

ft_dfs_delta <- list(basic_delta = basic_md, meta_delta = meta_md, 
                     taxa_delta = taxa_md, micom_delta = micom_md, 
                     pathway_delta = pathway_md, metabo_delta = metabo_md, 
                     all_delta = all_md)

ft_colors <- c(meta_lg = "#960018", taxa_lg = "#8B3A3A", 
               micom_lg = "#996759", pathway_lg = "#D8A39D", 
               metabo_lg  = "#C76374", all_lg = "#C13427", 
               basic_lg = "#582F2B")

ft_colors_delta <- c(meta_delta = "#960018", taxa_delta = "#8B3A3A", 
                     micom_delta = "#996759", pathway_delta = "#D8A39D", 
                     metabo_delta  = "#C76374", all_delta = "#C13427", 
                     basic_delta = "#582F2B")
all_ft_dfs <- list(
  lg = list(data = ft_dfs, colors = ft_colors, label = "long"),
  delta = list(data = ft_dfs_delta, colors = ft_colors_delta, label = "delta"))

plots <- list()

for (type in names(all_ft_dfs)) {
  ft_set <- all_ft_dfs[[type]]$data
  color_map <- all_ft_dfs[[type]]$colors
  label <- all_ft_dfs[[type]]$label  # "long" or "delta"
  
  for (ft_name in names(ft_set)) {
    omic <- ft_set[[ft_name]]
    
    omic_parsed <- omic %>%
      dplyr::mutate(Top_15_Feature_Importances = gsub("'", '"', Top_15_Feature_Importances)) %>%
      dplyr::mutate(parsed = map(Top_15_Feature_Importances, fromJSON))
    
    omic_long <- omic_parsed %>%
      unnest(parsed) %>%
      dplyr::select(-Top_15_Feature_Importances) %>%
      filter(Feature != "time")
    
    avg_importance <- omic_long %>%
      group_by(Feature) %>%
      summarise(avg_importance = mean(Importance, na.rm = TRUE), .groups = "drop") %>%
      mutate(Feature = reorder(Feature, avg_importance)) %>%
      head(10)
    
    plot <- ggplot(avg_importance, aes(x = Feature, y = avg_importance)) +
      geom_col(fill = color_map[[ft_name]], width = 0.7) +
      labs(title = paste("Avg. MERF Ft. Importance -", toupper(ft_name)),
           x = "Feature", y = "Average Importance") +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "none",
        plot.margin = ggplot2::margin(t = 1, r = 1, b = 1, l = 3, unit = "cm"),
        plot.title = element_text(hjust = 0.2, face = "bold", size = 16),
        axis.text.y = element_text(face = "bold", hjust = 1, size = 14),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.title.x = element_text(margin = ggplot2::margin(1, 0, 0, 0, unit = "cm"), size = 16),
        axis.title.y = element_text(margin = ggplot2::margin(0, 1, 0, 0, unit = "cm"), size = 16),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      ) +
      scale_y_continuous(limits = c(-0.4, 0.7), expand = expansion(mult = c(0, 0.3))) +
      coord_flip(clip = "off")
    
    plot_name <- paste0(ft_name, "_", label)
    plots[[plot_name]] <- plot
    
    # Optional: Save the plot
    # ggsave(filename = file.path("omic_out_dir", paste0("ft_imp_merf_", label, "/merf_avg_importance_", ft_name, ".png")),
    #        plot = plot, width = 9, height = 5, dpi = 300)
  }
}

######## OLD: 

ft_dfs <- list(basic_lg = basic, meta_lg = meta, taxa_lg = taxa, 
          micom_lg = micom, pathway_lg = pathway, 
          metabo_lg = metabo, all_lg = all)

ft_dfs_delta <- list(basic_delta = basic_md, meta_delta = meta_md, taxa_delta = taxa_md, 
               micom_delta = micom_md, pathway_delta = pathway_md, 
               metabo_delta = metabo_md, all_delta = all_md)

ft_colors <- c(meta_lg = "#960018", taxa_lg = "#8B3A3A", micom_lg = "#996759",     
  pathway_lg = "#D8A39D", metabo_lg  = "#C76374", all_lg = "#C13427", basic_lg = "#582F2B")

ft_colors_delta <- c(meta_delta = "#960018", taxa_delta = "#8B3A3A", micom_delta = "#996759",     
               pathway_delta = "#D8A39D", metabo_delta  = "#C76374", 
               all_delta = "#C13427", basic_delta = "#582F2B")

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
      plot.margin = ggplot2::margin(t = 1, r = 1, b = 1, l = 3, unit = "cm"),  # More left margin
      plot.title = element_text(hjust = 0.2, face = "bold", size = 16),
      axis.text.y = element_text(face = "bold", hjust = 1, size = 14),
      axis.text.x = element_text(face = "bold", size = 14),
      axis.title.x = element_text(margin = ggplot2::margin(1, 0, 0, 0, unit = "cm"), size = 16),
      axis.title.y = element_text(margin = ggplot2::margin(0, 1, 0, 0, unit = "cm"), size = 16),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)) +
    scale_y_continuous(limits = c(-0.4, 0.7), 
                       expand = expansion(mult = c(0, 0.3))) +  # adds space to the right
    coord_flip(clip = "off") 

  plots[[ft_name]] <- merf_long_ft_imp # store or print
  
  #ggsave(filename = file.path(omic_out_dir, paste0("ft_imp_merf_long/merf_avg_importance_april29_", ft_name, ".png")),
  #ggsave(filename = file.path(omic_out_dir, paste0("ft_imp_merf_delta/merf_avg_importance_april29_", ft_name, ".png")),     
  #        plot = merf_long_ft_imp, width = 9, height = 5, dpi = 300)
}


### GLMMLASSO Top features plotting loops ### you have to chose long or delta below

lg_lass_ft_imp <- list(basic_lg = gl_ftimp_long_basic, meta_lg = gl_ftimp_long_meta, 
                       taxa_lg = gl_ftimp_long_taxa, micom_lg = gl_ftimp_long_micom, 
                       pathway_lg = gl_ftimp_long_pathway, 
                       metabo_lg = gl_ftimp_long_metabo, all_lg = gl_ftimp_long_all)

imap_dfr(lg_lass_ft_imp, ~ 
           summarise(.x, min = min(Estimate, na.rm = TRUE), max = max(Estimate, na.rm = TRUE)) %>%
           mutate(name = .y)) %>% bind_rows() %>% dplyr::select(name, min, max)

delta_lass_ft_imp <- list(basic_delta = gl_ftimp_delta_basic, meta_delta = gl_ftimp_delta_meta, 
                       taxa_delta = gl_ftimp_delta_taxa, micom_delta = gl_ftimp_delta_micom, 
                       pathway_delta = gl_ftimp_delta_pathway, 
                       metabo_delta = gl_ftimp_delta_metabo, all_delta = gl_ftimp_delta_all)

imap_dfr(delta_lass_ft_imp, ~ 
           summarise(.x, min = min(Estimate, na.rm = TRUE), max = max(Estimate, na.rm = TRUE)) %>%
           mutate(name = .y)) %>% bind_rows() %>% dplyr::select(name, min, max)

ft_colors <- c(meta_lg = "#960018", taxa_lg = "#8B3A3A", micom_lg = "#996759",     
               pathway_lg = "#D8A39D", metabo_lg  = "#C76374", all_lg = "#C13427",
               basic_lg = "#582F2B")

ft_colors <- c(meta_delta = "#960018", taxa_delta = "#8B3A3A", micom_delta = "#996759",     
               pathway_delta = "#D8A39D", metabo_delta  = "#C76374", all_delta = "#C13427",
               basic_delta = "#582F2B")

lass_plots <- list()
for (ft_name in names(lg_lass_ft_imp)) {
  top_features <- lg_lass_ft_imp[[ft_name]]
  ymin <- -1.0
  ymax <- 2.0
  plot_feat <- top_features %>%
    dplyr::filter(!str_detect(tolower(Feature), "time")) %>%
    slice_head(n = 10) %>%
    mutate(Estimate_capped = pmin(pmax(Estimate, ymin), ymax))  # Truncate values
  
  features <- ggplot(plot_feat, aes(x = reorder(Feature, Estimate_capped), y = Estimate_capped)) +
    geom_col(fill = ft_colors[[ft_name]], width = 0.8) +
    ggtitle(paste("Top Features & Coefficients from glmmLasso Models", ft_name)) +
    xlab("Feature") + ylab("Coefficient Estimate") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 1, r = 1, b = 1, l = 3, unit = "cm"),
      plot.title = element_text(hjust = 0.2, face = "bold", size = 16),
      axis.text.y = element_text(face = "bold", hjust = 1, size = 14),
      axis.text.x = element_text(face = "bold", size = 14),
      axis.title.x = element_text(margin = ggplot2::margin(1, 0, 0, 0, unit = "cm"), size = 16),
      axis.title.y = element_text(margin = ggplot2::margin(0, 1, 0, 0, unit = "cm"), size = 16),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    scale_y_continuous(limits = c(ymin, ymax), expand = expansion(mult = c(0, 0.3))) +
    coord_flip(clip = "off")
  
  lass_plots[[ft_name]] <- features # store or print
  
  #ggsave(filename = file.path(omic_out_dir, paste0("ft_imp_glmmlasso/glm_avg_imp_delta_april29_", ft_name, ".png")),
  #       plot = features, width = 9, height = 5, dpi = 300)
}

### HEATMAPS OF SCORES #########################################################

# Extract numeric columns from each dataframe
common_cols <- c("bmi", "y_new_basic", "y_new_meta", "y_new_taxa", 
                 "y_new_micom", "y_new_pathway", "y_new_metabo", "y_new_all")

df_list_long <- list(merf_long = merf_long[common_cols],
  gl_lasso_long = gl_lasso_long[common_cols])

df_list_long_sorted <- imap(df_list_long, ~ {
  arrange(.x, desc(bmi))})

# Rename columns to indicate source
df_list_named <- map2(df_list_long_sorted, names(df_list_long_sorted), 
                      ~ setNames(.x, paste0(.x %>% names(), "_", .y)))
df_combined <- bind_cols(df_list_named) %>% 
               rename(`Actual BMI Score` = bmi_merf_long,
                      `Basic BMI Score Merf` = y_new_basic_merf_long, 
                      `Meta BMI Score Merf` = y_new_meta_merf_long,
                      `All Omic BMI Score Merf` = y_new_all_merf_long,
                      `Taxa BMI Score Merf` = y_new_taxa_merf_long,
                      `Metabolomics BMI Score Merf` = y_new_metabo_merf_long,
                      `Pathways BMI Score Merf` = y_new_pathway_merf_long,
                      `Micom BMI Score Merf` = y_new_micom_merf_long,
                      `Basic BMI Score glmmLasso` = y_new_basic_gl_lasso_long, 
                      `Meta BMI Score glmmLasso` = y_new_meta_gl_lasso_long,
                      `All Omic BMI Score glmmLasso` = y_new_all_gl_lasso_long,
                      `Taxa BMI Score glmmLasso` = y_new_taxa_gl_lasso_long,
                      `Metabolomics BMI Score glmmLasso` = y_new_metabo_gl_lasso_long,
                      `Pathways BMI Score glmmLasso` = y_new_pathway_gl_lasso_long,
                      `Micom BMI Score glmmLasso` = y_new_micom_gl_lasso_long) %>% 
  dplyr::select(-c("bmi_gl_lasso_long"))

corr_mat <- cor(df_combined, use = "pairwise.complete.obs")
col_order <- colnames(corr_mat)
col_order <- c(grep("Merf$", colnames(corr_mat), value = TRUE),                # Merf columns first
  setdiff(colnames(corr_mat), grep("Merf$", colnames(corr_mat), value = TRUE)))
corr_mat_ordered <- corr_mat[col_order, col_order] # Reorder 
corr_long <- reshape2::melt(corr_mat_ordered) # Melt 

var_levels <- levels(factor(corr_long$Var1)) 
actual_idx <- which(var_levels == "Actual BMI Score")
ggplot(corr_long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#043362", high = "#C96374", mid = "white", 
                       midpoint = 0.25, limit = c(-0.5, 1), name = "Correlation") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.margin = ggplot2::margin(1, 1, 1, 1, unit = "cm")) +
  labs(title = "Longitudinal BMI Scores Correlation Heatmap", x = "", y = "") +
  # Draw lines above and below the "Actual BMI Score"
  geom_hline(yintercept = actual_idx - 0.5, color = "grey70", size = 1) +
  geom_hline(yintercept = actual_idx + 0.5, color = "grey70", size = 1) +
  geom_vline(xintercept = actual_idx - 0.5, color = "grey70", size = 1) +
  geom_vline(xintercept = actual_idx + 0.5, color = "grey70", size = 1)

# Do the same for the delta value 
df_list_delta <- list(merf_delta = merf_delta[common_cols],
                      gl_lasso_delta = gl_lasso_delta[common_cols])


########### PAIRS

# Example: y_new_meta across dataframes
compare_df <- tibble(merf_long = merf_long$y_new_meta,
  gl_lasso_long = gl_lasso_long$y_new_meta)
pairs(compare_df, main = "Pairwise Scatter Plots for y_new_meta")

compare_df <- tibble(merf_delta = merf_delta$y_new_meta,
                     gl_lasso_delta = gl_lasso_delta$y_new_meta)
ggpairs(compare_df, title = "Pairwise Delta Correlations for y_new_meta")

#### Individual heatmaps
plot_heatmap <- function(df, title) {
  corr <- cor(df[common_cols], use = "pairwise.complete.obs")
  pheatmap(corr, 
           main = title,
           color = colorRampPalette(c("#043362", "white", "#C96374"))(100),
           fontsize = 10)
}

plot_heatmap(merf_long, "MERF Long")
plot_heatmap(merf_delta, "MERF Delta")
plot_heatmap(gl_lasso_long, "GL Lasso Long")
plot_heatmap(gl_lasso_delta, "GL Lasso Delta")

#### correlation diff

# Step 1: Compute correlations for each dataframe
cor_matrices_long <- map(df_list_long, ~ cor(.x, use = "complete.obs"))
cor_matrices_delta <- map(df_list_delta, ~ cor(.x, use = "complete.obs"))

# Fix: Proper binding of melted correlation matrices
cor_long_df <- list(
  merf_long = reshape2::melt(cor_matrices_long$merf_long) %>% 
              mutate(method = "Merf", type = "Longitudinal BMI Scores"),
  gl_lasso_long = reshape2::melt(cor_matrices_long$gl_lasso_long) %>% 
                  mutate(method = "glmmLasso", type = "Longitudinal BMI Scores"),
  merf_delta = reshape2::melt(cor_matrices_delta$merf_delta) %>% 
               mutate(method = "Merf", type = "BMI Change Scores"),
  gl_lasso_delta = reshape2::melt(cor_matrices_delta$gl_lasso_delta) %>% 
                   mutate(method = "glmmLasso", type = "BMI Change Scores")) %>% 
  bind_rows() %>%
    mutate(Var1 = fct_recode(Var1,
                        "Actual BMI" = "bmi",
                        "Basic Score" = "y_new_basic",
                        "Meta Score" = "y_new_meta",
                        "Taxa Score" = "y_new_taxa",
                        "Micom Score" = "y_new_micom",
                        "Pathway Score" = "y_new_pathway",
                        "Metabolomics Score" = "y_new_metabo",
                        "All Omic Score" = "y_new_all"),
      Var2 = fct_recode(Var2,
                        "Actual BMI" = "bmi",
                        "Basic Score" = "y_new_basic",
                        "Meta Score" = "y_new_meta",
                        "Taxa Score" = "y_new_taxa",
                        "Micom Score" = "y_new_micom",
                        "Pathway Score" = "y_new_pathway",
                        "Metabolomics Score" = "y_new_metabo",
                        "All Omic Score" = "y_new_all"))

ggplot(cor_long_df, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  facet_grid(method ~ type) +
  scale_fill_gradient2(low = "#a3c6d4", high = "#C96374", mid = "white", midpoint = 0.5,
                       limit = c(0, 1), name = "Correlation") +
  theme_minimal(base_size = 16) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15), 
        strip.text = element_text(face = "bold", size = 16)) +
  labs(title = "Correlation Heatmaps: Long vs Delta", x = "", y = "")

##### VENN DIAGRAM FOR GLMLASSO ##### ##### ##### ##### ##### ##### ##### ##### 
omic_types <- c("basic", "meta", "taxa", "micom", "pathway", "metabo", "all")
venn_plots <- list()  # To store plots

for (omic in omic_types) {
  long_df <- lg_lass_ft_imp[[paste0(omic, "_lg")]]
  delta_df <- delta_lass_ft_imp[[paste0(omic, "_delta")]]
  
  top_long <- long_df %>%
    arrange(desc(abs_estimate)) %>%
    slice_head(n = 20) %>%
    pull(Feature)
  
  top_delta <- delta_df %>%
    arrange(desc(abs_estimate)) %>%
    slice_head(n = 20) %>%
    pull(Feature)
  
  # Generate venn plot grob
  venn_grob <- venn.diagram(
    x = list(Longitudinal = unique(top_long), Delta = unique(top_delta)),
    filename = NULL,
    fill = c("#a3c6d4", "#C96374"),
    alpha = 0.6,
    cex = 1.5,
    cat.cex = 1.4,
    cat.fontface = "bold",
    print.mode = "percent",
    main = paste("Top 20 Feature Overlap -", toupper(omic)),
    main.cex = 2)
  
  # Store grob in list
  venn_plots[[omic]] <- venn_grob
  #png(file.path(omic_out_dir, paste0("venn_glmmLasso_", omic, ".png")), width = 800, height = 600)
  grid.draw(venn_grob)
  dev.off()
}

grid.newpage()
grid.draw(venn_plots$basic)
grid.draw(venn_plots$meta)
grid.draw(venn_plots$taxa)  # or $basic, $meta, etc.
grid.draw(venn_plots$meta)
grid.draw(venn_plots$micom)
grid.draw(venn_plots$metabo)
grid.draw(venn_plots$all)
# Create Venn grobs for each pair
grid.arrange(
  venn_plots$basic,
  venn_plots$meta,
  venn_plots$taxa,
  venn_plots$micom,
  venn_plots$metabo,
  venn_plots$all,
  ncol = 2  # or adjust as needed
)

#### VENN DIAGRAMM for MERF ################################################

pair_keys <- c("basic", "micom", "taxa", "meta", "metabo", "all", "pathway")

# Create a function to extract top 20 features by average importance
get_top_features <- function(df) {
  df %>%
    mutate(Top_15_Feature_Importances = gsub("'", '"', Top_15_Feature_Importances)) %>%
    mutate(parsed = map(Top_15_Feature_Importances, fromJSON)) %>%
    unnest(parsed) %>%
    filter(Feature != "time") %>%
    group_by(Feature) %>%
    summarise(avg_importance = mean(Importance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(avg_importance)) %>%
    slice_head(n = 20) %>%
    pull(Feature)
}

# Store overlapping features
overlap_results <- list()
for (key in pair_keys) {
  df_long <- ft_dfs[[paste0(key, "_lg")]]
  df_delta <- ft_dfs_delta[[paste0(key, "_delta")]]
  if (!is.null(df_long) && !is.null(df_delta)) {
    top_lg <- get_top_features(df_long)
    top_delta <- get_top_features(df_delta)
    overlap <- intersect(top_lg, top_delta)
    overlap_results[[key]] <- list(
      long = top_lg,
      delta = top_delta,
      overlap = overlap)
    cat("Pair:", key, "\n")
    cat("Overlap (", length(overlap), " features):\n", paste(overlap, collapse = ", "), "\n\n")
  }
}

grid.newpage()
# Loop over all key pairs
for (key in pair_keys) {
  overlap <- overlap_results[[key]]
  if (is.null(overlap)) next
  file_path <- file.path(venn_dir, paste0("venn_", key, ".png")) # save path
  png(filename = file_path, width = 800, height = 800, res = 150) # Open PNG 
  
  # Draw Venn
  grid.newpage()
  draw.pairwise.venn(
    area1 = length(overlap$long),
    area2 = length(overlap$delta),
    cross.area = length(overlap$overlap),
    category = c(paste0(key, "_lg"), paste0(key, "_delta")),
    fill = c("#D8A39D", ""),
    cat.cex = 1.2,
    cex = 1.5)
  
  grid.text(paste("MERF models Top 20 Feature Overlap -", toupper(key)),  # custom title
            y = unit(0.95, "npc"),
            gp = gpar(fontsize = 16, fontface = "bold"))
  dev.off()
}


#"#a9746e" , "#d7837f", "#5d2a3b", "#f1c27d" , "#a0522d", "#f4cccc", "#872657"

# Function to extract top 20 features from MERF-style data
get_top_features <- function(df) {
  df %>%
    mutate(Top_15_Feature_Importances = gsub("'", '"', Top_15_Feature_Importances)) %>%
    mutate(parsed = map(Top_15_Feature_Importances, fromJSON)) %>%
    unnest(parsed) %>%
    filter(Feature != "time") %>%
    group_by(Feature) %>%
    summarise(avg_importance = mean(Importance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(avg_importance)) %>%
    slice_head(n = 20) %>%
    pull(Feature)
}

# Define omics and assign colors
omic_types <- c("basic", "meta", "taxa", "micom", "pathway", "metabo", "all")
venn_colors <- c("#d2495e", "#e9919c", "#5d2a3b", "#f1c28d", "#c0522d", "#f1cabb", "#872657")

# Output directory
venn_out_dir <- "venn_outputs"  # Change this if needed
dir.create(venn_out_dir, showWarnings = FALSE)

# Loop through each omic type
venn_plots_4_ <- list() 
for (i in seq_along(omic_types)) {
  omic <- omic_types[i]
  color <- venn_colors[i %% length(venn_colors) + 1]  # rotate if more omics than colors
  
  # Extract top 20 from GLMM
  top_glmm_lg <- lg_lass_ft_imp[[paste0(omic, "_lg")]] %>%
    arrange(desc(abs_estimate)) %>%
    slice_head(n = 20) %>%
    pull(Feature)
  
  top_glmm_delta <- delta_lass_ft_imp[[paste0(omic, "_delta")]] %>%
    arrange(desc(abs_estimate)) %>%
    slice_head(n = 20) %>%
    pull(Feature)
  
  # Extract top 20 from MERF
  top_merf_lg <- get_top_features(ft_dfs[[paste0(omic, "_lg")]])
  top_merf_delta <- get_top_features(ft_dfs_delta[[paste0(omic, "_delta")]])
  
  # Prepare input
  venn_input <- list(
    GLMM_LG = unique(top_glmm_lg),
    GLMM_DELTA = unique(top_glmm_delta),
    MERF_LG = unique(top_merf_lg),
    MERF_DELTA = unique(top_merf_delta))
  
  # Generate Venn diagram
  venn_plot <- venn.diagram(
    x = venn_input,
    filename = NULL,
    fill = rep(color, 4),
    alpha = 0.5,
    cex = 1.3,
    cat.cex = 1.2,
    cat.fontface = "bold",
    main = paste("Top 20 Feature Overlap -", toupper(omic)),
    main.cex = 2,
    margin = 0.1)
  
  venn_plots_4_[[omic]] <- venn_plot
  
  # Save to file
  file_path <- file.path(venn_out_dir, paste0("venn_4set_", omic, ".png"))
  png(file_path, width = 1000, height = 800, res = 150)
  grid.newpage()
  grid.draw(venn_plot)
  dev.off()
}

grid.arrange(
  venn_plots_4_$basic,
  venn_plots_4_$meta,
  venn_plots_4_$taxa,
  venn_plots_4_$micom,
  venn_plots_4_$metabo,
  venn_plots_4_$all,
  ncol = 2  # or adjust as needed
)

grid.arrange(
  grobTree(venn_plots_4_$basic),
  grobTree(venn_plots_4_$meta),
  grobTree(venn_plots_4_$taxa),
  grobTree(venn_plots_4_$micom),
  grobTree(venn_plots_4_$metabo),
  grobTree(venn_plots_4_$all),
  ncol = 2, 
  padding = unit(2, "cm"))












