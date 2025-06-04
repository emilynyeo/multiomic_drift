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
    anova_out <- anova(model_1, model_2)  # This works with lmerMod objects
    #tidied_result <- tidy(anova_result)  # Tidy the ANOVA result using broom::tidy()
    anova_out$model_comparison <- paste(model_pair[1], "vs", model_pair[2])
    anova_results[[length(anova_results) + 1]] <- anova_out
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
    dplyr::select(Compared_Model, df, AIC, statistic, p.value, Significance)
  
  r2_annotated <- r2_df %>%
    left_join(anova_labels, by = c("Model" = "Compared_Model"))
  
  for_table <- r2_annotated %>% dplyr::select(-c(R2c)) %>% 
               rename("R²*" = R2m, 
                      "P-Value" = p.value,
                      "Chisq Statistic" = statistic, 
                      "Df" = df)
  
  for_table %>% gt() %>%
    tab_header(title = "Model Comparison (ANOVA)") %>%
    fmt_number(columns = where(is.numeric), decimals = 2) %>%
    fmt_number(columns = "Df", decimals = 0) %>% 
    sub_missing(missing_text = "-") %>%
    cols_width(1 ~ pct(20),`R²*` ~ pct(15), Df ~ pct(8),
      `Chisq Statistic` ~ pct(15), Significance ~ pct(15),
      `P-Value` ~ pct(15), AIC ~ pct(11), Significance ~ pct(11)) %>%
    cols_align(align = "center",
      columns = everything()) %>% 
    cols_align("left", columns = Model) %>%
    tab_options(table.width = pct(60),
                table.align = "center") %>%
    opt_table_font(font = list(
        google_font("Arial"),   # fallback to system font if unavailable
        default_fonts())) %>%
    tab_style(style = cell_text(weight = "bold"),
      locations = cells_column_labels())
    
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
merf_dfs_long <- list(basic_lg = basic, meta_lg = meta, 
               taxa_lg = taxa, micom_lg = micom, 
               pathway_lg = pathway, metabo_lg = metabo, 
               all_lg = all)

merf_dfs_delta <- list(basic_delta = basic_md, meta_delta = meta_md, 
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

# Combined feature sets
all_ft_dfs <- list(merf_lg = list(data = merf_dfs_long, colors = ft_colors, 
                             label = "long", model_label = "MERF",
                             x_ax_label = "Feature Importance"),
                   merf_delta = list(data = merf_dfs_delta, colors = ft_colors_delta, 
                                label = "delta", model_label = "MERF",
                                x_ax_label = "Feature Importance"))

merf_plots <- list()
for (type in names(all_ft_dfs)) {
  ft_set <- all_ft_dfs[[type]]$data
  color_map <- all_ft_dfs[[type]]$colors
  label <- all_ft_dfs[[type]]$label
  model_label <- all_ft_dfs[[type]]$model_label  # "MERF" or "GlmmLasso"
  x_ax_label <- all_ft_dfs[[type]]$x_ax_label
  
  output_dir <- file.path("omic_out_dir", paste0("ft_imp_", model_label, "_", label))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (ft_name in names(ft_set)) {
    omic <- ft_set[[ft_name]]
    
    omic_parsed <- omic %>%
      mutate(Top_15_Feature_Importances = gsub("'", '"', Top_15_Feature_Importances)) %>%
      mutate(parsed = map(Top_15_Feature_Importances, fromJSON))
    
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
      labs(title = paste(model_label, toupper(ft_name)),
        x = "Top 10 Model Features", 
        y = paste("Average ", x_ax_label)) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none",
        plot.margin = ggplot2::margin(t = 1, r = 1, b = 1, l = 3, unit = "cm"),
        plot.title = element_text(hjust = 0.2, face = "bold", size = 16),
        axis.text.y = element_text(face = "bold", hjust = 1, size = 13),
        axis.text.x = element_text(face = "bold", size = 13),
        axis.title.x = element_text(margin = ggplot2::margin(0.5, 0, 0, 0, unit = "cm"), size = 14),
        axis.title.y = element_text(margin = ggplot2::margin(0, 0.5, 0, 0, unit = "cm"), size = 14),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)) +
      scale_y_continuous(limits = c(-0.1, 0.6), expand = expansion(mult = c(0, 0.05))) +
      #scale_x_continuous(limits = c(-0.1, 0.6)) +
      coord_flip(clip = "off")
    
    # Save to named R object
    #object_name <- paste0(ft_set, "_", ft_name, "_", label, "_plot")
    object_name <- paste0("plot_", type, "_", ft_name, "_", label)
    assign(object_name, plot, envir = .GlobalEnv)
    merf_plots[[object_name]] <- plot
    png_name <- file.path(output_dir, paste0("feature_importance_", ft_name, "_", label, ".png"))
    #ggsave(filename = png_name, plot = plot, width = 9, height = 5, dpi = 300)
  }
}

merf_long_plots <- merf_plots[c(2, 3, 4, 5, 6, 7)]
combined_merf_long_plot <- wrap_plots(merf_long_plots, ncol = 2, nrow = 3) +
  plot_annotation(title = "MERF longitudinal BMI Prediction Features",
    theme = theme(plot.title = element_text(size = 18, face = "bold"))) &
  theme(plot.margin = ggplot2::margin(8, 8, 8, 8, unit = "pt"))

merf_delta_plots <- merf_plots[c(9, 13, 14)]
combined_merf_delta_plot <- wrap_plots(merf_delta_plots, nrow = 3) +
  plot_annotation(title = "MERF BMI Change Prediction Features",
                  theme = theme(plot.title = element_text(size = 18, face = "bold"))) &
  theme(plot.margin = ggplot2::margin(8, 40, 8, 40, unit = "pt"))

print(combined_merf_long_plot)
print(combined_merf_delta_plot)

### GLMMLASSO Top features plotting 

lass_dfs_long <- list(basic_lg = gl_ftimp_long_basic, meta_lg = gl_ftimp_long_meta, 
                      taxa_lg = gl_ftimp_long_taxa, micom_lg = gl_ftimp_long_micom, 
                      pathway_lg = gl_ftimp_long_pathway, metabo_lg = gl_ftimp_long_metabo, 
                      all_lg = gl_ftimp_long_all)

lass_dfs_delta <- list(basic_delta = gl_ftimp_delta_basic, meta_delta = gl_ftimp_delta_meta, 
                       taxa_delta = gl_ftimp_delta_taxa, micom_delta = gl_ftimp_delta_micom, 
                       pathway_delta = gl_ftimp_delta_pathway, metabo_delta = gl_ftimp_delta_metabo, 
                       all_delta = gl_ftimp_delta_all)

# Combined feature sets
all_lass_ft_dfs <- list(lass_lg = list(data = lass_dfs_long, colors = ft_colors, 
                                  label = "long", model_label = "GlmmLasso",
                                  x_ax_label = "Feature Co-efficient"),
                   lass_delta = list(data = lass_dfs_delta, colors = ft_colors_delta, 
                                     label = "delta", model_label = "GlmmLasso",
                                     x_ax_label = "Feature Co-efficient"))
imap_dfr(lass_dfs_long, ~ 
           summarise(.x, min = min(Estimate, na.rm = TRUE), max = max(Estimate, na.rm = TRUE)) %>%
           mutate(name = .y)) %>% bind_rows() %>% dplyr::select(name, min, max)

imap_dfr(lass_dfs_delta, ~ 
           summarise(.x, min = min(Estimate, na.rm = TRUE), max = max(Estimate, na.rm = TRUE)) %>%
           mutate(name = .y)) %>% bind_rows() %>% dplyr::select(name, min, max)

lass_plots <- list()
for (type in names(all_lass_ft_dfs)) {
  ft_set <- all_lass_ft_dfs[[type]]$data
  color_map <- all_lass_ft_dfs[[type]]$colors
  label <- all_lass_ft_dfs[[type]]$label
  model_label <- all_lass_ft_dfs[[type]]$model_label  # "MERF" or "GlmmLasso"
  x_ax_label <- all_lass_ft_dfs[[type]]$x_ax_label
  
  for (ft_name in names(ft_set)) {
    top_features <- ft_set[[ft_name]]
    ymin <- -1.0
    ymax <- 1.5
    plot_feat <- top_features %>%
      dplyr::filter(!str_detect(tolower(Feature), "time")) %>%
      slice_head(n = 10) %>%
      mutate(Estimate_capped = pmin(pmax(Estimate, ymin), ymax))  # Truncate values
    
    features <- ggplot(plot_feat, aes(x = reorder(Feature, Estimate_capped), y = Estimate_capped)) +
      geom_col(fill = color_map[[ft_name]], width = 0.7) +
      ggtitle(paste("Top Features & Coefficients from glmmLasso Models", ft_name)) +
      labs(title = paste(model_label, toupper(ft_name)),
           x = "Top 10 Model Features", 
           y = paste("Average ", x_ax_label)) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none",
            plot.margin = ggplot2::margin(t = 1, r = 1, b = 1, l = 3, unit = "cm"),
            plot.title = element_text(hjust = 0.2, face = "bold", size = 16),
            axis.text.y = element_text(face = "bold", hjust = 1, size = 13),
            axis.text.x = element_text(face = "bold", size = 13),
            axis.title.x = element_text(margin = ggplot2::margin(0.5, 0, 0, 0, unit = "cm"), size = 14),
            axis.title.y = element_text(margin = ggplot2::margin(0, 0.5, 0, 0, unit = "cm"), size = 14),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA)) +
      scale_y_continuous(limits = c(ymin, ymax), 
                         expand = expansion(mult = c(0, 0.1))) +
      coord_flip(clip = "off")
    
    object_name <- paste0("plot_", type, "_", ft_name, "_", label)
    assign(object_name, features, envir = .GlobalEnv)
    lass_plots[[ft_name]] <- features # store or print
    png_name <- file.path(output_dir, paste0("feature_importance_", ft_name, "_", label, ".png"))
    #ggsave(filename = png_name, plot = plot, width = 9, height = 5, dpi = 300)
  }
}

lass_long_plots <- lass_plots[c(2, 5, 6, 7)]
combined_lass_long_plot <- wrap_plots(lass_long_plots, ncol = 2) +
  plot_annotation(title = "GlmmLasso longitudinal BMI Prediction Features",
                  theme = theme(plot.title = element_text(size = 18, face = "bold"))) &
  theme(plot.margin = ggplot2::margin(20, 20, 20, 20, unit = "pt"))

lass_delta_plots <- lass_plots[c(9, 13, 14)]
combined_lass_delta_plot <- wrap_plots(lass_delta_plots, nrow = 3) +
  plot_annotation(title = "GlmmLasso BMI Change Prediction Features",
                  theme = theme(plot.title = element_text(size = 18, face = "bold"))) &
  theme(plot.margin = ggplot2::margin(20, 40, 20, 40, unit = "pt"))

print(combined_lass_long_plot)
print(combined_lass_delta_plot)

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
  scale_fill_gradient2(low = "white", high = "#E96374", 
                       limit = c(-0.5, 1), name = "Correlation") +
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

final_heat <- ggplot(cor_long_df, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  facet_grid(method ~ type) +
  scale_fill_gradient2(low = "#09236c", high = "#C96374", mid = "white", midpoint = 0.0,
                       limit = c(-0.2, 1), name = "Correlation") +
  theme_minimal(base_size = 16) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15), 
        strip.text = element_text(face = "bold", size = 16)) +
  labs(title = "Correlation Heatmaps: Long vs Delta", x = "", y = "")

final_heat

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
  top_glmm_lg <- get(paste0("gl_ftimp_long_", omic)) %>%
    arrange(desc(abs_estimate)) %>%
    slice_head(n = 20) %>%
    pull(Feature)
  
  top_glmm_delta <- get(paste0("gl_ftimp_delta_", omic)) %>%
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


# R2 anova plots
print(combined_plot)
print(combined_linear_plot)
# Co-efficients
print(combined_merf_long_plot)
print(combined_merf_delta_plot)
print(combined_lass_long_plot)
print(combined_lass_delta_plot)
# Heat map
final_heat
# Venn Diagram Overlap 
grid.arrange(
  grobTree(venn_plots_4_$basic),
  grobTree(venn_plots_4_$meta),
  grobTree(venn_plots_4_$taxa),
  grobTree(venn_plots_4_$micom),
  grobTree(venn_plots_4_$metabo),
  grobTree(venn_plots_4_$all),
  ncol = 2, 
  padding = unit(2, "cm"))












