# Feature important plots 

# Read in data / objects / packages (you might not need all of them)
load("data_for_plot.RData")

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

### MERF FEATURE IMPORTANCES ###

merf_dfs_long <- list(Basic = basic, Meta = meta, 
                      Taxa = taxa, Micom = micom, 
                      Pathway = pathway, Metabo = metabo, 
                      All = all)

merf_dfs_delta <- list(Basic = basic_md, Meta = meta_md, 
                       Taxa = taxa_md, Micom = micom_md, 
                       Pathway = pathway_md, Metabo = metabo_md, 
                       All = all_md)

ft_colors <- c(Meta = "#960018", Taxa = "#8B3A3A", 
               Micom = "#996759", Pathway = "#D8A39D", 
               Metabo  = "#C76374", All = "#C13427", 
               Basic = "#582F2B")

ft_colors_delta <- c(Meta = "#960018", Taxa = "#8B3A3A", 
                     Micom = "#996759", Pathway = "#D8A39D", 
                     Metabo  = "#C76374", All = "#C13427", 
                     Basic = "#582F2B")

# Combined feature sets
all_ft_dfs <- list(merf_lg = list(data = merf_dfs_long, 
                                  colors = ft_colors, 
                                  label = "long", 
                                  model_label = "MERF (Longitudinal BMI Prediction)",
                                  x_ax_label = "Feature Importance"),
                   merf_delta = list(data = merf_dfs_delta, 
                                     colors = ft_colors_delta, 
                                     label = "delta", 
                                     model_label = "MERF (Change in BMI (Delta) Prediction)",
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
      geom_col(fill = color_map[[ft_name]], width = 0.5) +
      labs(title = paste("", toupper(ft_name)),
           x = "Top 10 Model Features", 
           y = paste(" ", x_ax_label)) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none",
            plot.margin = ggplot2::margin(t = 1, r = 1, b = 1, l = 3, unit = "cm"),
            plot.title = element_text(hjust = 0, face = "bold", size = 16),
            axis.text.y = element_text(hjust = 1, face = "bold", size = 12),
            axis.text.x = element_text(face = "bold", size = 12),
            axis.title.x = element_text(margin = ggplot2::margin(0.5, 0, 0, 0, unit = "cm"), size = 14),
            axis.title.y = element_text(margin = ggplot2::margin(0, 0.5, 0, 0, unit = "cm"), size = 14),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA)) +
      scale_y_continuous(limits = c(-0.1, 0.43), expand = expansion(mult = c(0, 0.05))) +
      coord_flip(clip = "off")
    
    # Save to named R object
    #object_name <- paste0(ft_set, "_", ft_name, "_", label, "_plot")
    object_name <- paste0("plot_", type, "_", ft_name, "_", label)
    assign(object_name, plot, envir = .GlobalEnv)
    merf_plots[[object_name]] <- plot
    png_name <- file.path(output_dir, paste0("feature_importance_", ft_name, "_", label, ".png"))
    ggsave(filename = png_name, plot = plot, width = 9, height = 5, dpi = 300)
  }
}

merf_long_plots <- merf_plots[c(2, 5, 6, 7)]
combined_merf_long_plot <- wrap_plots(merf_long_plots, ncol = 2, nrow = 2) +
  plot_annotation(title = "MERF longitudinal BMI Prediction Features",
                  theme = theme(plot.title = element_text(size = 18, face = "bold"))) &
  theme(plot.margin = ggplot2::margin(8, 8, 8, 8, unit = "pt"))

merf_delta_plots <- merf_plots[c(2, 3, 5, 6)]
combined_merf_delta_plot <- wrap_plots(merf_delta_plots, nrow = 2, ncol = 2) +
  plot_annotation(title = "MERF BMI Change Prediction Features",
                  theme = theme(plot.title = element_text(size = 18, face = "bold"))) &
  theme(plot.margin = ggplot2::margin(8, 20, 8, 20, unit = "pt"))

print(combined_merf_delta_plot)
print(combined_merf_long_plot)

### GLMMLASSO FEATURE IMPORTANCES ###

lass_dfs_long <- list(Basic = gl_ftimp_long_basic, Meta = gl_ftimp_long_meta, 
                      Taxa = gl_ftimp_long_taxa, Micom = gl_ftimp_long_micom, 
                      Pathway = gl_ftimp_long_pathway, Metabo = gl_ftimp_long_metabo, 
                      All = gl_ftimp_long_all)

lass_dfs_delta <- list(Basic = gl_ftimp_delta_basic, Meta = gl_ftimp_delta_meta, 
                       Taxa = gl_ftimp_delta_taxa, Micom = gl_ftimp_delta_micom, 
                       Pathway = gl_ftimp_delta_pathway, Metabo = gl_ftimp_delta_metabo, 
                       All = gl_ftimp_delta_all)

# Combined feature sets
all_lass_ft_dfs <- list(lass_lg = list(data = lass_dfs_long, colors = ft_colors, 
                                       label = "long", 
                                       model_label = "GlmmLasso (Longitudinal BMI Prediction)",
                                       x_ax_label = "Feature Co-efficient"),
                        lass_delta = list(data = lass_dfs_delta, colors = ft_colors_delta, 
                                          label = "delta", 
                                          model_label = "GlmmLasso (Change in BMI (Delta) Prediction)",
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
    ymax <- 1.1
    plot_feat <- top_features %>%
      dplyr::filter(!str_detect(tolower(Feature), "time")) %>%
      slice_head(n = 10) %>%
      mutate(Estimate_capped = pmin(pmax(Estimate, ymin), ymax))  # Truncate values
    
    features <- ggplot(plot_feat, aes(x = reorder(Feature, Estimate_capped), 
                                      y = Estimate_capped)) +
      geom_col(fill = color_map[[ft_name]], width = 0.7) +
      ggtitle(paste("Top Features & Coefficients from glmmLasso Models", ft_name)) +
      labs(title = paste(model_label, toupper(ft_name)),
           x = "Top 10 Model Features", 
           y = paste("", x_ax_label)) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none",
            plot.margin = ggplot2::margin(t = 1, r = 1, b = 1, 
                                          l = 3, unit = "cm"),
            plot.title = element_text(hjust = 0, face = "bold", size = 14, 
                                      margin = ggplot2::margin(t = 0, r = 0, 
                                                               b = 20, l = -150, 
                                                               unit = "pt")),
            axis.text.y = element_text(hjust = 1, size = 13),
            axis.text.x = element_text(size = 13),
            axis.title.x = element_text(margin = ggplot2::margin(0.5, 0, 0, 0, 
                                                                 unit = "cm"), 
                                        size = 14),
            axis.title.y = element_text(margin = ggplot2::margin(0, 0.5, 0, 0, 
                                                                 unit = "cm"), 
                                        size = 14),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA)) +
      scale_y_continuous(limits = c(ymin, ymax), 
                         expand = expansion(mult = c(0, 0.1))) +
      coord_flip(clip = "off")
    
    object_name <- paste0("plot_", type, "_", ft_name, "_", label)
    assign(object_name, features, envir = .GlobalEnv)
    lass_plots[[ft_name]] <- features # store or print
    png_name <- file.path(output_dir, paste0("lass_ft_imp/feature_importance_", 
                                             ft_name, "_", label, ".png"))
    #ggsave(filename = png_name, plot = plot, width = 9, height = 5, dpi = 300)
  }
}