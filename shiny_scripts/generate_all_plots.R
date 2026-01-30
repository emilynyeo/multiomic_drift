# ============================================================================
# Script to generate all plots and tables from rendered_word_plots.Rmd
# ============================================================================
# This script extracts all R code from the R Markdown file and runs it
# to produce the same output files (plots, tables, etc.)
#
# Usage: source("shiny_scripts/generate_all_plots.R")
#        or run: Rscript shiny_scripts/generate_all_plots.R
#
# Output files will be saved to: paper_plots/
# ============================================================================

# Setup
rm(list = ls())

base_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions"
source(file.path(base_dir, "shiny_scripts/plotting_utils.R"))
source(file.path(base_dir, "shiny_scripts/plot_input_data.R"))

# Load required packages
library(ggpattern)  # Required for scale_pattern_manual and geom_col_pattern
library(patchwork)  # Required for combining plots

# Set output directories
omic_out_dir <- file.path(base_dir, "paper_plots_jan18")
venn_dir <- file.path(base_dir, "paper_plots_jan18")

# Create output subdirectories if they don't exist
output_dirs <- c("shap_02", "anova_tables_02", "aov_feats_combined_02", "venns_02", "lass_ft_imp_02", 
                 "20251110", "pred_vs_actual", "aov_feats_combined", "venns", "venn_outputs")
for (dir in output_dirs) {
  dir.create(file.path(omic_out_dir, dir), recursive = TRUE, showWarnings = FALSE)
}

# Color palettes
cluster_palette_contrast <- c(
  "#004E64", "#B0E0E6", "#4682B4", "#2F4F4F", "#556B2F", "#6B8E23",
  "#708090", "#A9CBB7", "#F0E68C", "#E6BE8A", "#FFF5E1", "#B5A642")

# Main feature type colors
ft_colors <- c(
  `Clinical+`     = "#fbbdd1",
  `Clinical`     = "#fbbdd1",
  Taxa         = "#8B3A3A",
  MICOM        = "#996759",
  Pathway      = "#D8A39D",
  Metabolomics = "#C76374",
  Combined     = "#b44241",
  `Combined*`     = "#b44241",
  GRS          = "#E18476",
  Basic        = "#582F2B")

title_map <- c(
  shap_long_all_no_age_sex   = "MERF Combined",
  shap_long_metabo   = "MERF Metabolomics",
  shap_long_grs      = "MERF GRS",
  shap_long_pathway  = "MERF Pathway",
  shap_long_meta_no_age_sex = "MERF Clinical")

shap_long_plot_list <- list()
heights_vec <- numeric(0)

# Filter shap_dfs to only entries that have a title in title_map
valid_names <- intersect(names(shap_dfs), names(title_map))
shap_dfs_filtered <- shap_dfs[valid_names]

for (name in names(shap_dfs_filtered)) {
  df <- shap_dfs[[name]]

  # Determine group for color mapping
  group_key <- case_when(
  grepl("meta_no_age_sex$", name) ~ "Clinical",
  grepl("meta$", name)            ~ "Clinical+",
  grepl("all_no_clin$", name)     ~ "Combined*",
  grepl("all_no_age_sex$", name)  ~ "Combined",
  grepl("taxa$", name)            ~ "Taxa",
  grepl("micom$", name)           ~ "MICOM",
  grepl("pathway$", name)         ~ "Pathway",
  grepl("metabo$", name)          ~ "Metabolomics",
  grepl("grs$", name)             ~ "GRS",
  grepl("basic$", name)           ~ "Basic",
  TRUE                            ~ "Combined")
  
  high_col <- ft_colors[[group_key]]

  # Count features & limit to top 10
  top_features <- df %>%
    filter(Feature != "time") %>%
    group_by(Feature) %>%
    summarise(max_abs_shap = max(abs(SHAP), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_abs_shap)) %>%
    slice_head(n = 10) %>%
    pull(Feature)
  
  df_top10 <- df %>% filter(Feature %in% top_features)

  # Calculate height
  n_features <- length(unique(df_top10$Feature))
  heights_vec <- c(heights_vec, n_features + 2)

  # Create plot
  plot_title <- title_map[[name]]
  
  p <- plot_shap_beeswarm(df_top10, plot_title, high_color = "#b44241") +
        combined_plot_theme(base_size = 18) +
        theme(plot.margin = ggplot2::margin(10, 10, 10, 10))
  
  shap_long_plot_list[[name]] <- p
}

# Normalize heights to sum to 1
heights_norm <- heights_vec / sum(heights_vec)
combined_shap_long <- wrap_plots(shap_long_plot_list, ncol = 1, 
                                 heights = heights_norm, guides = "collect") +
  plot_annotation(
    tag_levels = "A",        # Adds panel tags A., B., C., ...
    tag_prefix = "",         # Optional custom prefix
    tag_suffix = ".",        # Suffix dot after the letter
    theme = theme(legend.position = "bottom",
                  plot.tag = element_text(face = "bold", size = 35)))

# Save large PDF
ggsave(
  filename = file.path(omic_out_dir, "shap_02/all_shap_plots_combined_w_outlier.pdf"),
  plot = combined_shap_long,
  width = 18,
  height = sum(heights_vec) * 0.5,
  units = "in",
  limitsize = FALSE)

title_map_delta <- c(
  shap_delta_all_no_age_sex      = "MERF Combined",
  shap_delta_all_no_clin      = "MERF Combined*",
  shap_delta_metabo   = "MERF Metabolomics",
  shap_delta_pathway  = "MERF Pathway",
  shap_delta_taxa     = "MERF Taxonomy",
  shap_delta_meta_no_age_sex     = "MERF Clinical")

shap_delta_plot_list <- list()
heights_vec <- numeric(0)

# Filter shap_delta_dfs to only entries that have a title in title_map
valid_names_delta <- intersect(names(shap_delta_dfs), names(title_map_delta))
shap_dfs_filtered_delta <- shap_delta_dfs[valid_names_delta]

for (name in names(shap_dfs_filtered_delta)) {
  df <- shap_dfs_filtered_delta[[name]]

  # Determine group for color mapping
  group_key <- case_when(
    grepl("meta$", name)      ~ "Clinical+",
    grepl("taxa$", name)      ~ "Taxa",
    grepl("micom$", name)     ~ "MICOM",
    grepl("pathway$", name)   ~ "Pathway",
    grepl("metabo$", name)    ~ "Metabolomics",
    grepl("all$", name)       ~ "Combined",
    grepl("grs$", name)       ~ "GRS",
    grepl("basic$", name)     ~ "Basic",
    TRUE ~ "Combined")
  high_col <- ft_colors[[group_key]]

  # Count features & limit to top 10
  top_features <- df %>%
    filter(Feature != "time") %>%
    group_by(Feature) %>%
    summarise(max_abs_shap = max(abs(SHAP), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_abs_shap)) %>%
    slice_head(n = 10) %>%
    pull(Feature)
  df_top10 <- df %>% filter(Feature %in% top_features)

  # Calculate height
  n_features <- length(unique(df_top10$Feature))
  heights_vec <- c(heights_vec, n_features + 2)

  # Create plot
  plot_title <- title_map_delta[[name]]
  
  p <- plot_shap_beeswarm(df_top10, plot_title, high_color = "#b44241") +
        combined_plot_theme(base_size = 18) +
        theme(plot.margin = ggplot2::margin(10, 10, 10, 10))
  
  shap_delta_plot_list[[name]] <- p
}

# Normalize heights to sum to 1
heights_norm <- heights_vec / sum(heights_vec)
combined_shap_delta <- wrap_plots(shap_delta_plot_list,  # Combine vertically
                                 ncol = 1, heights = heights_norm, 
                                 guides = "collect") +
  plot_annotation(
    tag_levels = "A",        # Adds panel tags A., B., C., ...
    tag_prefix = "",         # Optional custom prefix
    tag_suffix = ".",        # Suffix dot after the letter
    theme = theme(legend.position = "bottom",
                  plot.tag = element_text(face = "bold", size = 35)))

# Save large PDF
ggsave(
  filename = file.path(omic_out_dir, "shap_02/all_delta_shap_plots_combined_outliers.pdf"),
  plot = combined_shap_delta,
  width = 18,
  height = sum(heights_vec) * 0.5,
  units = "in",
  limitsize = FALSE)

# Setup and modeling for each omic
omic_dfs <- list(
  merf_long = merf_long,
  gl_lasso_long = gl_lasso_long,
  merf_delta = merf_delta, 
  gl_lasso_delta = gl_lasso_delta)

omic_name_labels <- c(
  "merf_long"      = "MERF BMI",
  "gl_lasso_long"  = "GlmmLasso BMI",
  "merf_delta"     = "MERF ΔBMI",
  "gl_lasso_delta" = "GlmmLasso ΔBMI")

# Create containers for later use
all_r2_dfs <- list()
all_anova_labels <- list()

for (omic_name in names(omic_dfs)) {
  print(glue::glue("Running for: {omic_name}"))
  omic <- omic_dfs[[omic_name]]
  
  ## Fit LMER Models on scores 
  lmer_basic <- lmer(bmi ~ y_new_basic + (1|Cluster), data = omic, REML = FALSE)
  lmer_meta_b <- lmer(bmi ~ y_new_basic + y_new_meta_noas + (1|Cluster), data = omic, REML = FALSE)
  lmer_grs_b <- lmer(bmi ~ y_new_basic + y_new_grs + (1|Cluster), data = omic, REML = FALSE)
  lmer_micom_b <- lmer(bmi ~ y_new_basic + y_new_micom + (1|Cluster), data = omic, REML = FALSE)
  lmer_path_b <- lmer(bmi ~ y_new_basic + y_new_pathway + (1|Cluster), data = omic, REML = FALSE)
  lmer_tax_b <- lmer(bmi ~ y_new_basic + y_new_taxa + (1|Cluster), data = omic, REML = FALSE)
  lmer_metabo_b <- lmer(bmi ~ y_new_basic + y_new_metabo + (1|Cluster), data = omic, REML = FALSE)
  lmer_all_b <- lmer(bmi ~ y_new_basic + y_new_all_noas + (1|Cluster), data = omic, REML = FALSE)

  models <- list(
    Basic = lmer_basic,
    `Basic + Clinical` = lmer_meta_b,
    `Basic + GRS` = lmer_grs_b,
    `Basic + MICOM` = lmer_micom_b,
    `Basic + Pathway` = lmer_path_b,
    `Basic + Taxa` = lmer_tax_b,
    `Basic + Metabolomics` = lmer_metabo_b,
    `Basic + Combined` = lmer_all_b)

  # Extract RÂ² values
  r2_df <- do.call(rbind, lapply(names(models), function(name) {
    r2_vals <- r.squaredGLMM(models[[name]])
    data.frame(Model = name, R2m = r2_vals[1, "R2m"], R2c = r2_vals[1, "R2c"])
  }))

  assign(paste0(omic_name, "_r2_df"), r2_df)
  all_r2_dfs[[omic_name]] <- r2_df

  ##  ANOVA 
  comparisons <- list(
    c("lmer_basic", "lmer_meta_b"),
    c("lmer_basic", "lmer_grs_b"),
    c("lmer_basic", "lmer_tax_b"),
    c("lmer_basic", "lmer_micom_b"),
    c("lmer_basic", "lmer_path_b"),
    c("lmer_basic", "lmer_metabo_b"),
    c("lmer_basic", "lmer_all_b"))

  # Store all models in global environment (needed for `get()`)
  list2env(models, envir = .GlobalEnv)

  anova_results <- list()
  for (pair in comparisons) {
    model_1 <- get(pair[1])
    model_2 <- get(pair[2])
    anova_out <- anova(model_1, model_2)
    anova_out <- as.data.frame(anova_out)
    anova_out$model_comparison <- paste(pair[1], "vs", pair[2])
    anova_results[[length(anova_results) + 1]] <- anova_out
  }

  anova_table <- do.call(rbind, anova_results) %>%
    filter(!is.na(`Chisq`) & !is.na(`Pr(>Chisq)`)) %>%
    mutate(across(where(is.numeric), \(x) round(x, 3))) %>% 
    dplyr::select(-npar)
  rownames(anova_table) <- NULL

  # Annotate comparisons
  anova_labels <- anova_table %>%
    dplyr::mutate(Compared_Model = dplyr::recode(model_comparison,
        "lmer_basic vs lmer_meta_b" = "Basic + Clinical",
        "lmer_basic vs lmer_grs_b" = "Basic + GRS",
        "lmer_basic vs lmer_tax_b" = "Basic + Taxa",
        "lmer_basic vs lmer_micom_b" = "Basic + MICOM",
        "lmer_basic vs lmer_path_b" = "Basic + Pathway",
        "lmer_basic vs lmer_metabo_b" = "Basic + Metabolomics",
        "lmer_basic vs lmer_all_b" = "Basic + Combined"),
      Significance = case_when(
        `Pr(>Chisq)` < 0.001 ~ "***",
        `Pr(>Chisq)` < 0.01  ~ "**",
        `Pr(>Chisq)` < 0.05  ~ "*",
        TRUE ~ ""),
      `Pr(>Chisq)` = case_when(
        `Pr(>Chisq)` < 0.001 ~ "<0.001",
        `Pr(>Chisq)` < 0.01 ~ "<0.01",
        TRUE ~ as.character(`Pr(>Chisq)`)))

  all_anova_labels[[omic_name]] <- anova_labels
}

for (omic_name in names(omic_dfs)) {
  omic <- omic_dfs[[omic_name]]
  pred_data <- omic

  cluster_palette <- cluster_palette_contrast[1:12]

  gl_long <- pred_data %>%
    pivot_longer(cols = starts_with("y_new"), 
                 names_to = "predictor_type", 
                 values_to = "prediction") %>%
    mutate(Time = as.factor(Time))

  long_plot <- ggplot(gl_long, aes(x = bmi, y = prediction, color = predictor_type)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_smooth(method = "lm", se = FALSE, size = 1.25) +
    scale_color_manual(values = cluster_palette) +
    labs(title = glue("{omic_name_labels[[omic_name]]} Predicted BMI Values vs. Actual BMI"),
         x = "Actual Test Set BMI", y = "Omic Risk Score Predicted BMI", color = "Predictor") +
    theme_minimal()

  print(long_plot)
  # Save plot
  ggsave(
    filename = file.path(omic_out_dir, paste0("pred_vs_actual/", omic_name, "_prediction_plot_w_outliers.png")),
    plot = long_plot,
    width = 8, height = 6, dpi = 300
  )
}

for (omic_name in names(all_r2_dfs)) {
  r2_df <- all_r2_dfs[[omic_name]]
  anova_labels <- all_anova_labels[[omic_name]]
  
  # Merge with R2 table
  r2_annotated <- left_join(r2_df, anova_labels, 
                            by = c("Model" = "Compared_Model"))
  
  # Fill AIC for Basic model
  basic_model <- get("lmer_basic")
  basic_aic <- AIC(basic_model)
  
  r2_annotated <- r2_annotated %>%
    mutate(AIC = ifelse(Model == "Basic", round(basic_aic, 3), AIC))

  # Format table
  final_table <- r2_annotated %>%
    dplyr::select(-BIC, -logLik, -Df, -model_comparison, -`-2*log(L)`) %>%
    rename("RÂ²*" = R2m, "RÂ²c" = R2c, "Chi²" = Chisq, "P-Value" = `Pr(>Chisq)`) %>%
    gt() %>%
    tab_header(title = "Model Comparison (ANOVA)") %>%
    fmt_number(columns = where(is.numeric), decimals = 2) %>%
    sub_missing(missing_text = "-") %>%
    cols_align("center", columns = everything()) %>%
    cols_align("left", columns = Model) %>%
    cols_width(1 ~ pct(18),
               `RÂ²*` ~ pct(11), 
               `RÂ²c` ~ pct(11),
               `Chi²` ~ pct(11),
              Significance ~ pct(11),
              `P-Value` ~ pct(11), 
              AIC ~ pct(10)) %>%
    tab_options(table.width = pct(60), table.align = "center") %>%
    opt_table_font(font = list(google_font("Arial"), default_fonts())) %>%
    tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels())

  assign(paste0(omic_name, "_anova_table"), final_table)
  print(final_table)
}

delta_labels <- c(
  "Basic" = "Δ Basic",
  "Basic + Clinical" = "+ Δ Clinical",
  "Basic + GRS" = "+ Δ GRS",
  "Basic + Pathway" = "+ Δ Pathway",
  "Basic + Taxa" = "+ Δ Taxa",
  "Basic + MICOM" = "+ Δ MICOM",
  "Basic + Metabolomics" = "+ Δ Metabolomics",
  "Basic + Combined" = "+ Δ Combined")

long_labels <- c(
  "Basic" = " Basic",
  "Basic + Clinical" = "+ Clinical",
  "Basic + GRS" = "+ GRS",
  "Basic + Pathway" = "+ Pathway",
  "Basic + Taxa" = "+ Taxa",
  "Basic + MICOM" = "+ MICOM",
  "Basic + Metabolomics" = "+ Metabolomics",
  "Basic + Combined" = "+ Combined")

for (omic_name in names(all_r2_dfs)) {
  r2_df <- all_r2_dfs[[omic_name]]
  anova_labels <- all_anova_labels[[omic_name]]
  
  # Merge R2 with significance labels
  r2_annotated <- left_join(r2_df, anova_labels, by = c("Model" = "Compared_Model"))

  # Set custom model order
  custom_order <- c("Basic", "Basic + Clinical", "Basic + GRS", "Basic + Pathway",
                    "Basic + Taxa", "Basic + MICOM", 
                    "Basic + Metabolomics", "Basic + Combined")
  fill_vals <- setNames(sapply(custom_order, get_color), custom_order)
  r2_annotated$Model <- factor(r2_annotated$Model, levels = custom_order)

  # Fill in AIC for "Basic" model if missing
  basic_model <- get("lmer_basic")
  basic_aic <- AIC(basic_model)
  r2_annotated <- r2_annotated %>%
    mutate(AIC = ifelse(Model == "Basic", round(basic_aic, 3), AIC))

  # Extract R² for horizontal line
  basic_r2 <- r2_annotated$R2m[r2_annotated$Model == "Basic"]

  # Select correct label set
  label_set <- if (omic_name %in% c("gl_lasso_delta", "merf_delta")) {
    delta_labels
  } else if (omic_name %in% c("gl_lasso_long", "merf_long")) {
    long_labels
  } else {
    NULL  
  }
  if (omic_name %in% c("gl_lasso_delta", "merf_delta")) {
  y_lim <- c(0, 0.6)
  y_lim_r2c <- c(0, 0.75)
  } else {
  y_lim <- c(0, 0.3)
  y_lim_r2c <- c(0, 0.92)
  }

  # Bar plot
  sig_plot <- ggplot(r2_annotated, aes(x = Model, y = R2m, fill = Model)) +
    geom_col(width = 0.7, color = "white", size = 0.5) +
    geom_text(aes(y = R2m + 0.03, label = Significance), size = 8, color = "black", fontface = "bold") +
    scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = fill_vals) +
    scale_x_discrete(labels = label_set) +  
    #labs(title = glue("{omic_name_labels[[omic_name]]}"),
    labs(title = NULL,
         x = NULL,
         y = "Marginal R\u00B2") +
    combined_plot_theme(base_size = 24)

  print(sig_plot)
  assign(paste0(omic_name, "_sig_plot"), sig_plot)
  
  # Create new plot with both R2m and R2c bars
  r2_annotated_long <- r2_annotated %>%
    pivot_longer(cols = c(R2m, R2c), names_to = "R2_type", values_to = "R2_value") %>%
    mutate(R2_type = factor(R2_type, levels = c("R2m", "R2c"), labels = c("R²m", "R²c")))
  
  sig_plot_r2c <- ggplot(r2_annotated_long, aes(x = Model, y = R2_value, fill = Model, pattern = R2_type)) +
    ggpattern::geom_col_pattern(
      position = position_dodge(preserve = "single"),
      width = 0.7,
      color = "white", 
      size = 0.5,
      pattern_fill = "white",
      pattern_color = NA,
      pattern_density = 0.4,
      pattern_spacing = 0.02
    ) +
    scale_pattern_manual(values = c("R²m" = "none", "R²c" = "stripe")) +
    geom_text(aes(y = R2_value + 0.03, label = ifelse(R2_type == "R²m", Significance, ""), hjust = 1.2), 
              position = position_dodge(preserve = "single"),
              size = 8, color = "black", fontface = "bold") +
    scale_y_continuous(limits = y_lim_r2c, expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = fill_vals) +
    scale_x_discrete(labels = label_set) +
    labs(title = NULL,
         x = NULL,
         y = "R²") +
    combined_plot_theme(base_size = 24)
  
  print(sig_plot_r2c)
  assign(paste0(omic_name, "_sig_plot_r2c"), sig_plot_r2c)
}

# Assuming you saved the cleaned data frames *before* turning them into gt tables:
all_anova_clean_tables <- list(
  merf_long    = merf_long_anova_table,
  gl_lasso_long = gl_lasso_long_anova_table,
  merf_delta    = merf_delta_anova_table,
  gl_lasso_delta = gl_lasso_delta_anova_table)

raw_tables <- lapply(all_anova_clean_tables, function(gt_obj) gt_obj[["_data"]])

combined_anova_df <- bind_rows(
  lapply(names(raw_tables), function(name) {
    df <- raw_tables[[name]]
    df$Model_Type <- name
    df
    }),
  .id = "Source")

# Combined ANOVA table
final_anova_table_initial <- combined_anova_df %>% 
  dplyr::select(-c("Source")) %>% 
  mutate(Model_Type = dplyr::case_when(
    Model_Type == "merf_long" ~ "MERF Longitudinal BMI prediction",
    Model_Type == "gl_lasso_long" ~ "GlmmLasso BMI Prediction",
    Model_Type == "merf_delta" ~ "MERF ΔBMI Prediction)",
    Model_Type == "gl_lasso_delta" ~ "GlmmLasso ΔBMI Prediction)",
    TRUE ~ Model_Type)) %>% 
  gt(groupname_col = "Model_Type") %>%
    tab_header(title = "Additive Model Comparison's (ANOVA)") %>%
    fmt_number(columns = where(is.numeric), decimals = 3) %>%
    sub_missing(missing_text = "-") %>%
    cols_width(1 ~ pct(14),
               `RÂ²*` ~ pct(9), 
               `RÂ²c` ~ pct(9),
               `Chi²` ~ pct(9),
               Significance ~ pct(9),
               `P-Value` ~ pct(9), 
               AIC ~ pct(9), 
               Model_Type ~ pct(12)) %>%
    cols_align(align = "center", columns = everything()) %>% 
    cols_align("left", columns = Model) %>%
    tab_options(table.width = pct(90), table.align = "center") %>%
    opt_table_font(font = list(google_font("Arial"), default_fonts())) %>%
    tab_style(style = cell_text(weight = "bold"), 
              locations = cells_column_labels()) %>% 
    tab_style(style = list(
    cell_fill(color = "lightgray"),
    cell_text(weight = "bold",
              size = 12)),
  locations = cells_row_groups())

final_anova_table_initial
gtsave(final_anova_table_initial, 
       filename = paste0(omic_out_dir, "/anova_tables_02/anova_basic_table_w_outliers.png"))

# Function to create Word documents from tables
create_word_doc <- function(tables, filename) {
  doc <- read_docx()
  for (name in names(tables)) {
    table_data <- tables[[name]]
    
    ft <- flextable(table_data) %>%
      autofit() %>%
      theme_vanilla() %>%
      set_caption(caption = name) %>%
      fontsize(size = 10, part = "all") %>%
      align(align = "center", part = "all") %>%
      border_remove() %>%
      border_outer(border = fp_border(color = "black", width = 1)) %>%
      border_inner(border = fp_border(color = "gray", width = 0.5))

    doc <- doc %>%
      body_add_par(name, style = "heading 2") %>%
      body_add_flextable(ft) %>%
      body_add_par("", style = "Normal")
  }
  print(doc, target = paste0(omic_out_dir, "/anova_tables_02/", filename))
}

create_word_doc(raw_tables, "anova_basic_combined_w_outliers.docx")

omic_dfs <- list(
  merf_long = merf_long,
  gl_lasso_long = gl_lasso_long,
  merf_delta = merf_delta, 
  gl_lasso_delta = gl_lasso_delta)

omic_name_labels <- c(
  "merf_long"      = "MERF",
  "gl_lasso_long"  = "GlmmLasso",
  "merf_delta"     = "MERF",
  "gl_lasso_delta" = "GlmmLasso")

# Fit models and get R squared 
for (omic_name in names(omic_dfs)) {
  omic <- omic_dfs[[omic_name]]
  mod_dat <- omic
  
  # Fit models
  lmer_basic_b       <- lmer(bmi ~ y_new_basic + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_basic_meta_b  <- lmer(bmi ~ y_new_basic + y_new_meta + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_meta_b        <- lmer(bmi ~ y_new_meta + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_meta_grs      <- lmer(bmi ~ y_new_meta + y_new_grs + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_meta_micom    <- lmer(bmi ~ y_new_meta + y_new_micom + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_meta_pathway  <- lmer(bmi ~ y_new_meta + y_new_pathway + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_meta_taxa     <- lmer(bmi ~ y_new_meta + y_new_taxa + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_meta_metabo   <- lmer(bmi ~ y_new_meta + y_new_metabo + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_meta_all      <- lmer(bmi ~ y_new_meta + y_new_all_nclin + (1|Cluster), data = mod_dat, REML = FALSE)

  models <- list(
    `Clinical+`           = lmer_meta_b,
    `Clinical+ + GRS`     = lmer_meta_grs,
    `Clinical+ + MICOM`   = lmer_meta_micom,
    `Clinical+ + Pathway` = lmer_meta_pathway,
    `Clinical+ + Taxa`    = lmer_meta_taxa,
    `Clinical+ + Metabolomics`  = lmer_meta_metabo,
    `Clinical+ + Combined*`     = lmer_meta_all)

  r2_df_clin_plus <- do.call(rbind, lapply(names(models), function(name) {
    r2_vals_clin_plus <- r.squaredGLMM(models[[name]])
    data.frame(Model = name,
               R2m = r2_vals_clin_plus[1, "R2m"],
               R2c = r2_vals_clin_plus[1, "R2c"])
  }))
  assign(paste0(omic_name, "_r2_df"), r2_df_clin_plus)
}

for (omic_name in names(omic_dfs)) {
  # Get RÂ² results and model
  r2_df_clin_plus <- get(paste0(omic_name, "_r2_df"))
  omic <- omic_dfs[[omic_name]]
  mod_dat <- omic
  
  # Fit models for this specific omic_name
  lmer_meta_b        <- lmer(bmi ~ y_new_meta + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_meta_grs      <- lmer(bmi ~ y_new_meta + y_new_grs + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_meta_micom    <- lmer(bmi ~ y_new_meta + y_new_micom + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_meta_pathway  <- lmer(bmi ~ y_new_meta + y_new_pathway + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_meta_taxa     <- lmer(bmi ~ y_new_meta + y_new_taxa + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_meta_metabo   <- lmer(bmi ~ y_new_meta + y_new_metabo + (1|Cluster), data = mod_dat, REML = FALSE)
  lmer_meta_all      <- lmer(bmi ~ y_new_meta + y_new_all_nclin + (1|Cluster), data = mod_dat, REML = FALSE)

  lmer_models <- list(
    c("lmer_meta_b", "lmer_meta_b"),
    c("lmer_meta_b", "lmer_meta_grs"),
    c("lmer_meta_b", "lmer_meta_taxa"),
    c("lmer_meta_b", "lmer_meta_micom"),
    c("lmer_meta_b", "lmer_meta_pathway"),
    c("lmer_meta_b", "lmer_meta_metabo"),
    c("lmer_meta_b", "lmer_meta_all"))

  # Run ANOVA comparisons
  lmer_models_clin_p <- list()
  for (model_pair in lmer_models) {
    model_1 <- get(model_pair[1])
    model_2 <- get(model_pair[2])
    anova_out <- anova(model_1, model_2)
    anova_out$model_comparison <- paste(model_pair[1], "vs", model_pair[2])
    lmer_models_clin_p[[length(lmer_models_clin_p) + 1]] <- anova_out
  }

  anova_table <- do.call(rbind, lmer_models_clin_p) %>%
    filter(!is.na(`Chisq`) & !is.na(`Pr(>Chisq)`)) %>%
    mutate(across(where(is.numeric), round, 3)) %>%
    select(-npar)

  rownames(anova_table) <- NULL
  anova_table$model_1 <- gsub(" vs .*", "", anova_table$model_comparison)
  anova_table$model_2 <- gsub(".* vs ", "", anova_table$model_comparison)

  # Add significance & rename models
  anova_labels <- anova_table %>%
    mutate(
      Compared_Model = recode(model_comparison,
        "lmer_meta_b vs lmer_meta_b"     = "Clinical+",
        "lmer_meta_b vs lmer_meta_grs"   = "Clinical+ + GRS",
        "lmer_meta_b vs lmer_meta_taxa"  = "Clinical+ + Taxa",
        "lmer_meta_b vs lmer_meta_micom" = "Clinical+ + MICOM",
        "lmer_meta_b vs lmer_meta_pathway" = "Clinical+ + Pathway",
        "lmer_meta_b vs lmer_meta_metabo" = "Clinical+ + Metabolomics",
        "lmer_meta_b vs lmer_meta_all"     = "Clinical+ + Combined*"),
      Significance = case_when(
        `Pr(>Chisq)` < 0.001 ~ "***",
        `Pr(>Chisq)` < 0.01  ~ "**",
        `Pr(>Chisq)` < 0.05  ~ "*",
        TRUE                 ~ ""),
      `Pr(>Chisq)` = case_when(
        `Pr(>Chisq)` < 0.001 ~ "<0.001",
        `Pr(>Chisq)` < 0.01  ~ "<0.01",
        TRUE ~ as.character(`Pr(>Chisq)`))) %>%
    select(Compared_Model, AIC, Chisq, `Pr(>Chisq)`, Significance)

  custom_order <- c("Clinical+", 
                    "Clinical+ + GRS", 
                    "Clinical+ + Pathway", 
                    "Clinical+ + Taxa", 
                    "Clinical+ + MICOM", 
                    "Clinical+ + Metabolomics",
                    "Clinical+ + Combined*")
  
  clinical_aic <- AIC(lmer_meta_b)
  r2_annotated <- r2_df_clin_plus %>%
    left_join(anova_labels, by = c("Model" = "Compared_Model")) %>%
    mutate(Model = factor(Model, levels = custom_order)) %>%
    arrange(Model) %>%
    mutate(AIC = ifelse(Model == "Clinical+", 
                        round(clinical_aic, 3), AIC)) %>%
    rename("RÂ²" = R2m, "RÂ²c" = R2c, "Chi²" = Chisq, "P-Value" = `Pr(>Chisq)`) %>%
    mutate(across(c(`P-Value`, Significance), ~replace_na(., "-")))

  final_anova_table <- r2_annotated %>%
    gt() %>%
    tab_header(title = "Model Comparison (ANOVA)") %>%
    fmt_number(columns = where(is.numeric), decimals = 2) %>%
    sub_missing(missing_text = "-") %>%
    cols_align(align = "center", columns = everything()) %>%
    cols_width(Model ~ pct(18),
               `RÂ²` ~ pct(11),
               `RÂ²c` ~ pct(11),
               `Chi²` ~ pct(11),
               Significance ~ pct(11),
               `P-Value` ~ pct(11),
               AIC ~ pct(11)) %>%
    tab_options(table.width = pct(60), table.align = "center") %>%
    opt_table_font(font = list(google_font("Arial"), default_fonts())) %>%
    tab_style(style = cell_text(size = 14), locations = cells_body()) %>%
    tab_style(style = cell_text(size = 16, weight = "bold"), locations = cells_column_labels()) %>%
    tab_style(style = cell_text(size = 18, weight = "bold"), locations = cells_title())

  assign(paste0(omic_name, "_anova_table_clin_p"), final_anova_table)
  assign(paste0(omic_name, "_r2_annotated_clin_plus"), r2_annotated)
}

# Exact x-axis labels
delta_labels <- c(
  "Clinical+" = "Δ Clinical+",
  "Clinical+ + GRS" = "+ Δ GRS",
  "Clinical+ + Pathway" = "+ Δ Pathway",
  "Clinical+ + Taxa" = "+ Δ Taxa",
  "Clinical+ + MICOM" = "+ Δ MICOM",
  "Clinical+ + Metabolomics" = "+ Δ Metabolomics",
  "Clinical+ + Combined*" = "+ Δ Combined*")

long_labels <- c(
  "Clinical+" = " Clinical+",
  "Clinical+ + GRS" = "+ GRS",
  "Clinical+ + Pathway" = "+ Pathway",
  "Clinical+ + Taxa" = "+ Taxa",
  "Clinical+ + MICOM" = "+ MICOM",
  "Clinical+ + Metabolomics" = "+ Metabolomics",
  "Clinical+ + Combined*" = "+ Combined*")

# Fixed order used across all plots
custom_order <- c(
  "Clinical+",
  "Clinical+ + GRS",
  "Clinical+ + Pathway",
  "Clinical+ + Taxa",
  "Clinical+ + MICOM",
  "Clinical+ + Metabolomics",
  "Clinical+ + Combined*")

# Build an explicit fill palette mapped to these levels.
# Clinical+ baseline in gray; others pulled from your ft_colors.
fill_vals <- c(
  "Clinical+"                 = "#808080",
  "Clinical+ + GRS"           = ft_colors[["GRS"]],
  "Clinical+ + Pathway"       = ft_colors[["Pathway"]],
  "Clinical+ + Taxa"          = ft_colors[["Taxa"]],
  "Clinical+ + MICOM"         = ft_colors[["MICOM"]],
  "Clinical+ + Metabolomics"  = ft_colors[["Metabolomics"]],
  "Clinical+ + Combined*"     = ft_colors[["Combined*"]])

# Safety: ensure no missing colors for expected levels
missing_cols <- setdiff(custom_order, names(fill_vals))
if (length(missing_cols) > 0) {
  stop(paste("fill_vals missing names:", paste(missing_cols, collapse = ", ")))
}

for (omic_name in names(omic_dfs)) {
  r2_annotated <- get(paste0(omic_name, "_r2_annotated_clin_plus"))

  # Harmonize factor levels and clean significance labels
  r2_annotated <- r2_annotated %>%
    dplyr::mutate(
      Model = factor(Model, levels = custom_order),
      Significance = dplyr::na_if(Significance, "-"),
      Significance = dplyr::na_if(Significance, ""),
      Significance = dplyr::na_if(Significance, "–"),
      Significance = ifelse(is.na(Significance), "", Significance))

  # Axis labels per dataset family and y limits
  label_set <- if (omic_name %in% c("gl_lasso_delta", "merf_delta")) {
    delta_labels
  } else {
    long_labels
  }
  y_lim <- if (omic_name %in% c("gl_lasso_delta", "merf_delta")) c(0, 0.73) else c(0, 0.4)
  y_lim_r2c <- if (omic_name %in% c("gl_lasso_delta", "merf_delta")) c(0, 0.8) else c(0, 0.95)

  # Final plot: exact manual scale with named palette, no recycling
  sig_plot <- ggplot(r2_annotated, aes(x = Model, y = `RÂ²`, fill = Model)) +
    geom_col(width = 0.7, color = "white", size = 0.5) +
    geom_text(aes(y = `RÂ²` + 0.036, label = Significance),
              size = 10, color = "black", fontface = "bold") +
    scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = fill_vals, drop = FALSE) +
    scale_x_discrete(labels = label_set) +
    labs(title = NULL, x = NULL, y = "Marginal R\u00B2") +
    combined_plot_theme(base_size = 24)

  print(sig_plot)
  assign(paste0(omic_name, "_sig_plot_clin_plus"), sig_plot)
  
  # Create new plot with both R2m and R2c bars
  r2_annotated_long <- r2_annotated %>%
    pivot_longer(cols = c(`RÂ²`, `RÂ²c`), names_to = "R2_type", values_to = "R2_value") %>%
    mutate(R2_type = factor(R2_type, levels = c("RÂ²", "RÂ²c"), labels = c("R²m", "R²c")))
  
  sig_plot_r2c <- ggplot(r2_annotated_long, aes(x = Model, y = R2_value, fill = Model, pattern = R2_type)) +
    ggpattern::geom_col_pattern(
      position = position_dodge(preserve = "single"),
      width = 0.7,
      color = "white", 
      size = 0.5,
      pattern_fill = "white",
      pattern_color = NA,
      pattern_density = 0.4,
      pattern_spacing = 0.02) +
    scale_pattern_manual(values = c("R²m" = "none", "R²c" = "stripe")) +
    geom_text(aes(y = R2_value + 0.036, label = ifelse(R2_type == "R²m", Significance, ""), hjust = 1.2), 
              position = position_dodge(preserve = "single"),
              size = 10, color = "black", fontface = "bold") +
    scale_y_continuous(limits = y_lim_r2c, expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = fill_vals, drop = FALSE) +
    scale_x_discrete(labels = label_set) +
    labs(title = NULL, x = NULL, y = "R²") +
    combined_plot_theme(base_size = 24)
  
  print(sig_plot_r2c)
  assign(paste0(omic_name, "_sig_plot_clin_plus_r2c"), sig_plot_r2c)
}

# Assuming you saved the cleaned data frames *before* turning them into gt tables:
all_anova_clean_tables_clin_p <- list(
  merf_long    = merf_long_anova_table_clin_p,
  gl_lasso_long = gl_lasso_long_anova_table_clin_p,
  merf_delta    = merf_delta_anova_table_clin_p,
  gl_lasso_delta = gl_lasso_delta_anova_table_clin_p)

raw_tables_clin_p <- lapply(all_anova_clean_tables_clin_p, 
                            function(gt_obj) gt_obj[["_data"]])

combined_anova_df_clin_p <- bind_rows(
  lapply(names(raw_tables_clin_p), function(name) {
    df <- raw_tables_clin_p[[name]]
    df$Model_Type <- name
    df
    }),
  .id = "Source")

# Step 1: Prepare and clean data
combined_anova_df_clin_p <- combined_anova_df_clin_p %>%
  dplyr::mutate(Model_Type = case_when(
    Model_Type == "merf_long"      ~ "D. MERF BMI prediction",
    Model_Type == "gl_lasso_long"  ~ "B. GlmmLasso BMI Prediction",
    Model_Type == "merf_delta"     ~ "C. MERF ΔBMI Prediction",
    Model_Type == "gl_lasso_delta" ~ "A. GlmmLasso ΔBMI Prediction",
    TRUE ~ Model_Type))

table_data_clin_p <- combined_anova_df_clin_p %>%
  dplyr::select(-Source) %>%
  arrange(Model_Type) %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  dplyr::select(-Model_Type)

# Split by model type
tables_clin_p <- split(table_data_clin_p, combined_anova_df_clin_p$Model_Type)

# Combined ANOVA table
final_anova_table_clin_p2 <- combined_anova_df_clin_p %>% 
  dplyr::select(-c("Source")) %>% 
  mutate(Model_Type = dplyr::case_when(
    Model_Type == "merf_long" ~ "MERF BMI Prediction",
    Model_Type == "gl_lasso_long" ~ "GlmmLasso BMI Prediction",
    Model_Type == "merf_delta" ~ "MERF ΔBMI Prediction",
    Model_Type == "gl_lasso_delta" ~ "GlmmLasso ΔBMI Prediction",
    TRUE ~ Model_Type)) %>% 
  gt(groupname_col = "Model_Type") %>%
    tab_header(title = "Additive Model Comparison's (ANOVA)") %>%
    fmt_number(columns = where(is.numeric), decimals = 3) %>%
    sub_missing(missing_text = " ") %>%
    cols_width(1 ~ pct(14),
               `RÂ²` ~ pct(9), 
               `RÂ²c` ~ pct(9),
               `Chi²` ~ pct(9),
               Significance ~ pct(9),
               `P-Value` ~ pct(9), 
               AIC ~ pct(9), 
               Model_Type ~ pct(12)) %>%
    cols_align(align = "center",
      columns = everything()) %>% 
    cols_align("left", columns = Model) %>%
    tab_options(table.width = pct(90),
                table.align = "center") %>%
    opt_table_font(font = list(
        google_font("Arial"),   # fallback to system font if unavailable
        default_fonts())) %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels()) %>% 
  tab_style(style = list(
    cell_fill(color = "lightgray"),
    cell_text(weight = "bold",
              size = 12)),
  locations = cells_row_groups())

final_anova_table_clin_p2
gtsave(final_anova_table_clin_p2, 
       filename = paste0(omic_out_dir, "/anova_tables_02/anova_table_clin_w_outliers.png"))

create_word_doc(tables_clin_p, "anova_clinical_combined_w_outliers.docx")

# Function to add % RÂ² Increase column
add_r2_increase <- function(tbl) {
  base_r2 <- tbl$`RÂ²`[1]
  r2_increase <- c(NA, round(100 * (tbl$`RÂ²`[-1] - base_r2) / base_r2, 2))
  tibble::add_column(tbl, `% RÂ² Increase` = r2_increase, .after = "RÂ²")
}

# Create tables with RÂ² increase and save to Word document
tables_clin_p_plus <- lapply(tables_clin_p, add_r2_increase)
create_word_doc(tables_clin_p_plus, "anova_clinical_combined_w_outliers_increase.docx")

# Define plot names for later use
plot_names <- glue("{names(omic_dfs)}_sig_plot")
plot_names_clin_plus <- glue("{names(omic_dfs)}_sig_plot_clin_plus")

all_long_tables <- list(
  merf_long    = merf_long_anova_table,
  gl_lasso_long = gl_lasso_long_anova_table,
  merf_long_c    = merf_long_anova_table_clin_p,
  gl_lasso_long_c = gl_lasso_long_anova_table_clin_p)

all_delta_tables <- list(
  merf_delta    = merf_delta_anova_table,
  gl_lasso_delta = gl_lasso_delta_anova_table,
  merf_delta_c    = merf_delta_anova_table_clin_p,
  gl_lasso_delta_c = gl_lasso_delta_anova_table_clin_p)

raw_tables_long <- lapply(all_long_tables, function(gt_obj) gt_obj[["_data"]])
raw_tables_delta <- lapply(all_delta_tables, function(gt_obj) gt_obj[["_data"]])

round_df <- function(df, digits = 3) {
  df %>% mutate(across(where(is.numeric), ~ round(.x, digits)))
}

raw_tables_long <- lapply(all_long_tables, function(gt_obj) {
  round_df(gt_obj[["_data"]], digits = 3)
})

raw_tables_delta <- lapply(all_delta_tables, function(gt_obj) {
  round_df(gt_obj[["_data"]], digits = 3)
})

create_word_doc(raw_tables_long, "anova_long.docx")
create_word_doc(raw_tables_delta, "anova_delta.docx")

merf_dfs_long <- list(Basic = basic, 
                      `Clinical+` = meta, 
                      `Clinical` = meta_no_age_sex,
                      GRS = grs,
                      Taxa = taxa, 
                      MICOM = micom, 
                      Pathway = pathway, 
                      Metabolomics = metabo, 
                      `Combined` = all_no_age_sex,
                      `Combined*` = all_no_clin)

merf_dfs_delta <- list(Basic = basic_md, 
                       `Clinical+` = meta_md, 
                       `Clinical` = meta_no_age_sex_md,
                       GRS = grs_md, 
                       Taxa = taxa_md, 
                       MICOM = micom_md, 
                       Pathway = pathway_md, 
                       Metabolomics = metabo_md, 
                      `Combined` = all_no_age_sex_md,
                      `Combined*` = all_no_clin_md)

# Combined feature sets
all_ft_dfs <- list(merf_lg = list(data = merf_dfs_long, 
                                  colors = ft_colors, 
                                  label = "Long", 
                                  model_label = "MERF",
                                  x_ax_label = "Feature Importance"),
                   
                   merf_delta = list(data = merf_dfs_delta, 
                                     colors = ft_colors, 
                                     label = "Delta", 
                                     model_label = "MERF",
                                     x_ax_label = "Feature Importance"))

merf_plots <- list()
for (type in names(all_ft_dfs)) {
  ft_set <- all_ft_dfs[[type]]$data
  color_map <- all_ft_dfs[[type]]$colors
  label <- all_ft_dfs[[type]]$label
  model_label <- all_ft_dfs[[type]]$model_label  # "MERF" or "GlmmLasso"
  x_ax_label <- all_ft_dfs[[type]]$x_ax_label
  
  output_dir <- file.path(omic_out_dir, paste0("02_ft_imp_", model_label, "_", label))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (ft_name in names(ft_set)) {
    omic <- ft_set[[ft_name]]
    
    omic_parsed <- omic %>%
      mutate(Top_15_Feature_Importances = gsub("'", '"', Top_15_Feature_Importances)) %>%
      mutate(parsed = map(Top_15_Feature_Importances, fromJSON))
    
    omic_long <- omic_parsed %>%
      unnest(parsed) %>%
      dplyr::select(-Top_15_Feature_Importances) %>%
      dplyr::filter(Feature != "time") %>% 
      dplyr::filter(Model == "MSE Model")
    
    avg_importance <- omic_long %>%
      group_by(Feature) %>%
      summarise(avg_importance = mean(Importance, na.rm = TRUE), .groups = "drop") %>%
      slice_max(order_by = avg_importance, n = 10, with_ties = FALSE) %>%
      dplyr::filter(avg_importance != 0)
    
    avg_imp_name <- paste0("av_imp_", type, "_", ft_name, "_", label)
    assign(avg_imp_name, avg_importance, envir = .GlobalEnv)
    
    max_val <- max(avg_importance$avg_importance, na.rm = TRUE)
    x_limit <- ceiling(max_val / 0.05) * 0.05
 
    plot <- ggplot(avg_importance, aes(x = 0, y = reorder(Feature, avg_importance))) +
    geom_segment(aes(xend = avg_importance, 
                     yend = reorder(Feature, avg_importance)),
                     lineend = "round", linewidth = 12, 
                     color = color_map[[ft_name]]) +
    labs(title = paste(model_label, " ~ ",toupper(ft_name)),
       x = paste(" ", x_ax_label),
       y = NULL) +
    combined_plot_theme(base_size = 18) + 
    scale_x_continuous(limits = c(0, x_limit), 
                       expand = expansion(mult = c(0, 0.02))) 
    
    # Save to named R object and png
    object_name <- paste0("plot_", type, "_", ft_name, "_", label)
    assign(object_name, plot, envir = .GlobalEnv)
    merf_plots[[object_name]] <- plot
    png_name <- file.path(output_dir, paste0("feature_importance_", ft_name, "_", label, "_w_outliers.png"))
    ggsave(filename = png_name, plot = plot, width = 12, height = 5, dpi = 300)
  }
}

lass_dfs_long <- list(Basic = gl_ftimp_long_basic, 
                      `Clinical+` = gl_ftimp_long_meta, 
                      `Clinical` = gl_ftimp_long_meta_no_age_sex,
                      GRS = gl_ftimp_long_grs,
                      Taxa = gl_ftimp_long_taxa, 
                      MICOM = gl_ftimp_long_micom, 
                      Pathway = gl_ftimp_long_pathway, 
                      Metabolomics = gl_ftimp_long_metabo, 
                      `Combined` = gl_ftimp_long_all_no_age_sex,
                      `Combined*` = gl_ftimp_long_all_no_clin)

lass_dfs_delta <- list(Basic = gl_ftimp_delta_basic, 
                       `Clinical+` = gl_ftimp_delta_meta, 
                       `Clinical` = gl_ftimp_delta_meta_no_age_sex,
                       GRS = gl_ftimp_delta_grs,
                       Taxa = gl_ftimp_delta_taxa, 
                       MICOM = gl_ftimp_delta_micom, 
                       Pathway = gl_ftimp_delta_pathway, 
                       Metabolomics = gl_ftimp_delta_metabo, 
                       `Combined` = gl_ftimp_delta_all_no_age_sex,
                       `Combined*` = gl_ftimp_delta_all_no_clin)

# Combined feature sets
all_lass_ft_dfs <- list(lass_lg = list(data = lass_dfs_long, colors = ft_colors, 
                                  label = "long", 
                                  model_label = "GlmmLasso",
                                  x_ax_label = "Feature Co-efficient"),
                        
                   lass_delta = list(data = lass_dfs_delta, colors = ft_colors, 
                                     label = "delta", 
                                     model_label = "GlmmLasso",
                                     x_ax_label = "Feature Co-efficient"))
imap_dfr(lass_dfs_long, ~ 
           summarise(.x, min = min(Estimate, na.rm = TRUE), max = max(Estimate, na.rm = TRUE)) %>%
           mutate(name = .y)) %>% bind_rows() %>% dplyr::select(name, min, max)

imap_dfr(lass_dfs_delta, ~ 
           summarise(.x, min = min(Estimate, na.rm = TRUE), max = max(Estimate, na.rm = TRUE)) %>%
           mutate(name = .y)) %>% bind_rows() %>% dplyr::select(name, min, max)

lass_plots_lg <- list()
lass_plots_delta <- list()

for (type in names(all_lass_ft_dfs)) {
  ft_set <- all_lass_ft_dfs[[type]]$data
  color_map <- all_lass_ft_dfs[[type]]$colors
  label <- all_lass_ft_dfs[[type]]$label
  model_label <- all_lass_ft_dfs[[type]]$model_label
  x_ax_label <- all_lass_ft_dfs[[type]]$x_ax_label

  for (ft_name in names(ft_set)) {
    top_features <- ft_set[[ft_name]]
    ymin <- -1.0
    ymax <- 1.1
    plot_feat <- top_features %>%
      dplyr::filter(!str_detect(tolower(Feature), "time")) %>%
      slice_head(n = 10) %>%
      dplyr::filter(Estimate != 0) %>%   
      mutate(Estimate_capped = pmin(pmax(Estimate, ymin), ymax))

    features <- ggplot(plot_feat,
       aes(x = 0,
           y = reorder(Feature, Estimate_capped))) +
      geom_segment(aes(xend = Estimate_capped,
                     yend = reorder(Feature, Estimate_capped)),
                 lineend = "round", linewidth = 12,
                 color = color_map[[ft_name]]) +
      labs(title = paste(model_label, toupper(ft_name)),
           y = NULL,
           x = paste("", x_ax_label)) +
      combined_plot_theme(base_size = 18) + 
      scale_y_discrete(expand = expansion(mult = c(0.05, 0.05)))

    object_name <- paste0("plot_", type, "_", ft_name, "_", label)
    assign(object_name, features, envir = .GlobalEnv)

    # Save to appropriate list
    if (type == "lass_lg") {
      lass_plots_lg[[ft_name]] <- features
    } else if (type == "lass_delta") {
      lass_plots_delta[[ft_name]] <- features
    }

    # PNG saving
    png_name <- file.path(omic_out_dir, "lass_ft_imp_02", 
                          paste0("feature_importance_", ft_name, "_", label, "_w_outliers.png"))
    ggsave(filename = png_name, plot = features, width = 9, height = 5, dpi = 300)
  }
}

library(patchwork)
bl <- c(mget(plot_names[c(1, 2)]), mget(plot_names_clin_plus[c(1, 2)]))
bl <- standardize_plot_theme(bl, skip_indices = c(1, 2))

# Add titles above each plot
bl[[1]] <- bl[[1]] + combined_plot_theme()
bl[[2]] <- bl[[2]] + combined_plot_theme()
bl[[3]] <- bl[[3]] + combined_plot_theme() 
bl[[4]] <- bl[[4]] + combined_plot_theme() 

left_col <- (bl[[1]] / plot_spacer() / bl[[3]] +
  plot_layout(heights = c(0.2, 0.02, 0.2))
)
right_col <- (bl[[2]] / plot_spacer() / bl[[4]]) + 
  plot_layout(heights = c(0.2, 0.02, 0.2))

p <- (left_col | right_col) +
  plot_layout(widths = c(0.48, 0.48)) +
  plot_annotation(tag_levels = "A", tag_suffix = ".") &
  theme(plot.tag = element_text(face = "bold", size = 35),
        plot.tag.position = c(0.1, 0.96))

ggsave(file.path(omic_out_dir, "20251110/figure_1_temp_w_outliers.png"),
  p,
  width = 40, height = 30, dpi = 500, bg = "white")

library(patchwork)
plot_names_r2c <- glue("{names(omic_dfs)}_sig_plot_r2c")
plot_names_clin_plus_r2c <- glue("{names(omic_dfs)}_sig_plot_clin_plus_r2c")

# Filter out empty strings and ensure objects exist
plot_names_r2c <- plot_names_r2c[nchar(plot_names_r2c) > 0]
plot_names_clin_plus_r2c <- plot_names_clin_plus_r2c[nchar(plot_names_clin_plus_r2c) > 0]

# Get the plots, filtering out any that don't exist
bl_r2c <- mget(plot_names_r2c[c(1, 2)], ifnotfound = list(NULL))
bl_clin_r2c <- mget(plot_names_clin_plus_r2c[c(1, 2)], ifnotfound = list(NULL))
bl <- c(bl_r2c, bl_clin_r2c)
bl <- bl[!sapply(bl, is.null)]
bl <- standardize_plot_theme(bl, skip_indices = c(1, 2))

# Add titles above each plot
bl[[1]] <- bl[[1]] + combined_plot_theme()
bl[[2]] <- bl[[2]] + combined_plot_theme()
bl[[3]] <- bl[[3]] + combined_plot_theme() 
bl[[4]] <- bl[[4]] + combined_plot_theme() 

left_col <- (bl[[1]] / plot_spacer() / bl[[3]] +
  plot_layout(heights = c(0.2, 0.02, 0.2))
)
right_col <- (bl[[2]] / plot_spacer() / bl[[4]]) + 
  plot_layout(heights = c(0.2, 0.02, 0.2))

p <- (left_col | right_col) +
  plot_layout(widths = c(0.48, 0.48)) +
  plot_annotation(tag_levels = "A", tag_suffix = ".") &
  theme(plot.tag = element_text(face = "bold", size = 35),
        plot.tag.position = c(0.1, 0.96))

ggsave(file.path(omic_out_dir, "20251110/figure_1_r2m_r2c.png"),
  p,
  width = 40, height = 30, dpi = 500, bg = "white")

library(patchwork)
bl <- c(merf_plots[c(3, 4, 7, 8, 9)], lass_plots_lg[c(2, 4, 9)])
bl <- standardize_plot_theme(bl)

merf_clinical <- bl[[1]]
merf_grs <- bl[[2]]
merf_pathway <- bl[[3]]
merf_metabolomics <- bl[[4]]
merf_combined_ <- bl[[5]]
glm_clinical <- bl[[6]]
glm_grs <- bl[[7]]
glm_metabolomics <- bl[[8]]

left_col <- (merf_clinical / plot_spacer() / 
             merf_grs / plot_spacer() /
             merf_pathway / plot_spacer() /
             merf_metabolomics) + 
  plot_layout(heights = c(0.10, 0.02, 
                          0.05, 0.02, 
                          0.12, 0.02,
                          0.12, 0.02))

right_col <- (merf_combined_ / plot_spacer() /
              glm_clinical / plot_spacer() / 
              glm_grs / plot_spacer() / 
              glm_metabolomics) + 
  plot_layout(heights = c(0.12, 0.02, 
                          0.10, 0.02, 
                          0.05, 0.02,
                          0.12, 0.02))

p <- (left_col | right_col) +
  plot_layout(widths = c(0.48, 0.48)) +
  plot_annotation(tag_levels = "A", tag_suffix = ".") &
  theme(plot.tag = element_text(face = "bold", size = 35),
        plot.tag.position = c(0.1, 0.96))

ggsave(file.path(omic_out_dir, "20251110/figure_2_temp_w_outliers.png"),
  p, width = 40, height = 30, dpi = 500, bg = "white")

library(patchwork)
bl <- c(mget(plot_names[c(3, 4)]), mget(plot_names_clin_plus[c(3, 4)]))
bl <- standardize_plot_theme(bl, skip_indices = c(1, 2))

# Add titles above each plot
#bl[[1]] <- bl[[1]] + combined_plot_theme() 
#bl[[2]] <- bl[[2]] + combined_plot_theme() 
#bl[[3]] <- bl[[3]] + combined_plot_theme() 
#bl[[4]] <- bl[[4]] + combined_plot_theme() 

left_col <- (bl[[1]] / plot_spacer() / bl[[3]] +
  plot_layout(heights = c(0.2, 0.02, 0.2)))
right_col <- (bl[[2]] / plot_spacer() / bl[[4]]) + 
  plot_layout(heights = c(0.2, 0.02, 0.2))

p <- (left_col | right_col) +
  plot_layout(widths = c(0.48, 0.48)) +
  plot_annotation(tag_levels = "A", tag_suffix = ".") &
  theme(plot.tag = element_text(face = "bold", size = 35),
    plot.tag.position = c(0.1, 0.96))

ggsave(file.path(omic_out_dir, "20251110/figure_3_tump_w_outiers.png"),
  p,width = 40, height = 30, dpi = 500, bg = "white")

library(patchwork)
# Filter out empty strings and ensure objects exist (if not already done)
if(!exists("plot_names_r2c")) {
  plot_names_r2c <- glue("{names(omic_dfs)}_sig_plot_r2c")
  plot_names_clin_plus_r2c <- glue("{names(omic_dfs)}_sig_plot_clin_plus_r2c")
  plot_names_r2c <- plot_names_r2c[nchar(plot_names_r2c) > 0]
  plot_names_clin_plus_r2c <- plot_names_clin_plus_r2c[nchar(plot_names_clin_plus_r2c) > 0]
}

# Get the plots, filtering out any that don't exist
bl_r2c <- mget(plot_names_r2c[c(3, 4)], ifnotfound = list(NULL))
bl_clin_r2c <- mget(plot_names_clin_plus_r2c[c(3, 4)], ifnotfound = list(NULL))
bl <- c(bl_r2c, bl_clin_r2c)
bl <- bl[!sapply(bl, is.null)]
bl <- standardize_plot_theme(bl, skip_indices = c(1, 2))

# Add titles above each plot
#bl[[1]] <- bl[[1]] + combined_plot_theme() 
#bl[[2]] <- bl[[2]] + combined_plot_theme() 
#bl[[3]] <- bl[[3]] + combined_plot_theme() 
#bl[[4]] <- bl[[4]] + combined_plot_theme() 

left_col <- (bl[[1]] / plot_spacer() / bl[[3]] +
  plot_layout(heights = c(0.2, 0.02, 0.2)))
right_col <- (bl[[2]] / plot_spacer() / bl[[4]]) + 
  plot_layout(heights = c(0.2, 0.02, 0.2))

p <- (left_col | right_col) +
  plot_layout(widths = c(0.48, 0.48)) +
  plot_annotation(tag_levels = "A", tag_suffix = ".") &
  theme(plot.tag = element_text(face = "bold", size = 35),
    plot.tag.position = c(0.1, 0.96))

ggsave(file.path(omic_out_dir, "20251110/figure_3_r2m_r2c.png"),
  p,width = 40, height = 30, dpi = 500, bg = "white")

library(patchwork)
# merf = all but grs and micom, lass = c* met cb cb*
bl <- c(merf_plots[c(13,15,17,18, 19,20)], lass_plots_delta[c(3, 8, 9, 10)])
bl <- standardize_plot_theme(bl)

merf_clinical <- bl[[1]]
merf_taxa <- bl[[2]]
merf_pathway <- bl[[3]]
merf_metabolomics <- bl[[4]]
merf_combined_ <- bl[[5]]
merf_combined_m <- bl[[6]]

glm_clinical <- bl[[7]]
glm_metabolomics <- bl[[8]]
glm_combined_ <- bl[[9]]
glm_combined_m <- bl[[10]]

left_col <- (merf_clinical / plot_spacer() / 
             merf_taxa / plot_spacer() / 
             merf_pathway / plot_spacer() / 
             merf_metabolomics / plot_spacer() / 
             merf_combined_) + 
   plot_layout(heights = c(0.08, 0.02, 
                           0.12, 0.02, 
                           0.12, 0.02,
                           0.12, 0.02,
                           0.12))
right_col <- (merf_combined_m / plot_spacer() / 
              glm_clinical / plot_spacer() / 
              glm_metabolomics / plot_spacer() / 
              glm_combined_ / plot_spacer() / 
              glm_combined_m) + 
   plot_layout(heights = c(0.13, 0.02, 
                           0.04, 0.02, 
                           0.13, 0.02,
                           0.13, 0.02, 
                           0.13))

p <- (left_col | right_col) +
  plot_layout(widths = c(0.48, 0.48)) +
  plot_annotation(tag_levels = "A", tag_suffix = ".") &
  theme(plot.tag = element_text(face = "bold", size = 35),
    plot.tag.position = c(0.1, 0.96))

ggsave(file.path(omic_out_dir, "20251110/figure_4_temp_w_outlier.png"),
  p,width = 40, height = 40, dpi = 500, bg = "white")

library(patchwork)
bl <- c(mget(plot_names[c(1,2)]), merf_plots[c(2, 3, 6, 7, 8)],  lass_plots_lg[c(2, 7)])
bl <- standardize_plot_theme(bl, skip_indices = c(1,2))

left_col <- (bl[[1]] / plot_spacer() / bl[[3]] / plot_spacer() / 
             bl[[4]] / plot_spacer() /  bl[[6]] / plot_spacer() / bl[[8]]) +
             plot_layout(heights = c(0.2, 0.02, 0.1, 0.02, 0.04, 0.02, 0.2, 0.02, 0.2))

right_col <- (bl[[2]] / plot_spacer() / bl[[5]] / plot_spacer() / 
              bl[[7]] / plot_spacer() / bl[[9]]) +
              plot_layout(heights = c(0.2, 0.02, 0.2, 0.02, 0.2, 0.02, 0.2))

p <- (left_col | right_col) +
  plot_layout(widths = c(0.49, 0.49)) +
  plot_annotation(tag_levels = "A", tag_suffix = ".") &
  theme(plot.tag = element_text(face = "bold", size = 35),
        plot.tag.position = c(0.1, 0.96))

ggsave(file.path(omic_out_dir, "aov_feats_combined/basic_long_two_col_patchwork_spacers.png"),
       p, width = 40, height = 40, dpi = 500, bg = "white")

library(patchwork)
bd <- c(mget(plot_names[c(3,4)]), merf_plots[c(10, 12, 15, 16)],  lass_plots_delta[c(2, 7)])
bd <- standardize_plot_theme(bd, skip_indices = c(1,2))

# Left column plots
left_col <- (bd[[1]] / plot_spacer() / bd[[3]] / plot_spacer() /
             bd[[6]] / plot_spacer() / bd[[8]]) + 
             plot_layout(heights = c(0.22, 0.02, 0.22, 0.02, 0.22, 0.02, 0.22))

# Right column plots
right_col <- (bd[[2]] / plot_spacer() / bd[[4]] / plot_spacer() /
              bd[[5]] / plot_spacer() / bd[[7]]) + 
              plot_layout(heights = c(0.215, 0.02, 0.215, 0.02, 0.215, 0.02, 0.215))

# Combine columns
pd <- (left_col | right_col) +
  plot_layout(widths = c(0.49, 0.49)) +
  plot_annotation(tag_levels = "A", tag_suffix = ".") & 
  theme(plot.tag = element_text(face = "bold", size = 35),
            plot.tag.position = c(0.1, 0.96))

# Save
ggsave(paste0(omic_out_dir, "/aov_feats_combined/basic_delta_test_patchwork.png"),
       pd, width = 40, height = 40, dpi = 500, bg = 'white')

cl <- c(mget(plot_names_clin_plus[c(1,2)]), merf_plots[c(7, 8)],  lass_plots_lg[c(3,7)])
cl <- standardize_plot_theme(cl, skip_indices = c(1,2))

# Left column plots
left_col <- (cl[[1]] / plot_spacer() / cl[[3]] / plot_spacer() / cl[[6]] / plot_spacer()) + 
             plot_layout(heights = c(0.2, 0.02, 0.2, 0.02, 0.2, 0.22))

# Right column plots
right_col <- (cl[[2]] / plot_spacer() / cl[[4]] / plot_spacer() / cl[[5]] / plot_spacer()) + 
              plot_layout(heights = c(0.2, 0.02, 0.2, 0.02, 0.2, 0.22))

# Combine columns
cl <- (left_col | right_col) +
  plot_layout(widths = c(0.48, 0.48)) +
  plot_annotation(tag_levels = "A", tag_suffix = ".") & 
  theme(plot.tag = element_text(face = "bold", size = 35),
            plot.tag.position = c(0.1, 0.96))

ggsave(paste0(omic_out_dir, "/aov_feats_combined/clin_long_test_patchwork.png"), cl, 
       width = 40, height = 40, dpi = 500, bg = 'white')

cd <- c(mget(plot_names_clin_plus[c(3,4)]), merf_plots[c(12, 15, 16)],  lass_plots_delta[c(7)])
cd <- standardize_plot_theme(cd, skip_indices = c(1,2))

# Left column plots
left_col <- (cd[[1]] / plot_spacer() / cd[[3]] / plot_spacer() / cd[[6]] / plot_spacer()) + 
             plot_layout(heights = c(0.2, 0.02, 0.2, 0.02, 0.2, 0.22))

# Right column plots
right_col <- (cd[[2]] / plot_spacer() / cd[[4]] / plot_spacer() / cd[[5]] / plot_spacer()) + 
              plot_layout(heights = c(0.2, 0.02, 0.2, 0.02, 0.2, 0.22))

# Combine columns
cd <- (left_col | right_col) +
  plot_layout(widths = c(0.49, 0.49)) +
  plot_annotation(tag_levels = "A", tag_suffix = ".") & 
  theme(plot.tag = element_text(face = "bold", size = 35),
            plot.tag.position = c(0.1, 0.96))

ggsave(paste0(omic_out_dir, "/aov_feats_combined/clin_delta_test_patchwork.png"), cd, 
       width = 40, height = 40, dpi = 500, bg = 'white')

library(patchwork)
bl <- c(mget(plot_names[c(1,2)]), shap_long_plot_list[c(2, 3, 5, 7, 8)],  lass_plots_lg[c(2, 7)])
bl <- standardize_plot_theme(bl, skip_indices = c(1,2))

left_col <- (bl[[1]] / plot_spacer() / bl[[3]] / plot_spacer() / 
             bl[[4]] / plot_spacer() /  bl[[6]] / plot_spacer() / bl[[8]]) +
             plot_layout(heights = c(0.2, 0.02, 0.1, 0.02, 0.04, 0.02, 0.2, 0.02, 0.2))

right_col <- (bl[[2]] / plot_spacer() / bl[[5]] / plot_spacer() / 
              bl[[7]] / plot_spacer() / bl[[9]]) +
              plot_layout(heights = c(0.2, 0.02, 0.2, 0.02, 0.2, 0.02, 0.2))

p <- (left_col | right_col) +
  plot_layout(widths = c(0.49, 0.49)) +
  plot_annotation(tag_levels = "A", tag_suffix = ".") &
  theme(plot.tag = element_text(face = "bold", size = 35),
        plot.tag.position = c(0.1, 0.96))

ggsave(file.path(omic_out_dir, "aov_feats_combined/shap_basic_long_two_col_patchwork.png"),
       p, width = 40, height = 40, dpi = 550, bg = "white")

bd <- c(mget(plot_names[c(3,4)]), shap_delta_plot_list[c(2, 4, 7, 8)],  lass_plots_delta[c(2, 7, 8)])
bd <- standardize_plot_theme(bd, skip_indices = c(1,2))

# Left column plots
left_col <- (bd[[1]]/ plot_spacer()/ bd[[3]]/ plot_spacer()/ bd[[6]]/ plot_spacer()/
             bd[[8]]/ plot_spacer()) + 
             plot_layout(heights = c(0.22, 0.02, 0.2, 0.02, 0.2, 0.02, 0.25, 0.01))

# Right column plots
right_col <- (bd[[2]]/ plot_spacer()/ bd[[4]]/ plot_spacer()/ bd[[5]]/ plot_spacer()/ 
              bd[[7]]/ plot_spacer()/ bd[[9]]) + 
              plot_layout(heights = c(0.23, 0.02, 0.2, 0.02, 
                                      0.2, 0.02, 0.05, 0.015, 0.18))

# Combine columns
pd <- (left_col | right_col) +
  plot_layout(widths = c(0.49, 0.49)) +
  plot_annotation(tag_levels = "A", tag_suffix = ".") & 
  theme(plot.tag = element_text(face = "bold", size = 35),
        plot.tag.position = c(0.1, 0.96))

# Save
ggsave(paste0(omic_out_dir, "/aov_feats_combined/shap_basic_delta_test_patchwork.png"),
       pd, width = 40, height = 40, dpi = 500, bg = 'white')

cl <- c(mget(plot_names_clin_plus[c(1,2)]), shap_long_plot_list[c(7, 8)],  lass_plots_lg[c(3,7)])
cl <- standardize_plot_theme(cl, skip_indices = c(1,2))

# Left column plots
left_col <- (cl[[1]] / plot_spacer() / cl[[3]] / plot_spacer() / cl[[6]] / plot_spacer()) + 
             plot_layout(heights = c(0.2, 0.02, 0.2, 0.02, 0.2, 0.25))

# Right column plots
right_col <- (cl[[2]] / plot_spacer() / cl[[4]] / plot_spacer() / cl[[5]] / plot_spacer()) + 
              plot_layout(heights = c(0.2, 0.02, 0.2, 0.02, 0.2, 0.25))

# Combine columns
cl <- (left_col | right_col) +
  plot_layout(widths = c(0.49, 0.49)) +
  plot_annotation(tag_levels = "A", tag_suffix = ".") & 
  theme(plot.tag = element_text(face = "bold", size = 35),
            plot.tag.position = c(0.1, 0.96))

ggsave(paste0(omic_out_dir, "/aov_feats_combined/shap_clin_long_test_patchwork.png"), cl, 
       width = 40, height = 40, dpi = 500, bg = 'white')

cd <- c(mget(plot_names_clin_plus[c(3,4)]), shap_delta_plot_list[c(4 , 7, 8)],  lass_plots_delta[c(7)])
cd <- standardize_plot_theme(cd, skip_indices = c(1,2))

# Left column plots
left_col <- (cd[[1]] / plot_spacer() / cd[[3]] / plot_spacer() / cd[[5]] / plot_spacer()) + 
             plot_layout(heights = c(0.2, 0.02, 0.2, 0.02, 0.2, 0.25))

# Right column plots
right_col <- (cd[[2]] / plot_spacer() / cd[[4]] / plot_spacer() / cd[[6]] / plot_spacer()) + 
              plot_layout(heights = c(0.2, 0.02, 0.2, 0.02, 0.2, 0.25))

# Combine columns
cd <- (left_col | right_col) +
  plot_layout(widths = c(0.495, 0.495)) +
  plot_annotation(tag_levels = "A", tag_suffix = ".") & 
  theme(plot.tag = element_text(face = "bold", size = 35),
            plot.tag.position = c(0.1, 0.96))

ggsave(paste0(omic_out_dir, "/aov_feats_combined/shap_clin_delta_test_patchwork.png"), cd, 
       width = 40, height = 40, dpi = 500, bg = 'white')

common_cols <- c("bmi", "y_new_basic", "y_new_grs",
                 "y_new_meta", "y_new_taxa", 
                 "y_new_micom", "y_new_pathway", 
                 "y_new_metabo", "y_new_all")

df_list_long <- list(merf_long = merf_long[common_cols],
                 gl_lasso_long = gl_lasso_long[common_cols])
df_list_delta <- list(merf_delta = merf_delta[common_cols],
                  gl_lasso_delta = gl_lasso_delta[common_cols])

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
                        "GRS Score" = "y_new_grs",
                        "Taxa Score" = "y_new_taxa",
                        "MICOM Score" = "y_new_micom",
                        "Pathway Score" = "y_new_pathway",
                        "Metabolomics Score" = "y_new_metabo",
                        "Combined Omic Score" = "y_new_all"),
      Var2 = fct_recode(Var2,
                        "Actual BMI" = "bmi",
                        "Basic Score" = "y_new_basic",
                        "Meta Score" = "y_new_meta",
                        "GRS Score" = "y_new_grs",
                        "Taxa Score" = "y_new_taxa",
                        "MICOM Score" = "y_new_micom",
                        "Pathway Score" = "y_new_pathway",
                        "Metabolomics Score" = "y_new_metabo",
                        "Combined Omic Score" = "y_new_all"))

final_heat <- ggplot(cor_long_df, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  facet_grid(method ~ type) +
  scale_fill_gradient2(low = "#09236c", high = "#781720",
                       mid = "white", midpoint = 0.0,
                       limit = c(-0.3, 1), name = "Correlation") +
  theme_minimal(base_size = 16) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15), 
        strip.text = element_text(face = "bold", size = 16)) +
  labs(title = "Correlation Heatmaps: Long vs Delta", x = "", y = "")

png_name <- file.path(omic_out_dir, paste0("final_heat_map_scores", ".png"))
ggsave(filename = png_name, plot = final_heat, width = 14, height = 10, dpi = 300)

final_heat

rename_map <- c(
  "Actual BMI" = "bmi",
  "Basic Score" = "y_new_basic",
  "Clinical Score" = "y_new_meta",
  "GRS Score" = "y_new_grs",
  "Taxa Score" = "y_new_taxa",
  "MICOM Score" = "y_new_micom",
  "Pathway Score" = "y_new_pathway",
  "Metabolomics Score" = "y_new_metabo",
  "Combined Omic Score" = "y_new_all")

score_colors <- c(
  "Basic Score" = ft_colors["Basic"],
  "Clinical Score" = ft_colors["Clinical"],
  "GRS Score" = ft_colors["GRS"],
  "Taxa Score" = ft_colors["Taxa"],
  "MICOM Score" = ft_colors["MICOM"],
  "Pathway Score" = ft_colors["Pathway"],
  "Metabolomics Score" = ft_colors["Metabolomics"],
  "Combined Omic Score" = ft_colors["Combined"],
  "Actual BMI" = "grey50")


merf_delta_bar <- merf_delta %>% rename(!!!rename_map) %>% 
                  filter(Time == 6)
merf_long_bar <- merf_long %>% rename(!!!rename_map) %>% 
                  filter(Time == 6)
gl_lasso_delta_bar <- gl_lasso_delta %>% rename(!!!rename_map) %>% 
                  filter(Time == 1)
gl_lasso_long_bar <- gl_lasso_long %>% rename(!!!rename_map) %>% 
                  filter(Time == 2)

# Consolidated function for plotting score means
plot_score_means <- function(df, title = "Mean Scores with Error Bars", y_limits = NULL) {
  summary_df <- df %>%
    select(`Basic Score`, `Clinical Score`, `GRS Score`, `Taxa Score`,
           `MICOM Score`, `Pathway Score`, `Metabolomics Score`,  
           `Combined Omic Score`, `Actual BMI`) %>%
    pivot_longer(cols = everything(), 
                 names_to = "Score Type", 
                 values_to = "Value") %>%
    group_by(`Score Type`) %>%
    summarise(Mean = mean(Value, na.rm = TRUE),
              SE = sd(Value, na.rm = TRUE) / sqrt(n()),
              .groups = "drop")

  # Set desired order of bars
  score_levels <- c(
    "Actual BMI", "Basic Score", "Clinical Score",
    "GRS Score", "Taxa Score", "MICOM Score",
    "Pathway Score", "Metabolomics Score", "Combined Omic Score"
  )
  summary_df$`Score Type` <- factor(summary_df$`Score Type`, levels = score_levels)

  # Create base plot
  p <- ggplot(summary_df, aes(x = `Score Type`, y = Mean, fill = `Score Type`)) +
    geom_col() +
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
    scale_fill_manual(values = score_colors) +
    labs(x = NULL, y = "Mean Score", title = title) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  # Add y-axis limits if specified
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  return(p)
}

p10 <- plot_score_means(merf_delta_bar, title = "MERF Delta: Mean Scores Â± SE")
p20 <- plot_score_means(merf_long_bar, title = "MERF Long: Mean Scores Â± SE", y_limits = c(25, 34))
p30 <- plot_score_means(gl_lasso_delta_bar, title = "GLMM Delta: Mean Scores Â± SE")
p40 <- plot_score_means(gl_lasso_long_bar, title = "GLMM Long: Mean Scores Â± SE", y_limits = c(25, 34))

combined_raw_scores_plot <- (p10 | p20) / (p30 | p40)
ggsave(paste0(omic_out_dir, "score_means_grid.png"), 
       combined_raw_scores_plot, width = 12, height = 10, dpi = 300)


