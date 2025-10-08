# plotting_utils

library(pacman)

p_load(tools, reticulate, viridis, tidyplots, patchwork, jsonlite, maps, ggvenn, 
       caret, caretEnsemble, glmnet, xgboost, ggplot2, glmmLasso, corrplot,
       readr, plyr, dplyr, tidyr, purrr, tibble, stringr, psych, randomForest,  
       reshape2, scales, gridExtra, plotly, sf, tidyverse, naniar, VIM, gridExtra,
       sjPlot, htmltools, officer, flextable, webshot, apaTables, MuMIn, lme4, 
       glue, grid, rsq, pheatmap, GGally, VennDiagram, glmmTMB, broom.mixed, gt,
       patchwork, tidyverse, ggbeeswarm, scales, viridis, magrittr)

##### FUNCTIONS ###########################################################

'%ni%' <- Negate('%in%')

r2_general <-function(preds,actual){ 
  return(1- sum((preds - actual) ^ 2)/sum((actual - mean(actual))^2))
}

##########################################################################

q_squared <- function(actual, predicted) {
  ss_res <- sum((actual - predicted)^2)
  ss_tot <- sum((actual - mean(actual))^2)
  1 - (ss_res / ss_tot)
}

##########################################################################

geom_top_rounded_rect <- function(mapping = NULL, data = NULL,
                                  stat = "identity", position = "identity",
                                  radius = grid::unit(6, "pt"),
                                  ...,
                                  na.rm = FALSE,
                                  show.legend = NA,
                                  inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomTopRoundedRect,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      radius = radius,
      na.rm = na.rm,
      ...
    )
  )
}

GeomTopRoundedRect <- ggplot2::ggproto(
  "GeomTopRoundedRect", ggplot2::Geom,
  
  default_aes = ggplot2::aes(
    colour = NA, fill = "grey35", size = 0.5, linetype = 1, alpha = NA
  ),
  
  required_aes = c("xmin", "xmax", "ymin", "ymax"),
  
  draw_panel = function(self, data, panel_params, coord,
                        radius = grid::unit(6, "pt")) {
    
    coords <- coord$transform(data, panel_params)
    
    grobs <- lapply(1:length(coords$xmin), function(i) {
      
      gridGeometry::polyclipGrob(
        grid::roundrectGrob(
          coords$xmin[i], coords$ymax[i],
          width = (coords$xmax[i] - coords$xmin[i]),
          height = (coords$ymax[i] - coords$ymin[i]),
          r = radius,
          default.units = "native",
          just = c("left", "top")
        ),
        grid::rectGrob(
          coords$xmin[i], coords$ymax[i] - (coords$ymax[i] - coords$ymin[i]) / 2,
          width = (coords$xmax[i] - coords$xmin[i]),
          height = (coords$ymax[i] - coords$ymin[i]) / 2,
          default.units = "native",
          just = c("left", "top")
        ),
        op = "union",
        gp = grid::gpar(
          col = coords$colour[i],
          fill = alpha(coords$fill[i], coords$alpha[i]),
          lwd = coords$size[i] * .pt,
          lty = coords$linetype[i],
          lineend = "butt"
        )
      )
    })
    
    grobs <- do.call(grid::gList, grobs)
    ggplot2:::ggname("geom_top_rounded_rect", grid::grobTree(children = grobs))
  },
  draw_key = ggplot2::draw_key_polygon)

##########################################################################

geom_top_rounded_col <- function(mapping = NULL, data = NULL,
                                 position = ggplot2::position_stack(reverse = TRUE),
                                 radius = grid::unit(3, "pt"), ..., width = NULL,
                                 na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(
    data = data, mapping = mapping, stat = "identity",
    geom = GeomTopRoundedCol, position = position, show.legend = show.legend,
    inherit.aes = inherit.aes, params = list(
      width = width, radius = radius, na.rm = na.rm, ...
    )
  )
}

GeomTopRoundedCol <- ggproto(
  "GeomTopRoundedCol", GeomTopRoundedRect,
  required_aes = c("x", "y"),
  
  setup_params = function(data, params) {
    params$flipped_aes <- has_flipped_aes(data, params)
    params
  },
  non_missing_aes = c("xmin", "xmax", "ymin", "ymax"),
  setup_data = function(data, params) {
    data$width <- data$width %||%
      params$width %||% (resolution(data$x, FALSE) * 0.9)
    transform(data,
              ymin = pmin(y, 0), ymax = pmax(y, 0),
              xmin = x - width / 2, xmax = x + width / 2, width = NULL
    )
  },
  draw_panel = function(self, data, panel_params, coord, width = NULL, radius = grid::unit(3, "pt")) {
    ggproto_parent(GeomTopRoundedRect, self)$draw_panel(data, panel_params, coord, radius = radius)
  }
)


##########################################################################

# Renaming function
rename_features_in_json <- function(json_str, feature_rename_map) {
  # Fix the invalid JSON format by replacing single quotes with double quotes
  json_str <- gsub("'", "\"", json_str)
  
  # Parse the fixed JSON string into a list
  feature_list <- jsonlite::fromJSON(json_str, simplifyDataFrame = FALSE)
  
  # Rename each feature if it's in the mapping
  renamed_list <- lapply(feature_list, function(item) {
    original_name <- item$Feature
    if (original_name %in% names(feature_rename_map)) {
      item$Feature <- feature_rename_map[[original_name]]
    }
    item
  })
  # Convert back to JSON
  toJSON(renamed_list, auto_unbox = TRUE)
}

##### BeSwarm Plot Function ################################################### 

plot_shap_beeswarm <- function(shap_df, plot_title, high_color) {
  shap_df %>%
    filter(Feature != "time") %>%
    
    # Step 1: Filter features based on SHAP value proportion
    group_by(Feature) %>%
    mutate(prop_small_shap = mean(SHAP > -0.005 & SHAP < 0.005, na.rm = TRUE)) %>%
    filter(prop_small_shap <= 0.5) %>%
    
    # Step 2: Reorder and normalize
    mutate(
      Feature = fct_reorder(Feature, abs(SHAP), .fun = max),
      Norm_Feature_Value = scales::rescale(Feature_Value, 
                                           from = range(Feature_Value, 
                                                        na.rm = TRUE))) %>%  
    #Norm_Feature_Value = rescale(Feature_Value, to = c(0, 1))) %>%
    ungroup() %>%
    
    # Step 3: Plot
    ggplot(aes(x = SHAP, y = Feature, color = Norm_Feature_Value)) +
    geom_quasirandom(alpha = 0.5, width = 0.1, size = 4.5, groupOnX = FALSE) +
    
    scale_color_gradientn(
      colors = c("#004E64", high_color),
      values = scales::rescale(c(0, 1)),
      name = "Feature Value",
      breaks = c(0, 1),
      labels = c("Low", "High")) +
    
    labs(title = plot_title,
         x = "SHAP Value (Impact on Model Output)",
         y = NULL) +
    
    scale_y_discrete(expand = expansion(add = 0.75)) +
    theme_minimal(base_size = 18) +
    theme(legend.position = "right",
          legend.box.just = "center",
          legend.title = element_text(face = "bold", size = 25, 
                                      hjust = 0.5, vjust = 2,
                                      margin = ggplot2::margin(b = 15)),
          legend.text = element_text(size = 28),
          axis.text = element_text(face = "bold"), 
          plot.background = element_blank(),
          plot.title = element_text(hjust = 0, face = "bold", size = 28),
          plot.margin = ggplot2::margin(r = 10, l = 2, t = 3, b = 3),
          axis.text.y = element_text(face = "bold", size = 32),
          axis.text.x = element_text(face = "bold", size = 32),
          axis.title.x = element_text(size = 28, margin = ggplot2::margin(t = 14)),
          panel.grid.major.y = element_blank(), 
          panel.border     = element_blank(),  
          panel.background = element_blank(),
          panel.grid.minor = element_blank())
}

