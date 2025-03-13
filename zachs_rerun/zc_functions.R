# Useful functions 

### Function to load RData file into a new environment and return the environment
load_rdata <- function(file_path) {
  if (file_ext(file_path) == "RData") {
    env <- new.env()
    load(file_path, envir = env)
    message("RData file loaded successfully.")
    return(env)
  } else {
    stop("The file is not an RData file.")
  }
}


### Make new columns 
make_new_columns <- function(data, column_name) {
  data <- data %>%
    mutate(
      subject_id = str_split(column_name, "\\.", simplify = TRUE)[, 1],
      TIMEPOINT = str_split(column_name, "\\.", simplify = TRUE)[, 2]
    )
  return(data)
}


### Filter data
filter_data <- function(data, column, value) {
  data <- data %>%
    filter(column == value)
  return(data)
}


### Merge data
merge_data <- function(data1, data2, join_type, columnname) {
  data <- join_type(data1, data2, by = columnname)
  
  # Remove the duplicated columns and rename as necessary
  data <- data %>%
    dplyr::select(-matches(paste0(columnname, "\\.y$"))) %>%
    rename_with(~ gsub("\\.x$", "", .), ends_with(".x"))
  
  return(data)
}


### Remove columns 
remove_columns <- function(data, columns_to_remove = NULL, pattern = NULL) {
  # Remove specified columns if provided
  if (!is.null(columns_to_remove)) {
    data <- data %>% dplyr::select(-all_of(columns_to_remove))
  }
  
  # Remove columns matching the pattern if provided
  if (!is.null(pattern)) {
    data <- data %>% dplyr::select(-matches(pattern))
  }
  
  return(data)
}


### Extract columns 
extract_columns <- function(data, columns_to_extract = NULL, pattern = NULL) {
  # Case 1: Both columns_to_extract and pattern are specified
  if (!is.null(columns_to_extract) && !is.null(pattern)) {
    data <- data %>% 
      select(all_of(intersect(names(data), 
                              columns_to_extract)), 
             matches(pattern))
    
    # Case 2: Only pattern is specified
  } else if (is.null(columns_to_extract) && !is.null(pattern)) {
    data <- data %>% select(matches(pattern))
    
    # Case 3: Only columns_to_extract is specified
  } else if (!is.null(columns_to_extract) && is.null(pattern)) {
    data <- data %>% select(all_of(columns_to_extract))
  }
  
  # Return data with the selected columns
  return(data)
}


### Function to automatically save all data frames/matrices in an environment as CSV files
save_env_to_csv <- function(env, output_dir) {
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Loop through the objects in the environment
  for (obj_name in ls(envir = env)) {
    obj <- get(obj_name, envir = env)
    
    # Check if the object is a data frame or matrix
    if (is.data.frame(obj) || is.matrix(obj)) {
      file_path <- file.path(output_dir, paste0(obj_name, ".csv"))
      write.csv(obj, file = file_path)
      message(paste("Saved:", file_path))
    } else {
      message(paste("Skipping:", obj_name, "as it is not a data frame or matrix."))
    }
  }
}

### Function to split column into id and time col
make_new_columns <- function(data, column_name) {
  data <- data %>%
    mutate(
      subject_id = str_split(column_name, "\\.", simplify = TRUE)[, 1],
      TIMEPOINT = str_split(column_name, "\\.", simplify = TRUE)[, 2]
    )
  return(data)
}

### Function to deal with Taxa names 
rename_columns_species_to_domain <- function(dataframe) {
  # Define the order of taxonomic levels
  order <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__")
  
  # Get all column names
  columns <- colnames(dataframe)
  
  # Filter columns starting with 'd__'
  species_columns <- columns[grepl("^d__", columns)]
  
  # Loop through each species column
  for (column in species_columns) {
    # Split by '_{any single letter}__'
    split_column <- unlist(strsplit(column, "_[a-z]__"))
    split_column[1] <- sub("d__", "", split_column[1]) # Remove 'd__' prefix
    
    # Remove empty elements
    split_column <- split_column[split_column != ""]
    
    # Get the lowest level of taxonomy
    order_index <- length(split_column) - 1
    new_column_name <- paste0(order[order_index + 1], split_column[length(split_column)])
    
    # Rename the column in the data frame
    colnames(dataframe)[colnames(dataframe) == column] <- new_column_name
  }
  
  return(dataframe)
}


### Function to standardize and impute the data
process_data <- function(data, columns_to_remove, columns_to_standardize, impute_method = "medianImpute") {
  data_cleaned <- remove_columns(data, columns_to_remove)
  data_standardized <- preprocess_data(data_cleaned, columns_to_standardize, impute_method)
  return(data_standardized)
}

### Function to preprocess the data 
preprocess_data <- function(data, columns_to_standardize, imputation_method) {
  data_imputed <- predict(preProcess(data, method = c(imputation_method)), data)
  data_standardized <- predict(preProcess(data_imputed[, columns_to_standardize], method = c("center", "scale")), data_imputed)
  return(data_standardized)
}

### Function to train all models 
train_all_models <- function(data, target, train_control) {
  lasso_model <- train(
    x = as.data.frame(data[, -which(names(data) %in% c(target, "subject_id"))]),
    y = as.numeric(data[[target]]),
    method = "glmnet",
    trControl = train_control,
    tuneGrid = expand.grid(alpha = 1, lambda = seq(0.1, 1, 0.1))
  )
  
  ridge_model <- train(
    x = as.data.frame(data[, -which(names(data) %in% c(target, "subject_id"))]),
    y = as.numeric(data[[target]]),
    method = "glmnet",
    trControl = train_control,
    tuneGrid = expand.grid(alpha = 0, lambda = seq(0.1, 1, 0.1))
  )
  
  elastic_net_model <- train(
    x = as.data.frame(data[, -which(names(data) %in% c(target, "subject_id"))]),
    y = as.numeric(data[[target]]),
    method = "glmnet",
    trControl = train_control,
    tuneGrid = expand.grid(alpha = seq(0.1, 1, 0.1), lambda = seq(0.1, 1, 0.1))
  )
  
  rf_model <- train(
    x = as.data.frame(data[, -which(names(data) %in% c(target, "subject_id"))]),
    y = as.numeric(data[[target]]),
    method = "rf",
    trControl = train_control,
    tuneGrid = expand.grid(mtry = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)),
    ntree = 10000
  )
  
  xgb_model <- train(
    x = as.data.frame(data[, -which(names(data) %in% c(target, "subject_id"))]),
    y = as.numeric(data[[target]]),
    method = "xgbTree",
    trControl = train_control
  )
  
  caret.list <- caretList(
    x = as.data.frame(data[, -which(names(data) %in% c(target, "subject_id"))]),
    y = as.numeric(data[[target]]),
    trControl = train_control,
    methodList = c("glmnet", "rf", "xgbTree")
  )
  
  ens <- caretEnsemble(caret.list)
  
  return(list(
    lasso_model = lasso_model,
    ridge_model = ridge_model,
    elastic_net_model = elastic_net_model,
    rf_model = rf_model,
    xgb_model = xgb_model,
    ens = ens
  ))
}

### Function to extract importances 
extract_importance_df <- function(model, label) {
  importance <- varImp(model)$importance %>%
    as.data.frame() %>%
    rownames_to_column("Variable")
  colnames(importance)[2] <- label
  return(importance)
}

### Function to combine importances 
combine_importances <- function(model_list, labels) {
  importance_dfs <- map2(model_list, labels, extract_importance_df)
  reduce(importance_dfs, full_join, by = "Variable")
}

### Function to extract the best beta values 
extract_best_betas <- function(model_list, labels) {
  beta_dfs <- map2(model_list, labels, function(model, label) {
    best_lambda <- model$bestTune$lambda
    betas <- as.data.frame(as.matrix(coef(model$finalModel, s = best_lambda))) %>%
      rownames_to_column("Variable")
    colnames(betas)[2] <- label
    return(betas)
  })
  beta_combined <- reduce(beta_dfs, full_join, by = "Variable") %>%
    mutate(across(everything(), ~ replace_na(., 0)))
  return(beta_combined)
}

### Function to replace NA values 
replace_na <- function(x, value) {
  ifelse(is.na(x), value, x)
}

### Helper function to calculate model performance metrics for either training or testing data
calculate_metrics <- function(model, data, target_var, model_name, data_type) {
  predictions <- predict(model, data)
  actuals <- data[[target_var]]
  
  r2 <- caret::R2(predictions, actuals)
  mae <- caret::MAE(predictions, actuals)
  rmse <- caret::RMSE(predictions, actuals)
  
  return(data.frame(Model = model_name, DataType = data_type, R2 = r2, MAE = mae, RMSE = rmse))
}


### Function to run all the models 
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
  if (!dir.exists("drift_fs/csv/results/")) {
    dir.create("drift_fs/csv/results/", recursive = TRUE)
  }
  
  write.csv(feature_importance, paste0("drift_fs/csv/results/", 
                                       result_prefix, "_feature_importance.csv"), 
            row.names = FALSE)
  
  # Extract best betas
  beta_coefficients <- extract_best_betas(
    list(results$lasso_model, results$ridge_model, results$elastic_net_model),
    c("Lasso_Beta", "Ridge_Beta", "Enet_Beta")
  )
  write.csv(beta_coefficients, paste0("drift_fs/csv/results/", result_prefix, "_beta.csv"), row.names = FALSE)
  
  # Initialize an empty DataFrame to store performance metrics
  metrics_df <- data.frame(Model = character(), DataType = character(), R2 = numeric(), MAE = numeric(), RMSE = numeric(), stringsAsFactors = FALSE)
  
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
  write.csv(metrics_df, paste0("drift_fs/csv/results/", result_prefix, "_metrics.csv"), row.names = FALSE)
  
  # ensure the model directory exists
  if (!dir.exists("drift_fs/models/")) {
    dir.create("drift_fs/models/", recursive = TRUE)
  }
  # Save model results
  saveRDS(results, paste0("drift_fs/models/", result_prefix, "_results.rds"))
  
  return(results)
}

####Same fucntion as above but multiple models 
train_and_save_multiple_models <- function(datasets, target_vars, 
                                           train_control, result_prefixes, 
                                           test_size = 0.3) {
  
  # Check if lengths of datasets, target_vars, and result_prefixes are equal
  if (length(datasets) != length(target_vars) || length(datasets) != length(result_prefixes)) {
    stop("The number of datasets, target variables, and result prefixes must be the same.")
  }
  
  # Loop over each dataset, target variable, and result prefix
  for (i in seq_along(datasets)) {
    data <- datasets[[i]]
    target_var <- target_vars[i]
    result_prefix <- result_prefixes[i]
    
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
    
    # Ensure output directory exists
    if (!dir.exists("drift_fs/csv/results/feb20/")) {
      dir.create("drift_fs/csv/results/feb20/", recursive = TRUE)
    }
    
    write.csv(feature_importance, 
              paste0("drift_fs/csv/results/feb20/", 
                     result_prefix, "_feature_importance.csv"), 
              row.names = FALSE)
    
    # Extract best betas
    beta_coefficients <- extract_best_betas(
      list(results$lasso_model, results$ridge_model, results$elastic_net_model),
      c("Lasso_Beta", "Ridge_Beta", "Enet_Beta")
    )
    write.csv(beta_coefficients, paste0("drift_fs/csv/results/feb20/", 
                                        result_prefix, "_beta.csv"), row.names = FALSE)
    
    # Initialize an empty DataFrame to store performance metrics
    metrics_df <- data.frame(Model = character(), DataType = character(), 
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
    write.csv(metrics_df, paste0("drift_fs/csv/results/feb20/", 
                                 result_prefix, "_metrics.csv"), row.names = FALSE)
    
    # Ensure the model directory exists
    if (!dir.exists("drift_fs/models/feb20/")) {
      dir.create("drift_fs/models/feb20/", recursive = TRUE)
    }
    
    # Save model results
    saveRDS(results, paste0("drift_fs/models/feb20/", result_prefix, "_results.rds"))
  }
  
  return("Model training completed for all datasets.")
}

####

# Function to get top N important features for a given model
get_top_n_features <- function(feature_importance, model_importance_column, n = 20) {
  feature_importance %>%
    select(Variable, all_of(model_importance_column)) %>%
    arrange(desc(get(model_importance_column))) %>%
    head(n)
}

# Function to get top N important features for all models in the feature importance dataframe
get_top_n_features_all_models <- function(feature_importance, n = 20) {
  # Identify columns ending with '_Importance'
  importance_columns <- names(feature_importance)[grepl("_Importance$", names(feature_importance))]
  
  # Initialize a list to store results for each model
  top_features_list <- list()
  
  # Loop through each importance column and get top N features
  for (column in importance_columns) {
    model_name <- gsub("_Importance", "", column) # Extract model name
    top_features_list[[model_name]] <- get_top_n_features(feature_importance, column, n)
  }
  
  return(top_features_list)
}

##### 
train_multiple_models <- function(datasets, target_vars, 
                                  train_control, result_prefixes, 
                                  test_size = 0.3) {
  
  # Check if lengths of datasets, target_vars, and result_prefixes are equal
  if (length(datasets) != length(target_vars) || length(datasets) != length(result_prefixes)) {
    stop("The number of datasets, target variables, and result prefixes must be the same.")
  }
  
  # Initialize an empty list to store the results
  all_results <- list()
  
  # Loop over each dataset, target variable, and result prefix
  for (i in seq_along(datasets)) {
    data <- datasets[[i]]
    target_var <- target_vars[i]
    result_prefix <- result_prefixes[i]
    
    # Split data into training and testing sets
    set.seed(123) # Ensure reproducibility
    train_indices <- sample(seq_len(nrow(data)), size = (1 - test_size) * nrow(data))
    train_data <- data[train_indices, ]
    test_data <- data[-train_indices, ]
    
    # Train models
    results <- train_all_models(train_data, target_var, train_control)
    
    # Save results with a unique name based on result_prefix
    all_results[[result_prefix]] <- results
  }
  
  # Return the list of all results
  return(all_results)
}


#####

# Function to plot Venn diagram for top features
plot_venn_diagram <- function(feature_sets, colors = NULL, 
                              figurename = "venn_diagram.png", 
                              output_dir = "drift_fs/figures/") {
  if (is.null(colors)) {
    colors <- viridis(length(feature_sets))
  }
  
  venn.plot <- venn.diagram(
    x = feature_sets,
    category.names = names(feature_sets),
    filename = NULL, # Plot directly to the object
    output = TRUE,
    col = "transparent",
    fill = colors,
    alpha = 0.3,
    cex = 1.5,
    cat.cex = 1.2,
    cat.pos = 0,
    margin = 0.1
  )
  
  ggsave(file.path(output_dir, figurename), plot = venn.plot, dpi = 600, bg = "white")
}

# Function to plot feature importance or beta values
plot_importance_or_beta <- function(data, value_column, plot_title, y_label, figurename, palette = "viridis") {
  long_format <- data %>%
    pivot_longer(
      cols = ends_with(value_column),
      names_to = "Model",
      values_to = value_column
    )
  
  ggplot(long_format, aes(y = reorder(Variable, !!sym(value_column)), x = !!sym(value_column), fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
    scale_fill_viridis_d(option = palette) +
    labs(
      title = plot_title,
      x = paste(value_column, "Score"),
      y = y_label
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 5), # Adjust text size
      axis.title.y = element_text(vjust = 1)
    ) +
    scale_y_discrete(expand = expansion(mult = c(0.1, 0.2))) # Add padding between features
  
  ggsave(file.path("drift_fs/figures/", figurename), dpi = 600, bg = "white")
}

# Function to plot model performance metrics
plot_performance_metrics <- function(metrics, dataset_name) {
  metrics_long <- metrics %>%
    pivot_longer(
      cols = c("R2", "MAE", "RMSE"),
      names_to = "Metric",
      values_to = "Value"
    )
  
  # Reorder factors for better readability
  metrics_long$DataType <- factor(metrics_long$DataType, levels = c("Train", "Test"))
  metrics_long$Metric <- factor(metrics_long$Metric, levels = c("R2", "MAE", "RMSE"))
  
  ggplot(metrics_long, aes(x = Model, y = Value, fill = DataType)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    labs(
      title = paste("Model Performance Metrics (R², MAE, RMSE) -", dataset_name),
      x = "Model",
      y = "Metric Value",
      fill = "Dataset"
    ) +
    scale_fill_viridis(discrete = TRUE) +
    facet_wrap(~Metric, scales = "free_y", nrow = 3) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 12),
      legend.position = "top",
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
  
  ggsave(file.path("drift_fs/figures/", paste0(dataset_name, "_performance_metrics.png")), dpi = 600, bg = "white")
}


### Function to process data and generate plots for a dataset
process_and_plot_data <- function(data_list, dataset_name, n = 20) {
  beta <- data_list$beta
  feature_importance <- data_list$feature_importance
  metrics <- data_list$metrics
  
  # Get top N features for all models
  top_20_features <- get_top_n_features_all_models(feature_importance, n)
  
  # Create feature sets for Venn diagram
  feature_sets <- lapply(top_20_features, function(df) df$Variable)
  
  # Plot Venn diagram
  plot_venn_diagram(feature_sets, figurename = paste0(dataset_name, "_venn_diagram.png"))
  
  # Extract top features and sort by total importance
  all_top_features <- unique(unlist(feature_sets))
  filtered_feature_importance <- feature_importance %>%
    filter(Variable %in% all_top_features) %>%
    mutate(Total_Importance = rowSums(select(., ends_with("_Importance")), na.rm = TRUE)) %>%
    arrange(desc(Total_Importance))
  
  # Remove total importance column
  filtered_feature_importance <- select(filtered_feature_importance, -Total_Importance)
  
  # Plot Feature Importances
  plot_importance_or_beta(filtered_feature_importance, "Importance", "Top 20 Feature Importances", "Features", paste0(dataset_name, "_feature_importance.png"))
  
  # Plot model performance metrics
  plot_performance_metrics(metrics, dataset_name)
}


### Function to add points and recalculate regression
generate_plot <- function(x, y, modified_x = NULL, modified_y = NULL) {
  # Create the data frame
  plot_data <- data.frame(x = x, y = y)
  
  # Add modified points if provided
  if (!is.null(modified_x) && !is.null(modified_y)) {
    modified_data <- data.frame(x = modified_x, y = modified_y)
    plot_data <- rbind(plot_data, modified_data)
  }
  
  # Fit the linear regression model
  model <- lm(y ~ x, data = plot_data)
  r_squared <- summary(model)$r.squared
  
  # Generate the base plot
  p <- ggplot(plot_data, aes(x = x, y = y)) +
    geom_point(aes(color = ifelse(x %in% modified_x, "Modified", "Original")), size = 5) +
    scale_color_manual(values = c("Original" = "#1C7C54", "Modified" = "red")) +
    theme_minimal() +
    xlim(0, max(plot_data$x)) +
    ylim(min(plot_data$y) - 2, max(plot_data$y) + 2) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.ticks.length = unit(0.3, "cm"), # Adjust tick size
      axis.title.x = element_text(size = 22, face = "bold"), # Larger x-axis title
      axis.title.y = element_text(size = 22, face = "bold") # Larger y-axis title
    )
  
  # Add regression line and R² annotation
  p1 <- p +
    geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray", alpha = 0.5) +
    annotate("text",
             x = max(plot_data$x) - 2, y = mean(plot_data$y),
             label = paste("R² =", round(r_squared, 2)),
             color = "black", size = 8
    )
  
  # Return both plots in a list
  return(list(base_plot = p, regression_plot = p1))
}

### Define a function to create plots
create_plots <- function(data_list, max_r2, titles) {
  plots <- lapply(seq_along(data_list), function(i) {
    ggplot(data_list[[i]], aes(x = R2, y = Model, fill = Model)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(
        title = titles[i],
        x = "R²",
        y = "Model",
        fill = "Model"
      ) +
      coord_cartesian(xlim = c(0, max_r2)) +
      #coord_cartesian(xlim(0, 0.5)) +
      theme_minimal() +
      scale_fill_viridis_d() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15, angle = 45),
        legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 20)
      )
  })
  # Combine the plots vertically
  Reduce(`/`, plots)
}

### Updated function to create and save plots
create_feature_plot <- function(features, title, save_path) {
  # Prepare data for plotting, excluding the `Model` column from pivoting
  features_long <- features %>%
    pivot_longer(
      cols = c(-Variable), # Exclude `Variable` and `Model`
      names_to = "Model",
      values_to = "Importance"
    )
  
  # Create the plot
  feature_plot <- ggplot(features_long, aes(x = Importance, y = reorder(Variable, Importance), fill = Model)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      title = title,
      x = "Importance",
      y = "Feature",
      fill = "Model"
    ) +
    scale_fill_viridis_d() +
    xlim(0, 100) +  # Set x-axis range from 0 to 100
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 13),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      title = element_text(size = 15),
      legend.position = "top"
    )
  
  # Print and save the plot
  print(feature_plot)
  pdf(file = save_path, width = 10, height = 10)
  print(feature_plot)
  dev.off()
}


### Function to extract top features
extract_top_features <- function(dataset, top_n = 15) {
  # Extract columns ending with '_Importance'
  importance_columns <- names(dataset$feature_importance) %>%
    grep("_Importance$", ., value = TRUE)
  
  # Check if importance columns exist
  if (length(importance_columns) == 0) {
    stop("No columns ending with '_Importance' found in feature_importance.")
  }
  
  # get a column with the total importance
  dataset$feature_importance <- dataset$feature_importance %>%
    mutate(Total_Importance = rowSums(select(., ends_with("_Importance")), na.rm = TRUE))
  
  # sort by total importance from high to low
  dataset$feature_importance <- dataset$feature_importance %>%
    arrange(desc(Total_Importance))
  
  # get the top n features
  top_features <- dataset$feature_importance %>%
    select(Variable, Total_Importance) %>%
    head(10)
  
  # remove the total importance column
  dataset$feature_importance <- select(dataset$feature_importance, -Total_Importance) %>% head(10)
  
  return(dataset$feature_importance)
}


### Function to get top 20 features for each model
get_top_features <- function(dataset, model_name, importance_col) {
  dataset$feature_importance %>%
    select(Variable, all_of(importance_col)) %>%
    arrange(desc(get(importance_col))) %>%
    head(20)
}

### Extract top 20 features for each model and dataset
get_features_for_dataset <- function(dataset_name) {
  dataset <- datasets[[dataset_name]]
  map2(
    model_names, model_names_map_to,
    ~ get_top_features(dataset, .x, .y)
  ) %>%
    set_names(model_names)
}

### Function to get top models based on R²
get_top_models <- function(dataset, top_n = 3) {
  dataset$metrics %>%
    filter(DataType == "Test") %>%
    arrange(desc(R2)) %>%
    slice_head(n = top_n) %>%
    pull(Model)
}


# Function to extract metrics and calculate max R²
extract_metrics <- function(dataset_name, datasets) {
  metrics <- datasets[[dataset_name]]$metrics %>%
    filter(DataType == "Test") %>%
    select(Model, R2)
  max_r2 <- max(metrics$R2, na.rm = TRUE)
  list(metrics = metrics, max_r2 = max_r2)
}


#### Get LASSO features
get_lasso_features <- function(lasso_model) {
  # Ensure the model is of class 'train' and that 'finalModel' is available
  if (!inherits(lasso_model, "train")) {
    stop("The provided object is not a valid 'train' model object.")
  }
  
  # Extract the best lambda value based on cross-validation
  best_lambda <- lasso_model$bestTune$lambda
  
  # Extract coefficients for the final model at the best lambda
  best_coefs <- coef(lasso_model$finalModel, s = best_lambda)
  
  # Get selected features (non-zero coefficients)
  selected_features <- rownames(best_coefs)[which(best_coefs != 0)]
  
  # Remove intercept from selected features
  selected_features <- selected_features[selected_features != "(Intercept)"]
  
  # Return the outcomes as a list
  result <- list(
    best_lambda = best_lambda,
    selected_features = selected_features,
    coefficients = best_coefs
  )
  
  return(result)
}

#### Extract best model from LASSO ###

get_lasso_features <- function(lasso_model, data, outcome_col) {
  # Ensure the model is of class 'train' and that 'finalModel' is available
  if (!inherits(lasso_model, "train")) {
    stop("The provided object is not a valid 'train' model object.")
  }
  
  # Extract the best lambda value based on cross-validation
  best_lambda <- lasso_model$bestTune$lambda
  
  # Extract coefficients for the final model at the best lambda
  best_coefs <- coef(lasso_model$finalModel, s = best_lambda)
  
  # Get selected features (non-zero coefficients)
  selected_features <- rownames(best_coefs)[which(best_coefs != 0)]
  
  # Remove intercept from selected features
  selected_features <- selected_features[selected_features != "(Intercept)"]
  
  # Extract the suffix from the model name (after the last underscore)
  model_name <- deparse(substitute(lasso_model))
  model_suffix <- sub(".*_(.*)", "\\1", model_name)
  
  # Create a new data frame with the selected features
  train_selected <- data[, c(outcome_col, selected_features)]
  
  # Create the formula for the final model using the selected features
  final_formula <- as.formula(paste(outcome_col, "~", paste(selected_features, collapse = " + ")))
  
  # Fit the final model (use lm, or modify depending on your desired model type)
  final_model <- lm(final_formula, data = train_selected)
  
  # Store the final model with a dynamic name
  assign(paste0("final_model_", model_suffix), final_model, envir = .GlobalEnv)
  
  # Return the outcomes as a list, including the model and selected features
  result <- list(
    best_lambda = best_lambda,
    selected_features = selected_features,
    coefficients = best_coefs,
    final_model = final_model
  )
  
  return(result)
}

### Get predictions 
best_lasso_predictions <- function(lasso_model, data, outcome_col, test_data) { 
  # Ensure the model is of class 'train' and that 'finalModel' is available
  if (!inherits(lasso_model, "train")) {
    stop("The provided object is not a valid 'train' model object.")
  }
  
  # Extract the best lambda value based on cross-validation
  best_lambda <- lasso_model$bestTune$lambda
  
  # Extract coefficients for the final model at the best lambda
  best_coefs <- coef(lasso_model$finalModel, s = best_lambda)
  
  # Get selected features (non-zero coefficients)
  selected_features <- rownames(best_coefs)[which(best_coefs != 0)]
  
  # Remove intercept from selected features
  selected_features <- selected_features[selected_features != "(Intercept)"]
  
  # Extract the suffix from the model name (after the last underscore)
  model_name <- deparse(substitute(lasso_model))
  model_suffix <- sub(".*_(.*)", "\\1", model_name)
  
  # Create a new data frame with the selected features
  train_selected <- data[, c(outcome_col, selected_features)]
  
  # Create the formula for the final model using the selected features
  final_formula <- as.formula(paste(outcome_col, "~", paste(selected_features, collapse = " + ")))
  
  # Fit the final model (use lm, or modify depending on your desired model type)
  final_model <- lm(final_formula, data = train_selected)
  
  # Store the final model with a dynamic name
  assign(paste0("final_model_", model_suffix), final_model, envir = .GlobalEnv)
  
  # Extract summary information of the model
  mymodsum <- summary(final_model)
  mod_coef_df <- as.data.frame(mymodsum[["coefficients"]])  # Extract coefficients from model summary
  
  # prediction on reserved validation samples
  test <- test_data[complete.cases(test_data), ]  # remove rows with missing values
  
  # Grab only the columns that we have coefficients for (starting at index 2 removes the intercept)
  test_vars <- test %>% dplyr::select(rownames(mod_coef_df)[2:nrow(mod_coef_df)])
  test_vars$Intercept <- 1  # add a column for the intercept
  test_vars <- test_vars %>% dplyr::select(Intercept, everything())
  
  # Predict the risk scores without covariates
  pred_risk_scores <- as.matrix(test_vars) %*% as.matrix(mod_coef_df$Estimate)
  
  # Combine the predicted and actual values into a data frame
  pred_df <- as.data.frame(cbind(
    test$subject_id,
    test$time,
    test$bmi_bL_6m,
    scale(pred_risk_scores)  # Scale the predictions if needed
  ))
  
  # Rename the columns for clarity
  colnames(pred_df) <- c("subject_id", "time", "actual", "predicted")
  
  # Return the prediction data frame as the result
  return(pred_df)
}

####
best_lasso_predictions <- function(lasso_model, data, outcome_col, test_data, s = NULL) {  
  # Ensure the model is of class 'train' and that 'finalModel' is available
  if (!inherits(lasso_model, "train")) {
    stop("The provided object is not a valid 'train' model object.")
  }
  
  # Extract the best lambda value based on cross-validation if s is NULL
  best_lambda <- if (is.null(s)) {
    lasso_model$bestTune$lambda
  } else {
    s  # Use the provided s if not NULL
  }
  
  # Extract coefficients for the final model at the selected lambda (s)
  best_coefs <- coef(lasso_model$finalModel, s = best_lambda)
  
  # Get selected features (non-zero coefficients)
  selected_features <- rownames(best_coefs)[which(best_coefs != 0)]
  
  # Remove intercept from selected features
  selected_features <- selected_features[selected_features != "(Intercept)"]
  
  # Extract the suffix from the model name (after the last underscore)
  model_name <- deparse(substitute(lasso_model))
  model_suffix <- sub(".*_(.*)", "\\1", model_name)
  
  # Create a new data frame with the selected features
  train_selected <- data[, c(outcome_col, selected_features)]
  
  # Create the formula for the final model using the selected features
  final_formula <- as.formula(paste(outcome_col, "~", paste(selected_features, collapse = " + ")))
  
  # Fit the final model (use lm, or modify depending on your desired model type)
  final_model <- lm(final_formula, data = train_selected)
  
  # Store the final model with a dynamic name
  assign(paste0("final_model_", model_suffix), final_model, envir = .GlobalEnv)
  
  # Extract summary information of the model
  mymodsum <- summary(final_model)
  mod_coef_df <- as.data.frame(mymodsum[["coefficients"]])  # Extract coefficients from model summary
  
  # Prediction on reserved validation samples
  test <- test_data[complete.cases(test_data), ]  # Remove rows with missing values
  
  # Grab only the columns that we have coefficients for (starting at index 2 removes the intercept)
  test_vars <- test %>% dplyr::select(rownames(mod_coef_df)[2:nrow(mod_coef_df)])
  test_vars$Intercept <- 1  # Add a column for the intercept
  test_vars <- test_vars %>% dplyr::select(Intercept, everything())
  
  # Predict the risk scores without covariates
  pred_risk_scores <- as.matrix(test_vars) %*% as.matrix(mod_coef_df$Estimate)
  
  # Combine the predicted and actual values into a data frame
  pred_df <- as.data.frame(cbind(
    test$subject_id,
    test$time,
    test$bmi_bL_6m,
    scale(pred_risk_scores)  # Scale the predictions if needed
  ))
  
  # Rename the columns for clarity
  colnames(pred_df) <- c("subject_id", "time", "actual", "predicted")
  
  # Return the prediction data frame as the result
  return(pred_df)
}

### Get RF predictions
best_rf_predictions <- function(rf_model, data, outcome_col, test_data) {  
  # Ensure the model is of class 'train' and that 'finalModel' is available
  if (!inherits(rf_model, "train")) {
    stop("The provided object is not a valid 'train' model object.")
  }
  
  # Extract the best mtry value from the model
  best_mtry <- rf_model$bestTune$mtry
  
  # Extract the Random Forest model (final model)
  rf_final_model <- rf_model$finalModel
  
  # Ensure that the test data is formatted similarly to the training data (excluding outcome column)
  test_data_prepared <- test_data[, -which(names(test_data) %in% c(outcome_col, "subject_id"))]
  
  # Predict using the Random Forest model
  predictions <- predict(rf_final_model, newdata = test_data_prepared)
  
  # Combine the predictions with actual values from the test set
  pred_df <- data.frame(
    subject_id = test_data$subject_id, 
    time = test_data$time, 
    actual = test_data[[outcome_col]], 
    predicted = predictions
  )
  
  # Return the prediction data frame
  return(pred_df)
}
