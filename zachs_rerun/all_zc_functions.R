# all_zc_functions.R

# Function to load RData file into a new environment and return the environment
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

# Function to automatically save all data frames/matrices in an environment as CSV files
save_env_to_csv <- function(env, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (obj_name in ls(envir = env)) {
    obj <- get(obj_name, envir = env)
    
    if (is.data.frame(obj) || is.matrix(obj)) {
      file_path <- file.path(output_dir, paste0(obj_name, ".csv"))
      write.csv(obj, file = file_path)
      message(paste("Saved:", file_path))
    } else {
      message(paste("Skipping:", obj_name, "as it is not a data frame or matrix."))
    }
  }
}

make_new_columns <- function(data, column_name) {
  data <- data %>%
    mutate(
      subject_id = str_split(column_name, "\\.", simplify = TRUE)[, 1],
      TIMEPOINT = str_split(column_name, "\\.", simplify = TRUE)[, 2]
    )
  return(data)
}

filter_data <- function(data, column, value) {
  data <- data %>%
    filter(column == value)
  return(data)
}

merge_data <- function(data1, data2, join_type, columnname) {
  data <- join_type(data1, data2, by = columnname)
  
  data <- data %>%
    select(-matches(paste0(columnname, "\\.y$"))) %>%
    rename_with(~ gsub("\\.x$", "", .), ends_with(".x"))
  
  return(data)
}

remove_columns <- function(data, columns_to_remove = NULL, pattern = NULL) {
  if (!is.null(columns_to_remove)) {
    data <- data %>% select(-all_of(columns_to_remove))
  }
  
  if (!is.null(pattern)) {
    data <- data %>% select(-matches(pattern))
  }
  
  return(data)
}

extract_columns <- function(data, columns_to_extract = NULL, pattern = NULL) {
  if (!is.null(columns_to_extract) && !is.null(pattern)) {
    data <- data %>% 
      select(all_of(intersect(names(data), 
                              columns_to_extract)), 
             matches(pattern))
  } else if (is.null(columns_to_extract) && !is.null(pattern)) {
    data <- data %>% select(matches(pattern))
  } else if (!is.null(columns_to_extract) && is.null(pattern)) {
    data <- data %>% select(all_of(columns_to_extract))
  }
  
  return(data)
}

rename_columns_species_to_domain <- function(dataframe) {
  order <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__")
  columns <- colnames(dataframe)
  species_columns <- columns[grepl("^d__", columns)]
  
  for (column in species_columns) {
    split_column <- unlist(strsplit(column, "_[a-z]__"))
    split_column[1] <- sub("d__", "", split_column[1])
    split_column <- split_column[split_column != ""]
    order_index <- length(split_column) - 1
    new_column_name <- paste0(order[order_index + 1], split_column[length(split_column)])
    colnames(dataframe)[colnames(dataframe) == column] <- new_column_name
  }
  
  return(dataframe)
}

process_data <- function(data, columns_to_remove, columns_to_standardize, impute_method = "medianImpute") {
  data_cleaned <- remove_columns(data, columns_to_remove)
  data_standardized <- preprocess_data(data_cleaned, columns_to_standardize, impute_method)
  return(data_standardized)
}

preprocess_data <- function(data, columns_to_standardize, imputation_method) {
  data_imputed <- predict(preProcess(data, method = c(imputation_method)), data)
  data_standardized <- predict(preProcess(data_imputed[, columns_to_standardize], method = c("center", "scale")), data_imputed)
  return(data_standardized)
}

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

# Add any other functions