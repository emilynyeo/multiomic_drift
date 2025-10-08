# Functions needed for processing 
library(pacman)
pacman::p_load(knitr, data.table, dplyr, tidyr, tableone, kableExtra, readxl,
               readr, car, RColorBrewer, gridExtra, mlbench, earth, ggplot2, missForest,
               AppliedPredictiveModeling, caret, reshape2, corrplot, stringr,
               summarytools, grid, mice, plyr, mlmRev, cowplot, ape, e1071,
               jtools, broom, patchwork, phyloseq, microbiome, glmnet, ISLR,
               MicrobiomeStat, ANCOMBC, ape, vegan, zCompositions, janitor, MASS,
               RColorBrewer, DT, ggpubr, microbiomeutilities, compositions, VIM,
               naniar, tibble, table1)

## Functions

#########################################################################

# Function to process the names BL
process_names_bl <- function(names) {
  bl_names <- grep("\\.BL$", names, value = TRUE) # ending with "BL"
  extracted_numbers_BL <- sub(".*-(\\d+)\\.BL$", 
                              "\\1", bl_names) # Numbers after "-", before "."
  cleaned_numbers_bl <- sub("^0+", "", extracted_numbers_BL) # Remove leading zeros
  cleaned_numbers_bl} # Return the cleaned numbers

#########################################################################

# Function to process the names ALL
process_names_all <- function(names) {
  cleaned_numbers_all <- character(length(names)) #vector to store cleaned numbers
  all_names <- grep("\\..*$", names, value = TRUE) #names matching the pattern
  for (i in seq_along(names)) { # Extract numbers and remove leading zeros
    if (names[i] %in% all_names) {
      extracted_number <- sub(".*-(\\d+)\\..*$", "\\1", names[i])
      cleaned_number <- sub("^0+", "", extracted_number)
      cleaned_numbers_all[i] <- cleaned_number
    } else {
      cleaned_numbers_all[i] <- names[i] # If name don't match, keep the original
    }
  }
  return(cleaned_numbers_all)  # Return the cleaned numbers
}

#########################################################################

# Plot Density
my_plt_density <- function(df, cols, xlim_vals, plot_title) {
  path_long <- pivot_longer(df, cols = cols, names_to = "variable", values_to = "value")
  
  # Create a color palette based on the number of variables
  num_vars <- length(unique(path_long$variable))
  color_palette <- rep("skyblue", num_vars)  # Skyblue for all variables (or change as needed)
  color_palette_lines <- rep("black", num_vars)  # Black for all variables (or change as needed)
  
  ggplot(path_long, aes(x = value, fill = variable, color = variable)) + 
    geom_density(alpha = 0.5) + 
    theme_minimal() + 
    xlim(xlim_vals) + 
    labs(title = plot_title, x = "Value", y = "Density") + 
    scale_fill_manual(values = color_palette) + 
    scale_color_manual(values = color_palette_lines) +
    guides(fill = FALSE, color = FALSE)  # Hide the color and fill legends
}

#########################################################################

# Check normality and skewness 
check_normality_and_skewness <- function(df, cols) {
  # Shapiro-Wilk test for normality
  normality_summary <- table(ifelse(apply(df[, cols], 2, 
                                          function(x) shapiro.test(x)$p.value) > 0.05, 
                                    "Normal", "Non-Normal"))
  
  # Skewness calculation
  skew_summary <- data.frame(Variable = names(df)[cols],
                             Skewness = apply(df[, cols], 2, skewness, na.rm = TRUE))
  
  # Summary of skewness (skewed or not)
  skew_summary_count <- table(ifelse(abs(skew_summary$Skewness) > 1, "Skewed", "Not Skewed"))
  
  # Automatically print normality and skewness results
  print("Normality Summary:")
  print(normality_summary)
  
  print("Skewness Summary:")
  print(skew_summary_count)
  
  # Return the results as a list
  return(list(
    normality_summary = normality_summary,
    skew_summary_count = skew_summary_count
  ))
}
# Make unique
make_unique_names <- function(names) {
  name_counts <- table(names)  # Create a table of name counts
  occurrence <- integer(length(names)) # Vector to track occurrences
  
  # Loop through each name and assign a unique suffix if necessary
  for (i in seq_along(names)) {
    if (name_counts[names[i]] > 1) {
      occurrence[i] <- sum(names[1:i] == names[i])  # Count occurrences up to current index
      names[i] <- paste0(names[i], "_", occurrence[i])  # Append the occurrence number
    }
  }
  
  # Ensure unique names
  for (i in seq_along(names)) {
    while (sum(names[i] == names) > 1) {
      occurrence[i] <- occurrence[i] + 1
      names[i] <- sub("_\\d+$", "", names[i])  # Remove existing suffix
      names[i] <- paste0(names[i], "_", occurrence[i])  # Append new occurrence number
    }
  }
  return(names)
}

#########################################################################

# fancy_processing
fancy_process <- function(data, min_sample_threshold = 0.1, 
                          var_threshold = 0.1, 
                          correlation_threshold = 0.95, 
                          center = TRUE, scale = TRUE) {
  
  # Check if the data is a data.frame or matrix
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input data must be a data frame or matrix.")
  }
  
  # Initialize the removed_info list
  removed_info <- list()
  
  # Step 0: Determine the most suitable transformation based on data distribution
  cat("\nStep 0: Determining the most suitable transformation...\n")
  
  # Function to check normality with Shapiro-Wilk test
  normality_test <- function(x) {
    if (length(x) < 3) return(TRUE)  # Handle small data cases
    p_value <- shapiro.test(x)$p.value
    return(p_value > 0.05)  # p-value > 0.05 means data is normally distributed
  }
  
  # Isolate numeric columns only
  numeric_data <- data[sapply(data, is.numeric)]
  
  # Step 0.1: Test if the numeric data is normal (Shapiro-Wilk test) and skewness
  normal_distrib <- apply(numeric_data, 2, normality_test)  # Apply to all numeric columns
  skew_vals <- apply(numeric_data, 2, function(x) skewness(x, na.rm = TRUE))
  
  # Check if the majority of the data columns are normal
  if (all(normal_distrib)) {
    transformation <- "zscore"
    cat("Data appears to be normally distributed. Using Z-score normalization.\n")
  } else {
    # If data is not normally distributed, check for skewness
    skew_columns <- which(abs(skew_vals) > 0.5)  # Skew threshold (positive or negative skew)
    
    if (length(skew_columns) > 0) {
      if (any(skew_vals[skew_columns] > 0)) {
        # If positively skewed
        transformation <- "log10"
        cat("Data is positively skewed. Using Log10 transformation.\n")
      } else {
        # If negatively skewed or zero/negative values present
        transformation <- "arcsinh"
        cat("Data is negatively skewed or has zero/negative values. Using Arcsinh transformation.\n")
      }
    } else {
      # If data isn't skewed, use Box-Cox for stabilization
      transformation <- "boxcox"
      cat("Data is not skewed. Using Box-Cox transformation.\n")
    }
  }
  
  # Step 0.2: Detect if the data is relative abundance (compositional data)
  is_relative_abundance <- all(abs(rowSums(numeric_data) - 1) < 0.01)  # Adjust threshold as needed
  if (is_relative_abundance) {
    transformation <- "clr"
    cat("Data appears to be relative abundance data. Using CLR transformation.\n")
  }
  
  # Step 1: Normalize the numeric data using selected transformation
  cat("\nStep 1: Normalizing the data...\n")
  if (transformation == "log10") {
    numeric_data <- log10(numeric_data + 1)
    cat("Applied log10 transformation.\n")
  } else if (transformation == "clr") {
    gm <- apply(numeric_data, 1, function(x) exp(sum(log(x[x > 0])) / sum(x > 0)))  # Geometric mean/row
    numeric_data <- t(apply(numeric_data, 1, function(x) log(x / gm)))
    cat("Applied CLR transformation.\n")
  } else if (transformation == "zscore") {
    numeric_data <- scale(numeric_data)
    cat("Applied Z-score normalization.\n")
  } else if (transformation == "boxcox") {
    numeric_data <- apply(numeric_data, 2, function(x) ifelse(all(x > 0), boxcox(x ~ 1)$y, x))
    cat("Applied Box-Cox transformation.\n")
  } else if (transformation == "arcsinh") {
    numeric_data <- asinh(numeric_data)
    cat("Applied Arcsinh transformation.\n")
  }
  
  removed_info$step1 <- list(removed_rows = 0, removed_columns = 0)
  
  # Step 2: Centering and Scaling (standardize data)
  cat("\nStep 2: Centering and Scaling the data...\n")
  if (center) {
    numeric_data <- scale(numeric_data, center = TRUE, scale = FALSE)  # Centering only
    cat("Applied centering (subtracting mean).\n")
  }
  if (scale) {
    numeric_data <- scale(numeric_data, center = FALSE, scale = TRUE)  # Scaling only
    cat("Applied scaling (dividing by standard deviation).\n")
  }
  
  removed_info$step2 <- list(removed_rows = 0, removed_columns = 0)
  
  # Step 3: Remove variables not present in a minimum % threshold of samples
  cat("\nStep 3: Removing variables not present in a minimum percentage of samples...\n")
  min_samples <- floor(min_sample_threshold * nrow(numeric_data))
  cols_before_step3 <- ncol(numeric_data)
  numeric_data <- numeric_data[, colSums(!is.na(numeric_data)) >= min_samples]
  cols_removed_step3 <- cols_before_step3 - ncol(numeric_data)
  cat(cols_removed_step3, "columns removed based on the minimum sample threshold of", 
      min_sample_threshold * 100, "%.\n")
  removed_info$step3 <- list(removed_rows = 0, removed_columns = cols_removed_step3)
  
  # Step 4: Remove variables with variance below the variance threshold
  cat("\nStep 4: Removing variables with variance below", var_threshold, "...\n")
  vars_before_step4 <- ncol(numeric_data)
  numeric_data <- numeric_data[, apply(numeric_data, 2, var, na.rm = TRUE) >= var_threshold]
  vars_removed_step4 <- vars_before_step4 - ncol(numeric_data)
  cat(vars_removed_step4, "columns removed due to variance below the threshold.\n")
  removed_info$step4 <- list(removed_rows = 0, removed_columns = vars_removed_step4)
  
  # Step 5: Remove highly collinear features (Pearson's > 0.95)
  cat("\nStep 5: Removing highly collinear features (Pearson's > 0.95)...\n")
  non_zero_variance_data <- numeric_data[, apply(numeric_data, 2, var, na.rm = TRUE) > 0]
  if (ncol(non_zero_variance_data) > 1) {
    correlation_matrix <- cor(non_zero_variance_data, use = "pairwise.complete.obs")
    upper_tri <- correlation_matrix[upper.tri(correlation_matrix)]
    highly_correlated <- which(abs(upper_tri) > correlation_threshold, arr.ind = TRUE)
    if (length(highly_correlated) > 0) {
      high_corr_idx <- unique(c(highly_correlated[, 1], highly_correlated[, 2]))
      cols_before_step5 <- ncol(numeric_data)
      numeric_data <- numeric_data[, -high_corr_idx]
      cols_removed_step5 <- cols_before_step5 - ncol(numeric_data)
      cat(cols_removed_step5, "columns removed due to high correlation ( > 0.95).\n")
      removed_info$step5 <- list(removed_rows = 0, removed_columns = cols_removed_step5)
    } else {
      cat("No highly correlated features found (  > 0.95).\n")
      removed_info$step5 <- list(removed_rows = 0, removed_columns = 0)
    }
  } else {
    cat("Not enough columns left for correlation analysis.\n")
    removed_info$step5 <- list(removed_rows = 0, removed_columns = 0)
  }
  
  # Replace numeric data back into the original data frame, keeping non-numeric columns unchanged
  non_numeric_data <- data[sapply(data, Negate(is.numeric))]
  final_data <- cbind(non_numeric_data, numeric_data)
  
  # Return the processed data
  return(final_data)
}

#########################################################################

# Define a function to print taxa removal summary
print_taxa_summary <- function(before, after, step) {
  removed <- before - after
  print(paste0("# Taxa removed in step ", step, ": ", removed))
  print(paste0("# Taxa remaining after step ", step, ": ", after))
}

#########################################################################

## Check the satus function 
check_status <- function(df, id_track = "Sample ID") {
  if (!id_track %in% names(df)) {
    stop(paste("Column", id_track, "not found in the dataframe."))
  }
  n_samples <- df %>% pull(all_of(id_track)) %>% unique() %>% length()
  n_features <- ncol(df) - 1  # subtract 1 for the ID column
  cat("Tracking ID column:", id_track, "\n")
  cat("Unique IDs:", n_samples, "\n")
  cat("Number of Features:", n_features, "\n")
}