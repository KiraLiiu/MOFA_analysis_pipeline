# ====
# Functions for MOFA Analysis Pipeline (Standalone)
# ====

# Function to reorder and fill missing columns with NA
reorder_and_fill_na <- function(df, desired_order) {
  missing_cols <- setdiff(desired_order, colnames(df))
  na_df <- matrix(NA, nrow = nrow(df), ncol = length(missing_cols)) %>% as.data.frame()
  colnames(na_df) <- missing_cols
  df <- cbind(df, na_df)
  df <- df[, desired_order, drop = FALSE]
  return(df)
}

# Function for selecting high variance features
select_high_variance_features <- function(data_matrix, variance_threshold = 0.8) {
  feature_variance <- apply(data_matrix, 1, var, na.rm = TRUE)
  sorted_var <- sort(feature_variance, decreasing = TRUE)
  cumulative_var <- cumsum(sorted_var) / sum(sorted_var)
  num_features <- which(cumulative_var >= variance_threshold)[1]
  top_features <- names(sorted_var)[1:num_features]
  selected_features <- data_matrix[top_features, , drop = FALSE]
  return(as.data.frame(selected_features))
}

# MOFA analysis function with per-dataset feature selection
run_mofa_analysis <- function(data_list,
                              view_names = NULL,
                              num_factors = 10,
                              scale_views = FALSE,
                              convergence_mode = "medium",
                              seed = 42,
                              max_iter = 1000,
                              save_model = TRUE,
                              outfile = "MOFA_model.hdf5",
                              rds_file = "MOFA_model.rds",
                              drop_factor_threshold = -1,
                              groups = NULL,
                              use_gpu = FALSE,
                              feature_selection = NULL) {
  start_time <- Sys.time()
  if (!is.list(data_list)) stop("data_list must be a list of data matrices")
  if (!is.null(feature_selection)) {
    cat("Performing feature selection...\n")
    for (i in seq_along(data_list)) {
      view_name = names(data_list)[i]
      if (!is.null(feature_selection[[view_name]])) {
        threshold <- feature_selection[[view_name]]
      } else if (!is.null(feature_selection$default)) {
        threshold <- feature_selection$default
      } else {
        threshold <- 0.8
      }
      data_list[[i]] <- select_high_variance_features(data_list[[i]], variance_threshold = threshold)
      cat(sprintf("View %s: %d features selected (variance threshold = %.2f)\n", view_name, nrow(data_list[[i]]), threshold))
    }
  }
  for (i in seq_along(data_list)) {
    data_list[[i]] <- as.matrix(data_list[[i]])
  }
  if (!is.null(view_names) && length(view_names) == length(data_list)) {
    names(data_list) <- view_names
  } else if (is.null(names(data_list))) {
    names(data_list) <- paste0("view", seq_along(data_list))
  }
  cat("Data dimensions after feature selection:\n")
  for (i in seq_along(data_list)) {
    cat(sprintf("%s: %d features x %d samples\n", names(data_list)[i], nrow(data_list[[i]]), ncol(data_list[[i]])))
  }
  cat("Creating MOFA object...\n")
  MOFAobject <- MOFA2::create_mofa(data_list, groups = groups)
  cat("Generating data overview plot...\n")
  data_overview_plot <- MOFA2::plot_data_overview(MOFAobject)
  if (save_model && dir.exists(dirname(outfile))) {
    plot_file <- file.path(dirname(outfile), "data_overview_plot.pdf")
    pdf(plot_file, width = 10, height = 7)
    print(data_overview_plot)
    dev.off()
    cat(sprintf("Data overview plot saved to: %s\n", plot_file))
  }
  cat("Setting data options...\n")
  data_opts <- MOFA2::get_default_data_options(MOFAobject)
  data_opts$scale_views <- scale_views
  cat("Setting model options...\n")
  model_opts <- MOFA2::get_default_model_options(MOFAobject)
  model_opts$num_factors <- num_factors
  cat("Setting training options...\n")
  train_opts <- MOFA2::get_default_training_options(MOFAobject)
  train_opts$convergence_mode <- convergence_mode
  train_opts$seed <- seed
  train_opts$maxiter <- max_iter
  train_opts$drop_factor_threshold <- drop_factor_threshold
  train_opts$gpu_mode <- use_gpu
  cat("Preparing MOFA object...\n")
  MOFAobject <- MOFA2::prepare_mofa(MOFAobject,
                                    data_options = data_opts,
                                    model_options = model_opts,
                                    training_options = train_opts)
  cat("Training MOFA model...\n")
  if (save_model && !is.null(outfile)) {
    MOFAobject <- MOFA2::run_mofa(MOFAobject, outfile = outfile)
  } else {
    MOFAobject <- MOFA2::run_mofa(MOFAobject)
  }
  if (save_model && !is.null(rds_file)) {
    cat(sprintf("Saving model to %s...\n", rds_file))
    saveRDS(MOFAobject, rds_file)
  }
  end_time <- Sys.time()
  execution_time <- difftime(end_time, start_time, units = "mins")
  cat(sprintf("MOFA analysis completed in %.2f minutes\n", as.numeric(execution_time)))
  return(MOFAobject)
}

# function to extract top features for each factor
extract_top_features <- function(data, factor_column, value_column, factors = 1:9, output_dir = "./results/Top_Features", top_n = 1000) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  results_list <- list()
  for (factor_id in factors) {
    factor_name <- paste0("Factor", factor_id)
    factor_weights <- subset(data, data[[factor_column]] == factor_name)
    if (nrow(factor_weights) == 0) {
      message(paste("No data found for", factor_name))
      next
    }
    all_features_path <- file.path(output_dir, paste0(factor_name, "_all_features.xlsx"))
    write.xlsx(factor_weights, all_features_path)
    factor_top_features <- factor_weights[order(-abs(factor_weights[[value_column]])), ][1:min(top_n, nrow(factor_weights)), ]
    top_features_path <- file.path(output_dir, paste0(factor_name, "_top_features.xlsx"))
    write.xlsx(factor_top_features, top_features_path)
    results_list[[factor_name]] <- list(
      all_features = factor_weights,
      top_features = factor_top_features
    )
    message(paste("Saved results for", factor_name, "- All:", all_features_path, "Top:", top_features_path))
  }
  message("Feature extraction completed for all selected factors.")
  return(results_list)
} 