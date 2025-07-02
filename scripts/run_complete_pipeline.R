# =============================================================================
# MOFA + MultiGSEA Complete Pipeline Runner
# =============================================================================
# This script runs the complete MOFA + MultiGSEA pipeline sequentially
# with comprehensive error handling and logging for robust execution.
# 
# Author: Multi-omics Analysis Team
# Date: 2025-01-27
# Version: 1.0.0
# 
# Input: Preprocessed multi-omics data
# Output: Complete MOFA and MultiGSEA analysis results
# =============================================================================

# Load required libraries
library(MOFA2)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(openxlsx)
library(data.table)
library(multiGSEA)

# Source standalone functions (outside ref)
source("functions/MOFA_analysis_function.R")
source("functions/Complex_Heatmap_for_MultiGSEA_Results.R")

# =============================================================================
# CONFIGURATION
# =============================================================================

# Pipeline configuration
pipeline_config <- list(
  # MOFA parameters
  num_factors = 3,
  convergence_mode = "slow",
  scale_views = FALSE,
  seed = 42,
  max_iter = 1000,
  
  # Feature selection
  variance_threshold = 1.0,
  
  # MultiGSEA parameters
  pathway_databases = c("kegg", "reactome"),
  organism = "mmusculus",
  
  # Output settings
  save_plots = TRUE,
  plot_formats = c("svg", "pdf"),
  
  # Step control
  run_step_01 = TRUE,
  run_step_02 = TRUE,
  run_step_03 = TRUE,
  run_step_04 = TRUE,
  run_step_05 = TRUE
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Function to log messages with timestamps
log_message <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s: %s\n", timestamp, level, message))
}

# Function to check if file exists
check_file_exists <- function(file_path, step_name) {
  if (!file.exists(file_path)) {
    log_message(sprintf("Required file not found: %s", file_path), "ERROR")
    log_message(sprintf("Step %s cannot proceed", step_name), "ERROR")
    return(FALSE)
  }
  return(TRUE)
}

# Function to create directory if it doesn't exist
ensure_directory <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    log_message(sprintf("Created directory: %s", dir_path))
  }
}

# Function to run a step with error handling
run_pipeline_step <- function(step_name, step_function, required_files = NULL) {
  log_message(sprintf("Starting %s", step_name), "INFO")
  
  # Check required files
  if (!is.null(required_files)) {
    for (file_path in required_files) {
      if (!check_file_exists(file_path, step_name)) {
        return(FALSE)
      }
    }
  }
  
  # Run the step
  tryCatch({
    start_time <- Sys.time()
    result <- step_function()
    end_time <- Sys.time()
    execution_time <- difftime(end_time, start_time, units = "mins")
    
    log_message(sprintf("%s completed successfully in %.2f minutes", 
                       step_name, as.numeric(execution_time)), "SUCCESS")
    return(TRUE)
  }, error = function(e) {
    log_message(sprintf("Error in %s: %s", step_name, e$message), "ERROR")
    return(FALSE)
  })
}

# =============================================================================
# STEP 01: DATA PREPARATION
# =============================================================================

step_01_data_preparation <- function() {
  log_message("Loading normalized omics datasets...")
  
  # Load metabolomics data
  metabolomics_normalized <- read.csv("data/processed/metabolomics_normalized/metabolomics_normalized.csv", 
                                     check.names = FALSE, row.names = 1)
  
  # Load proteomics data
  proteomics_normalized <- read.csv("data/processed/proteomics_normalized/proteomics_normalized.csv", 
                                   check.names = FALSE, row.names = 1)
  
  # Load transcriptomics data
  transcriptomics_normalized <- read.csv("data/processed/RNAseq_DESeq2_results/Results_without_annotation/vst_normalized_counts.csv", 
                                        check.names = FALSE, row.names = 1)
  
  # Filter transcriptomics data
  log_message("Filtering transcriptomics data for significant genes...")
  transcriptomics_summary <- read.csv("data/processed/RNAseq_DESeq2_results/Results_without_annotation/DESeq2_HALI vs Control_results.csv", 
                                     check.names = FALSE)
  
  colnames(transcriptomics_summary)[1] <- "gene_id"
  colnames(transcriptomics_normalized)[1] <- "gene_id"
  
  transcriptomics_normalized <- transcriptomics_normalized %>% 
    merge(transcriptomics_summary, by = "gene_id") %>%
    filter(padj < 0.05) %>%
    select(-c(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj))
  
  rownames(transcriptomics_normalized) <- transcriptomics_normalized$gene_id
  transcriptomics_normalized <- transcriptomics_normalized %>% select(-gene_id)
  
  # Feature selection
  log_message("Performing feature selection...")
  source("scripts/functions/MOFA_analysis_function.R")
  
  transcriptomics_top_features <- select_high_variance_features(transcriptomics_normalized, 
                                                               variance_threshold = pipeline_config$variance_threshold)
  proteomics_top_features <- select_high_variance_features(proteomics_normalized, 
                                                          variance_threshold = pipeline_config$variance_threshold)
  metabolomics_top_features <- select_high_variance_features(metabolomics_normalized, 
                                                            variance_threshold = pipeline_config$variance_threshold)
  
  # Create MOFA data structure
  mofa_data <- list(
    rna = as.matrix(transcriptomics_top_features),
    protein = as.matrix(proteomics_top_features),
    metabolite = as.matrix(metabolomics_top_features)
  )
  
  # Create output directory
  output_dir <- "results/mofa_analysis"
  ensure_directory(output_dir)
  
  # Save prepared data
  saveRDS(mofa_data, file.path(output_dir, "mofa_data_prepared.rds"))
  saveRDS(transcriptomics_top_features, file.path(output_dir, "transcriptomics_selected_features.rds"))
  saveRDS(proteomics_top_features, file.path(output_dir, "proteomics_selected_features.rds"))
  saveRDS(metabolomics_top_features, file.path(output_dir, "metabolomics_selected_features.rds"))
  
  # Create and save MOFA object
  MOFAobject <- create_mofa(mofa_data)
  saveRDS(MOFAobject, file.path(output_dir, "MOFA_object_created.rds"))
  
  # Generate data overview plot
  data_overview_plot <- plot_data_overview(MOFAobject)
  ggsave(file.path(output_dir, "data_overview_plot.pdf"), 
         plot = data_overview_plot, width = 10, height = 7)
  
  log_message("Data preparation completed successfully")
  return(TRUE)
}

# =============================================================================
# STEP 02: MOFA TRAINING
# =============================================================================

step_02_mofa_training <- function() {
  log_message("Loading prepared data...")
  
  # Load prepared data
  mofa_data <- readRDS("results/mofa_analysis/mofa_data_prepared.rds")
  
  # Create output directory
  output_dir <- "results/mofa_analysis"
  hdf5_file <- file.path(output_dir, "MOFA_multi-omics.hdf5")
  rds_file <- file.path(output_dir, "MOFA_multi-omics.rds")
  
  # Run MOFA analysis
  log_message("Training MOFA model...")
  source("scripts/functions/MOFA_analysis_function.R")
  
  MOFAobject <- run_mofa_analysis(
    data_list = mofa_data,
    num_factors = pipeline_config$num_factors,
    scale_views = pipeline_config$scale_views,
    convergence_mode = pipeline_config$convergence_mode,
    seed = pipeline_config$seed,
    max_iter = pipeline_config$max_iter,
    outfile = hdf5_file,
    rds_file = rds_file
  )
  
  log_message("MOFA training completed successfully")
  return(TRUE)
}

# =============================================================================
# STEP 03: MOFA EXPLANATION
# =============================================================================

step_03_mofa_explanation <- function() {
  log_message("Loading trained MOFA model...")
  
  # Load trained model
  MOFAobject <- load_model("results/mofa_analysis/MOFA_multi-omics.hdf5")
  
  # Create output directories
  ensure_directory("results/MOFA_multiomics_results")
  ensure_directory("results/figures")
  
  # Generate factor correlation plot
  log_message("Generating factor correlation plot...")
  factor_cor_plot <- plot_factor_cor(MOFAobject)
  ggsave("results/figures/mofa_factor_correlation.svg", 
         plot = factor_cor_plot, width = 8, height = 6)
  
  # Generate variance explained plot
  log_message("Generating variance explained plot...")
  variance_plot <- plot_variance_explained(MOFAobject, max_r2 = 15)
  ggsave("results/figures/mofa_variance_explained.svg", 
         plot = variance_plot, width = 8, height = 6)
  
  # Load metadata and correlate with factors
  log_message("Correlating factors with covariates...")
  MOFA_metadata <- read.xlsx("data/metadata/sample_info.xlsx", sheet = 2)
  samples_metadata(MOFAobject) <- MOFA_metadata
  
  covariate_plot <- correlate_factors_with_covariates(MOFAobject, 
                                                     covariates = c("HALI", "Control"), 
                                                     factors = "all",
                                                     plot = "r")
  ggsave("results/figures/mofa_factor_correlation_with_covariates.svg", 
         plot = covariate_plot, width = 8, height = 6)
  
  # Extract top features
  log_message("Extracting top features for each factor...")
  weights_df <- get_weights(MOFAobject, as.data.frame = TRUE)
  
  source("scripts/functions/MOFA_analysis_function.R")
  mofa_results <- extract_top_features(weights_df, 
                                      factor_column = "factor", 
                                      value_column = "value", 
                                      factors = 1:pipeline_config$num_factors, 
                                      output_dir = "results/MOFA_multiomics_results", 
                                      top_n = 1000)
  
  log_message("MOFA explanation completed successfully")
  return(TRUE)
}

# =============================================================================
# STEP 04: MULTIGSEA ANALYSIS
# =============================================================================

step_04_multigsea_analysis <- function() {
  log_message("Loading MOFA factor results...")
  
  # Load factor results
  features_factor1 <- read.xlsx("results/MOFA_multiomics_results/Factor1_all_features.xlsx")
  
  # Separate features by omics type
  rna_features <- features_factor1 %>% filter(view == "rna")
  protein_features <- features_factor1 %>% filter(view == "protein")
  metabolite_features <- features_factor1 %>% filter(view == "metabolite")
  
  # Load original omics data
  rna_results <- read.csv("data/processed/RNAseq_DESeq2_results/Results_without_annotation/DESeq2_HALI vs Control_results.csv", 
                          check.names = FALSE)
  protein_results <- read.xlsx("data/raw/Protein_DEPs.xlsx", sheet = 1)
  metabolite_results <- read.xlsx("data/raw/Metabolomics_Clean.xlsx", sheet = 1)
  
  # Prepare data for MultiGSEA
  log_message("Preparing data for MultiGSEA analysis...")
  
  # Transcriptomics data
  df_transcriptome <- rna_features %>%
    left_join(rna_results, by = c("feature" = "gene_id")) %>%
    mutate(feature_cleaned = str_replace(feature, "\\.[0-9]+$", "")) %>%
    dplyr::select(c(feature_cleaned, value, pvalue, padj)) %>%
    drop_na()
  
  # Proteomics data
  df_proteome <- protein_features %>%
    left_join(protein_results, by = c("feature" = "Gene.Name")) %>%
    dplyr::select(c(feature, value, pvalue)) %>%
    mutate(padj = pvalue) %>%
    drop_na()
  
  # Metabolomics data
  metabolite_annotation <- read.csv("data/processed/metabolomics_normalized/name_map.csv")
  df_metabolome <- metabolite_features %>%
    left_join(metabolite_results, by = c("feature" = "Metabolites.Name")) %>%
    left_join(metabolite_annotation, by = c("feature" = "Query")) %>%
    dplyr::select(c(feature, HMDB, value, p.value)) %>%
    dplyr::rename(pvalue = p.value) %>%
    mutate(padj = pvalue) %>%
    drop_na()
  
  # Process duplicates in metabolomics
  df_metabolome <- df_metabolome %>%
    group_by(HMDB) %>%
    summarise(
      features_combined = paste(unique(feature), collapse = "; "),
      mean_value = mean(value, na.rm = TRUE),
      min_pvalue = min(pvalue, na.rm = TRUE),
      min_padj = min(padj, na.rm = TRUE)
    ) %>%
    ungroup()
  
  # Initialize MultiGSEA
  log_message("Initializing MultiGSEA analysis...")
  omics_data <- initOmicsDataStructure(layer = c("transcriptome", "proteome", "metabolome"))
  
  omics_data$transcriptome <- df_transcriptome$value
  names(omics_data$transcriptome) <- df_transcriptome$feature_cleaned
  
  omics_data$proteome <- df_proteome$value
  names(omics_data$proteome) <- df_proteome$feature
  
  omics_data$metabolome <- df_metabolome$mean_value
  names(omics_data$metabolome) <- df_metabolome$HMDB
  
  # Get pathways
  log_message("Loading pathway databases...")
  pathways <- getMultiOmicsFeatures(
    dbs = pipeline_config$pathway_databases, 
    layer = names(omics_data),
    returnTranscriptome = "ENSEMBL",
    returnProteome = "SYMBOL",
    returnMetabolome = "HMDB",
    organism = pipeline_config$organism,
    useLocal = TRUE
  )
  
  # Run MultiGSEA
  log_message("Running MultiGSEA analysis...")
  enrichment_scores <- multiGSEA(pathways, omics_data)
  
  # Extract results
  df_enrichment_results <- extractPvalues(enrichmentScores = enrichment_scores,
                                         pathwayNames = names(pathways[[1]]))
  
  df_enrichment_results$combined_pval <- combinePvalues(df_enrichment_results)
  df_enrichment_results$combined_padj <- p.adjust(df_enrichment_results$combined_pval, method = "BH")
  
  df_enrichment_results <- cbind(data.frame(pathway = names(pathways[[1]])), df_enrichment_results)
  
  # Clean results
  enrichment_results_cleaned <- df_enrichment_results %>%
    mutate(
      source = str_extract(pathway, "(?<=\\().+?(?=\\))"),
      pathway = str_replace(pathway, "^\\(.+?\\)\\s*", "")
    ) %>%
    dplyr::select(
      pathway, source,
      transcriptome_pval, proteome_pval, metabolome_pval, combined_pval,
      transcriptome_padj, proteome_padj, metabolome_padj, combined_padj
    )
  
  # Save results
  ensure_directory("results/tables")
  write.csv(enrichment_results_cleaned,
            file = "results/tables/MOFA_Factor1_multiGSEA_results.csv",
            row.names = FALSE)
  write.xlsx(enrichment_results_cleaned,
             file = "results/tables/MOFA_Factor1_multiGSEA_results.xlsx")
  
  log_message("MultiGSEA analysis completed successfully")
  return(TRUE)
}

# =============================================================================
# STEP 05: VISUALIZATION
# =============================================================================

step_05_visualization <- function() {
  log_message("Generating final visualizations...")
  
  # Load trained model
  MOFAobject <- load_model("results/mofa_analysis/MOFA_multi-omics.hdf5")
  
  # Create output directory
  ensure_directory("results/figures")
  
  # Generate MOFA visualizations
  log_message("Generating MOFA visualizations...")
  
  # Data overview
  data_overview_plot <- plot_data_overview(MOFAobject)
  ggsave("results/figures/mofa_data_overview.svg", 
         plot = data_overview_plot, width = 10, height = 7)
  
  # Factor correlation
  factor_cor_plot <- plot_factor_cor(MOFAobject)
  ggsave("results/figures/mofa_factor_correlation.svg", 
         plot = factor_cor_plot, width = 8, height = 6)
  
  # Variance explained
  variance_plot <- plot_variance_explained(MOFAobject, max_r2 = 15)
  ggsave("results/figures/mofa_variance_explained.svg", 
         plot = variance_plot, width = 8, height = 6)
  
  # Load metadata for covariate correlation
  MOFA_metadata <- read.xlsx("data/metadata/sample_info.xlsx", sheet = 2)
  samples_metadata(MOFAobject) <- MOFA_metadata
  
  covariate_plot <- correlate_factors_with_covariates(MOFAobject, 
                                                     covariates = c("HALI", "Control"), 
                                                     factors = "all",
                                                     plot = "r")
  ggsave("results/figures/mofa_factor_correlation_with_covariates.svg", 
         plot = covariate_plot, width = 8, height = 6)
  
  # Generate MultiGSEA heatmaps
  log_message("Generating MultiGSEA heatmaps...")
  source("scripts/functions/Complex_Heatmap_for_MultiGSEA_Results.R")
  
  # Load MultiGSEA results
  enrichment_results <- read.csv("results/tables/MOFA_Factor1_multiGSEA_results.csv")
  
  # Create heatmaps for different pathway sources
  ensure_directory("results/figures/MOFA_Factor1_MultiGSEA_heatmaps")
  
  # KEGG pathways heatmap
  kegg_heatmap <- generate_pathway_heatmaps(enrichment_results, 
                                           selected_sources = "KEGG", 
                                           show_source_annotation = TRUE,
                                           save_output = TRUE,
                                           output_format = "svg",
                                           plot_width = 28,
                                           plot_height = 15,
                                           output_dir = "results/figures/MOFA_Factor1_MultiGSEA_heatmaps")
  
  # Reactome pathways heatmap
  reactome_heatmap <- generate_pathway_heatmaps(enrichment_results, 
                                               selected_sources = "REACTOME", 
                                               show_source_annotation = FALSE,
                                               save_output = TRUE,
                                               output_format = "svg",
                                               plot_width = 35,
                                               plot_height = 15,
                                               output_dir = "results/figures/MOFA_Factor1_MultiGSEA_heatmaps")
  
  log_message("Visualization completed successfully")
  return(TRUE)
}

# =============================================================================
# MAIN PIPELINE EXECUTION
# =============================================================================

main_pipeline <- function() {
  log_message("Starting MOFA + MultiGSEA Pipeline", "INFO")
  log_message("Configuration:", "INFO")
  log_message(sprintf("  Number of factors: %d", pipeline_config$num_factors))
  log_message(sprintf("  Convergence mode: %s", pipeline_config$convergence_mode))
  log_message(sprintf("  Pathway databases: %s", paste(pipeline_config$pathway_databases, collapse = ", ")))
  
  # Track overall success
  pipeline_success <- TRUE
  
  # Step 01: Data Preparation
  if (pipeline_config$run_step_01) {
    required_files_01 <- c(
      "data/processed/metabolomics_normalized/metabolomics_normalized.csv",
      "data/processed/proteomics_normalized/proteomics_normalized.csv",
      "data/processed/RNAseq_DESeq2_results/Results_without_annotation/vst_normalized_counts.csv",
      "data/processed/RNAseq_DESeq2_results/Results_without_annotation/DESeq2_HALI vs Control_results.csv"
    )
    
    if (!run_pipeline_step("Step 01: Data Preparation", step_01_data_preparation, required_files_01)) {
      pipeline_success <- FALSE
    }
  }
  
  # Step 02: MOFA Training
  if (pipeline_config$run_step_02 && pipeline_success) {
    required_files_02 <- c("results/mofa_analysis/mofa_data_prepared.rds")
    
    if (!run_pipeline_step("Step 02: MOFA Training", step_02_mofa_training, required_files_02)) {
      pipeline_success <- FALSE
    }
  }
  
  # Step 03: MOFA Explanation
  if (pipeline_config$run_step_03 && pipeline_success) {
    required_files_03 <- c(
      "results/mofa_analysis/MOFA_multi-omics.hdf5",
      "data/metadata/sample_info.xlsx"
    )
    
    if (!run_pipeline_step("Step 03: MOFA Explanation", step_03_mofa_explanation, required_files_03)) {
      pipeline_success <- FALSE
    }
  }
  
  # Step 04: MultiGSEA Analysis
  if (pipeline_config$run_step_04 && pipeline_success) {
    required_files_04 <- c(
      "results/MOFA_multiomics_results/Factor1_all_features.xlsx",
      "data/raw/Protein_DEPs.xlsx",
      "data/raw/Metabolomics_Clean.xlsx",
      "data/processed/metabolomics_normalized/name_map.csv"
    )
    
    if (!run_pipeline_step("Step 04: MultiGSEA Analysis", step_04_multigsea_analysis, required_files_04)) {
      pipeline_success <- FALSE
    }
  }
  
  # Step 05: Visualization
  if (pipeline_config$run_step_05 && pipeline_success) {
    required_files_05 <- c(
      "results/mofa_analysis/MOFA_multi-omics.hdf5",
      "results/tables/MOFA_Factor1_multiGSEA_results.csv"
    )
    
    if (!run_pipeline_step("Step 05: Visualization", step_05_visualization, required_files_05)) {
      pipeline_success <- FALSE
    }
  }
  
  # Final summary
  if (pipeline_success) {
    log_message("Pipeline completed successfully!", "SUCCESS")
    log_message("Results are available in the 'results' directory", "INFO")
  } else {
    log_message("Pipeline failed. Check error messages above.", "ERROR")
  }
  
  return(pipeline_success)
}

# =============================================================================
# EXECUTE PIPELINE
# =============================================================================

# Run the pipeline
if (interactive()) {
  # If running interactively, ask for confirmation
  cat("This will run the complete MOFA + MultiGSEA pipeline.\n")
  cat("Make sure all required data files are in place.\n")
  response <- readline("Do you want to continue? (y/n): ")
  
  if (tolower(response) %in% c("y", "yes")) {
    main_pipeline()
  } else {
    cat("Pipeline execution cancelled.\n")
  }
} else {
  # If running non-interactively (e.g., from command line), run directly
  main_pipeline()
} 