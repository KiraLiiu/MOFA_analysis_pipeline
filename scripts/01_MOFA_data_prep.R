# =============================================================================
# MOFA + MultiGSEA Pipeline - Step 01: Data Preparation
# =============================================================================
# This script prepares multi-omics data for MOFA analysis by performing
# quality control, feature selection, and data normalization across
# transcriptomics, proteomics, and metabolomics datasets.
# 
# Author: Kira Liu
# Date: 2025-01-27
# 
# Input: Normalized omics datasets from preprocessing pipeline
# Output: MOFA-ready data structure and quality control plots
# =============================================================================

# Load required libraries
library(MOFA2)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(openxlsx)
library(data.table)

# =============================================================================
# Source feature selection and visualization utilities
# =============================================================================
source("functions/Omics_Feature_Selection_Utils.R")

# =============================================================================
# SET OUTPUT DIRECTORIES
# =============================================================================
# Create organized output structure
output_base <- "results/mofa_analysis"
dirs_to_create <- c(
  file.path(output_base, "data_preparation"),
  file.path(output_base, "data_preparation/elbow_plots"),
  file.path(output_base, "data_preparation/distribution_analysis"),
  file.path(output_base, "data_preparation/feature_selection")
)

for (dir in dirs_to_create) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

# =============================================================================
# 1. LOAD NORMALIZED OMICS DATASETS (WITHOUT OUTLIERS)
# =============================================================================
cat("Loading normalized omics datasets (without outliers)...\n")
metabolomics_normalized <- read.csv("./data/processed/metabolomics/outlier_removed/Limma_results/normalized/metabolomics_normalized.csv", check.names = FALSE, row.names = 1)
proteomics_normalized <- read.csv("./data/processed/proteomics/outlier_removed/Limma_results/normalized/proteomics_normalized.csv", check.names = FALSE, row.names = 1)
transcriptomics_normalized <- read.csv("./data/processed/transcriptomics/DESeq2_results/raw_results/vst_normalized_counts.csv", check.names = FALSE, row.names = 1)

# =============================================================================
# 2. DATA QUALITY CHECK AND DIMENSIONS
# =============================================================================
cat("Checking data dimensions and quality...\n")
cat("Original data dimensions:\n")
cat(sprintf("Metabolomics: %d features x %d samples\n", nrow(metabolomics_normalized), ncol(metabolomics_normalized)))
cat(sprintf("Proteomics: %d features x %d samples\n", nrow(proteomics_normalized), ncol(proteomics_normalized)))
cat(sprintf("Transcriptomics: %d features x %d samples\n", nrow(transcriptomics_normalized), ncol(transcriptomics_normalized)))
cat("\nMissing values check:\n")
cat(sprintf("Metabolomics missing values: %d\n", sum(is.na(metabolomics_normalized))))
cat(sprintf("Proteomics missing values: %d\n", sum(is.na(proteomics_normalized))))
cat(sprintf("Transcriptomics missing values: %d\n", sum(is.na(transcriptomics_normalized))))

# =============================================================================
# 3. CUMULATIVE VARIANCE ANALYSIS (ELBOW PLOTS)
# =============================================================================
cat("\n", strrep("=", 60), "\n", sep = "")
cat("CUMULATIVE VARIANCE ANALYSIS (ELBOW PLOTS)\n")
cat(strrep("=", 60), "\n", sep = "")
cat("Generating cumulative variance curves to help select HVF thresholds...\n")

elbow_dir <- file.path(output_base, "data_preparation/elbow_plots")
transcriptomics_elbow <- plot_cumulative_variance(transcriptomics_normalized, "Transcriptomics", save_plot = TRUE, output_dir = elbow_dir)
proteomics_elbow <- plot_cumulative_variance(proteomics_normalized, "Proteomics", save_plot = TRUE, output_dir = elbow_dir)
metabolomics_elbow <- plot_cumulative_variance(metabolomics_normalized, "Metabolomics", save_plot = TRUE, output_dir = elbow_dir)
cat("\nElbow plots generated! Check the plots in:", elbow_dir, "\n")
cat("Use these plots to decide on appropriate variance thresholds for HVF selection.\n")

# =============================================================================
# 4. HIGH VARIANCE FEATURE SELECTION (FOR ALL THREE OMICS LAYERS)
# =============================================================================
cat("\n", strrep("=", 60), "\n", sep = "")
cat("HIGH VARIANCE FEATURE SELECTION\n")
cat(strrep("=", 60), "\n", sep = "")
cat("Based on the elbow plots above, you can adjust these thresholds:\n")
cat("Current thresholds: Transcriptomics=0.75, Proteomics=0.90, Metabolomics=0.95\n")
cat("To change thresholds, modify the variance_threshold values below.\n\n")

cat("Performing high variance feature selection for all omics layers...\n")
transcriptomics_hvf <- select_high_variance_features(transcriptomics_normalized, variance_threshold = 0.75)
proteomics_hvf <- select_high_variance_features(proteomics_normalized, variance_threshold = 0.90)
metabolomics_hvf <- select_high_variance_features(metabolomics_normalized, variance_threshold = 0.95)

cat("High variance features selected:\n")
cat(sprintf("Transcriptomics: %d features (from %d, %.1f%%)\n", 
            nrow(transcriptomics_hvf), nrow(transcriptomics_normalized), 
            nrow(transcriptomics_hvf)/nrow(transcriptomics_normalized)*100))
cat(sprintf("Proteomics: %d features (from %d, %.1f%%)\n", 
            nrow(proteomics_hvf), nrow(proteomics_normalized), 
            nrow(proteomics_hvf)/nrow(proteomics_normalized)*100))
cat(sprintf("Metabolomics: %d features (from %d, %.1f%%)\n", 
            nrow(metabolomics_hvf), nrow(metabolomics_normalized), 
            nrow(metabolomics_hvf)/nrow(metabolomics_normalized)*100))

# =============================================================================
# 5. ANALYZE MEAN AND STANDARD DEVIATION DISTRIBUTIONS (POST-HVF)
# =============================================================================
cat("\nAnalyzing mean and standard deviation distributions after HVF selection...\n")

dist_dir <- file.path(output_base, "data_preparation/distribution_analysis")
metabolomics_stats <- analyze_distributions(metabolomics_hvf, "Metabolomics", dist_dir)
proteomics_stats <- analyze_distributions(proteomics_hvf, "Proteomics", dist_dir)
transcriptomics_stats <- analyze_distributions(transcriptomics_hvf, "Transcriptomics", dist_dir)

# Save distribution statistics
write.csv(metabolomics_stats, file.path(dist_dir, "01_Metabolomics_distributions_post_HVF.csv"), row.names = FALSE)
write.csv(proteomics_stats, file.path(dist_dir, "02_Proteomics_distributions_post_HVF.csv"), row.names = FALSE)
write.csv(transcriptomics_stats, file.path(dist_dir, "03_Transcriptomics_distributions_post_HVF.csv"), row.names = FALSE)
cat("\nDistribution analysis completed! Check the plots in:", dist_dir, "\n")
cat("Manual assessment: If distributions are very different across layers, consider z-score normalization.\n")

# =============================================================================
# 6. DECIDE ON Z-SCORE NORMALIZATION (MANUAL ASSESSMENT)
# =============================================================================
cat("\n", strrep("=", 60), "\n", sep = "")
cat("Z-SCORE NORMALIZATION DECISION\n")
cat(strrep("=", 60), "\n", sep = "")
cat("Based on the distribution analysis above, decide if z-score normalization is needed.\n")
cat("If the mean and SD distributions are very different across omics layers,\n")
cat("you may want to apply z-score normalization to make them more comparable.\n\n")

apply_zscore <- TRUE  # Change to TRUE if z-score is needed
if (apply_zscore) {
  cat("Applying z-score normalization to all omics layers...\n")
  apply_zscore_normalization <- function(data) {
    data_scaled <- t(scale(t(data)))
    return(as.data.frame(data_scaled))
  }
  transcriptomics_final <- apply_zscore_normalization(transcriptomics_hvf)
  proteomics_final <- apply_zscore_normalization(proteomics_hvf)
  metabolomics_final <- apply_zscore_normalization(metabolomics_hvf)
  cat("Z-score normalization applied.\n")
} else {
  cat("No z-score normalization applied. Using HVF-selected data as is.\n")
  transcriptomics_final <- transcriptomics_hvf
  proteomics_final <- proteomics_hvf
  metabolomics_final <- metabolomics_hvf
}

# =============================================================================
# 7. CREATE MOFA DATA STRUCTURE
# =============================================================================
# =============================================================================
# 补全所有 sample，缺的列补 NA，顺序对齐
# =============================================================================
cat("\nAligning samples across all omics layers (filling missing samples with NA)...\n")
all_samples <- Reduce(union, list(
  colnames(transcriptomics_final),
  colnames(proteomics_final),
  colnames(metabolomics_final)
))

align_samples <- function(mat, all_samples) {
  aligned <- matrix(NA, nrow = nrow(mat), ncol = length(all_samples),
                    dimnames = list(rownames(mat), all_samples))
  idx <- match(colnames(mat), all_samples)
  aligned[, idx[!is.na(idx)]] <- as.matrix(mat[, !is.na(idx), drop = FALSE])
  return(as.data.frame(aligned))
}

transcriptomics_final <- align_samples(transcriptomics_final, all_samples)
proteomics_final     <- align_samples(proteomics_final, all_samples)
metabolomics_final   <- align_samples(metabolomics_final, all_samples)

cat("Sample alignment done. All omics layers now have the same samples and order.\n")

cat("\nCreating MOFA data structure...\n")
mofa_data <- list(
  rna = as.matrix(transcriptomics_final),
  protein = as.matrix(proteomics_final),
  metabolite = as.matrix(metabolomics_final)
)

# =============================================================================
# 8. CREATE MOFA OBJECT AND VISUALIZE DATA STRUCTURE
# =============================================================================
cat("Creating MOFA object and generating data overview...\n")
MOFAobject <- create_mofa(mofa_data)
data_overview_plot <- plot_data_overview(MOFAobject)
ggsave(file.path(output_base, "01_Data_Overview_Plot.pdf"), plot = data_overview_plot, width = 10, height = 7)

# =============================================================================
# 9. SAVE PREPARED DATA
# =============================================================================
cat("Saving prepared data...\n")
feature_dir <- file.path(output_base, "data_preparation/feature_selection")

# Save MOFA data structure
saveRDS(mofa_data, file.path(output_base, "01_MOFA_Data_Structure.rds"))

# Save individual feature matrices
write.csv(transcriptomics_final, file.path(feature_dir, "01_Transcriptomics_Selected_Features.csv"))
write.csv(proteomics_final, file.path(feature_dir, "02_Proteomics_Selected_Features.csv"))
write.csv(metabolomics_final, file.path(feature_dir, "03_Metabolomics_Selected_Features.csv"))

# Save MOFA object
saveRDS(MOFAobject, file.path(output_base, "01_MOFA_Object_Created.rds"))

# Save feature selection summary
feature_summary <- data.frame(
  Omics_Layer = c("Transcriptomics", "Proteomics", "Metabolomics"),
  Original_Features = c(nrow(transcriptomics_normalized), nrow(proteomics_normalized), nrow(metabolomics_normalized)),
  Selected_Features = c(nrow(transcriptomics_final), nrow(proteomics_final), nrow(metabolomics_final)),
  Selection_Percentage = c(
    round(nrow(transcriptomics_final)/nrow(transcriptomics_normalized)*100, 1),
    round(nrow(proteomics_final)/nrow(proteomics_normalized)*100, 1),
    round(nrow(metabolomics_final)/nrow(metabolomics_normalized)*100, 1)
  ),
  Samples = c(ncol(transcriptomics_final), ncol(proteomics_final), ncol(metabolomics_final)),
  Z_Score_Normalization = c(apply_zscore, apply_zscore, apply_zscore)
)
write.csv(feature_summary, file.path(feature_dir, "04_Feature_Selection_Summary.csv"), row.names = FALSE)

# =============================================================================
# 10. SUMMARY REPORT
# =============================================================================
cat("\n", strrep("=", 50), "\n", sep = "")
cat("DATA PREPARATION SUMMARY\n")
cat(strrep("=", 50), "\n", sep = "")
cat(sprintf("Transcriptomics: %d features x %d samples\n", nrow(transcriptomics_final), ncol(transcriptomics_final)))
cat(sprintf("Proteomics: %d features x %d samples\n", nrow(proteomics_final), ncol(proteomics_final)))
cat(sprintf("Metabolomics: %d features x %d samples\n", nrow(metabolomics_final), ncol(metabolomics_final)))
cat(sprintf("Total features: %d\n", sum(nrow(transcriptomics_final), nrow(proteomics_final), nrow(metabolomics_final))))
cat(sprintf("Z-score normalization: %s\n", ifelse(apply_zscore, "Applied", "Not applied")))
cat(sprintf("Output directory: %s\n", output_base))
cat("Data preparation completed successfully!\n")
cat(strrep("=", 50), "\n", sep = "")

rm(list = setdiff(ls(), c("mofa_data", "MOFAobject", "output_base")))
gc()
cat("Step 01: Data preparation completed!\n")
cat("Next step: Run 02_MOFA_training.R\n") 
