# =============================================================================
# MOFA + MultiGSEA Pipeline - Step 03: MOFA Explanation
# =============================================================================
# This script analyzes and interprets the trained MOFA model to extract
# factor information, variance explained, and feature weights for
# downstream analysis and visualization.
# 
# Author: Kira Liu
# Date: 2025-01-27
# 
# Input: Trained MOFA model from Step 02
# Output: Factor tables, variance explained data, and summary reports
# Note: Visualization plots are generated in Step 05
# =============================================================================

library(MOFA2)
library(tidyverse)
library(openxlsx)
library(data.table)
library(ggplot2)

# Source MOFA analysis functions
source("functions/MOFA_analysis_function.R")

# =============================================================================
# Set output directories
# =============================================================================
output_base <- "results/mofa_analysis"
explanation_dir <- file.path(output_base, "model_explanation")
tables_dir <- file.path(explanation_dir, "tables")

for (dir in c(explanation_dir, tables_dir)) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# =============================================================================
# 1. Load trained MOFA model
# =============================================================================
MOFAobject <- load_model(file.path(output_base, "model_training/02_MOFA_Model_Trained.hdf5"))

# =============================================================================
# 2. Prepare model explanation data (no plot saving here)
# =============================================================================
# All plots are generated in 05_Visualization.R. Only data for downstream visualization is prepared here.
plot_factor_cor(MOFAobject)
plot_variance_explained(MOFAobject, max_r2 = 15)
plot_variance_explained(MOFAobject, plot_total = T)[[2]]

MOFA_metadata <- read.xlsx("data/metadata/sample_info.xlsx", sheet = 2)
samples_metadata(MOFAobject) <- MOFA_metadata
correlate_factors_with_covariates(MOFAobject, covariates = c("HALI", "Control"), factors = "all", plot = "r")

# =============================================================================
# 3. Extract top features for each factor
# =============================================================================
weights_df <- get_weights(MOFAobject, as.data.frame = TRUE)
extract_top_features(weights_df, factor_column = "factor", value_column = "value", factors = 1:3, output_dir = tables_dir, top_n = 1000)

# =============================================================================
# 4. Save additional model information
# =============================================================================
factor_values <- get_factors(MOFAobject, as.data.frame = TRUE)
write.csv(factor_values, file.path(tables_dir, "03_Factor_Values_All_Samples.csv"), row.names = FALSE)

variance_explained <- get_variance_explained(MOFAobject, as.data.frame = TRUE)
write.csv(variance_explained, file.path(tables_dir, "03_Variance_Explained_By_Factors.csv"), row.names = FALSE)

dimensions <- get_dimensions(MOFAobject)
dimensions_df <- data.frame(
  Dimension = c("Number_of_Factors", "Number_of_Views", "Number_of_Samples"),
  Value = c(dimensions$K, dimensions$M, dimensions$N)
)
write.csv(dimensions_df, file.path(tables_dir, "03_Model_Dimensions.csv"), row.names = FALSE)

# =============================================================================
# 5. Create summary report
# =============================================================================
summary_report <- data.frame(
  Analysis_Step = c("Model_Loading", "Factor_Correlation", "Variance_Explained", 
                    "Covariate_Correlation", "Top_Features_Extraction", "Additional_Info_Saved"),
  Status = rep("Completed", 6),
  Output_Files = c(
    "Model loaded successfully",
    "Factor correlation data (for 05_Visualization.R)",
    "Variance explained data (for 05_Visualization.R)",
    "Factor-covariate correlation data (for 05_Visualization.R)",
    "Top features CSV files in tables directory",
    "Model information CSV files in tables directory"
  )
)
write.csv(summary_report, file.path(explanation_dir, "03_Analysis_Summary_Report.csv"), row.names = FALSE)

message("MOFA explanation completed. Tables saved in: ", tables_dir)
message("Summary report saved in: ", explanation_dir)
message("All visualization is handled in 05_Visualization.R.")
message("Next step: Run 04_MultiGSEA.R") 
