# =============================================================================
# MOFA + MultiGSEA Pipeline - Step 02: MOFA Training
# =============================================================================
# This script trains the MOFA model using prepared multi-omics data
# to identify shared factors across transcriptomics, proteomics, and
# metabolomics datasets.
# 
# Author: Kira Liu
# Date: 2025-01-27
# 
# Input: MOFA-ready data structure from Step 01
# Output: Trained MOFA model (HDF5 and RDS formats)
# =============================================================================

library(MOFA2)
library(tidyverse)
library(openxlsx)
library(data.table)

# Source standalone MOFA functions
source("functions/MOFA_analysis_function.R")

# =============================================================================
# SET OUTPUT DIRECTORIES
# =============================================================================
output_base <- "results/mofa_analysis"
training_dir <- file.path(output_base, "model_training")

if (!dir.exists(training_dir)) {
  dir.create(training_dir, recursive = TRUE)
}

# =============================================================================
# 1. LOAD PREPARED DATA
# =============================================================================
cat("Loading prepared data from Step 01...\n")
mofa_data <- readRDS(file.path(output_base, "01_MOFA_Data_Structure.rds"))

# =============================================================================
# 2. SET MOFA PARAMETERS
# =============================================================================
num_factors <- 3
convergence_mode <- "slow"
scale_views <- FALSE
seed <- 42
max_iter <- 1000

# Define output files with organized naming
hdf5_file <- file.path(training_dir, "02_MOFA_Model_Trained.hdf5")
rds_file <- file.path(training_dir, "02_MOFA_Model_Trained.rds")

# Save training parameters
training_params <- data.frame(
  Parameter = c("Number_of_Factors", "Convergence_Mode", "Scale_Views", "Seed", "Max_Iterations"),
  Value = c(num_factors, convergence_mode, scale_views, seed, max_iter)
)
write.csv(training_params, file.path(training_dir, "02_Training_Parameters.csv"), row.names = FALSE)

# =============================================================================
# 3. RUN MOFA ANALYSIS
# =============================================================================
cat("Running MOFA analysis...\n")
cat(sprintf("Parameters: %d factors, %s convergence, max %d iterations\n", 
            num_factors, convergence_mode, max_iter))

MOFAobject <- run_mofa_analysis(
  data_list = mofa_data,
  num_factors = num_factors,
  scale_views = scale_views,
  convergence_mode = convergence_mode,
  seed = seed,
  max_iter = max_iter,
  outfile = hdf5_file,
  rds_file = rds_file
)

# =============================================================================
# 4. SAVE TRAINING SUMMARY
# =============================================================================
cat("Saving training summary...\n")

# Get model information
model_info <- data.frame(
  Metric = c("Number_of_Factors", "Number_of_Views", "Number_of_Samples", "Converged", "ELBO_Value"),
  Value = c(
    get_dimensions(MOFAobject)$K,
    get_dimensions(MOFAobject)$M,
    get_dimensions(MOFAobject)$N,
    is_converged(MOFAobject),
    round(get_elbo(MOFAobject), 2)
  )
)
write.csv(model_info, file.path(training_dir, "02_Model_Information.csv"), row.names = FALSE)

cat("\nMOFA training completed!\n")
cat(sprintf("Model saved as: %s and %s\n", hdf5_file, rds_file))
cat(sprintf("Training outputs saved in: %s\n", training_dir))
cat("Next step: Run 03_MOFA_explanation.R\n") 

# =============================================================================
# clean up
# =============================================================================
rm(list = setdiff(ls(), c("MOFAobject", "output_base")))
gc()
