# =============================================================================
# Bioinformatics Project Setup Script
# Multi-omics Analysis Pipeline
# =============================================================================
# This script sets up the directory structure for a multi-omics analysis project
# Author: Multi-omics Analysis Pipeline
# Date: 2024

library(fs)
library(here)

# =============================================================================
# 1. Define project structure
# =============================================================================
project_structure <- list(
  # Core analysis directories
  "scripts" = list(
    "multi-omics" = list(
      "01_MOFA_data_prep.R",
      "02_MOFA_training.R", 
      "03_MOFA_explaination.R",
      "04_MultiGSEA.R",
      "05_Visualization.R",
      "run_complete_pipeline.R",
      "OUTPUT_STRUCTURE_OPTIMIZATION.md",
      "MOFA_MultiGSEA_Pipeline_Overview.md"
    ),
    "preprocessing" = list(
      "Metabolomics" = list(
        "01_Data_Processing_Metabolomics.R",
        "02_Data_Processing_Metabolomics_RemoveOutliers.R",
        "README.md"
      ),
      "Proteomics" = list(
        "01_Data_Processing_Proteomics.R",
        "02_Data_Processing_Proteomics_RemoveOutliers.R",
        "README.md"
      ),
      "Transcriptomics" = list(
        "01_Data_Processing_Transcriptomics.R",
        "02_Gene_Annotation.R",
        "03_Duplicate_Gene_Handler.R",
        "README.md"
      ),
      "README.md"
    )
  ),
  
  # Data directories
  "data" = list(
    "raw" = list(
      "metabolomics" = list(),
      "proteomics" = list(),
      "transcriptomics" = list()
    ),
    "processed" = list(
      "metabolomics_normalized" = list(),
      "proteomics_normalized" = list(),
      "transcriptomics_normalized" = list()
    ),
    "metadata" = list(
      "sample_info.xlsx"
    )
  ),
  
  # Results directories
  "results" = list(
    "mofa_analysis" = list(
      "data_preparation" = list(
        "elbow_plots" = list(),
        "distribution_analysis" = list(),
        "feature_selection" = list()
      ),
      "model_training" = list(),
      "model_explanation" = list(
        "figures" = list(),
        "tables" = list()
      ),
      "multigsea_analysis" = list(
        "tables" = list(),
        "selected_pathways" = list()
      ),
      "visualization" = list(
        "mofa_figures" = list(),
        "multigsea_figures" = list(
          "pathway_heatmaps" = list(),
          "complete_pathways" = list()
        )
      )
    ),
    "tables" = list(),
    "figures" = list()
  ),
  
  # Functions and utilities
  "functions" = list(
    "MOFA_analysis_function.R",
    "Complex_Heatmap_for_MultiGSEA_Results.R",
    "Data_Quality_Check_Functions.R",
    "Duplicate_Gene_Handler_Functions.R",
    "Gene_Annotation_Functions.R"
  ),
  
  # Configuration and documentation
  "config" = list(
    "config.yaml"
  ),
  "doc" = list(
    "reference" = list()
  ),
  
  # Project files
  "README.md",
  ".gitignore",
  ".Rprofile"
)

# =============================================================================
# 2. Create directory structure
# =============================================================================
create_project_structure <- function(structure, base_path = ".") {
  cat("Creating project directory structure...\n")
  
  for (item in names(structure)) {
    if (is.list(structure[[item]])) {
      # Create directory
      dir_path <- file.path(base_path, item)
      if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
        cat("Created directory:", dir_path, "\n")
      }
      
      # Recursively create subdirectories
      create_project_structure(structure[[item]], dir_path)
    } else {
      # Create file
      file_path <- file.path(base_path, item)
      if (!file.exists(file_path)) {
        file.create(file_path)
        cat("Created file:", file_path, "\n")
      }
    }
  }
}

# =============================================================================
# 3. Create configuration file
# =============================================================================
create_config_file <- function() {
  config_content <- '
# Multi-omics Analysis Configuration
# =================================

# Data paths
data:
  raw:
    metabolomics: "data/raw/metabolomics/"
    proteomics: "data/raw/proteomics/"
    transcriptomics: "data/raw/transcriptomics/"
  processed:
    metabolomics: "data/processed/metabolomics_normalized/"
    proteomics: "data/processed/proteomics_normalized/"
    transcriptomics: "data/processed/transcriptomics_normalized/"
  metadata: "data/metadata/sample_info.xlsx"

# Analysis parameters
analysis:
  mofa:
    num_factors: 3
    convergence_mode: "slow"
    scale_views: false
    seed: 42
  multigsea:
    databases: ["kegg", "reactome"]
    organism: "mmusculus"
    pvalue_threshold: 0.05

# Output settings
output:
  base_dir: "results/mofa_analysis/"
  formats: ["pdf", "svg"]
  plot_width: 35
  plot_height: 15
'
  
  config_path <- "config/config.yaml"
  writeLines(config_content, config_path)
  cat("Created configuration file:", config_path, "\n")
}

# =============================================================================
# 4. Create README file
# =============================================================================
create_readme_file <- function() {
  readme_content <- '
# Multi-omics Analysis Pipeline

## Overview
This project implements a comprehensive multi-omics analysis pipeline using MOFA (Multi-Omics Factor Analysis) and MultiGSEA for pathway enrichment analysis.

## Project Structure
```
├── scripts/           # Analysis scripts
├── data/             # Data files
├── results/          # Analysis results
├── functions/        # Custom functions
├── config/           # Configuration files
└── doc/              # Documentation
```

## Quick Start
1. Place your data files in the appropriate `data/` subdirectories
2. Update `config/config.yaml` with your parameters
3. Run the analysis pipeline:
   ```r
   source("scripts/multi-omics/run_complete_pipeline.R")
   ```

## Analysis Steps
1. **Data Preparation** (`01_MOFA_data_prep.R`)
2. **MOFA Training** (`02_MOFA_training.R`)
3. **Model Explanation** (`03_MOFA_explaination.R`)
4. **MultiGSEA Analysis** (`04_MultiGSEA.R`)
5. **Visualization** (`05_Visualization.R`)

## Requirements
- R >= 4.0
- Required packages: MOFA2, multiGSEA, tidyverse, pheatmap

## Documentation
- See `scripts/multi-omics/MOFA_MultiGSEA_Pipeline_Overview.md` for detailed usage
- See `scripts/multi-omics/OUTPUT_STRUCTURE_OPTIMIZATION.md` for file structure details
'
  
  writeLines(readme_content, "README.md")
  cat("Created README file: README.md\n")
}

# =============================================================================
# 5. Create .gitignore file
# =============================================================================
create_gitignore_file <- function() {
  gitignore_content <- '
# R specific
.Rhistory
.RData
.Ruserdata
*.Rproj

# Data files (large files)
data/raw/
data/processed/
*.csv
*.xlsx
*.hdf5
*.rds

# Results (can be regenerated)
results/
*.pdf
*.png
*.svg

# Temporary files
*.tmp
*.temp
.DS_Store

# Logs
*.log

# Test outputs
test_output/
'
  
  writeLines(gitignore_content, ".gitignore")
  cat("Created .gitignore file: .gitignore\n")
}

# =============================================================================
# 6. Create .Rprofile file
# =============================================================================
create_rprofile_file <- function() {
  rprofile_content <- '
# Multi-omics Analysis Project R Profile
# =====================================

# Set working directory to project root
if (file.exists("config/config.yaml")) {
  cat("Multi-omics Analysis Project loaded\n")
}

# Load common libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(MOFA2)
  library(multiGSEA)
})

# Set default options
options(stringsAsFactors = FALSE)
options(scipen = 999)
'
  
  writeLines(rprofile_content, ".Rprofile")
  cat("Created .Rprofile file: .Rprofile\n")
}

# =============================================================================
# 7. Main setup function
# =============================================================================
setup_project <- function() {
  cat("Setting up Multi-omics Analysis Project...\n")
  cat("==========================================\n\n")
  
  # Create directory structure
  create_project_structure(project_structure)
  
  # Create configuration files
  create_config_file()
  create_readme_file()
  create_gitignore_file()
  create_rprofile_file()
  
  cat("\n==========================================\n")
  cat("Project setup completed successfully!\n")
  cat("Next steps:\n")
  cat("1. Add your data files to the data/ directories\n")
  cat("2. Update config/config.yaml with your parameters\n")
  cat("3. Run the analysis pipeline\n")
  cat("4. Check scripts/multi-omics/ for detailed documentation\n")
}

# =============================================================================
# 8. Run setup
# =============================================================================
if (interactive()) {
  setup_project()
} else {
  cat("Run setup_project() to create the project structure\n")
} 