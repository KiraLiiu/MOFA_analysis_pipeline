# =============================================================================
# MOFA + MultiGSEA Pipeline - Step 05: Visualization
# =============================================================================
# This script generates publication-ready visualizations for MOFA and
# MultiGSEA results including factor plots, variance explained plots,
# and pathway enrichment heatmaps.
# 
# Author: Kira Liu
# Date: 2025-01-27
# 
# Input: MOFA model and MultiGSEA results from previous steps
# Output: Publication-ready figures in multiple formats (PDF, SVG, PNG)
# =============================================================================

library(MOFA2)
library(tidyverse)
library(openxlsx)
library(data.table)
library(ggplot2)

# =============================================================================
# SET OUTPUT DIRECTORIES
# =============================================================================
output_base <- "results/mofa_analysis"
visualization_dir <- file.path(output_base, "visualization")
mofa_figures_dir <- file.path(visualization_dir, "mofa_figures")
multigsea_figures_dir <- file.path(visualization_dir, "multigsea_figures")
heatmap_dir <- file.path(multigsea_figures_dir, "pathway_heatmaps")

dirs_to_create <- c(visualization_dir, mofa_figures_dir, multigsea_figures_dir, heatmap_dir)
for (dir in dirs_to_create) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

# =============================================================================
# 1. LOAD DATA
# =============================================================================
cat("Loading MOFA model and MultiGSEA results...\n")
MOFAobject <- load_model(file.path(output_base, "model_training/02_MOFA_Model_Trained.hdf5"))

# =============================================================================
# 2. GENERATE MOFA VISUALIZATIONS
# =============================================================================
cat("Generating MOFA visualizations...\n")

# Data overview plot
data_overview_plot <- plot_data_overview(MOFAobject)
ggsave(file.path(mofa_figures_dir, "05_MOFA_Data_Overview.svg"), plot = data_overview_plot, width = 10, height = 7)
ggsave(file.path(mofa_figures_dir, "05_MOFA_Data_Overview.png"), plot = data_overview_plot, width = 10, height = 7, dpi = 300)

# Factor correlation plot
pdf(file.path(mofa_figures_dir, "05_MOFA_Factor_Correlation.pdf"), width = 8, height = 6)
plot_factor_cor(MOFAobject)
dev.off()

svg(file.path(mofa_figures_dir, "05_MOFA_Factor_Correlation.svg"), width = 8, height = 6)
plot_factor_cor(MOFAobject)
dev.off()

# Variance explained plot
pdf(file.path(mofa_figures_dir, "05_MOFA_Variance_Explained.pdf"), width = 8, height = 6)
plot_variance_explained(MOFAobject, max_r2 = 15)
dev.off()

svg(file.path(mofa_figures_dir, "05_MOFA_Variance_Explained.svg"), width = 8, height = 6)
plot_variance_explained(MOFAobject, max_r2 = 15)
dev.off()

pdf(file.path(mofa_figures_dir, "05_MOFA_Variance_Explained_Total.pdf"), width = 8, height = 6)
plot_variance_explained(MOFAobject, plot_total = T)[[2]]
dev.off()

svg(file.path(mofa_figures_dir, "05_MOFA_Variance_Explained_Total.svg"), width = 8, height = 6)
plot_variance_explained(MOFAobject, plot_total = T)[[2]]
dev.off()

# Factor-covariate correlation plot
MOFA_metadata <- read.xlsx("data/metadata/sample_info.xlsx", sheet = 2)
samples_metadata(MOFAobject) <- MOFA_metadata

pdf(file.path(mofa_figures_dir, "05_MOFA_Factor_Covariate_Correlation.pdf"), width = 8, height = 6)
correlate_factors_with_covariates(MOFAobject, covariates = c("HALI", "Control"), factors = "all", plot = "r")
dev.off()

svg(file.path(mofa_figures_dir, "05_MOFA_Factor_Covariate_Correlation.svg"), width = 8, height = 6)
correlate_factors_with_covariates(MOFAobject, covariates = c("HALI", "Control"), factors = "all", plot = "r")
dev.off()

# =============================================================================
# 3. GENERATE MULTIGSEA HEATMAPS
# =============================================================================
cat("Generating MultiGSEA heatmaps...\n")

enrichment_results <- read.csv(file.path(output_base, "multigsea_analysis/tables/04_MultiGSEA_Complete_Results_selected.csv"))

# Source standalone heatmap functions
source("functions/Complex_Heatmap_for_MultiGSEA_Results.R")

# Generate all heatmaps (KEGG and Reactome in both PDF and SVG formats)
sources <- c("KEGG", "REACTOME")
formats <- c("pdf", "svg")

for (source in sources) {
  for (format in formats) {
    cat(paste("Generating", source, "heatmap in", format, "format...\n"))
    
    generate_pathway_heatmaps(
      data = enrichment_results,
      selected_sources = source,
      show_source_annotation = (source == "KEGG"),  # Only show annotation for KEGG
      save_output = TRUE,
      output_format = format,
      plot_width = 35,
      plot_height = 15,
      output_dir = heatmap_dir
    )
  }
}
cat("All heatmaps generated successfully!\n")

# =============================================================================
# 4. CREATE VISUALIZATION SUMMARY
# =============================================================================
cat("Creating visualization summary...\n")

# List all generated files
mofa_files <- list.files(mofa_figures_dir, full.names = TRUE)
multigsea_files <- list.files(multigsea_figures_dir, full.names = TRUE, recursive = TRUE)

visualization_summary <- data.frame(
  Category = c("MOFA_Figures", "MultiGSEA_Figures", "Total_Files"),
  Count = c(length(mofa_files), length(multigsea_files), length(mofa_files) + length(multigsea_files)),
  Location = c(mofa_figures_dir, multigsea_figures_dir, "Multiple directories")
)
write.csv(visualization_summary, file.path(visualization_dir, "05_Visualization_Summary.csv"), row.names = FALSE)

# Create detailed file list
file_details <- data.frame(
  File_Path = c(mofa_files, multigsea_files),
  File_Type = c(rep("MOFA", length(mofa_files)), rep("MultiGSEA", length(multigsea_files))),
  File_Size_KB = round(file.size(c(mofa_files, multigsea_files)) / 1024, 2)
)
write.csv(file_details, file.path(visualization_dir, "05_Generated_Files_List.csv"), row.names = FALSE)

# =============================================================================
# 5. CREATE ANALYSIS REPORT
# =============================================================================
cat("Creating analysis report...\n")

analysis_report <- data.frame(
  Analysis_Step = c("MOFA_Visualizations", "MultiGSEA_Heatmaps", "File_Organization", "Summary_Generation"),
  Status = c("Completed", "Completed", "Completed", "Completed"),
  Details = c(
    paste0("Generated ", length(mofa_files), " MOFA visualization files"),
    paste0("Generated ", length(multigsea_files), " MultiGSEA heatmap files"),
    "Organized files into structured directories",
    "Created visualization summary and file lists"
  )
)
write.csv(analysis_report, file.path(visualization_dir, "05_Visualization_Analysis_Report.csv"), row.names = FALSE)

cat("\nVisualization completed!\n")
cat(sprintf("MOFA figures saved in: %s\n", mofa_figures_dir))
cat(sprintf("MultiGSEA figures saved in: %s\n", multigsea_figures_dir))
cat(sprintf("Visualization summary saved in: %s\n", visualization_dir))
cat("Pipeline completed successfully!\n") 
