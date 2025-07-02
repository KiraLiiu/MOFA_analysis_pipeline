# =============================================================================
# MOFA + MultiGSEA Pipeline - Step 04: MultiGSEA Analysis
# =============================================================================
# This script performs pathway enrichment analysis on MOFA factors using
# multi-omics data to identify biological pathways enriched across
# transcriptomics, proteomics, and metabolomics layers.
# 
# Author: Kira Liu
# Date: 2025-01-27
# 
# Input: MOFA factor features from Step 03 and omics statistical results
# Output: Pathway enrichment results and summary statistics
# Note: Visualization plots are generated in Step 05
# =============================================================================

library(multiGSEA)
library(tidyverse)
library(openxlsx)
library(data.table)

# =============================================================================
# Set output directories
# =============================================================================
output_base <- "results/mofa_analysis"
multigsea_dir <- file.path(output_base, "multigsea_analysis")
tables_dir <- file.path(multigsea_dir, "tables")

for (dir in c(multigsea_dir, tables_dir)) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# =============================================================================
# 1. Load top features for Factor 1
# =============================================================================
features_factor1 <- read.xlsx(file.path(output_base, "model_explanation/tables/Factor1_all_features.xlsx"))
# features_factor2 <- read.xlsx(file.path(output_base, "model_explanation/tables/Factor2_all_features.xlsx"))

# =============================================================================
# 2. Prepare omics data for MultiGSEA
# =============================================================================
# All plots are generated in 05_Visualization.R. Only data for downstream visualization is prepared here.
rna_features <- features_factor1 %>% filter(view == "rna")
protein_features <- features_factor1 %>% filter(view == "protein")
metabolite_features <- features_factor1 %>% filter(view == "metabolite")


rna_results <- read.csv("data/processed/transcriptomics/final_results_remove_duplicates/HALI vs Control_final_cleaned.csv", check.names = FALSE)
protein_results <- read.xlsx("data/raw/proteomics/Protein_DEPs.xlsx", sheet = 1)
metabolite_results <- read.xlsx("data/raw/metabolomics/Metabolomics_Clean.xlsx", sheet = 1)
metabolite_annotation <- read.csv("data/raw/metabolomics/name_map.csv")

# Merge and prepare for MultiGSEA
df_transcriptome <- rna_features %>%
  left_join(rna_results, by = c("feature" = "gene_id")) %>%
  dplyr::select(c(feature, value, pvalue, padj)) %>%
  drop_na()

any(duplicated(df_transcriptome$feature))

df_proteome <- protein_features %>%
  left_join(protein_results, by = c("feature" = "Gene.Name")) %>%
  dplyr::select(c(feature, value, pvalue)) %>%
  mutate(padj = pvalue) %>%
  drop_na()

any(duplicated(df_proteome$feature))

df_metabolome <- metabolite_features %>%
  left_join(metabolite_results, by = c("feature" = "Metabolites.Name")) %>%
  left_join(metabolite_annotation, by = c("feature" = "Query")) %>%
  dplyr::select(c(feature, HMDB, value, p.value)) %>%
  dplyr::rename(pvalue = p.value) %>%
  mutate(padj = pvalue) %>%
  drop_na()

any(duplicated(df_metabolome$feature))
any(duplicated(df_metabolome$HMDB))

df_metabolome <- df_metabolome %>%
  group_by(HMDB) %>%
  summarise(
    features_combined = paste(unique(feature), collapse = "; "),
    mean_value = mean(value, na.rm = TRUE),
    min_pvalue = min(pvalue, na.rm = TRUE),
    min_padj = min(padj, na.rm = TRUE)
  ) %>%
  ungroup()

any(duplicated(df_metabolome$HMDB))

write.csv(df_transcriptome, file.path(tables_dir, "04_Transcriptome_Data_for_MultiGSEA.csv"), row.names = FALSE)
write.csv(df_proteome, file.path(tables_dir, "04_Proteome_Data_for_MultiGSEA.csv"), row.names = FALSE)
write.csv(df_metabolome, file.path(tables_dir, "04_Metabolome_Data_for_MultiGSEA.csv"), row.names = FALSE)

# =============================================================================
# 3. Run MultiGSEA analysis
# =============================================================================
omics_data <- initOmicsDataStructure(layer = c("transcriptome", "proteome", "metabolome"))
omics_data$transcriptome <- df_transcriptome$value
names(omics_data$transcriptome) <- df_transcriptome$feature
omics_data$proteome <- df_proteome$value
names(omics_data$proteome) <- df_proteome$feature
omics_data$metabolome <- df_metabolome$mean_value
names(omics_data$metabolome) <- df_metabolome$HMDB

pathways <- getMultiOmicsFeatures(
  dbs = c("kegg", "reactome"),
  layer = names(omics_data),
  returnTranscriptome = "ENSEMBL",
  returnProteome = "SYMBOL",
  returnMetabolome = "HMDB",
  organism = "mmusculus",
  useLocal = TRUE
)

enrichment_scores <- multiGSEA(pathways, omics_data)
df_enrichment_results <- extractPvalues(enrichmentScores = enrichment_scores, pathwayNames = names(pathways[[1]]))
df_enrichment_results$combined_pval <- combinePvalues(df_enrichment_results)
df_enrichment_results$combined_padj <- p.adjust(df_enrichment_results$combined_pval, method = "BH")
df_enrichment_results <- cbind(data.frame(pathway = names(pathways[[1]])), df_enrichment_results)
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

# =============================================================================
# 4. Save enrichment results
# =============================================================================
write.csv(enrichment_results_cleaned, file.path(tables_dir, "04_MultiGSEA_Complete_Results.csv"), row.names = FALSE)
write.xlsx(enrichment_results_cleaned, file.path(tables_dir, "04_MultiGSEA_Complete_Results.xlsx"))

significant_results <- enrichment_results_cleaned %>%
  filter(combined_padj < 0.05) %>%
  arrange(combined_padj)
write.csv(significant_results, file.path(tables_dir, "04_MultiGSEA_Significant_Results.csv"), row.names = FALSE)
write.xlsx(significant_results, file.path(tables_dir, "04_MultiGSEA_Significant_Results.xlsx"))

kegg_results <- enrichment_results_cleaned %>%
  filter(source == "KEGG") %>%
  arrange(combined_padj) %>%
  head(20)
write.csv(kegg_results, file.path(tables_dir, "04_MultiGSEA_Top_KEGG_Pathways.csv"), row.names = FALSE)

reactome_results <- enrichment_results_cleaned %>%
  filter(source == "Reactome") %>%
  arrange(combined_padj) %>%
  head(20)
write.csv(reactome_results, file.path(tables_dir, "04_MultiGSEA_Top_Reactome_Pathways.csv"), row.names = FALSE)

# =============================================================================
# 5. Create summary statistics
# =============================================================================
summary_stats <- data.frame(
  Metric = c("Total_Pathways_Analyzed", "Significant_Pathways_padj_0.05", 
             "KEGG_Pathways", "Reactome_Pathways", "Transcriptome_Features", 
             "Proteome_Features", "Metabolome_Features"),
  Value = c(
    nrow(enrichment_results_cleaned),
    nrow(significant_results),
    sum(enrichment_results_cleaned$source == "KEGG", na.rm = TRUE),
    sum(enrichment_results_cleaned$source == "Reactome", na.rm = TRUE),
    length(omics_data$transcriptome),
    length(omics_data$proteome),
    length(omics_data$metabolome)
  )
)
write.csv(summary_stats, file.path(tables_dir, "04_MultiGSEA_Summary_Statistics.csv"), row.names = FALSE)

# =============================================================================
# 6. Create analysis report
# =============================================================================
analysis_report <- data.frame(
  Analysis_Step = c("Feature_Loading", "Data_Preparation", "MultiGSEA_Analysis", 
                    "Results_Processing", "Output_Generation"),
  Status = rep("Completed", 5),
  Details = c(
    paste0("Loaded ", nrow(features_factor1), " features for Factor 1"),
    "Prepared transcriptome, proteome, and metabolome data",
    "Performed pathway enrichment analysis using KEGG and Reactome",
    "Processed and filtered results",
    "Generated complete and significant results files"
  )
)
write.csv(analysis_report, file.path(multigsea_dir, "04_MultiGSEA_Analysis_Report.csv"), row.names = FALSE)

message("MultiGSEA analysis completed. Complete results saved in: ", tables_dir)
message("Significant results saved in: ", tables_dir)
message("Analysis report saved in: ", multigsea_dir)
message("All visualization is handled in 05_Visualization.R.")
message("Next step: Run 05_Visualization.R") 

# =============================================================================
# clean up
rm(list = ls())
gc()
