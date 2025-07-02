# MOFA + MultiGSEA Pipeline Overview

## Pipeline Description
This pipeline performs Multi-Omics Factor Analysis (MOFA) followed by Multi-omics Gene Set Enrichment Analysis (MultiGSEA) to identify shared biological patterns across transcriptomics, proteomics, and metabolomics data.

## Pipeline Steps

### Step 01: Data Preparation (`01_mofa_data_preparing.R`)
**Purpose**: Prepare and preprocess multi-omics data for MOFA analysis

**Key Functions**:
- Load normalized omics datasets (transcriptomics, proteomics, metabolomics)
- Filter transcriptomics data for significant genes (padj < 0.05)
- Perform quality checks and report data dimensions
- Apply high variance feature selection to reduce dimensionality
- Create MOFA data structure and initial MOFA object
- Generate data overview plots
- Save prepared data for subsequent steps

**Outputs**:
- `mofa_data_prepared.rds`: Prepared data list
- `transcriptomics_selected_features.csv`: Filtered transcriptomics data
- `proteomics_selected_features.csv`: Filtered proteomics data
- `metabolomics_selected_features.csv`: Filtered metabolomics data
- `MOFA_object_created.rds`: Initial MOFA object
- `data_overview_plot.pdf`: Data structure visualization

### Step 02: MOFA Training (`02_MOFA_training.R`)
**Purpose**: Train the MOFA model to identify shared factors across omics layers

**Key Functions**:
- Load prepared data from Step 01
- Configure MOFA parameters (number of factors, convergence mode, etc.)
- Train the MOFA model
- Save trained model as HDF5 and RDS files

**Parameters**:
- `num_factors`: Number of factors to extract (default: 3-5)
- `convergence_mode`: "fast", "medium", or "slow"
- `scale_views`: Whether to scale views to unit variance
- `seed`: Random seed for reproducibility

**Outputs**:
- `MOFA_multi-omics.hdf5`: Trained MOFA model (HDF5 format)
- `MOFA_multi-omics.rds`: Trained MOFA model (R format)

### Step 03: MOFA Explanation (`03_MOFA_explaination.R`)
**Purpose**: Analyze and interpret the trained MOFA model

**Key Functions**:
- Load trained MOFA model
- Generate factor correlation plots
- Calculate and visualize variance explained by factors
- Correlate factors with sample metadata/covariates
- Extract factor weights and feature importance
- Generate factor plots and weight visualizations
- Extract top features for each factor

**Outputs**:
- Factor correlation plots
- Variance explained plots
- Factor-covariate correlation plots
- Top feature lists for each factor
- `Factor1_all_features.xlsx`, `Factor1_top_features.xlsx`, etc.

### Step 04: MultiGSEA Analysis (`04_MultiGSEA.R`)
**Purpose**: Perform pathway enrichment analysis on MOFA factors using multi-omics data

**Key Functions**:
- Load MOFA factor results and original omics data
- Prepare data for MultiGSEA (merge features with statistical results)
- Initialize multi-omics data structure
- Load pathway databases (KEGG, Reactome)
- Perform multi-omics GSEA analysis
- Calculate combined p-values across omics layers
- Generate enrichment results

**Databases Used**:
- KEGG pathways
- Reactome pathways
- Custom gene sets (optional)

**Outputs**:
- `MOFA_Factor1_multiGSEA_results.csv`: MultiGSEA enrichment results
- `MOFA_Factor1_multiGSEA_results.xlsx`: MultiGSEA enrichment results (Excel)

### Step 05: Visualization (`05_Visualization.R`)
**Purpose**: Create publication-ready visualizations of MOFA and MultiGSEA results

**Key Functions**:
- Generate MOFA result visualizations (factor plots, variance explained, etc.)
- Create complex heatmaps for MultiGSEA results
- Customize plot aesthetics and save in multiple formats
- Generate pathway-specific heatmaps

**Output Formats**:
- SVG (vector graphics for publications)
- PDF (high-quality print format)
- PNG (raster format for presentations)

## Required Functions

### MOFA Analysis Functions (`functions/MOFA_analysis_function.R`)
- `run_mofa_analysis()`: Main MOFA analysis function
- `select_high_variance_features()`: Feature selection based on variance
- `extract_top_features()`: Extract top features for each factor
- `reorder_and_fill_na()`: Data preprocessing utility

### MultiGSEA Visualization Functions (`functions/Complex_Heatmap_for_MultiGSEA_Results.R`)
- `generate_pathway_heatmaps()`: Main heatmap generation function
- `create_selected_sources_heatmap()`: Create heatmaps for selected pathway sources
- `create_source_separated_heatmap()`: Create individual source heatmaps
- `prepare_heatmap_data()`: Data preparation for heatmaps

## Data Requirements

### Input Data Structure
```
data/
├── processed/
│   ├── metabolomics_normalized/
│   │   └── metabolomics_normalized.csv
│   ├── proteomics_normalized/
│   │   └── proteomics_normalized.csv
│   └── RNAseq_DESeq2_results/
│       └── Results_without_annotation/
│           ├── vst_normalized_counts.csv
│           └── DESeq2_HALI vs Control_results.csv
├── raw/
│   ├── Metabolomics_Clean.xlsx
│   └── Protein_DEPs.xlsx
└── metadata/
    └── sample_info.xlsx
```

### Output Directory Structure
```
results/
├── mofa_analysis/
│   ├── MOFA_multi-omics.hdf5
│   ├── MOFA_multi-omics.rds
│   └── data_overview_plot.pdf
├── MOFA_multiomics_results/
│   ├── Factor1_all_features.xlsx
│   ├── Factor1_top_features.xlsx
│   └── ...
├── tables/
│   ├── MOFA_Factor1_multiGSEA_results.csv
│   └── MOFA_Factor1_multiGSEA_results.xlsx
└── figures/
    ├── mofa_data_overview.svg
    ├── mofa_factor_correlation.svg
    ├── mofa_variance_explained.svg
    └── MOFA_Factor1_MultiGSEA_heatmaps/
        ├── kegg_reactome_pathways_heatmap.svg
        └── ...
```

## Key Parameters

### MOFA Parameters
- **Number of Factors**: 3-5 (adjust based on data complexity)
- **Convergence Mode**: "slow" for final analysis, "fast" for exploration
- **Feature Selection**: Variance threshold = 1.0 (keep all features)
- **Scaling**: scale_views = FALSE (data already normalized)

### MultiGSEA Parameters
- **Pathway Databases**: KEGG, Reactome
- **Organism**: mmusculus (mouse)
- **P-value Combination**: Fisher's method
- **Multiple Testing Correction**: Benjamini-Hochberg (BH)

## Usage Instructions

1. **Prepare Environment**: Ensure all required R packages are installed
2. **Run Step 01**: Execute `01_mofa_data_preparing.R` to prepare data
3. **Run Step 02**: Execute `02_MOFA_training.R` to train MOFA model
4. **Run Step 03**: Execute `03_MOFA_explaination.R` to analyze results
5. **Run Step 04**: Execute `04_MultiGSEA.R` to perform pathway enrichment
6. **Run Step 05**: Execute `05_Visualization.R` to generate plots

## Troubleshooting

### Common Issues
1. **Memory Issues**: Reduce number of features or use faster convergence mode
2. **Convergence Problems**: Increase max iterations or adjust convergence mode
3. **Missing Data**: Check data quality in Step 01
4. **Pathway Mapping Issues**: Verify gene/metabolite annotations

### Performance Tips
- Use "fast" convergence mode for initial exploration
- Reduce feature set size for faster computation
- Use GPU mode if available (requires cupy installation)
- Monitor memory usage during training

## Citation
When using this pipeline, please cite:
- MOFA2: Argelaguet, R. et al. (2020). Multi-Omics Factor Analysis—a framework for unsupervised integration of multi-omics data sets. Molecular Systems Biology, 16(6), e9192.
- MultiGSEA: Please cite the multiGSEA package appropriately 