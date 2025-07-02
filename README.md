# MOFA + MultiGSEA Analysis Pipeline

## Overview
An implementation of Multi-Omics Factor Analysis (MOFA) and MultiGSEA pathway enrichment analysis for integrating transcriptomics, proteomics, and metabolomics data.

## Directory Structure (Minimal)
```
├── config/                  # Configuration files (config.yaml)
├── data/
│   ├── raw/                # Raw data (optional, not uploaded)
│   ├── processed/          # All normalized data (all omics together, e.g. metabolomics.csv, proteomics.csv, transcriptomics.csv)
│   └── metadata/           # Sample metadata (e.g. sample_info.xlsx)
├── results/                # All output results (no subfolders by default)
├── scripts/                # Analysis scripts
├── functions/              # Custom R functions
├── 0_setup/                # Project setup scripts
├── .gitignore
├── .gitkeep
└── README.md
```

## Quick Start
1. **Configure**
   Edit `config/config.yaml` with your parameters.
2. **Add Data**
   Place your preprocessed data (all omics) in `data/processed/`:
   - `data/processed/metabolomics.csv`
   - `data/processed/proteomics.csv`
   - `data/processed/transcriptomics.csv`
   - `data/metadata/sample_info.xlsx`
3. **Run Analysis**
   ```r
   source("scripts/run_complete_pipeline.R")
   run_pipeline()
   ```

## Output
- All results are saved in the `results/` directory (no subfolders by default, just a `.gitkeep` file).

## .gitignore Policy
- All raw data and large results are not uploaded by default; only the directory structure is kept (using `.gitkeep`).
- Example data, scripts, and configuration files can be uploaded.
- See the `.gitignore` file for details.

## Requirements
- R >= 4.0
- Packages: MOFA2, multiGSEA, tidyverse, yaml, here

## Citation
- MOFA: Multi‐Omics Factor Analysis—a framework for unsupervised integration of multi‐omics data sets Mol Syst Biol. (2018) 14: e8124 https://doi.org/10.15252/msb.20178124
- MultiGSEA: Canzler, S., Hackermüller, J. multiGSEA: a GSEA-based pathway enrichment analysis for multi-omics data. BMC Bioinformatics 21, 561 (2020). https://doi.org/10.1186/s12859-020-03910-x

---
*This is a minimal version of the MOFA analysis pipeline. For a more comprehensive multi-omics project, see the main repository.*
 