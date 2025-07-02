# MOFA + MultiGSEA Analysis Pipeline (Standalone)

## Overview
A minimal, standalone implementation of Multi-Omics Factor Analysis (MOFA) and MultiGSEA pathway enrichment analysis for integrating transcriptomics, proteomics, and metabolomics data.

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
1. **Setup**
   ```r
   source("0_setup/01_project_setup.R")
   setup_project()
   ```
2. **Configure**
   Edit `config/config.yaml` with your parameters.
3. **Add Data**
   Place your preprocessed data (all omics) in `data/processed/`:
   - `data/processed/metabolomics.csv`
   - `data/processed/proteomics.csv`
   - `data/processed/transcriptomics.csv`
   - `data/metadata/sample_info.xlsx`
4. **Run Analysis**
   ```r
   source("scripts/run_complete_pipeline.R")
   run_pipeline()
   ```

## Output
- All results are saved in the `results/` directory (no subfolders by default, just a `.gitkeep` file).

## .gitignore Policy
- All raw data和大体量结果默认不上传，保留目录结构（通过 .gitkeep）。
- 可上传小型示例数据和分析脚本、配置文件。
- 详见 `.gitignore` 文件注释。

## Requirements
- R >= 4.0
- Packages: MOFA2, multiGSEA, tidyverse, yaml, here

## Citation
- MOFA2: Argelaguet, R. et al. (2020). Multi-Omics Factor Analysis. Molecular Systems Biology, 16(6), e9192.
- MultiGSEA: Cite the multiGSEA package appropriately

---
*This is a minimal standalone version of the MOFA analysis pipeline. For a more comprehensive multi-omics project, see the main repository.*
 