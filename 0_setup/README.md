# Setup Directory

This directory contains scripts for setting up and initializing the multi-omics analysis project.

## Files

- `01_project_setup.R` - Main project setup script that creates directory structure and configuration files

## Usage

To set up the project structure and configuration:

```r
# Run the setup script
source("setup/01_project_setup.R")
setup_project()
```

## What it does

1. Creates the complete directory structure for the project
2. Generates configuration files (`config/config.yaml`)
3. Creates initial README and documentation files
4. Sets up `.gitignore` and `.Rprofile` files
5. Initializes all necessary directories for data, results, and analysis

## After Setup

Once the setup is complete:

1. Place your raw data files in the appropriate `data/raw/` subdirectories
2. Update `config/config.yaml` with your specific parameters
3. Run the analysis pipeline using the main script 