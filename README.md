# miic_helper

A toolkit for preprocessing single-cell genomics data and automating MIIC (Multivariate Information-based Inductive Causation) causal network analysis workflows.

## Purpose

This repository bridges single-cell RNA-seq analysis tools (Seurat/AnnData) with MIIC causal inference by providing:

- **Smart Feature Selection**: Mutual information-based gene selection to identify relevant features for causal analysis
- **Format Conversion**:  conversion between Seurat objects, AnnData objects, and MIIC-compatible formats
- **Automated Workflows**: End-to-end pipelines from raw single-cell data to causal network inference
- **Web Integration**: Automated submission to the MIIC web server at miic.curie.fr

## Core Functionality

### Data Processing & Feature Selection
- Filter Seurat objects by cell type, timepoint, and custom metadata
- Compute mutual information scores between genes and variables of interest
- Smart gene selection using MI scores with grouping and ranking strategies
- Handle both discrete (categorical) and continuous variables

### Multi-Format Support
- **Seurat Objects**: Native R single-cell analysis integration
- **AnnData Objects**: Python scanpy/anndata compatibility  
- **CSV Export**: Generate MIIC-ready CSV matrices and state files

### MIIC Integration
- Generate properly formatted input files (.csv + .st.txt) for MIIC analysis
- Define variable types, state orders, contextual variables, and causal consequences
- Automated web server submission with result tracking
- Local MIIC execution support

## Key Components

### R Functions (`src/preprocess_sureat_object.R`)
- `Filter_seurat_object()` - Apply metadata-based filters
- `wrap_MI_compute()` - Compute mutual information scores  
- `handle_MI_Selection()` - Smart feature selection with visualization
- `Generate_miic_files()` - Export MIIC-compatible datasets

### Python Modules
- `process_anndata_object.py` - AnnData creation, manipulation, and I/O
- `compute_MI_wrapper.py` - Python interface to R-based MI computation
- `miic_web_automation.py` - Selenium-based web server automation

### Command Line Tools
- `run_miic_select.R` - Standalone MI feature selection script

## Quick Start

See `tuto/Automatic_preprocess_template.Rmd` for a complete workflow example.

## Author

Ali Chemkhi (ali.chemkhi@curie.fr)