# Mapping Immune Resistance in Melanoma (Updated Analysis)

This repository contains a comprehensive suite of R scripts for analyzing immune resistance mechanisms in melanoma across multiple modalities: Single-Cell RNA-seq (scRNA-seq), bulk genomics (TCGA, MSK-IMPACT), CRISPR dependency screens (DepMap), and clinical trial data.

The codebase reproduces key molecular patterns associated with immune checkpoint blockade (ICB) resistance and provides tools for comparative analysis.

## Requirements

The analysis pipeline relies on standard Bioconductor and CRAN packages. Dependencies are managed via `Code/Utilities.R`.

**Key Libraries:** `SingleCellExperiment`, `scran`, `scater`, `ComplexHeatmap`, `DESeq2`, `limma`, `ggplot2`, `dplyr`, `pheatmap`.

## Code Structure and Functionality

The `Code/` directory contains the following primary analysis scripts. These are designed to be run individually as needed.

### 1. Single-Cell Analysis
*   **[Molecular_pattern.R](Code/Molecular_pattern.R)**: The core script for reproducing the molecular patterns of resistance.
    *   Generates OncoPrints of genomic alterations (Mutation, CNA) and integrates them with clinical metadata.
    *   Performs single-cell analysis: Quality Control (QC), Dimension Reduction (PCA, UMAP), Clustering, and Cell Type Annotation.
    *   Conducts Differential Expression Analysis (DEA) to identify resistance-associated gene signatures (e.g., PD1-Resistant vs Naive).
    *   Visualizes results using Volcano plots and Heatmaps.
*   **[Immuno_resistance.R](Code/Immuno_resistance.R)**: A dedicated pipeline for in-depth immune resistance analysis.
    *   Focuses on detailed scRNA-seq processing: normalization, feature selection, and manifold learning.
    *   Includes sub-clustering of malignant cell populations to identify intra-tumoral heterogeneity.
    *   Computes and visualizes marker genes for specific cell states.

### 2. Bulk Genomics & Multi-Omics
*   **[TCGA_melanoma.R](Code/TCGA_melanoma.R)**: Comprehensive analysis of The Cancer Genome Atlas (TCGA) Skin Cutaneous Melanoma (SKCM) data.
    *   **Genomics**: Mutation and Copy Number Variation (CNV) landscape.
    *   **Transcriptomics**: RNA-seq differential expression analysis (e.g., M0 vs M1 metastatic stages).
    *   **Proteomics**: RPPA protein expression profiling.
    *   **Interaction**: Mutation co-occurrence analysis and OncoPrint visualization.
*   **[MSK_melanoma.R](Code/MSK_melanoma.R)**: Targeted sequencing analysis using MSK-IMPACT data.
    *   Focuses on clinically relevant driver mutations and their frequencies.
    *   Analyzes Tumor Mutational Burden (TMB) and Overall Survival (OS).
    *   Generates OncoPrints for top altered genes.

### 3. Functional Screens & Clinical Data
*   **[Gene Dependency.R](Code/Gene%20Dependency.R)**: Analysis of CRISPR-Cas9 gene dependency screens from the DepMap project.
    *   Identifies essential genes in melanoma cell lines.
    *   Compares gene dependency scores between skin cancer lines and other lineages.
    *   Highlights potential therapeutic targets.
*   **[Clinical_trial.R](Code/Clinical_trial.R)**: Analytics for melanoma clinical trials.
    *   Mines clinical trial data to identify common drug interventions and combinations.
    *   Visualizes trial landscapes (monotherapy vs combination therapy).

### 4. Utilities
*   **[Utilities.R](Code/Utilities.R)**: Helper script for installing and loading required R packages.

## Data Requirement

The scripts expect a `Data/` directory in the project root containing the necessary datasets:
*   `Data/scData/`: Single-cell RNA-seq objects (e.g., `Mel.all.data.QC.rds`).
*   `Data/molecular_pattern/`: Genomic alteration data (`Alteration_data.xlsx`) and raw scRNA-seq counts.
*   `Data/skcm_tcga_pan_can_atlas_2018/`: TCGA datasets.
*   `Data/mel_mskimpact_2020/`: MSK-IMPACT datasets.
*   `Data/depmap/`: DepMap CRISPR gene effect and model metadata.
*   `Data/clinical_trials/`: Clinical trial records.

## Usage

To run an analysis, open the specific R script in RStudio or run via command line:

```R
# Example: Run the core molecular pattern analysis
source("Code/Molecular_pattern.R")
```

Outputs (Figures and Tables) are saved to the `Images/` and `Data/Results/` directories respectively.
