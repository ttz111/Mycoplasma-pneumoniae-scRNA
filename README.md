# Mycoplasma pneumoniae Pneumonia scRNA-Seq Analysis Pipeline

This repository contains the publication-grade bioinformatics analysis code and processed Seurat objects for the single-cell RNA-sequencing study evaluating lung immune cell remodeling in a murine *Mycoplasma pneumoniae* pneumonia model.

## Repository Structure & File Descriptions

- **`Complete_QC_and_Harmony_Pipeline.R`**: The comprehensive main analysis script covering:
  - Adaptive quality control using **miQC** and mitochondrial filtering.
  - Standard Seurat preprocessing workflow (`NormalizeData`, `FindVariableFeatures`, `ScaleData`).
  - Principal Component Analysis (PCA) and Elbow plot determination.
  - Batch effect correction across biological replicates using **Harmony**.
  - Graph-based cell clustering (`FindNeighbors`, `FindClusters`) and non-linear dimensionality reduction (**UMAP** & **t-SNE**).
  - Automated generation and export of publication-quality plots and summary metrics tables.
- **`simple_analysis.R`**: A streamlined, lightweight preliminary script for quick exploration and baseline visualization.
- **`alldata_merged.rds`**: The merged Seurat object containing raw or preprocessed multi-sample data prior to final cell-type subsetting.
- **`alldata_6_major_cell_types.rds`**: The fully processed and annotated Seurat object containing 6 major cell lineages identified in the study.
- **`con_4_seurat.rds` & `con_5_seurat.rds`**: Individual processed Seurat objects for control biological replicates.
- **`mol_7_seurat.rds` & `mol_8_seurat.rds`**: Individual processed Seurat objects for *M. pneumoniae*-infected treatment replicates.
- **`.gitattributes`**: Git Large File Storage (LFS) configuration tracking binary `.rds` files.

## System Requirements & Dependencies

The analyses were performed under an Ubuntu Linux environment using **R (version 4.x)**. The pipeline relies on the following core Bioconductor and CRAN packages:
* `Seurat`
* `harmony`
* `miQC`
* `flexmix`
* `dplyr`
* `ggplot2`
* `patchwork`
* `scales`
* `RColorBrewer`

## Data Availability

The raw sequencing reads and processed count matrices associated with this study have been deposited in the **NCBI Gene Expression Omnibus (GEO)** under accession number **GSE342902**.

## Usage

1. Clone the repository and ensure Git LFS is installed to correctly download the `.rds` objects:
   ```bash
   git lfs pull
