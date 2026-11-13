#' Title: Single-Cell RNA-Seq Data Loading and Initial Merging Pipeline
#' Description: Script for loading 10X Genomics matrices, creating Seurat objects, 
#'              merging replicates, and calculating fundamental QC metrics.
#' Dataset Accession: GEO GSE342902

# 1. Load Required Libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(patchwork)
  library(RColorBrewer)
})

# 2. Set Working Directory and Input Paths
# Note: Users should adjust the base directory according to their local environment.
base_dir <- "/analysis/yyy"
setwd(base_dir)

# Define sample identifiers
sample_ids <- c("con_4", "con_5", "mol_7", "mol_8")

# 3. Automated Seurat Object Creation and Initial Filtering
# Parameters: min.cells = 3, min.features = 200 (standard thresholds)
seurat_list <- lapply(sample_ids, function(sample_id) {
  data_path <- file.path(base_dir, sample_id, "filtered_cell_gene_matrix")
  
  # Load 10X data
  counts <- Read10X(data.dir = data_path)
  
  # Create Seurat object
  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 3,       # Filter features expressed in fewer than 3 cells
    min.features = 200   # Filter low-quality cells with < 200 features
  )
  
  # Save individual processed object
  saveRDS(obj, file = paste0(sample_id, "_seurat.rds"))
  return(obj)
})

# Assign names to the list for easy indexing
names(seurat_list) <- sample_ids

# 4. Merge All Replicates into a Single Object
# add.cell.ids prepends sample prefixes to cell barcodes to prevent duplication
alldata <- merge(
  x = seurat_list[["con_4"]],
  y = list(seurat_list[["con_5"]], seurat_list[["mol_7"]], seurat_list[["mol_8"]]),
  add.cell.ids = sample_ids,
  project = "MPP_merged_atlas"
)

# 5. Metadata Annotation and QC Metric Calculation
# Define experimental conditions (Control vs. Disease/Model)
alldata$group <- ifelse(grepl("^con", alldata$orig.ident), "Control", "Model")

# Calculate mitochondrial gene percentage (using lowercase "^mt-" for mouse data conformity)
alldata[["percent.mt"]] <- PercentageFeatureSet(alldata, pattern = "^mt-")

# 6. Summary Statistics Check
cat("=== Cell Count per Sample ===\n")
print(table(alldata$orig.ident))

cat("\n=== Cell Count per Group ===\n")
print(table(alldata$group))

# 7. Save Final Merged Raw Object
output_file <- "alldata_merged.rds"
saveRDS(alldata, file = output_file)
cat(paste0("\n[SUCCESS] Merged Seurat object successfully saved to: ", file.path(base_dir, output_file), "\n"))

