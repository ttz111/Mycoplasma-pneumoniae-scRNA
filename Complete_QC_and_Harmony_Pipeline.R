#' Title: Complete Single-Cell RNA-Seq Analysis Pipeline (miQC, Harmony, and Clustering)
#' Description: End-to-end R script covering adaptive quality control (miQC), 
#'              standard Seurat workflow, Harmony batch integration, dimensionality reduction, 
#'              clustering, and comprehensive metrics summarization.
#' Dataset Accession: GEO GSE342902

# 1. Load Required Libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(miQC)
  library(flexmix)
  library(scales)
  library(harmony)
  library(RColorBrewer)
})

# 2. Load Merged Seurat Object
# Ensure "alldata_merged.rds" is present in your working directory
alldata <- readRDS("alldata_merged.rds")

# ==========================================
# 3. Mitochondrial Percentage & miQC Adaptive QC
# ==========================================
alldata[["percent.mt"]] <- PercentageFeatureSet(alldata, pattern = "^mt-")

# Extract metadata for control and treatment groups separately
df_ctrl <- FetchData(alldata, vars = c("nFeature_RNA", "nCount_RNA", "percent.mt", "orig.ident")) %>%
  filter(grepl("^con", orig.ident))
df_ctrl$cell_id <- rownames(df_ctrl)

df_treat <- FetchData(alldata, vars = c("nFeature_RNA", "nCount_RNA", "percent.mt", "orig.ident")) %>%
  filter(grepl("^mol", orig.ident))
df_treat$cell_id <- rownames(df_treat)

# Robust miQC fitting function with error protection
fit_miqc_safe <- function(df) {
  set.seed(123)
  model <- flexmix(percent.mt ~ nFeature_RNA, data = df, k = 2)
  params <- parameters(model)
  post <- posterior(model)
  
  if (ncol(post) < 2) {
    df$prob_compromised <- 0 
  } else {
    mean_mt_1 <- sum(df$percent.mt * post[, 1]) / sum(post[, 1])
    mean_mt_2 <- sum(df$percent.mt * post[, 2]) / sum(post[, 2])
    comp_idx <- if (mean_mt_1 > mean_mt_2) 1 else 2
    df$prob_compromised <- post[, comp_idx]
  }
  return(list(df = df, params = params))
}

fit_c <- fit_miqc_safe(df_ctrl)
df_ctrl <- fit_c$df
params_ctrl <- fit_c$params

fit_t <- fit_miqc_safe(df_treat)
df_treat <- fit_t$df
params_treat <- fit_t$params

# Define Selection Thresholds and Assign QC Status
f_low <- 500
f_high <- 10000

assign_qc <- function(df) {
  df %>% mutate(
    Selection = case_when(
      nFeature_RNA < f_low & prob_compromised >= 0.5 ~ "Features & Mitochondrial Filter",
      nFeature_RNA < f_low ~ "Features Filter",
      nFeature_RNA > f_high ~ "Features Filter",
      prob_compromised >= 0.5 ~ "Mitochondrial Filter",
      TRUE ~ "Pass"
    )
  )
}

df_ctrl <- assign_qc(df_ctrl)
df_treat <- assign_qc(df_treat)

df_qc_all <- rbind(df_ctrl, df_treat)
alldata$Selection <- df_qc_all$Selection[match(colnames(alldata), df_qc_all$cell_id)]
alldata$prob_compromised <- df_qc_all$prob_compromised[match(colnames(alldata), df_qc_all$cell_id)]

cat("--- QC Selection Summary Table Across Samples ---\n")
print(table(alldata$orig.ident, alldata$Selection))

# ==========================================
# 4. Cell Filtering & Standard Seurat Workflow
# ==========================================
alldata_filtered <- subset(alldata, subset = Selection == "Pass")

cat(sprintf("\nTotal cells before filtering: %d\n", ncol(alldata)))
cat(sprintf("Total high-quality cells retained ('Pass'): %d\n", ncol(alldata_filtered)))

# Normalization & Variable Features
alldata_filtered <- NormalizeData(alldata_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
alldata_filtered <- FindVariableFeatures(alldata_filtered, selection.method = "vst", nfeatures = 2000)

# Plot and save labeled variable features
top10 <- head(VariableFeatures(alldata_filtered), 10)
p_vfg <- VariableFeaturePlot(alldata_filtered, cols = c("black", "red"))
p_vfg_labeled <- LabelPoints(plot = p_vfg, points = top10, repel = TRUE)
ggsave("Variable_Feature_Plot_Labeled.png", plot = p_vfg_labeled, width = 7, height = 5, dpi = 300)

# Scale data (regressing out percent.mt)
all.genes <- rownames(alldata_filtered)
alldata_filtered <- ScaleData(alldata_filtered, features = all.genes, vars.to.regress = c("percent.mt"))

# PCA & Elbow Plot
alldata_filtered <- RunPCA(alldata_filtered, npcs = 30, verbose = FALSE)
p_elbow <- ElbowPlot(alldata_filtered, ndims = 30)
ggsave("PCA_ElbowPlot.png", plot = p_elbow, width = 6, height = 5, dpi = 300)

pcs_use <- 1:20

# ==========================================
# 5. Batch Effect Correction (Harmony) & Clustering
# ==========================================
pca_embeddings <- Embeddings(alldata_filtered, reduction = "pca")
harmony_embeddings <- RunHarmony(
  data_mat = pca_embeddings,
  meta_data = alldata_filtered@meta.data,
  vars_use = "orig.ident"
)

alldata_filtered[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_embeddings,
  key = "harmony_",
  assay = DefaultAssay(alldata_filtered)
)

reduction_use <- "harmony"

alldata_filtered <- FindNeighbors(alldata_filtered, reduction = reduction_use, dims = pcs_use)
alldata_filtered <- FindClusters(alldata_filtered, resolution = 0.2)

# Non-linear Dimensionality Reduction (UMAP & t-SNE)
alldata_filtered <- RunUMAP(alldata_filtered, reduction = reduction_use, dims = pcs_use, reduction.name = "umap")
alldata_filtered <- RunTSNE(alldata_filtered, reduction = reduction_use, dims = pcs_use, reduction.name = "tsne")

# ==========================================
# 6. Visualization & Figure Exports
# ==========================================
p_umap_cluster <- DimPlot(alldata_filtered, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("UMAP by Cluster") + theme_minimal()

p_umap_group <- DimPlot(alldata_filtered, reduction = "umap", group.by = "group", pt.size = 0.5) + 
  ggtitle("UMAP by Group") + theme_minimal()

p_pca_group <- DimPlot(alldata_filtered, reduction = "pca", group.by = "group", pt.size = 0.5) + 
  ggtitle("PCA by Group") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) + theme_minimal()

ggsave("UMAP_by_Cluster.png", plot = p_umap_cluster, width = 6.5, height = 5.5, dpi = 300)
ggsave("UMAP_by_Group.png", plot = p_umap_group, width = 6.5, height = 5.5, dpi = 300)
ggsave("PCA_by_Group.png", plot = p_pca_group, width = 6.5, height = 5.5, dpi = 300)

# ==========================================
# 7. Summary Statistics & Final Export
# ==========================================
metadata_df <- alldata_filtered@meta.data

qc_summary_formatted <- metadata_df %>%
  group_by(Sample_ID = orig.ident) %>%
  summarize(
    Pass_Cells = n(),                            
    Median_Genes_per_Cell = median(nFeature_RNA),   
    Mean_Genes_per_Cell = mean(nFeature_RNA),       
    Median_UMIs_per_Cell = median(nCount_RNA),      
    Median_percent.mt = median(percent.mt)          
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 1)))

cat("\n=== Final Processed QC Summary Table ===\n")
print(qc_summary_formatted)
write.csv(qc_summary_formatted, file = "Processed_QC_Metrics_Summary.csv", row.names = FALSE)

# Save final processed Seurat object
output_rds <- "alldata_filtered_processed.rds"
saveRDS(alldata_filtered, file = output_rds)

cat(sprintf("\n[SUCCESS] Entire pipeline completed successfully! Final object saved to: %s\n", output_rds))