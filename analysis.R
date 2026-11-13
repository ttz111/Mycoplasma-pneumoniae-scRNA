library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(RColorBrewer)
# con_4
setwd("/home/c001/analysis/yyy")
counts <- Read10X(data.dir = "/home/c001/analysis/yyy/con_4/filtered_cell_gene_matrix")

con_4_seurat <- CreateSeuratObject(
  counts = counts,
  project = "con_4",
  min.cells = 3,      # 过滤掉表达少于3个细胞的基因
  min.features = 200  # 过滤掉表达基因数少于200的低质量细胞
)
saveRDS(con_4_seurat, file = "con_4_seurat.rds")
# con_5
counts <- Read10X(data.dir = "/home/c001/analysis/yyy/con_5/filtered_cell_gene_matrix")

con_5_seurat <- CreateSeuratObject(
  counts = counts,
  project = "con_5",
  min.cells = 3,      # 过滤掉表达少于3个细胞的基因
  min.features = 200  # 过滤掉表达基因数少于200的低质量细胞
)
saveRDS(con_5_seurat, file = "con_5_seurat.rds")
# mol_7
counts <- Read10X(data.dir = "/home/c001/analysis/yyy/mol_7/filtered_cell_gene_matrix")

mol_7_seurat <- CreateSeuratObject(
  counts = counts,
  project = "mol_7",
  min.cells = 3,      # 过滤掉表达少于3个细胞的基因
  min.features = 200  # 过滤掉表达基因数少于200的低质量细胞
)
saveRDS(mol_7_seurat, file = "mol_7_seurat.rds")

# mol_8
counts <- Read10X(data.dir = "/home/c001/analysis/yyy/mol_8/filtered_cell_gene_matrix")

mol_8_seurat <- CreateSeuratObject(
  counts = counts,
  project = "mol_8",
  min.cells = 3,      # 过滤掉表达少于3个细胞的基因
  min.features = 200  # 过滤掉表达基因数少于200的低质量细胞
)
saveRDS(mol_8_seurat, file = "mol_8_seurat.rds")

###all cell
con_4_seurat <- readRDS("con_4_seurat.rds")
con_5_seurat <- readRDS("con_5_seurat.rds")
mol_7_seurat <- readRDS("mol_7_seurat.rds")
mol_8_seurat <- readRDS("mol_8_seurat.rds")

# 2. 将 4 个对象合并为一个整体的 Seurat 对象
# add.cell.ids 会在 cell barcode 前加上样本名前缀，防止不同样本间 barcode 重名
alldata <- merge(
  x = con_4_seurat,
  y = list(con_5_seurat, mol_7_seurat, mol_8_seurat),
  add.cell.ids = c("con_4", "con_5", "mol_7", "mol_8"),
  project = "all"
)
alldata$group <- ifelse(grepl("^con", alldata$orig.ident), "con", "mol")

alldata[["percent.mt"]] <- PercentageFeatureSet(alldata, pattern = "^mt-")

table(alldata$orig.ident) # 查看每个样本包含的细胞数
table(alldata$group)      # 查看各组细胞数

saveRDS(alldata, file = "alldata_merged.rds")

alldata <- readRDS("/home/c001/analysis/yyy/alldata_merged.rds")
##质量控制
cols_miqc <- scale_color_gradientn(
  colors = rev(brewer.pal(n = 9, name = "Blues")[3:9]), # 移除极白颜色
  limits = c(0, 1.00),
  breaks = c(0, 0.25, 0.50, 0.75, 1.00),
  na.value = "black" # 异常值显示为黑色
)

cols_selection <- c(
  "Features & Mitochondrial Filter" = "#FF3300", # 红色/粉红
  "Features Filter" = "#FF9933",                 # 橙色
  "Mitochondrial Filter" = "#9966CC",             # 紫色
  "Pass" = "#33CC99"                            # 绿色
)

# 3. 定义通用绘图主题
my_theme <- theme_minimal() + 
  theme(
    text = element_text(color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    panel.grid.minor = element_blank(), # 隐藏细网格线
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )
##未过滤的原始细胞数
raw_counts <- table(alldata$orig.ident)

alldata_filtered <- subset(
  alldata, 
  subset = nFeature_RNA > 200 & 
    nFeature_RNA < 6000 & 
    percent.mt < 10 & 
    nCount_RNA > 500
)

filtered_counts <- table(alldata_filtered$orig.ident)

qc_summary <- data.frame(
  Sample = names(raw_counts),
  Raw_Cells = as.numeric(raw_counts),
  Pass_Cells = as.numeric(filtered_counts[names(raw_counts)]),
  Filtered_Out = as.numeric(raw_counts) - as.numeric(filtered_counts[names(raw_counts)]),
  Pass_Rate = paste0(round(as.numeric(filtered_counts[names(raw_counts)]) / as.numeric(raw_counts) * 100, 2), "%")
)

qc_summary
VlnPlot(
  alldata_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "orig.ident",
  ncol = 3,
  pt.size = 0
) & 
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
  )
ggsave("QC_VlnPlot.png", width = 11, height = 4.5, dpi = 300)






alldata_filtered <- NormalizeData(
  alldata_filtered, 
  normalization.method = "LogNormalize", 
  scale.factor = 10000
)

alldata_filtered <- FindVariableFeatures(
  alldata_filtered, 
  selection.method = "vst", 
  nfeatures = 2000
)

alldata_filtered <- ScaleData(
  alldata_filtered, 
  vars.to.regress = "percent.mt"
)

alldata_filtered <- RunPCA(alldata_filtered, npcs = 30, verbose = FALSE)

alldata_filtered <- FindNeighbors(alldata_filtered, dims = 1:20)
alldata_filtered <- FindClusters(alldata_filtered, resolution = 0.2)

alldata_filtered <- RunUMAP(alldata_filtered, dims = 1:20)


p1 <- DimPlot(
  alldata_filtered, 
  reduction = "umap", 
  group.by = "seurat_clusters", 
  label = TRUE, 
  pt.size = 0.5
) + labs(title = "UMAP Clusters")

ggsave("all_cell.png", plot = p1, width = 8, height = 6, dpi = 300)

head(alldata_filtered)
p2 <- DimPlot(
  alldata_filtered, 
  reduction = "umap", 
  group.by = "group", 
  pt.size = 0.5
)

ggsave("group_cell.png", plot = p2,width = 8, height = 6, dpi = 300)


###

library(miQC)
library(flexmix)
library(patchwork)
library(dplyr)
library(scales)
library(Seurat)


# ==========================================
# 1. 线粒体比例计算与数据提取
# ==========================================
alldata[["percent.mt"]] <- PercentageFeatureSet(alldata, pattern = "^mt-")

df_ctrl <- FetchData(alldata, vars = c("nFeature_RNA", "nCount_RNA", "percent.mt", "orig.ident")) %>%
  filter(grepl("^con", orig.ident))
df_ctrl$cell_id <- rownames(df_ctrl)

df_treat <- FetchData(alldata, vars = c("nFeature_RNA", "nCount_RNA", "percent.mt", "orig.ident")) %>%
  filter(grepl("^mol", orig.ident))
df_treat$cell_id <- rownames(df_treat)

# ==========================================
# 2. 修复后的 miQC 安全拟合函数
# ==========================================
fit_miqc_safe <- function(df) {
  set.seed(123)
  model <- flexmix(percent.mt ~ nFeature_RNA, data = df, k = 2)
  params <- parameters(model)
  post <- posterior(model)
  
  # 防护机制：如果 flexmix 只拟合出了 1 个成分（ncol(post) < 2）
  if (ncol(post) < 2) {
    df$prob_compromised <- 0 # 无法区分受损细胞，默认标记为低受损概率
  } else {
    # 使用后验概率加权求线粒体均值，完美避开下标越界与 NA 问题
    mean_mt_1 <- sum(df$percent.mt * post[, 1]) / sum(post[, 1])
    mean_mt_2 <- sum(df$percent.mt * post[, 2]) / sum(post[, 2])
    
    # 线粒体均值更高的成分确定为受损组 (Compromised)
    comp_idx <- if (mean_mt_1 > mean_mt_2) 1 else 2
    df$prob_compromised <- post[, comp_idx]
  }
  
  return(list(df = df, params = params))
}

# 执行拟合
fit_c <- fit_miqc_safe(df_ctrl)
df_ctrl <- fit_c$df
params_ctrl <- fit_c$params

fit_t <- fit_miqc_safe(df_treat)
df_treat <- fit_t$df
params_treat <- fit_t$params

# ==========================================
# 3. 划分 Selection 阈值并写回 alldata
# ==========================================
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

# 验证各样本中的 Selection 分布
cat("--- 各样本过滤分类结果统计 ---\n")
print(table(alldata$orig.ident, alldata$Selection))

# ==========================================
# 4. 绘制四张质控散点图 (a, b, c, d)
# ==========================================
# 图 a: Control 组 miQC 概率图
int1_c <- params_ctrl[1, 1]; slope1_c <- params_ctrl[2, 1]
int2_c   <- if(ncol(params_ctrl) >= 2) params_ctrl[1, 2] else NA
slope2_c <- if(ncol(params_ctrl) >= 2) params_ctrl[2, 2] else NA

p_a <- ggplot(df_ctrl, aes(x = nFeature_RNA, y = percent.mt, color = prob_compromised)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_gradientn(
    colors = c("#0B1D3A", "#1D4E89", "#3A86C8", "#72B4EB"),
    limits = c(0, 1.00),
    breaks = c(0.25, 0.50, 0.75, 1.00)
  ) +
  geom_abline(intercept = int1_c, slope = slope1_c, color = "black", linewidth = 1.2) +
  {if(!is.na(int2_c)) geom_abline(intercept = int2_c, slope = slope2_c, color = "black", linewidth = 1.2)} +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(
    title = "Control: miQC Probability",
    x = "Number of detected genes", 
    y = "Mitochondrial RNA (%)", 
    color = "Probability\ncompromised"
  ) +
  theme_minimal()

# 图 b: Treatment 组 miQC 概率图
int1_t <- params_treat[1, 1]; slope1_t <- params_treat[2, 1]
int2_t   <- if(ncol(params_treat) >= 2) params_treat[1, 2] else NA
slope2_t <- if(ncol(params_treat) >= 2) params_treat[2, 2] else NA

p_b <- ggplot(df_treat, aes(x = nFeature_RNA, y = percent.mt, color = prob_compromised)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_gradientn(
    colors = c("#0B1D3A", "#1D4E89", "#3A86C8", "#72B4EB"),
    limits = c(0, 1.00),
    breaks = c(0.25, 0.50, 0.75, 1.00)
  ) +
  geom_abline(intercept = int1_t, slope = slope1_t, color = "black", linewidth = 1.2) +
  {if(!is.na(int2_t)) geom_abline(intercept = int2_t, slope = slope2_t, color = "black", linewidth = 1.2)} +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(
    title = "Treatment: miQC Probability",
    x = "Number of detected genes", 
    y = "Mitochondrial RNA (%)", 
    color = "Probability\ncompromised"
  ) +
  theme_minimal()

scatter_colors <- c(
  "Pass"                             = "#33CC99", 
  "Features Filter"                  = "#FF9933", 
  "Mitochondrial Filter"             = "#9966CC", 
  "Features & Mitochondrial Filter"  = "#FF3300"
)

# 图 c: Control 组 UMI vs Gene
p_c <- ggplot(df_ctrl, aes(x = nCount_RNA, y = nFeature_RNA, color = Selection)) +
  geom_point(size = 1.2, alpha = 0.6) +
  scale_color_manual(values = scatter_colors) +
  scale_x_log10(breaks = c(1e3, 1e4, 1e5, 1e6), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = c(300, 1000, 3000, 10000, 30000)) +
  geom_smooth(data = filter(df_ctrl, Selection == "Pass"), method = "lm", color = "black", linewidth = 0.8, se = FALSE) +
  geom_hline(yintercept = f_low, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_hline(yintercept = f_high, linetype = "dashed", color = "blue", linewidth = 0.8) +
  labs(
    title = "Control: Selection Result",
    x = "Total UMI counts per cell", 
    y = "Number expressed genes",
    color = "Selection"
  ) +
  theme_minimal()

# 图 d: Treatment 组 UMI vs Gene
p_d <- ggplot(df_treat, aes(x = nCount_RNA, y = nFeature_RNA, color = Selection)) +
  geom_point(size = 1.2, alpha = 0.6) +
  scale_color_manual(values = scatter_colors) +
  scale_x_log10(breaks = c(1e3, 1e4, 1e5, 1e6), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = c(300, 1000, 3000, 10000, 30000)) +
  geom_smooth(data = filter(df_treat, Selection == "Pass"), method = "lm", color = "black", linewidth = 0.8, se = FALSE) +
  geom_hline(yintercept = f_low, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_hline(yintercept = f_high, linetype = "dashed", color = "blue", linewidth = 0.8) +
  labs(
    title = "Treatment: Selection Result",
    x = "Total UMI counts per cell", 
    y = "Number expressed genes",
    color = "Selection"
  ) +
  theme_minimal()

# ==========================================
# 5. 绘制样本百分比堆叠柱状图 (p_prop)
# ==========================================
df_stack <- alldata@meta.data %>%
  filter(!is.na(Selection)) %>%
  group_by(orig.ident, Selection) %>%
  summarise(count = n(), .groups = "drop")

df_stack$Selection <- factor(df_stack$Selection, levels = c(
  "Features & Mitochondrial Filter",
  "Features Filter",
  "Mitochondrial Filter",
  "Pass"
))

p_prop <- ggplot(df_stack, aes(x = orig.ident, y = count, fill = Selection)) +
  geom_bar(stat = "identity", position = "fill", width = 0.65) +
  scale_fill_manual(values = c(
    "Pass"                             = "#1EA362",
    "Mitochondrial Filter"             = "#9B51E0",
    "Features Filter"                  = "#E67E22",
    "Features & Mitochondrial Filter"  = "#EB5757"
  )) +
  scale_y_continuous(labels = scales::percent_format(scale = 100), expand = c(0, 0)) +
  labs(x = "Sample", y = "Proportion of Cells", fill = "Selection") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    legend.position = "bottom"
  )

# ==========================================
# 6. 显示与导出图像
# ==========================================
p_quad <- (p_a | p_b) / (p_c | p_d) + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 16, fontface = "bold"))

print(p_quad)
print(p_prop)

ggsave("QC_Scatter_Combined.png", p_quad, width = 11, height = 9, dpi = 300)
ggsave("QC_Proportion_Stack.png", p_prop, width = 8, height = 4.5, dpi = 300)

# 1. 过滤：只保留 Pass 状态的合格细胞，创建清洗后的新对象
alldata_filtered <- subset(alldata, subset = Selection == "Pass")

# 2. 检查过滤前后细胞数对比
cat("过滤前细胞数:", ncol(alldata), "\n")
cat("过滤后细胞数:", ncol(alldata_filtered), "\n")

# 3. 基于 alldata_filtered 进行后续 Standard Workflow（标准化、高变基因、PCA、UMAP）
alldata_filtered <- NormalizeData(alldata_filtered)
alldata_filtered <- FindVariableFeatures(alldata_filtered, selection.method = "vst", nfeatures = 2000)

# 4. 绘制高变基因图（此时用的是清洗后的高质量数据）
top10 <- head(VariableFeatures(alldata_filtered), 10)
p_vfg <- VariableFeaturePlot(alldata_filtered, cols = c("black", "red"))
p_vfg_labeled <- LabelPoints(plot = p_vfg, points = top10, repel = TRUE)
print(p_vfg_labeled)
ggsave("p_vfg_labeled.png", p_vfg_labeled, width = 7, height = 5, dpi = 300)
library(harmony) 

# ==========================================
# 1. 过滤细胞，生成 alldata_filtered
# ==========================================
alldata_filtered <- subset(alldata, subset = Selection == "Pass")

cat("过滤前细胞数:", ncol(alldata), "\n")
cat("过滤后细胞数:", ncol(alldata_filtered), "\n")

# ==========================================
# 2. 数据标准化 (NormalizeData)
# ==========================================
# LogNormalize 将表达量转换为 Log1p(CPM/100)
alldata_filtered <- NormalizeData(
  alldata_filtered, 
  normalization.method = "LogNormalize", 
  scale.factor = 10000
)

# ==========================================
# 3. 寻找高变基因 (FindVariableFeatures)
# ==========================================
alldata_filtered <- FindVariableFeatures(
  alldata_filtered, 
  selection.method = "vst", 
  nfeatures = 2000
)

# ==========================================
# 4. 数据归一化 (ScaleData)
# ==========================================
# 默认对 2000 个高变基因进行 Scale（也可传入 all.genes 对所有基因 Scale）
all.genes <- rownames(alldata_filtered)
alldata_filtered <- ScaleData(
  alldata_filtered, 
  features = all.genes,
  vars.to.regress = c("percent.mt") # 可选：去除线粒体基因比例带来的技术噪音影响
)

# ==========================================
# 5. 主成分分析 (RunPCA)
# ==========================================
alldata_filtered <- RunPCA(
  alldata_filtered, 
  features = VariableFeatures(object = alldata_filtered),
  verbose = FALSE
)

# 查看 PCA 碎石图（ElbowPlot）确定保留的主成分维度 (dims)
p_elbow <- ElbowPlot(alldata_filtered, ndims = 30)
print(p_elbow)

# 假设根据 ElbowPlot 选择前 20 个主成分维度
pcs_use <- 1:20

# ==========================================
# 6. 批次效应矫正 (Harmony, 强烈推荐)
# ==========================================
# 单细胞数据多样本（如 con_4, con_5, mol_7, mol_8）通常存在批次差异，使用 Harmony 矫正 orig.ident
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
# ==========================================
# 7. UMAP / t-SNE 降维与细胞聚类
# ==========================================
# 注：若未运行 Harmony，以下 reduction 参数请替换为 "pca"
reduction_use <- "harmony" 

# 构建细胞近邻图
alldata_filtered <- FindNeighbors(
  alldata_filtered, 
  reduction = reduction_use, 
  dims = pcs_use
)

# 寻找细胞亚群 (resolution 决定聚类粒度，通常在 0.2 - 1.2 之间)
alldata_filtered <- FindClusters(
  alldata_filtered, 
  resolution = 0.2
)

# 运行 UMAP 与 t-SNE 降维
alldata_filtered <- RunUMAP(
  alldata_filtered, 
  reduction = reduction_use, 
  dims = pcs_use,
  reduction.name = "umap"
)

alldata_filtered <- RunTSNE(
  alldata_filtered, 
  reduction = reduction_use, 
  dims = pcs_use,
  reduction.name = "tsne"
)

# ==========================================
# 8. 可视化检查
# ==========================================
# 按 Cluster 和 Sample 分别查看 UMAP
p_umap_cluster <- DimPlot(alldata_filtered, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("UMAP by Cluster")

p_umap_sample <- DimPlot(alldata_filtered, reduction = "umap", group.by = "group", pt.size = 0.5) + 
  ggtitle("UMAP by group")

p_pca_sample <- DimPlot(alldata_filtered, reduction = "pca", group.by = "group", pt.size = 0.5) + 
  ggtitle("PCA by Group") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("UMAP_by_Cluster.png", plot = p_umap_cluster, width = 6.5, height = 5.5, dpi = 300)
ggsave("UMAP_by_Group.png", plot = p_umap_sample, width = 6.5, height = 5.5, dpi = 300)
ggsave("PCA_by_Group.png", plot = p_pca_sample, width = 6.5, height = 5.5, dpi = 300)

# 保存最终对象
saveRDS(alldata_filtered, file = "alldata_filtered_processed.rds")

# 使用过滤后的对象统计最终质控概览
metadata_df <- alldata_filtered@meta.data

qc_summary_formatted <- metadata_df %>%
  group_by(Sample_ID = orig.ident) %>%
  summarize(
    Pass_Cells = n(),                               # 过滤后保留的细胞数
    Median_Genes_per_Cell = median(nFeature_RNA),   # 基因数中位数
    Mean_Genes_per_Cell = mean(nFeature_RNA),       # 基因数平均值
    Median_UMIs_per_Cell = median(nCount_RNA),      # UMI数/Reads中位数
    Median_percent.mt = median(percent.mt)          # 线粒体比例中位数
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 1)))

print(qc_summary_formatted)
