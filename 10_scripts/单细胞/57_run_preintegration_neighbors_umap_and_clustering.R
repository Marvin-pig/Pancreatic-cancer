# =========================
# 57_run_preintegration_neighbors_umap_and_clustering.R
# 目的：
# 1. 基于 GSE155698 tumor tissue post-QC 对象
# 2. 先不做 integration，直接跑标准 Seurat neighbors / UMAP / clustering
# 3. 用 cluster-by-sample 和 markers 判断是否需要 Harmony
# 4. 兼容 Seurat v5 layer 结构与 marker 提取问题
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(Seurat)
  library(ggplot2)
})

setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "06_scRNA"

registry_dir <- file.path(root_dir, "00_registry")
proc_dir     <- file.path(root_dir, "02_processed_data")
result_dir   <- file.path(root_dir, "03_results")
figure_dir   <- file.path(root_dir, "04_figures")
log_dir      <- file.path(root_dir, "06_logs")

dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Step 1. 读取 post-QC 对象
# =========================================================
obj_file <- file.path(proc_dir, "55_gse155698_tumor_tissue_post_qc.rds")
if (!file.exists(obj_file)) {
  stop("未找到 55_gse155698_tumor_tissue_post_qc.rds")
}

obj <- readRDS(obj_file)

meta <- obj@meta.data
sample_col <- if ("orig.sample" %in% colnames(meta)) {
  "orig.sample"
} else if ("orig.ident" %in% colnames(meta)) {
  "orig.ident"
} else {
  NA_character_
}

if (is.na(sample_col)) {
  stop("meta.data 中未找到 orig.sample 或 orig.ident")
}

# =========================================================
# Step 2. 标准化 / HVG / PCA（省内存版）
# 关键修正：
# 1. 只对 HVG 做 ScaleData，避免 16GB 内存上限
# 2. PCA 也只用 HVG
# =========================================================
obj <- NormalizeData(
  obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

hvg_features <- VariableFeatures(obj)

hvg_tab <- data.frame(
  rank = seq_along(hvg_features),
  gene = hvg_features,
  stringsAsFactors = FALSE
)

fwrite(
  hvg_tab,
  file.path(result_dir, "57_hvg_top2000.tsv"),
  sep = "\t",
  na = "NA"
)

gc()

obj <- ScaleData(
  obj,
  features = hvg_features,
  verbose = FALSE
)

obj <- RunPCA(
  obj,
  features = hvg_features,
  npcs = 30,
  verbose = FALSE
)

# 导出 PCA embedding
pca_embed <- Embeddings(obj, reduction = "pca") %>%
  as.data.frame()
pca_embed$cell_id <- rownames(pca_embed)
pca_embed$sample_name <- obj@meta.data[[sample_col]]

fwrite(
  pca_embed,
  file.path(result_dir, "57_pca_embeddings.tsv.gz"),
  sep = "\t",
  na = "NA"
)

# PCA 按样本均值摘要
pca_summary <- pca_embed %>%
  group_by(sample_name) %>%
  summarise(
    n_cells = n(),
    mean_PC1 = mean(PC_1, na.rm = TRUE),
    mean_PC2 = mean(PC_2, na.rm = TRUE),
    mean_PC3 = mean(PC_3, na.rm = TRUE),
    sd_PC1 = sd(PC_1, na.rm = TRUE),
    sd_PC2 = sd(PC_2, na.rm = TRUE),
    sd_PC3 = sd(PC_3, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(
  pca_summary,
  file.path(result_dir, "57_pca_batch_check_summary.tsv"),
  sep = "\t",
  na = "NA"
)

pdf(file.path(figure_dir, "57_pca_by_sample.pdf"), width = 8, height = 6)
print(DimPlot(obj, reduction = "pca", group.by = sample_col))
dev.off()

pdf(file.path(figure_dir, "57_elbowplot.pdf"), width = 6, height = 4)
print(ElbowPlot(obj, ndims = 30))
dev.off()

# =========================================================
# Step 3. 不整合版 neighbors / clustering / UMAP
# =========================================================
use_dims <- 1:20
resolution_use <- 0.5

obj <- FindNeighbors(
  obj,
  dims = use_dims,
  verbose = FALSE
)

obj <- FindClusters(
  obj,
  resolution = resolution_use,
  verbose = FALSE
)

obj <- RunUMAP(
  obj,
  dims = use_dims,
  verbose = FALSE
)

saveRDS(
  obj,
  file.path(proc_dir, "57_gse155698_tumor_preintegration_umap_clustered.rds")
)

# =========================================================
# Step 4. cluster × sample 组成表
# =========================================================
meta2 <- obj@meta.data
meta2$cell_id <- rownames(meta2)
meta2$sample_name <- meta2[[sample_col]]
meta2$cluster_id <- as.character(Idents(obj))

cluster_sample_tab <- meta2 %>%
  count(cluster_id, sample_name, name = "n_cells") %>%
  group_by(cluster_id) %>%
  mutate(
    cluster_total = sum(n_cells),
    pct_in_cluster = n_cells / cluster_total
  ) %>%
  ungroup()

fwrite(
  cluster_sample_tab,
  file.path(result_dir, "57_cluster_by_sample_composition.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 5. cluster 主导样本摘要
# =========================================================
cluster_dominance <- cluster_sample_tab %>%
  group_by(cluster_id) %>%
  arrange(desc(pct_in_cluster), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    dominance_flag = ifelse(
      pct_in_cluster >= 0.70,
      "HIGH_SAMPLE_DOMINANCE",
      "ACCEPTABLE_OR_MIXED"
    )
  )

fwrite(
  cluster_dominance,
  file.path(result_dir, "57_cluster_sample_dominance_summary.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 6. UMAP 图
# =========================================================
pdf(file.path(figure_dir, "57_umap_by_cluster.pdf"), width = 7, height = 6)
print(DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE))
dev.off()

pdf(file.path(figure_dir, "57_umap_by_sample.pdf"), width = 8, height = 6)
print(DimPlot(obj, reduction = "umap", group.by = sample_col))
dev.off()

# =========================================================
# Step 7. marker 初筛（Seurat v5 兼容版）
# 关键修正：
# 1. 先 JoinLayers()
# 2. markers 为空时不要继续 group_by(cluster)
# 3. 同时兼容 avg_log2FC / avg_logFC
# =========================================================
obj <- JoinLayers(obj)

markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)

if (is.null(markers) || nrow(markers) == 0) {
  
  warning("FindAllMarkers() returned empty results after JoinLayers().")
  
  markers <- data.frame(
    cluster = character(0),
    gene = character(0),
    avg_log2FC = numeric(0),
    pct.1 = numeric(0),
    pct.2 = numeric(0),
    p_val = numeric(0),
    p_val_adj = numeric(0),
    stringsAsFactors = FALSE
  )
  
  top10_markers <- data.frame(
    cluster = character(0),
    gene = character(0),
    avg_log2FC = numeric(0),
    stringsAsFactors = FALSE
  )
  
} else {
  
  fc_col <- if ("avg_log2FC" %in% colnames(markers)) {
    "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(markers)) {
    "avg_logFC"
  } else {
    stop("markers 结果中未找到 avg_log2FC 或 avg_logFC 列。")
  }
  
  top10_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(
      order_by = .data[[fc_col]],
      n = 10,
      with_ties = FALSE
    ) %>%
    ungroup()
}

fwrite(
  markers,
  file.path(result_dir, "57_allcluster_markers.tsv.gz"),
  sep = "\t",
  na = "NA"
)

fwrite(
  top10_markers,
  file.path(result_dir, "57_top10_markers_per_cluster.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 8. integration 决策登记
# =========================================================
n_high_dominance <- sum(
  cluster_dominance$dominance_flag == "HIGH_SAMPLE_DOMINANCE",
  na.rm = TRUE
)

decision_registry <- data.frame(
  item = c(
    "dataset",
    "analysis_scope",
    "preintegration_dims",
    "preintegration_resolution",
    "n_clusters",
    "n_high_sample_dominance_clusters",
    "current_decision_rule",
    "recommended_next_action"
  ),
  value = c(
    "GSE155698",
    "tumor tissue only",
    paste(use_dims, collapse = ","),
    resolution_use,
    length(unique(meta2$cluster_id)),
    n_high_dominance,
    "If many clusters are dominated by single samples and markers are not biologically coherent, move to Harmony",
    ifelse(
      n_high_dominance >= 3,
      "consider_harmony_integration",
      "proceed_to_manual_cluster_annotation_preliminary"
    )
  ),
  stringsAsFactors = FALSE
)

fwrite(
  decision_registry,
  file.path(registry_dir, "57_integration_decision_registry.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 9. 阶段摘要
# =========================================================
stage_summary <- data.frame(
  item = c(
    "stage57_completed",
    "hvg_written",
    "pca_embeddings_written",
    "cluster_composition_written",
    "cluster_dominance_written",
    "markers_written",
    "umap_cluster_written",
    "umap_sample_written",
    "clustered_object_written",
    "next_stage"
  ),
  value = c(
    TRUE,
    file.exists(file.path(result_dir, "57_hvg_top2000.tsv")),
    file.exists(file.path(result_dir, "57_pca_embeddings.tsv.gz")),
    file.exists(file.path(result_dir, "57_cluster_by_sample_composition.tsv")),
    file.exists(file.path(result_dir, "57_cluster_sample_dominance_summary.tsv")),
    file.exists(file.path(result_dir, "57_allcluster_markers.tsv.gz")),
    file.exists(file.path(figure_dir, "57_umap_by_cluster.pdf")),
    file.exists(file.path(figure_dir, "57_umap_by_sample.pdf")),
    file.exists(file.path(proc_dir, "57_gse155698_tumor_preintegration_umap_clustered.rds")),
    "58_integration_or_preliminary_annotation_decision"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  stage_summary,
  file.path(registry_dir, "57_stage_summary.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 10. 日志
# =========================================================
writeLines(
  c(
    "57 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Run pre-integration clustering and UMAP to judge whether Harmony is necessary.",
    "",
    "[Important fixes in this script]",
    "1. ScaleData uses HVGs only to avoid memory explosion.",
    "2. JoinLayers() is applied before FindAllMarkers().",
    "3. Marker extraction is made robust to empty outputs and Seurat version differences.",
    "",
    "[Next step]",
    "Inspect cluster/sample dominance and marker coherence before deciding Harmony."
  ),
  file.path(log_dir, "57_run_preintegration_neighbors_umap_and_clustering_notes.txt")
)

cat("57_run_preintegration_neighbors_umap_and_clustering.R finished successfully.\n")
cat("Generated files:\n")
cat("- 03_results/57_hvg_top2000.tsv\n")
cat("- 03_results/57_pca_embeddings.tsv.gz\n")
cat("- 03_results/57_pca_batch_check_summary.tsv\n")
cat("- 03_results/57_cluster_by_sample_composition.tsv\n")
cat("- 03_results/57_cluster_sample_dominance_summary.tsv\n")
cat("- 03_results/57_allcluster_markers.tsv.gz\n")
cat("- 03_results/57_top10_markers_per_cluster.tsv\n")
cat("- 04_figures/57_pca_by_sample.pdf\n")
cat("- 04_figures/57_elbowplot.pdf\n")
cat("- 04_figures/57_umap_by_cluster.pdf\n")
cat("- 04_figures/57_umap_by_sample.pdf\n")
cat("- 00_registry/57_integration_decision_registry.tsv\n")
cat("- 00_registry/57_stage_summary.tsv\n")
cat("- 02_processed_data/57_gse155698_tumor_preintegration_umap_clustered.rds\n")