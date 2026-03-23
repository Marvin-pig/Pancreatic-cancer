# =========================
# 56_normalize_hvg_and_batch_check_gse155698_tumor.R
# 目的：
# 1. 读取 55 阶段 post-QC Seurat 对象
# 2. 完成标准化、HVG 筛选、ScaleData、PCA
# 3. 输出按样本分层的 PCA / batch check 结果
# 4. 判断后续是否需要 integration
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
sample_col <- if ("orig.sample" %in% colnames(meta)) "orig.sample" else
  if ("orig.ident" %in% colnames(meta)) "orig.ident" else NA_character_

if (is.na(sample_col)) {
  stop("meta.data 中未找到 orig.sample 或 orig.ident")
}

# =========================================================
# =========================================================
# Step 2. 标准化 + HVG
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

hvg <- VariableFeatures(obj)

hvg_tab <- data.frame(
  rank = seq_along(hvg),
  gene = hvg,
  stringsAsFactors = FALSE
)

fwrite(
  hvg_tab,
  file.path(result_dir, "56_hvg_top2000.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 3. ScaleData + PCA
# 只 scale HVG，避免 16GB 内存爆掉
# =========================================================
gc()

obj <- ScaleData(
  obj,
  features = hvg,
  verbose = FALSE
)

gc()

obj <- RunPCA(
  obj,
  features = hvg,
  npcs = 30,
  verbose = FALSE
)
# =========================================================
# Step 4. PCA embeddings 导出
# =========================================================
pca_embed <- Embeddings(obj, reduction = "pca") %>%
  as.data.frame()
pca_embed$cell_id <- rownames(pca_embed)
pca_embed$sample_name <- obj@meta.data[[sample_col]]

fwrite(
  pca_embed,
  file.path(result_dir, "56_pca_embeddings.tsv.gz"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 5. PCA 图
# =========================================================
pdf(file.path(figure_dir, "56_pca_by_sample.pdf"), width = 8, height = 6)
print(DimPlot(obj, reduction = "pca", group.by = sample_col))
dev.off()

pdf(file.path(figure_dir, "56_elbowplot.pdf"), width = 6, height = 4)
print(ElbowPlot(obj, ndims = 30))
dev.off()

# =========================================================
# Step 6. 样本分层 PCA 摘要
# =========================================================
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
  file.path(result_dir, "56_pca_batch_check_summary.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 7. 记录策略判断
# =========================================================
integration_registry <- data.frame(
  item = c(
    "dataset",
    "analysis_scope",
    "n_cells_post_qc",
    "normalization_method",
    "hvg_n",
    "pca_npcs",
    "current_decision",
    "decision_note"
  ),
  value = c(
    "GSE155698",
    "tumor tissue only",
    ncol(obj),
    "LogNormalize",
    2000,
    30,
    "batch_check_first",
    "Inspect 56_pca_by_sample.pdf and 56_pca_batch_check_summary.tsv before deciding Harmony/RPCA integration"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  integration_registry,
  file.path(registry_dir, "56_batch_check_decision_registry.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 8. 阶段摘要
# =========================================================
stage_summary <- data.frame(
  item = c(
    "stage56_completed",
    "normalized_object_written",
    "hvg_written",
    "pca_embeddings_written",
    "pca_summary_written",
    "next_stage"
  ),
  value = c(
    TRUE,
    file.exists(file.path(proc_dir, "56_gse155698_tumor_tissue_post_qc_norm_pca.rds")),
    file.exists(file.path(result_dir, "56_hvg_top2000.tsv")),
    file.exists(file.path(result_dir, "56_pca_embeddings.tsv.gz")),
    file.exists(file.path(result_dir, "56_pca_batch_check_summary.tsv")),
    "57_neighbors_umap_preintegration_or_integration_decision"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  stage_summary,
  file.path(registry_dir, "56_stage_summary.tsv"),
  sep = "\t", na = "NA"
)

writeLines(
  c(
    "56 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Normalize post-QC tumor tissue object, compute HVGs, run PCA, and inspect sample-driven structure.",
    "",
    "[Next step]",
    "Decide whether to proceed with standard Seurat clustering or apply batch integration."
  ),
  file.path(log_dir, "56_normalize_hvg_and_batch_check_notes.txt")
)

cat("56_normalize_hvg_and_batch_check_gse155698_tumor.R finished successfully.\n")