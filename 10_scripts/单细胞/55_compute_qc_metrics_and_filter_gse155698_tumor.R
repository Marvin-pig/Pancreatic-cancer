# =========================
# 55_compute_qc_metrics_and_filter_gse155698_tumor.R
# 目的：
# 1. 读取 54 阶段最终 raw Seurat 对象
# 2. 计算标准 QC 指标
# 3. 输出过滤前细胞级/样本级审计结果
# 4. 固定过滤规则并生成 post-QC 对象
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
# Step 1. 读取最终 raw 对象
# =========================================================
obj_file <- "/mnt/data/54_gse155698_tumor_tissue_raw_seurat_final.rds"
if (!file.exists(obj_file)) {
  obj_file <- file.path(proc_dir, "54_gse155698_tumor_tissue_raw_seurat_final.rds")
}
if (!file.exists(obj_file)) {
  stop("未找到 54_gse155698_tumor_tissue_raw_seurat_final.rds")
}

obj <- readRDS(obj_file)

# =========================================================
# Step 2. 计算 QC 指标
# =========================================================
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RPL|^RPS")
obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^HB[ABDEGMQZ]")

meta <- obj@meta.data
meta$cell_id <- rownames(meta)

sample_col <- if ("orig.sample" %in% colnames(meta)) "orig.sample" else
  if ("orig.ident" %in% colnames(meta)) "orig.ident" else NA_character_

if (is.na(sample_col)) {
  stop("meta.data 中未找到 orig.sample 或 orig.ident")
}

meta$sample_name <- meta[[sample_col]]

# =========================================================
# Step 3. 输出细胞级 QC 表
# =========================================================
qc_cell <- meta %>%
  dplyr::select(
    cell_id,
    sample_name,
    nFeature_RNA,
    nCount_RNA,
    percent.mt,
    percent.ribo,
    percent.hb
  )

fwrite(
  qc_cell,
  file.path(result_dir, "55_qc_metrics_per_cell.tsv.gz"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 4. 输出样本级 QC 摘要
# =========================================================
qc_sample <- qc_cell %>%
  group_by(sample_name) %>%
  summarise(
    n_cells = n(),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    median_percent_mt = median(percent.mt, na.rm = TRUE),
    p95_percent_mt = quantile(percent.mt, 0.95, na.rm = TRUE),
    median_percent_ribo = median(percent.ribo, na.rm = TRUE),
    median_percent_hb = median(percent.hb, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(
  qc_sample,
  file.path(result_dir, "55_qc_metrics_by_sample_summary.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 5. 冻结过滤规则
# 先用较稳妥的一版，后续若分布异常再微调
# =========================================================
filter_rule <- data.frame(
  rule_item = c(
    "min_features",
    "max_features",
    "min_counts",
    "max_percent_mt",
    "analysis_scope",
    "rule_note"
  ),
  value = c(
    300,
    7500,
    500,
    20,
    "GSE155698 tumor tissue only",
    "Initial global QC rule; revise only if sample-level diagnostics show obvious distortion"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  filter_rule,
  file.path(registry_dir, "55_filter_rule_registry.tsv"),
  sep = "\t", na = "NA"
)

min_features <- 300
max_features <- 7500
min_counts   <- 500
max_percent_mt <- 20

# =========================================================
# Step 6. 生成过滤前 QC 图
# =========================================================
pdf(file.path(figure_dir, "55_qc_violin_nFeature.pdf"), width = 14, height = 6)
print(VlnPlot(obj, features = "nFeature_RNA", group.by = sample_col, pt.size = 0))
dev.off()

pdf(file.path(figure_dir, "55_qc_violin_nCount.pdf"), width = 14, height = 6)
print(VlnPlot(obj, features = "nCount_RNA", group.by = sample_col, pt.size = 0))
dev.off()

pdf(file.path(figure_dir, "55_qc_violin_percent_mt.pdf"), width = 14, height = 6)
print(VlnPlot(obj, features = "percent.mt", group.by = sample_col, pt.size = 0))
dev.off()

# =========================================================
# Step 7. 执行过滤
# =========================================================
obj_post <- subset(
  obj,
  subset =
    nFeature_RNA >= min_features &
    nFeature_RNA <= max_features &
    nCount_RNA >= min_counts &
    percent.mt <= max_percent_mt
)

# =========================================================
# Step 8. 输出过滤后摘要
# =========================================================
post_meta <- obj_post@meta.data
post_meta$sample_name <- post_meta[[sample_col]]

post_summary <- post_meta %>%
  mutate(cell_id = rownames(post_meta)) %>%
  group_by(sample_name) %>%
  summarise(
    n_cells_post_qc = n(),
    median_nFeature_RNA_post_qc = median(nFeature_RNA, na.rm = TRUE),
    median_nCount_RNA_post_qc = median(nCount_RNA, na.rm = TRUE),
    median_percent_mt_post_qc = median(percent.mt, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(
  post_summary,
  file.path(result_dir, "55_post_qc_summary_by_sample.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 9. 保存 post-QC 对象
# =========================================================
saveRDS(
  obj_post,
  file.path(proc_dir, "55_gse155698_tumor_tissue_post_qc.rds")
)

# =========================================================
# Step 10. 状态摘要
# =========================================================
stage_summary <- data.frame(
  item = c(
    "stage55_completed",
    "raw_object_loaded",
    "qc_metrics_per_cell_written",
    "qc_metrics_by_sample_written",
    "filter_rule_written",
    "post_qc_object_written",
    "n_cells_pre_qc",
    "n_cells_post_qc",
    "next_stage"
  ),
  value = c(
    TRUE,
    TRUE,
    file.exists(file.path(result_dir, "55_qc_metrics_per_cell.tsv.gz")),
    file.exists(file.path(result_dir, "55_qc_metrics_by_sample_summary.tsv")),
    file.exists(file.path(registry_dir, "55_filter_rule_registry.tsv")),
    file.exists(file.path(proc_dir, "55_gse155698_tumor_tissue_post_qc.rds")),
    ncol(obj),
    ncol(obj_post),
    "56_normalization_hvg_and_batch_check"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  stage_summary,
  file.path(registry_dir, "55_stage_summary.tsv"),
  sep = "\t", na = "NA"
)

writeLines(
  c(
    "55 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Compute QC metrics and filter GSE155698 tumor tissue raw Seurat object.",
    "",
    "[Next step]",
    "Run normalization, HVG selection, and pre-integration batch inspection."
  ),
  file.path(log_dir, "55_compute_qc_metrics_and_filter_notes.txt")
)

cat("55_compute_qc_metrics_and_filter_gse155698_tumor.R finished successfully.\n")