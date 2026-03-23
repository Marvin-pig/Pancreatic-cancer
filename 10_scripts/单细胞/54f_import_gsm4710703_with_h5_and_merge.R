# =========================
# 54f_import_gsm4710703_with_h5_and_merge.R
# 目的：
# 1. 使用 Read10X_h5 单独导入 GSM4710703_PDAC_TISSUE_14
# 2. 若已有 15 样本对象，则 merge 进去
# 3. 生成完整的 16 样本 tumor tissue 原始对象
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(Seurat)
  library(hdf5r)
})

setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "06_scRNA"

registry_dir <- file.path(root_dir, "00_registry")
proc_dir     <- file.path(root_dir, "02_processed_data")
table_dir    <- file.path(root_dir, "05_tables")
log_dir      <- file.path(root_dir, "06_logs")

dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

sample_name <- "GSM4710703_PDAC_TISSUE_14"
sample_root <- file.path(root_dir, "01_raw_data/GSE155698/gsm_extracted", sample_name)

if (!dir.exists(sample_root)) {
  stop("未找到样本目录: ", sample_root)
}

# ---------------------------------------------------------
# Step 1. 递归查找 h5 文件
# ---------------------------------------------------------
all_files <- list.files(sample_root, recursive = TRUE, full.names = TRUE)

h5_files <- all_files[grepl("\\.h5$", all_files, ignore.case = TRUE)]

if (length(h5_files) == 0) {
  stop("未找到 .h5 文件，请检查样本目录。")
}

# 若有多个 h5，优先选 filtered_feature_bc_matrix.h5
target_h5 <- h5_files[grepl("filtered_feature_bc_matrix\\.h5$", h5_files, ignore.case = TRUE)]
if (length(target_h5) == 0) {
  target_h5 <- h5_files[1]
} else {
  target_h5 <- target_h5[1]
}

h5_plan <- data.frame(
  sample_name = sample_name,
  h5_file = target_h5,
  stringsAsFactors = FALSE
)

fwrite(
  h5_plan,
  file.path(registry_dir, "54f_gsm4710703_h5_plan.tsv"),
  sep = "\t", na = "NA"
)

# ---------------------------------------------------------
# Step 2. 导入 sample14
# ---------------------------------------------------------
mat <- Read10X_h5(target_h5)

seu14 <- CreateSeuratObject(
  counts = mat,
  project = "GSE155698_PDACCohort",
  min.cells = 3,
  min.features = 200
)

seu14$orig.ident <- sample_name
seu14$sample_name <- sample_name
seu14$sample_scope <- "tumor_tissue"

saveRDS(
  seu14,
  file.path(proc_dir, "54f_gsm4710703_single_sample_seurat.rds")
)

qc14 <- data.frame(
  sample_name = sample_name,
  import_method = "Read10X_h5",
  n_cells = ncol(seu14),
  n_genes = nrow(seu14),
  median_nFeature_RNA = median(seu14$nFeature_RNA),
  median_nCount_RNA = median(seu14$nCount_RNA),
  stringsAsFactors = FALSE
)

fwrite(
  qc14,
  file.path(table_dir, "54f_gsm4710703_pre_qc_summary.tsv"),
  sep = "\t", na = "NA"
)

# ---------------------------------------------------------
# Step 3. merge 到已有对象
# 优先尝试 54c，再尝试 54b
# ---------------------------------------------------------
base_rds_candidates <- c(
  file.path(proc_dir, "54c_gse155698_tumor_tissue_raw_seurat.rds"),
  file.path(proc_dir, "54b_gse155698_tumor_tissue_raw_seurat.rds")
)

base_rds <- base_rds_candidates[file.exists(base_rds_candidates)][1]

if (is.na(base_rds)) {
  stop("未找到已有的 15 样本对象（54b/54c），请先确认基础对象已生成。")
}

base_obj <- readRDS(base_rds)

if (sample_name %in% unique(base_obj$sample_name)) {
  merge_status <- "sample_already_present"
  merged_obj <- base_obj
} else {
  merged_obj <- merge(base_obj, y = seu14)
  merge_status <- "merged_success"
}

saveRDS(
  merged_obj,
  file.path(proc_dir, "54f_gse155698_tumor_tissue_raw_seurat_16samples.rds")
)

merge_tab <- data.frame(
  item = c("sample_name", "base_object", "merge_status"),
  value = c(sample_name, base_rds, merge_status),
  stringsAsFactors = FALSE
)

fwrite(
  merge_tab,
  file.path(registry_dir, "54f_gsm4710703_merge_status.tsv"),
  sep = "\t", na = "NA"
)

writeLines(
  c(
    "54f completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Import GSM4710703_PDAC_TISSUE_14 via h5 and merge into tumor tissue object.",
    "",
    "[Next step]",
    "Proceed to QC metrics and filtering with full 16 tumor tissue samples."
  ),
  file.path(log_dir, "54f_import_gsm4710703_with_h5_and_merge_notes.txt")
)

cat("54f_import_gsm4710703_with_h5_and_merge.R finished successfully.\n")