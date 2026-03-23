# =========================
# 54d_rescue_gsm4710703_and_merge.R
# 目的：
# 1. 递归定位 GSM4710703_PDAC_TISSUE_14 的真实 10X 三件套
# 2. 单独导入该失败样本
# 3. 若已有 54c 对象，则直接 merge
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(Seurat)
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
sample_root <- file.path(
  root_dir,
  "01_raw_data/GSE155698/gsm_extracted",
  sample_name
)

if (!dir.exists(sample_root)) {
  stop("未找到样本目录: ", sample_root)
}

# ---------------------------------------------------------
# Step 1. 递归寻找三件套目录
# ---------------------------------------------------------
all_files <- list.files(sample_root, recursive = TRUE, full.names = TRUE)

file_tab <- data.frame(
  full_path = all_files,
  file_name = basename(all_files),
  dir_name  = dirname(all_files),
  stringsAsFactors = FALSE
)

dir_summary <- file_tab %>%
  group_by(dir_name) %>%
  summarise(
    has_matrix = any(str_detect(tolower(file_name), "^matrix\\.mtx(\\.gz)?$")),
    has_barcodes = any(str_detect(tolower(file_name), "^barcodes\\.tsv(\\.gz)?$")),
    has_features = any(str_detect(tolower(file_name), "^(features|genes)\\.tsv(\\.gz)?$")),
    .groups = "drop"
  ) %>%
  mutate(import_ready = has_matrix & has_barcodes & has_features)

fwrite(
  dir_summary,
  file.path(registry_dir, "54d_gsm4710703_candidate_dirs.tsv"),
  sep = "\t", na = "NA"
)

ready_dirs <- dir_summary %>% filter(import_ready)

if (nrow(ready_dirs) == 0) {
  stop("未找到可导入的三件套目录，请检查 54d_gsm4710703_candidate_dirs.tsv")
}

# 若有多个，优先选带 filtered_feature_bc_matrix 的
target_dir <- ready_dirs %>%
  mutate(priority = case_when(
    str_detect(dir_name, "filtered_feature_bc_matrix") ~ 1,
    str_detect(dir_name, "filtered_gene_bc_matrices") ~ 2,
    TRUE ~ 9
  )) %>%
  arrange(priority, dir_name) %>%
  slice(1) %>%
  pull(dir_name)

rescue_plan <- data.frame(
  sample_name = sample_name,
  chosen_dir = target_dir,
  stringsAsFactors = FALSE
)

fwrite(
  rescue_plan,
  file.path(registry_dir, "54d_gsm4710703_rescue_plan.tsv"),
  sep = "\t", na = "NA"
)

# ---------------------------------------------------------
# Step 2. 单独导入 sample14
# ---------------------------------------------------------
mat <- Read10X(data.dir = target_dir)

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
  file.path(proc_dir, "54d_gsm4710703_single_sample_seurat.rds")
)

qc14 <- data.frame(
  sample_name = sample_name,
  n_cells = ncol(seu14),
  n_genes = nrow(seu14),
  median_nFeature_RNA = median(seu14$nFeature_RNA),
  median_nCount_RNA = median(seu14$nCount_RNA),
  stringsAsFactors = FALSE
)

fwrite(
  qc14,
  file.path(table_dir, "54d_gsm4710703_pre_qc_summary.tsv"),
  sep = "\t", na = "NA"
)

# ---------------------------------------------------------
# Step 3. 若已有 54c merged object，则合并进去
# ---------------------------------------------------------
base_rds <- file.path(proc_dir, "54c_gse155698_tumor_tissue_raw_seurat.rds")
merge_status <- "base_object_not_found"

if (file.exists(base_rds)) {
  base_obj <- readRDS(base_rds)
  
  if (!(sample_name %in% unique(base_obj$sample_name))) {
    merged_obj <- merge(base_obj, y = seu14)
    saveRDS(
      merged_obj,
      file.path(proc_dir, "54d_gse155698_tumor_tissue_raw_seurat_merged.rds")
    )
    merge_status <- "merged_success"
  } else {
    merge_status <- "sample_already_present_in_base_object"
  }
}

merge_tab <- data.frame(
  item = c("sample_name", "target_dir", "merge_status"),
  value = c(sample_name, target_dir, merge_status),
  stringsAsFactors = FALSE
)

fwrite(
  merge_tab,
  file.path(registry_dir, "54d_gsm4710703_merge_status.tsv"),
  sep = "\t", na = "NA"
)

writeLines(
  c(
    "54d completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Rescue GSM4710703_PDAC_TISSUE_14 by recursive detection of 10X triplet.",
    "",
    "[Next step]",
    "If merge_status = merged_success, proceed to 55 QC workflow."
  ),
  file.path(log_dir, "54d_rescue_gsm4710703_notes.txt")
)

cat("54d_rescue_gsm4710703_and_merge.R finished successfully.\n")