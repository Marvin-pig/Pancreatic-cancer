# =========================
# 49_prepare_scRNA_download_and_metadata_audit.R
# 目的：
# 1. 为已固定的 PDAC 单细胞队列建立下载准备表
# 2. 建立文件级与元数据级审计框架
# 3. 为后续 Seurat 预处理前的数据核查做准备
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "06_scRNA"

registry_dir  <- file.path(root_dir, "00_registry")
raw_dir       <- file.path(root_dir, "01_raw_data")
proc_dir      <- file.path(root_dir, "02_processed_data")
result_dir    <- file.path(root_dir, "03_results")
table_dir     <- file.path(root_dir, "05_tables")
log_dir       <- file.path(root_dir, "06_logs")

dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Step 1. 读取当前已固定的数据集 registry
# =========================================================
registry_file <- file.path(registry_dir, "47_scRNA_dataset_registry.tsv")

if (!file.exists(registry_file)) {
  stop("缺少 47_scRNA_dataset_registry.tsv，请先完成 48。")
}

registry <- fread(registry_file)

# =========================================================
# Step 2. 生成下载准备表
# =========================================================
download_plan <- registry %>%
  mutate(
    local_dataset_dir = file.path(raw_dir, accession),
    download_status = "TODO",
    download_priority = case_when(
      include_as_main == "YES" ~ "HIGH",
      include_as_validation == "YES" ~ "HIGH",
      TRUE ~ "MEDIUM"
    ),
    expected_files = case_when(
      accession == "GSE155698" ~ "matrix / counts table / metadata / supplementary annotations if available",
      accession == "CRA001160" ~ "raw matrix / processed matrix / metadata / sample annotation",
      accession == "GSE205013" ~ "filtered matrices / metadata / treatment grouping information"
    ),
    expected_format = case_when(
      accession == "GSE155698" ~ "GEO supplementary files",
      accession == "CRA001160" ~ "NGDC/GSA archive files",
      accession == "GSE205013" ~ "GEO supplementary files"
    ),
    metadata_audit_needed = "YES",
    sample_level_manifest_needed = "YES",
    ready_for_seurat = "NO"
  )

fwrite(
  download_plan,
  file.path(registry_dir, "49_scRNA_download_plan.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 3. 生成文件级审计模板
# =========================================================
file_audit_template <- data.frame(
  accession = c("GSE155698", "CRA001160", "GSE205013"),
  local_dataset_dir = file.path(raw_dir, c("GSE155698", "CRA001160", "GSE205013")),
  file_name = NA_character_,
  file_type = NA_character_,
  file_role = NA_character_,
  file_exists = NA_character_,
  readable = NA_character_,
  contains_counts = NA_character_,
  contains_barcodes = NA_character_,
  contains_features = NA_character_,
  contains_metadata = NA_character_,
  contains_cell_annotation = NA_character_,
  note = NA_character_,
  stringsAsFactors = FALSE
)

fwrite(
  file_audit_template,
  file.path(registry_dir, "49_scRNA_file_audit_template.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 4. 生成元数据审计模板
# =========================================================
metadata_audit_template <- data.frame(
  accession = c("GSE155698", "CRA001160", "GSE205013"),
  metadata_file = NA_character_,
  sample_id_column = NA_character_,
  patient_id_column = NA_character_,
  tissue_column = NA_character_,
  treatment_column = NA_character_,
  diagnosis_column = NA_character_,
  response_column = NA_character_,
  tumor_normal_column = NA_character_,
  primary_metastatic_column = NA_character_,
  cell_barcode_column = NA_character_,
  celltype_annotation_column = NA_character_,
  metadata_complete_for_seurat = "NO",
  note = NA_character_,
  stringsAsFactors = FALSE
)

fwrite(
  metadata_audit_template,
  file.path(registry_dir, "49_scRNA_metadata_audit_template.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 5. 生成样本清单模板
# =========================================================
sample_manifest_template <- data.frame(
  accession = c("GSE155698", "CRA001160", "GSE205013"),
  sample_id = NA_character_,
  patient_id = NA_character_,
  tissue_type = NA_character_,
  disease_status = NA_character_,
  treatment_status = NA_character_,
  primary_or_metastatic = NA_character_,
  expected_use = NA_character_,
  include_for_main_analysis = NA_character_,
  exclude_reason = NA_character_,
  note = NA_character_,
  stringsAsFactors = FALSE
)

fwrite(
  sample_manifest_template,
  file.path(registry_dir, "49_scRNA_sample_manifest_template.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 6. 更新阶段状态
# =========================================================
status_summary <- data.frame(
  item = c(
    "registry_loaded",
    "candidate_dataset_n",
    "download_plan_written",
    "file_audit_template_written",
    "metadata_audit_template_written",
    "sample_manifest_template_written",
    "ready_for_manual_download_and_audit"
  ),
  value = c(
    TRUE,
    nrow(registry),
    file.exists(file.path(registry_dir, "49_scRNA_download_plan.tsv")),
    file.exists(file.path(registry_dir, "49_scRNA_file_audit_template.tsv")),
    file.exists(file.path(registry_dir, "49_scRNA_metadata_audit_template.tsv")),
    file.exists(file.path(registry_dir, "49_scRNA_sample_manifest_template.tsv")),
    TRUE
  ),
  stringsAsFactors = FALSE
)

fwrite(
  status_summary,
  file.path(registry_dir, "49_scRNA_status_summary.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 7. 输出日志
# =========================================================
writeLines(
  c(
    "49 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Prepare scRNA download plan and metadata/file audit framework before Seurat preprocessing.",
    "",
    "[Main dataset]",
    "GSE155698",
    "",
    "[Validation dataset]",
    "CRA001160",
    "",
    "[Supplementary dataset]",
    "GSE205013",
    "",
    "[Next step]",
    "Start 50: manually download files and fill file/metadata/sample audit templates."
  ),
  file.path(log_dir, "49_prepare_scRNA_download_and_metadata_audit_notes.txt")
)

cat("49_prepare_scRNA_download_and_metadata_audit.R finished successfully.\n")
cat("Generated files:\n")
cat("- 00_registry/49_scRNA_download_plan.tsv\n")
cat("- 00_registry/49_scRNA_file_audit_template.tsv\n")
cat("- 00_registry/49_scRNA_metadata_audit_template.tsv\n")
cat("- 00_registry/49_scRNA_sample_manifest_template.tsv\n")
cat("- 00_registry/49_scRNA_status_summary.tsv\n")
cat("- 06_logs/49_prepare_scRNA_download_and_metadata_audit_notes.txt\n")
cat("Next step:\n")
cat("- Start 50: download files and fill audit templates before Seurat preprocessing\n")