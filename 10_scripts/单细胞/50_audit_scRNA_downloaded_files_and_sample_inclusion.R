# =========================
# 50_audit_scRNA_downloaded_files_and_sample_inclusion.R
# 目的：
# 1. 扫描已下载的 scRNA 数据文件
# 2. 自动生成文件存在性审计结果
# 3. 为 sample-level 纳入/排除建立正式审计表
# 4. 为后续 Seurat 预处理前的数据接收确认做准备
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tools)
})

setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "06_scRNA"

registry_dir <- file.path(root_dir, "00_registry")
raw_dir      <- file.path(root_dir, "01_raw_data")
proc_dir     <- file.path(root_dir, "02_processed_data")
result_dir   <- file.path(root_dir, "03_results")
log_dir      <- file.path(root_dir, "06_logs")

dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Step 1. 读取 49 步产物
# =========================================================
download_plan_file <- file.path(registry_dir, "49_scRNA_download_plan.tsv")

if (!file.exists(download_plan_file)) {
  stop("缺少 49_scRNA_download_plan.tsv，请先完成 49。")
}

download_plan <- fread(download_plan_file)

# =========================================================
# Step 2. 扫描各 accession 目录下的实际文件
# =========================================================
safe_list_files <- function(path) {
  if (!dir.exists(path)) return(character(0))
  list.files(path, recursive = TRUE, full.names = TRUE, all.files = FALSE)
}

all_scan <- lapply(seq_len(nrow(download_plan)), function(i) {
  acc  <- download_plan$accession[i]
  ddir <- file.path(raw_dir, acc)
  
  files <- safe_list_files(ddir)
  
  if (length(files) == 0) {
    return(data.frame(
      accession = acc,
      local_dataset_dir = ddir,
      file_name = NA_character_,
      relative_path = NA_character_,
      ext = NA_character_,
      size_bytes = NA_real_,
      file_exists = "NO",
      likely_counts = "NO",
      likely_barcodes = "NO",
      likely_features = "NO",
      likely_metadata = "NO",
      likely_annotation = "NO",
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    accession = acc,
    local_dataset_dir = ddir,
    file_name = basename(files),
    relative_path = gsub(paste0("^", normalizePath(ddir, winslash = "/"), "/?"), "", normalizePath(files, winslash = "/")),
    ext = tolower(file_ext(files)),
    size_bytes = file.info(files)$size,
    file_exists = "YES",
    likely_counts = ifelse(str_detect(tolower(basename(files)), "matrix|count|counts|expr|expression"), "YES", "NO"),
    likely_barcodes = ifelse(str_detect(tolower(basename(files)), "barcode"), "YES", "NO"),
    likely_features = ifelse(str_detect(tolower(basename(files)), "feature|genes|gene"), "YES", "NO"),
    likely_metadata = ifelse(str_detect(tolower(basename(files)), "meta|annotation|sample|clinical"), "YES", "NO"),
    likely_annotation = ifelse(str_detect(tolower(basename(files)), "celltype|cell_type|annotation|label"), "YES", "NO"),
    stringsAsFactors = FALSE
  )
})

file_scan <- bind_rows(all_scan)

fwrite(
  file_scan,
  file.path(result_dir, "50_scRNA_downloaded_file_scan.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 3. 生成文件审计结果表
# =========================================================
file_audit_result <- file_scan %>%
  group_by(accession, local_dataset_dir) %>%
  summarise(
    n_files = sum(file_exists == "YES", na.rm = TRUE),
    has_any_file = ifelse(any(file_exists == "YES"), "YES", "NO"),
    has_likely_counts = ifelse(any(likely_counts == "YES"), "YES", "NO"),
    has_likely_barcodes = ifelse(any(likely_barcodes == "YES"), "YES", "NO"),
    has_likely_features = ifelse(any(likely_features == "YES"), "YES", "NO"),
    has_likely_metadata = ifelse(any(likely_metadata == "YES"), "YES", "NO"),
    has_likely_annotation = ifelse(any(likely_annotation == "YES"), "YES", "NO"),
    prelim_ready_for_manual_review = ifelse(any(file_exists == "YES"), "YES", "NO"),
    .groups = "drop"
  )

fwrite(
  file_audit_result,
  file.path(registry_dir, "50_scRNA_file_audit_result.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 4. 生成待人工填写的样本纳入审计表
# =========================================================
sample_inclusion_audit <- data.frame(
  accession = rep(download_plan$accession, each = 1),
  sample_id = NA_character_,
  patient_id = NA_character_,
  tissue_type = NA_character_,
  tumor_normal = NA_character_,
  primary_metastatic = NA_character_,
  treatment_status = NA_character_,
  file_source = NA_character_,
  include_for_main_analysis = NA_character_,
  include_for_validation_analysis = NA_character_,
  exclude_reason = NA_character_,
  note = NA_character_,
  stringsAsFactors = FALSE
)

fwrite(
  sample_inclusion_audit,
  file.path(registry_dir, "50_scRNA_sample_inclusion_audit.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 5. 更新状态摘要
# =========================================================
status_summary <- data.frame(
  item = c(
    "download_plan_loaded",
    "file_scan_written",
    "file_audit_result_written",
    "sample_inclusion_audit_written",
    "ready_for_manual_file_review",
    "ready_for_sample_level_inclusion_decision"
  ),
  value = c(
    TRUE,
    file.exists(file.path(result_dir, "50_scRNA_downloaded_file_scan.tsv")),
    file.exists(file.path(registry_dir, "50_scRNA_file_audit_result.tsv")),
    file.exists(file.path(registry_dir, "50_scRNA_sample_inclusion_audit.tsv")),
    TRUE,
    TRUE
  ),
  stringsAsFactors = FALSE
)

fwrite(
  status_summary,
  file.path(registry_dir, "50_scRNA_status_summary.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 6. 日志
# =========================================================
writeLines(
  c(
    "50 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Audit downloaded scRNA files and prepare sample-level inclusion review.",
    "",
    "[Key outputs]",
    "- 03_results/50_scRNA_downloaded_file_scan.tsv",
    "- 00_registry/50_scRNA_file_audit_result.tsv",
    "- 00_registry/50_scRNA_sample_inclusion_audit.tsv",
    "- 00_registry/50_scRNA_status_summary.tsv",
    "",
    "[Next step]",
    "Manually inspect downloaded files, identify true counts/metadata files, and fill sample inclusion audit before Seurat preprocessing."
  ),
  file.path(log_dir, "50_audit_scRNA_downloaded_files_and_sample_inclusion_notes.txt")
)

cat("50_audit_scRNA_downloaded_files_and_sample_inclusion.R finished successfully.\n")
cat("Generated files:\n")
cat("- 03_results/50_scRNA_downloaded_file_scan.tsv\n")
cat("- 00_registry/50_scRNA_file_audit_result.tsv\n")
cat("- 00_registry/50_scRNA_sample_inclusion_audit.tsv\n")
cat("- 00_registry/50_scRNA_status_summary.tsv\n")
cat("- 06_logs/50_audit_scRNA_downloaded_files_and_sample_inclusion_notes.txt\n")
cat("Next step:\n")
cat("- Manually review true counts/metadata files and fill sample-level inclusion decisions before Seurat preprocessing\n")