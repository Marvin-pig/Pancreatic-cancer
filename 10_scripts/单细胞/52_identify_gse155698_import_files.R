# =========================
# 52_identify_gse155698_import_files.R
# 目的：
# 1. 解压 GSE155698_RAW.tar
# 2. 识别可用于 Seurat 导入的真实文件
# 3. 生成 GSE155698 专用导入候选表
# 4. 更新当前单细胞策略：CRA001160 暂列 metadata_only
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
result_dir   <- file.path(root_dir, "03_results")
log_dir      <- file.path(root_dir, "06_logs")

dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Step 1. 仅处理主队列 GSE155698
# =========================================================
acc <- "GSE155698"
dataset_dir <- file.path(raw_dir, acc)
tar_file <- file.path(dataset_dir, "GSE155698_RAW.tar")
extract_dir <- file.path(dataset_dir, "extracted")

if (!file.exists(tar_file)) {
  stop("未找到 GSE155698_RAW.tar，请确认文件路径。")
}

dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

untar_status <- FALSE
untar_msg <- NA_character_

tryCatch({
  utils::untar(tar_file, exdir = extract_dir)
  untar_status <- TRUE
  untar_msg <- "untar_success"
}, error = function(e) {
  untar_msg <<- conditionMessage(e)
})

archive_audit <- data.frame(
  accession = acc,
  tar_file = tar_file,
  tar_exists = file.exists(tar_file),
  extract_dir = extract_dir,
  untar_success = untar_status,
  untar_message = untar_msg,
  stringsAsFactors = FALSE
)

fwrite(
  archive_audit,
  file.path(registry_dir, "52_gse155698_archive_audit.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 2. 扫描解压后文件
# =========================================================
scan_files <- function(path) {
  if (!dir.exists(path)) return(character(0))
  list.files(path, recursive = TRUE, full.names = TRUE, all.files = FALSE)
}

files <- scan_files(extract_dir)

inventory <- data.frame(
  full_path = files,
  file_name = basename(files),
  ext = tolower(file_ext(files)),
  size_bytes = file.info(files)$size,
  stringsAsFactors = FALSE
)

if (nrow(inventory) == 0) {
  inventory <- data.frame(
    full_path = character(0),
    file_name = character(0),
    ext = character(0),
    size_bytes = numeric(0),
    stringsAsFactors = FALSE
  )
}

inventory <- inventory %>%
  mutate(
    likely_matrix = ifelse(str_detect(tolower(file_name), "matrix|count|counts|mtx|expr"), "YES", "NO"),
    likely_barcodes = ifelse(str_detect(tolower(file_name), "barcode"), "YES", "NO"),
    likely_features = ifelse(str_detect(tolower(file_name), "feature|features|gene|genes"), "YES", "NO"),
    likely_metadata = ifelse(str_detect(tolower(file_name), "meta|annot|sample|clinical|csv|tsv|txt|xlsx"), "YES", "NO"),
    likely_object = ifelse(str_detect(tolower(file_name), "rds|rda|h5|h5ad|h5seurat"), "YES", "NO")
  )

fwrite(
  inventory,
  file.path(result_dir, "52_gse155698_file_inventory.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 3. 生成导入候选表
# =========================================================
import_candidates <- inventory %>%
  transmute(
    accession = acc,
    file_name = file_name,
    full_path = full_path,
    candidate_role = case_when(
      likely_matrix == "YES" ~ "matrix_candidate",
      likely_barcodes == "YES" ~ "barcodes_candidate",
      likely_features == "YES" ~ "features_candidate",
      likely_object == "YES" ~ "direct_object_candidate",
      likely_metadata == "YES" ~ "metadata_candidate",
      TRUE ~ "other"
    )
  )

fwrite(
  import_candidates,
  file.path(registry_dir, "52_gse155698_import_candidates.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 4. 更新当前策略
# =========================================================
strategy_update <- data.frame(
  item = c(
    "main_scRNA_dataset",
    "supplementary_scRNA_dataset",
    "cra001160_current_role",
    "current_focus",
    "next_step"
  ),
  value = c(
    "GSE155698",
    "GSE205013",
    "metadata_only",
    "identify_true_import_files_for_GSE155698",
    "write Seurat import precheck script after confirming matrix/barcodes/features or direct object"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  strategy_update,
  file.path(registry_dir, "52_scRNA_strategy_update.tsv"),
  sep = "\t",
  na = "NA"
)

writeLines(
  c(
    "52 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Main decision]",
    "Use GSE155698 as main scRNA dataset.",
    "Use GSE205013 as supplementary extension dataset.",
    "CRA001160 is currently metadata_only due to storage constraints.",
    "",
    "[Next step]",
    "Inspect 52_gse155698_import_candidates.tsv and identify the true matrix/barcodes/features or direct object."
  ),
  file.path(log_dir, "52_identify_gse155698_import_files_notes.txt")
)

cat("52_identify_gse155698_import_files.R finished successfully.\n")
cat("Generated files:\n")
cat("- 00_registry/52_gse155698_archive_audit.tsv\n")
cat("- 03_results/52_gse155698_file_inventory.tsv\n")
cat("- 00_registry/52_gse155698_import_candidates.tsv\n")
cat("- 00_registry/52_scRNA_strategy_update.tsv\n")
cat("- 06_logs/52_identify_gse155698_import_files_notes.txt\n")