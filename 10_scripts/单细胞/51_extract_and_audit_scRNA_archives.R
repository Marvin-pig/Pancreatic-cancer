# =========================
# 51_extract_and_audit_scRNA_archives.R
# 目的：
# 1. 解压已手动下载的 scRNA tar 文件
# 2. 重新扫描解压后的真实文件
# 3. 判断各数据集是否已具备进入 Seurat 的最低文件条件
# 4. 记录 CRA001160 当前仅 metadata 可用的状态
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
# Step 1. 定义待检查数据集
# =========================================================
dataset_tab <- data.frame(
  accession = c("GSE155698", "GSE205013", "CRA001160"),
  raw_subdir = c("GSE155698", "GSE205013", "CRA001160"),
  archive_file = c("GSE155698_RAW.tar", "GSE205013_RAW.tar", NA),
  metadata_file = c(NA, NA, "CRA001160.xlsx"),
  stringsAsFactors = FALSE
)

# =========================================================
# Step 2. 解压 tar 文件
# =========================================================
extract_status <- list()

for (i in seq_len(nrow(dataset_tab))) {
  acc <- dataset_tab$accession[i]
  subdir <- file.path(raw_dir, dataset_tab$raw_subdir[i])
  tar_name <- dataset_tab$archive_file[i]
  
  archive_path <- if (!is.na(tar_name)) file.path(subdir, tar_name) else NA_character_
  extract_dir  <- file.path(subdir, "extracted")
  
  archive_exists <- !is.na(archive_path) && file.exists(archive_path)
  
  extracted_ok <- FALSE
  extract_msg <- NA_character_
  
  if (archive_exists) {
    dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
    
    tryCatch({
      utils::untar(archive_path, exdir = extract_dir)
      extracted_ok <- TRUE
      extract_msg <- "untar_success"
    }, error = function(e) {
      extract_msg <<- conditionMessage(e)
    })
  } else {
    extract_msg <- "archive_not_found_or_not_applicable"
  }
  
  extract_status[[i]] <- data.frame(
    accession = acc,
    archive_path = ifelse(is.na(archive_path), NA, archive_path),
    archive_exists = archive_exists,
    extract_dir = extract_dir,
    extracted_ok = extracted_ok,
    extract_message = extract_msg,
    stringsAsFactors = FALSE
  )
}

extract_status_tab <- bind_rows(extract_status)

fwrite(
  extract_status_tab,
  file.path(registry_dir, "51_scRNA_archive_audit.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 3. 扫描解压后的文件
# =========================================================
safe_list_files <- function(path) {
  if (!dir.exists(path)) return(character(0))
  list.files(path, recursive = TRUE, full.names = TRUE, all.files = FALSE)
}

scan_rows <- list()

for (i in seq_len(nrow(dataset_tab))) {
  acc <- dataset_tab$accession[i]
  subdir <- file.path(raw_dir, dataset_tab$raw_subdir[i])
  extract_dir <- file.path(subdir, "extracted")
  metadata_file <- dataset_tab$metadata_file[i]
  
  files_extracted <- safe_list_files(extract_dir)
  
  if (length(files_extracted) > 0) {
    tmp <- data.frame(
      accession = acc,
      file_source = "extracted",
      full_path = files_extracted,
      file_name = basename(files_extracted),
      ext = tolower(file_ext(files_extracted)),
      size_bytes = file.info(files_extracted)$size,
      stringsAsFactors = FALSE
    )
    scan_rows[[length(scan_rows) + 1]] <- tmp
  }
  
  if (!is.na(metadata_file)) {
    meta_path <- file.path(subdir, metadata_file)
    if (file.exists(meta_path)) {
      tmp2 <- data.frame(
        accession = acc,
        file_source = "metadata_only",
        full_path = meta_path,
        file_name = basename(meta_path),
        ext = tolower(file_ext(meta_path)),
        size_bytes = file.info(meta_path)$size,
        stringsAsFactors = FALSE
      )
      scan_rows[[length(scan_rows) + 1]] <- tmp2
    }
  }
}

file_scan <- bind_rows(scan_rows)

if (nrow(file_scan) == 0) {
  file_scan <- data.frame(
    accession = character(0),
    file_source = character(0),
    full_path = character(0),
    file_name = character(0),
    ext = character(0),
    size_bytes = numeric(0),
    stringsAsFactors = FALSE
  )
}

file_scan <- file_scan %>%
  mutate(
    likely_matrix = ifelse(str_detect(tolower(file_name), "matrix|count|counts|mtx|expr"), "YES", "NO"),
    likely_barcodes = ifelse(str_detect(tolower(file_name), "barcode"), "YES", "NO"),
    likely_features = ifelse(str_detect(tolower(file_name), "feature|features|gene|genes"), "YES", "NO"),
    likely_metadata = ifelse(str_detect(tolower(file_name), "meta|annot|sample|clinical|xlsx|csv|tsv"), "YES", "NO"),
    likely_seurat_object = ifelse(str_detect(tolower(file_name), "rds|rda|h5seurat|h5ad"), "YES", "NO")
  )

fwrite(
  file_scan,
  file.path(result_dir, "51_scRNA_extracted_file_scan.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 4. 数据集就绪状态判断
# =========================================================
dataset_readiness <- file_scan %>%
  group_by(accession) %>%
  summarise(
    n_files = n(),
    has_likely_matrix = ifelse(any(likely_matrix == "YES"), "YES", "NO"),
    has_likely_barcodes = ifelse(any(likely_barcodes == "YES"), "YES", "NO"),
    has_likely_features = ifelse(any(likely_features == "YES"), "YES", "NO"),
    has_likely_metadata = ifelse(any(likely_metadata == "YES"), "YES", "NO"),
    has_likely_seurat_object = ifelse(any(likely_seurat_object == "YES"), "YES", "NO"),
    preliminary_ready_for_import = ifelse(
      any(likely_seurat_object == "YES") |
        (any(likely_matrix == "YES") & any(likely_metadata == "YES")),
      "YES", "NO"
    ),
    .groups = "drop"
  )

# 保证 CRA001160 即使只有 metadata 也写进去
if (!"CRA001160" %in% dataset_readiness$accession) {
  dataset_readiness <- bind_rows(
    dataset_readiness,
    data.frame(
      accession = "CRA001160",
      n_files = 0,
      has_likely_matrix = "NO",
      has_likely_barcodes = "NO",
      has_likely_features = "NO",
      has_likely_metadata = ifelse(file.exists(file.path(raw_dir, "CRA001160", "CRA001160.xlsx")), "YES", "NO"),
      has_likely_seurat_object = "NO",
      preliminary_ready_for_import = "NO",
      stringsAsFactors = FALSE
    )
  )
}

fwrite(
  dataset_readiness,
  file.path(registry_dir, "51_scRNA_dataset_readiness.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 5. 状态摘要
# =========================================================
status_summary <- data.frame(
  item = c(
    "gse155698_archive_present",
    "gse205013_archive_present",
    "cra001160_metadata_present",
    "file_scan_written",
    "dataset_readiness_written"
  ),
  value = c(
    file.exists(file.path(raw_dir, "GSE155698", "GSE155698_RAW.tar")),
    file.exists(file.path(raw_dir, "GSE205013", "GSE205013_RAW.tar")),
    file.exists(file.path(raw_dir, "CRA001160", "CRA001160.xlsx")),
    file.exists(file.path(result_dir, "51_scRNA_extracted_file_scan.tsv")),
    file.exists(file.path(registry_dir, "51_scRNA_dataset_readiness.tsv"))
  ),
  stringsAsFactors = FALSE
)

fwrite(
  status_summary,
  file.path(registry_dir, "51_scRNA_status_summary.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 6. 日志
# =========================================================
writeLines(
  c(
    "51 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Current practical strategy]",
    "Use GSE155698 as main dataset.",
    "Use GSE205013 as supplementary extension dataset.",
    "Record CRA001160 as metadata-only because raw files were not downloaded locally.",
    "",
    "[Next step]",
    "Identify true counts/metadata file paths and start Seurat import precheck for GSE155698."
  ),
  file.path(log_dir, "51_extract_and_audit_scRNA_archives_notes.txt")
)

cat("51_extract_and_audit_scRNA_archives.R finished successfully.\n")
cat("Generated files:\n")
cat("- 00_registry/51_scRNA_archive_audit.tsv\n")
cat("- 03_results/51_scRNA_extracted_file_scan.tsv\n")
cat("- 00_registry/51_scRNA_dataset_readiness.tsv\n")
cat("- 00_registry/51_scRNA_status_summary.tsv\n")
cat("- 06_logs/51_extract_and_audit_scRNA_archives_notes.txt\n")
cat("Next step:\n")
cat("- Identify true counts files for GSE155698 and prepare Seurat import precheck\n")