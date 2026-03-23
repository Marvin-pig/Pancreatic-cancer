# =========================
# 53_unpack_gse155698_gsm_archives_and_identify_import_ready_samples.R
# 目的：
# 1. 对 GSE155698 外层 tar 解压后的 GSM 子包继续二级解压
# 2. 识别真正可用于 Seurat 导入的样本级文件
# 3. 生成样本分层与主分析纳入建议
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
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

acc <- "GSE155698"
dataset_dir <- file.path(raw_dir, acc)
outer_extract_dir <- file.path(dataset_dir, "extracted")
gsm_extract_root  <- file.path(dataset_dir, "gsm_extracted")

dir.create(gsm_extract_root, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Step 1. 找到所有 GSM 子包
# =========================================================
gsm_archives <- list.files(
  outer_extract_dir,
  pattern = "^GSM.*\\.tar\\.gz$",
  full.names = TRUE
)

if (length(gsm_archives) == 0) {
  stop("未在外层解压目录中找到 GSM*.tar.gz 子包。")
}

sample_manifest <- data.frame(
  gsm_archive = basename(gsm_archives),
  gsm_id = str_extract(basename(gsm_archives), "GSM\\d+"),
  sample_label = str_remove(basename(gsm_archives), "\\.tar\\.gz$"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    sample_scope = case_when(
      str_detect(sample_label, "PDAC_TISSUE") ~ "tumor_tissue",
      str_detect(sample_label, "AdjNorm_TISSUE") ~ "adjacent_normal_tissue",
      str_detect(sample_label, "PDAC_PBMC") ~ "pdac_pbmc",
      str_detect(sample_label, "Healthy_PBMC") ~ "healthy_pbmc",
      TRUE ~ "other"
    ),
    include_for_main_tissue_analysis = ifelse(sample_scope == "tumor_tissue", "YES", "NO")
  )

fwrite(
  sample_manifest,
  file.path(registry_dir, "53_gse155698_sample_archive_manifest.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 2. 二级解压每个 GSM 子包
# =========================================================
unpack_status <- lapply(seq_along(gsm_archives), function(i) {
  f <- gsm_archives[i]
  sample_name <- str_remove(basename(f), "\\.tar\\.gz$")
  target_dir <- file.path(gsm_extract_root, sample_name)
  
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  
  ok <- FALSE
  msg <- NA_character_
  
  tryCatch({
    utils::untar(f, exdir = target_dir)
    ok <- TRUE
    msg <- "untar_success"
  }, error = function(e) {
    msg <<- conditionMessage(e)
  })
  
  data.frame(
    gsm_archive = basename(f),
    sample_name = sample_name,
    target_dir = target_dir,
    unpack_success = ok,
    unpack_message = msg,
    stringsAsFactors = FALSE
  )
})

unpack_status <- bind_rows(unpack_status)

fwrite(
  unpack_status,
  file.path(registry_dir, "53_gse155698_gsm_unpack_status.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 3. 扫描二级解压结果
# =========================================================
safe_list_files <- function(path) {
  if (!dir.exists(path)) return(character(0))
  list.files(path, recursive = TRUE, full.names = TRUE, all.files = FALSE)
}

all_files <- lapply(seq_len(nrow(unpack_status)), function(i) {
  td <- unpack_status$target_dir[i]
  sn <- unpack_status$sample_name[i]
  fs <- safe_list_files(td)
  
  if (length(fs) == 0) {
    return(data.frame(
      sample_name = sn,
      full_path = NA_character_,
      file_name = NA_character_,
      ext = NA_character_,
      size_bytes = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    sample_name = sn,
    full_path = fs,
    file_name = basename(fs),
    ext = tolower(file_ext(fs)),
    size_bytes = file.info(fs)$size,
    stringsAsFactors = FALSE
  )
})

nested_inventory <- bind_rows(all_files) %>%
  mutate(
    likely_matrix = ifelse(str_detect(tolower(file_name), "matrix|count|counts|mtx"), "YES", "NO"),
    likely_barcodes = ifelse(str_detect(tolower(file_name), "barcode"), "YES", "NO"),
    likely_features = ifelse(str_detect(tolower(file_name), "feature|features|gene|genes"), "YES", "NO"),
    likely_metadata = ifelse(str_detect(tolower(file_name), "meta|annot|sample|clinical|csv|tsv|txt|xlsx"), "YES", "NO"),
    likely_object = ifelse(str_detect(tolower(file_name), "rds|rda|h5|h5ad|h5seurat"), "YES", "NO")
  )

fwrite(
  nested_inventory,
  file.path(result_dir, "53_gse155698_nested_file_inventory.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 4. 判断哪些样本已具备导入条件
# =========================================================
import_ready <- nested_inventory %>%
  group_by(sample_name) %>%
  summarise(
    has_matrix = ifelse(any(likely_matrix == "YES", na.rm = TRUE), "YES", "NO"),
    has_barcodes = ifelse(any(likely_barcodes == "YES", na.rm = TRUE), "YES", "NO"),
    has_features = ifelse(any(likely_features == "YES", na.rm = TRUE), "YES", "NO"),
    has_metadata = ifelse(any(likely_metadata == "YES", na.rm = TRUE), "YES", "NO"),
    has_direct_object = ifelse(any(likely_object == "YES", na.rm = TRUE), "YES", "NO"),
    import_ready = ifelse(
      any(likely_object == "YES", na.rm = TRUE) |
        (
          any(likely_matrix == "YES", na.rm = TRUE) &
            any(likely_barcodes == "YES", na.rm = TRUE) &
            any(likely_features == "YES", na.rm = TRUE)
        ),
      "YES", "NO"
    ),
    .groups = "drop"
  ) %>%
  left_join(
    sample_manifest %>% select(sample_label, sample_scope, include_for_main_tissue_analysis),
    by = c("sample_name" = "sample_label")
  )

fwrite(
  import_ready,
  file.path(registry_dir, "53_gse155698_import_ready_samples.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 5. 生成分析范围建议
# =========================================================
analysis_scope <- data.frame(
  item = c(
    "primary_analysis_scope",
    "adjacent_normal_scope",
    "pbmc_scope",
    "current_recommended_import_order"
  ),
  value = c(
    "PDAC_TISSUE_*",
    "AdjNorm_TISSUE_*",
    "Exclude PDAC_PBMC_* and Healthy_PBMC_* from current tumor-tissue main analysis",
    "First import tumor_tissue samples that are import_ready == YES"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  analysis_scope,
  file.path(registry_dir, "53_gse155698_analysis_scope.tsv"),
  sep = "\t", na = "NA"
)

writeLines(
  c(
    "53 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Key point]",
    "GSE155698 outer archive contains nested GSM-level tar.gz files.",
    "Need second-level unpack before Seurat import.",
    "",
    "[Current recommendation]",
    "Use PDAC_TISSUE_* as main scRNA analysis scope.",
    "Do not mix PBMC into tumor tissue main analysis."
  ),
  file.path(log_dir, "53_unpack_gse155698_gsm_archives_notes.txt")
)

cat("53_unpack_gse155698_gsm_archives_and_identify_import_ready_samples.R finished successfully.\n")