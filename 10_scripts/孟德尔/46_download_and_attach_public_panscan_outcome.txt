# =========================
# 46_download_and_attach_public_panscan_outcome.R
# 目的：
# 1. 直接下载 PanScan 公共 analysis 文件
# 2. 自动检查表头和预览
# 3. 自动挑选最可能可用的 pancreatic cancer summary stats 文件
# 4. 尝试标准化为项目主 outcome 文件
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "05_MR"

raw_outcome_dir <- file.path(root_dir, "01_raw_data/outcome/PanScan_PanC4")
registry_dir    <- file.path(root_dir, "00_registry")
log_dir         <- file.path(root_dir, "06_logs")

dir.create(raw_outcome_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Step 1. 下载公共分析文件
# =========================================================
public_urls <- c(
  "https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs000206/analyses/phs000206.pha002874.txt.gz",
  "https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs000206/analyses/phs000206.pha002889.txt.gz",
  "https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs000206/analyses/phs000206.pha012458.txt"
)

download_status_list <- list()

for (i in seq_along(public_urls)) {
  url_i <- public_urls[i]
  dest_i <- file.path(raw_outcome_dir, basename(url_i))
  
  ok <- FALSE
  msg <- ""
  
  if (!file.exists(dest_i)) {
    tryCatch({
      download.file(url_i, destfile = dest_i, mode = "wb", quiet = FALSE)
      ok <- file.exists(dest_i)
    }, error = function(e) {
      msg <<- conditionMessage(e)
    })
  } else {
    ok <- TRUE
    msg <- "already_exists"
  }
  
  download_status_list[[i]] <- data.frame(
    url = url_i,
    dest = dest_i,
    exists = file.exists(dest_i),
    ok = ok,
    message = msg,
    size_bytes = ifelse(file.exists(dest_i), file.info(dest_i)$size, NA_real_),
    stringsAsFactors = FALSE
  )
}

download_status <- bind_rows(download_status_list)

fwrite(
  download_status,
  file.path(registry_dir, "46_download_status.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 2. 扫描表头与预览
# =========================================================
downloaded_files <- download_status$dest[file.exists(download_status$dest)]

if (length(downloaded_files) == 0) {
  stop("没有任何 PanScan 公共文件成功下载。请先检查网络或 URL。")
}

safe_read_header <- function(path) {
  tryCatch(names(fread(path, nrows = 0)), error = function(e) character(0))
}

safe_read_preview <- function(path, n = 10) {
  tryCatch(fread(path, nrows = n), error = function(e) NULL)
}

file_check_list <- list()

for (f in downloaded_files) {
  hdr <- safe_read_header(f)
  prv <- safe_read_preview(f, 10)
  
  # 导出 header
  fwrite(
    data.frame(column_name = hdr, stringsAsFactors = FALSE),
    file.path(registry_dir, paste0(basename(f), "_header.tsv")),
    sep = "\t"
  )
  
  # 导出 preview
  if (!is.null(prv)) {
    fwrite(
      prv,
      file.path(registry_dir, paste0(basename(f), "_preview.tsv")),
      sep = "\t"
    )
  }
  
  file_check_list[[basename(f)]] <- data.frame(
    file = basename(f),
    full_path = f,
    n_header_cols = length(hdr),
    has_ID = "ID" %in% hdr,
    has_Chr = "Chr" %in% hdr,
    has_Position = "Position" %in% hdr,
    has_Allele1 = "Allele1" %in% hdr,
    has_Allele2 = "Allele2" %in% hdr,
    has_Freq1 = "Freq1" %in% hdr,
    has_Effect = "Effect" %in% hdr,
    has_StdErr = "StdErr" %in% hdr,
    has_Pvalue = "P-value" %in% hdr,
    has_hg19_chr = "hg19_chr" %in% hdr,
    has_hg19_pos = "hg19_pos" %in% hdr,
    stringsAsFactors = FALSE
  )
}

file_check_tab <- bind_rows(file_check_list)

# 给一个简单评分：越接近 summary stats 模板分越高
file_check_tab <- file_check_tab %>%
  mutate(
    score = has_ID + has_Chr + has_Position + has_Allele1 + has_Allele2 +
      has_Freq1 + has_Effect + has_StdErr + has_Pvalue
  ) %>%
  arrange(desc(score), file)

fwrite(
  file_check_tab,
  file.path(registry_dir, "46_public_panscan_file_check.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 3. 自动选择最可能的主 outcome 文件
# =========================================================
best_file <- file_check_tab %>%
  slice(1)

if (nrow(best_file) == 0) {
  stop("未能识别任何候选 PanScan 文件。")
}

selected_file <- best_file$full_path[1]

selection_summary <- data.frame(
  item = c("selected_file", "selection_score"),
  value = c(selected_file, best_file$score[1]),
  stringsAsFactors = FALSE
)

fwrite(
  selection_summary,
  file.path(registry_dir, "46_selected_primary_outcome_file.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 4. 尝试标准化成 outcome 模板
# 模板要求：
# SNP, effect_allele, other_allele, beta, se, pval, eaf, outcome
# =========================================================
x <- fread(selected_file)

# 先记录原始列名
fwrite(
  data.frame(original_colnames = names(x), stringsAsFactors = FALSE),
  file.path(registry_dir, "46_primary_outcome_original_colnames.tsv"),
  sep = "\t"
)

# 自动映射（基于 PanScan 公共 summary 常见列）
mapped <- TRUE
missing_after_map <- character(0)

# 只在目标列不存在时才映射
if ("ID" %in% names(x) && !"SNP" %in% names(x)) x[, SNP := as.character(ID)]
if ("Allele2" %in% names(x) && !"effect_allele" %in% names(x)) x[, effect_allele := toupper(as.character(Allele2))]
if ("Allele1" %in% names(x) && !"other_allele" %in% names(x)) x[, other_allele := toupper(as.character(Allele1))]
if ("Effect" %in% names(x) && !"beta" %in% names(x)) x[, beta := as.numeric(Effect)]
if ("StdErr" %in% names(x) && !"se" %in% names(x)) x[, se := as.numeric(StdErr)]
if ("P-value" %in% names(x) && !"pval" %in% names(x)) x[, pval := as.numeric(`P-value`)]
if ("Freq1" %in% names(x) && !"eaf" %in% names(x)) x[, eaf := as.numeric(Freq1)]

# outcome 名称
if (!"outcome" %in% names(x)) x[, outcome := "Pancreatic cancer (PanScan/PanC4 public analysis)"]

required_cols <- c("SNP", "effect_allele", "other_allele", "beta", "se", "pval", "eaf", "outcome")
missing_after_map <- setdiff(required_cols, names(x))

attach_check <- data.frame(
  item = c(
    "selected_file_exists",
    "n_rows_raw",
    "n_cols_raw",
    "required_columns_ok",
    "missing_columns_after_mapping"
  ),
  value = c(
    file.exists(selected_file),
    nrow(x),
    ncol(x),
    length(missing_after_map) == 0,
    paste(missing_after_map, collapse = ",")
  ),
  stringsAsFactors = FALSE
)

fwrite(
  attach_check,
  file.path(registry_dir, "46_primary_outcome_attachment_check.tsv"),
  sep = "\t", na = "NA"
)

# 预览映射结果
preview_cols <- intersect(required_cols, names(x))
fwrite(
  x[1:min(20, nrow(x)), ..preview_cols],
  file.path(registry_dir, "46_primary_outcome_head_preview.tsv"),
  sep = "\t", na = "NA"
)

if (length(missing_after_map) > 0) {
  writeLines(
    c(
      "46 stopped after download + inspection because some required columns are still missing.",
      paste0("Created at: ", Sys.time()),
      paste0("Selected file: ", selected_file),
      paste0("Missing columns after mapping: ", paste(missing_after_map, collapse = ", "))
    ),
    file.path(log_dir, "46_primary_outcome_attachment_notes.txt")
  )
  
  cat("下载与检查完成，但标准化未完全成功。\n")
  cat("请把这些文件发给我：\n")
  cat("- 46_public_panscan_file_check.tsv\n")
  cat("- 46_selected_primary_outcome_file.tsv\n")
  cat("- 46_primary_outcome_original_colnames.tsv\n")
  cat("- 46_primary_outcome_head_preview.tsv\n")
  quit(save = "no")
}
