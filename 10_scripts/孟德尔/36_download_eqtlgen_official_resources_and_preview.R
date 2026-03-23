# =========================
# 36_download_eqtlgen_official_resources_and_preview.R
# 目的：
# 1. 下载 eQTLGen 官方 cis-eQTL资源
# 2. 生成本地下载状态表
# 3. 自动读取表头/前几行，检查 batch1 基因是否能匹配
# 4. 为后续 37 拆分 gene-level exposure 文件做准备
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ---------- 路径 ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "05_MR"

dir.create(file.path(root_dir, "00_registry"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "01_raw_data/exposure/eQTLGen"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "02_processed_data/batch1_queries"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "06_logs"), recursive = TRUE, showWarnings = FALSE)

batch1_file <- file.path(root_dir, "02_processed_data/batch1_queries/29_batch1_priority_genes.tsv")
if (!file.exists(batch1_file)) {
  stop("未找到 batch1 基因文件：", batch1_file)
}
batch1 <- fread(batch1_file)

if (!"gene" %in% colnames(batch1)) {
  stop("29_batch1_priority_genes.tsv 缺少 gene 列。")
}
batch1_genes <- unique(as.character(batch1$gene))

# ---------- eQTLGen 官方下载链接 ----------
# 来自 eQTLGen phase I cis-eQTL 页面
eqtlgen_dir <- file.path(root_dir, "01_raw_data/exposure/eQTLGen")

download_plan <- data.frame(
  item = c(
    "significant_cis_eqtls",
    "allele_frequencies",
    "readme_cis"
  ),
  url = c(
    "https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz",
    "https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/eQTLGen_AF.txt.gz",
    "https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/README_cis"
  ),
  dest = c(
    file.path(eqtlgen_dir, "36_eqtlgen_significant_cis_eqtls.txt.gz"),
    file.path(eqtlgen_dir, "36_eqtlgen_allele_frequencies.txt.gz"),
    file.path(eqtlgen_dir, "36_eqtlgen_README_cis.txt")
  ),
  stringsAsFactors = FALSE
)

# 可选：完整 full summary 很大，先默认不下
download_full_summary <- FALSE
full_summary_url <- "https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
full_summary_dest <- file.path(eqtlgen_dir, "36_eqtlgen_full_cis_summary.txt.gz")

if (download_full_summary) {
  download_plan <- bind_rows(
    download_plan,
    data.frame(
      item = "full_cis_summary",
      url = full_summary_url,
      dest = full_summary_dest,
      stringsAsFactors = FALSE
    )
  )
}

safe_download <- function(url, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  ok <- FALSE
  msg <- ""
  if (!file.exists(dest)) {
    tryCatch({
      download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
      ok <- file.exists(dest)
    }, error = function(e) {
      msg <<- conditionMessage(e)
    })
  } else {
    ok <- TRUE
    msg <- "already_exists"
  }
  list(ok = ok, msg = msg)
}

download_status <- lapply(seq_len(nrow(download_plan)), function(i) {
  res <- safe_download(download_plan$url[i], download_plan$dest[i])
  data.frame(
    item = download_plan$item[i],
    url = download_plan$url[i],
    dest = download_plan$dest[i],
    exists = file.exists(download_plan$dest[i]),
    ok = isTRUE(res$ok),
    message = res$msg,
    size_bytes = ifelse(file.exists(download_plan$dest[i]),
                        file.info(download_plan$dest[i])$size, NA),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

fwrite(
  download_status,
  file.path(root_dir, "00_registry/36_eqtlgen_download_status.tsv"),
  sep = "\t"
)

# ---------- 读取 significant cis-eQTL 文件表头/预览 ----------
sig_file <- file.path(eqtlgen_dir, "36_eqtlgen_significant_cis_eqtls.txt.gz")

preview_ok <- FALSE
header_names <- character(0)
preview_df <- NULL
gene_col_detected <- NA_character_

if (file.exists(sig_file)) {
  # 读表头
  header_names <- tryCatch(
    names(fread(sig_file, nrows = 0)),
    error = function(e) character(0)
  )
  
  # 读前几行
  preview_df <- tryCatch(
    fread(sig_file, nrows = 10),
    error = function(e) NULL
  )
  
  preview_ok <- !is.null(preview_df)
  
  # 自动识别基因列
  candidate_gene_cols <- header_names[
    str_detect(header_names, regex("gene|symbol|probe|gene_id|genesymbol", ignore_case = TRUE))
  ]
  if (length(candidate_gene_cols) > 0) {
    gene_col_detected <- candidate_gene_cols[1]
  }
}

# 保存表头
header_tab <- data.frame(
  column_name = header_names,
  stringsAsFactors = FALSE
)
fwrite(
  header_tab,
  file.path(root_dir, "00_registry/36_eqtlgen_significant_header.tsv"),
  sep = "\t"
)

# 保存预览
if (!is.null(preview_df)) {
  fwrite(
    preview_df,
    file.path(root_dir, "00_registry/36_eqtlgen_significant_preview.tsv"),
    sep = "\t"
  )
}

# ---------- batch1 基因命中检查 ----------
batch1_hit_tab <- data.frame(
  gene = batch1_genes,
  matched_in_preview = FALSE,
  stringsAsFactors = FALSE
)

if (!is.na(gene_col_detected) && !is.null(preview_df) && gene_col_detected %in% names(preview_df)) {
  vals <- unique(as.character(preview_df[[gene_col_detected]]))
  batch1_hit_tab$matched_in_preview <- batch1_hit_tab$gene %in% vals
}

fwrite(
  batch1_hit_tab,
  file.path(root_dir, "00_registry/36_batch1_gene_hit_in_preview.tsv"),
  sep = "\t"
)

# ---------- 总状态 ----------
status_summary <- data.frame(
  item = c(
    "significant_file_downloaded",
    "allele_freq_file_downloaded",
    "readme_downloaded",
    "significant_preview_readable",
    "gene_col_detected",
    "recommended_next_action"
  ),
  value = c(
    file.exists(file.path(eqtlgen_dir, "36_eqtlgen_significant_cis_eqtls.txt.gz")),
    file.exists(file.path(eqtlgen_dir, "36_eqtlgen_allele_frequencies.txt.gz")),
    file.exists(file.path(eqtlgen_dir, "36_eqtlgen_README_cis.txt")),
    preview_ok,
    ifelse(is.na(gene_col_detected), "NA", gene_col_detected),
    ifelse(
      preview_ok,
      "Proceed to 37: inspect column schema and split official eQTLGen file into gene-level candidate tables",
      "Check download integrity / file format before proceeding"
    )
  ),
  stringsAsFactors = FALSE
)

fwrite(
  status_summary,
  file.path(root_dir, "00_registry/36_eqtlgen_status_summary.tsv"),
  sep = "\t"
)

# ---------- 日志 ----------
notes <- c(
  "=== Code 36 Summary ===",
  paste0("Created at: ", Sys.time()),
  "",
  "[Downloaded resources]",
  paste(download_status$item, collapse = ", "),
  "",
  "[Significant file preview readable]",
  as.character(preview_ok),
  "",
  "[Detected gene column]",
  ifelse(is.na(gene_col_detected), "NA", gene_col_detected),
  "",
  "[Next step]",
  status_summary$value[status_summary$item == "recommended_next_action"]
)

writeLines(
  notes,
  file.path(root_dir, "06_logs/36_eqtlgen_download_and_preview_notes.txt")
)

cat("36_download_eqtlgen_official_resources_and_preview.R finished successfully.\n")
cat("Generated files:\n")
cat("- 00_registry/36_eqtlgen_download_status.tsv\n")
cat("- 00_registry/36_eqtlgen_significant_header.tsv\n")
cat("- 00_registry/36_eqtlgen_significant_preview.tsv (if readable)\n")
cat("- 00_registry/36_batch1_gene_hit_in_preview.tsv\n")
cat("- 00_registry/36_eqtlgen_status_summary.tsv\n")
cat("- 06_logs/36_eqtlgen_download_and_preview_notes.txt\n")