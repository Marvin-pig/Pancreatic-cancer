# =========================
# 35_download_exposure_support_and_prepare_eqtlgen_manifest.R
# 目的：
# 1. 下载 GTEx V8 Whole_Blood / Pancreas 的 SMR BESD 压缩包（支持性资源）
# 2. 生成 batch1 eQTLGen 手动下载清单（主 exposure 仍按你的既定方案保留为 eQTLGen）
# 3. 为后续 36 正式接入 exposure 做准备
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(jsonlite)
})

# ---------- 路径 ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "05_MR"

dir.create(file.path(root_dir, "00_registry"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "01_raw_data/exposure/eQTLGen"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "01_raw_data/exposure/GTEx_pancreas"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "01_raw_data/exposure/GTEx_whole_blood"), recursive = TRUE, showWarnings = FALSE)
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

batch1 <- batch1 %>%
  mutate(gene = as.character(gene)) %>%
  distinct(gene, .keep_all = TRUE)

# ---------- 35A. 自动下载 GTEx V8 支持性资源 ----------
# 这两个链接来自 Yang Lab 的 GTEx V8 cis-eQTL summary 页面：
# Whole_Blood.zip 与 Pancreas.zip 当前公开列在该页面上
# 若服务器临时不可达，代码会保留失败状态，不影响后续清单生成

download_plan <- data.frame(
  dataset = c("GTEx_V8_Whole_Blood", "GTEx_V8_Pancreas"),
  url = c(
    "https://yanglab.westlake.edu.cn/data/SMR/GTEx_V8_cis_eqtl_summary/Whole_Blood.zip",
    "https://yanglab.westlake.edu.cn/data/SMR/GTEx_V8_cis_eqtl_summary/Pancreas.zip"
  ),
  dest = c(
    file.path(root_dir, "01_raw_data/exposure/GTEx_whole_blood/Whole_Blood.zip"),
    file.path(root_dir, "01_raw_data/exposure/GTEx_pancreas/Pancreas.zip")
  ),
  downloaded = FALSE,
  size_bytes = NA_real_,
  stringsAsFactors = FALSE
)

safe_download <- function(url, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  ok <- FALSE
  msg <- NA_character_
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

download_log <- lapply(seq_len(nrow(download_plan)), function(i) {
  res <- safe_download(download_plan$url[i], download_plan$dest[i])
  download_plan$downloaded[i] <<- isTRUE(res$ok)
  if (file.exists(download_plan$dest[i])) {
    download_plan$size_bytes[i] <<- file.info(download_plan$dest[i])$size
  }
  data.frame(
    dataset = download_plan$dataset[i],
    ok = res$ok,
    message = ifelse(is.na(res$msg), "", res$msg),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

fwrite(
  download_plan,
  file.path(root_dir, "00_registry/35_support_download_status.tsv"),
  sep = "\t"
)

fwrite(
  download_log,
  file.path(root_dir, "00_registry/35_support_download_log.tsv"),
  sep = "\t"
)

# ---------- 35B. 生成 eQTLGen 主 exposure 手动下载/整理清单 ----------
# 保持你的既定主线：primary exposure 仍是 eQTLGen
# 这里不擅自改为 GTEx whole blood，以避免破坏你前面 30/31 的判定逻辑

eqtlgen_manifest <- batch1 %>%
  transmute(
    gene,
    eqtlgen_expected_tsv = file.path(root_dir, "01_raw_data/exposure/eQTLGen", paste0(gene, "_cis_eqtl.tsv")),
    gtex_pancreas_expected_zip = file.path(root_dir, "01_raw_data/exposure/GTEx_pancreas", "Pancreas.zip"),
    gtex_whole_blood_expected_zip = file.path(root_dir, "01_raw_data/exposure/GTEx_whole_blood", "Whole_Blood.zip"),
    eqtlgen_present = file.exists(eqtlgen_expected_tsv),
    priority = if ("run_priority" %in% colnames(batch1)) batch1$run_priority else "UNKNOWN"
  ) %>%
  arrange(desc(priority), gene)

fwrite(
  eqtlgen_manifest,
  file.path(root_dir, "00_registry/35_eqtlgen_batch1_manifest.tsv"),
  sep = "\t"
)

# ---------- 35C. 生成后续 outcome 提取所需的 SNP 占位文件 ----------
# 现在先预留，真正内容要等 36 从 exposure 文件中提取 instruments 后再填
writeLines(
  c("SNP", "# 36 会在 exposure 文件接入后自动填充"),
  file.path(root_dir, "02_processed_data/batch1_queries/35_bbj_outcome_snp_placeholder.tsv")
)

# ---------- 35D. 总状态 ----------
status_tab <- data.frame(
  item = c(
    "GTEx_whole_blood_zip_ready",
    "GTEx_pancreas_zip_ready",
    "any_support_resource_ready",
    "any_eqtlgen_primary_exposure_ready",
    "recommended_next_action"
  ),
  value = c(
    file.exists(file.path(root_dir, "01_raw_data/exposure/GTEx_whole_blood/Whole_Blood.zip")),
    file.exists(file.path(root_dir, "01_raw_data/exposure/GTEx_pancreas/Pancreas.zip")),
    any(file.exists(download_plan$dest)),
    any(eqtlgen_manifest$eqtlgen_present),
    ifelse(
      any(eqtlgen_manifest$eqtlgen_present),
      "Proceed to 36: standardize and audit exposure files, then extract BBJ outcomes",
      "Download or prepare at least one batch1 eQTLGen gene-level cis-eQTL TSV, then proceed to 36"
    )
  ),
  stringsAsFactors = FALSE
)

fwrite(
  status_tab,
  file.path(root_dir, "00_registry/35_status_summary.tsv"),
  sep = "\t"
)

# ---------- 35E. 日志 ----------
notes <- c(
  "=== Code 35 Summary ===",
  paste0("Created at: ", Sys.time()),
  "",
  "[What this script did]",
  "1. Attempted to download GTEx V8 Whole_Blood.zip",
  "2. Attempted to download GTEx V8 Pancreas.zip",
  "3. Generated batch1 eQTLGen manifest",
  "4. Created BBJ outcome SNP placeholder file",
  "",
  "[Important rule]",
  "Primary exposure remains eQTLGen (do not silently replace with GTEx).",
  "GTEx files are supporting resources only at this stage.",
  "",
  "[Next step]",
  status_tab$value[status_tab$item == "recommended_next_action"]
)

writeLines(
  notes,
  file.path(root_dir, "06_logs/35_download_and_manifest_notes.txt")
)

cat("35_download_exposure_support_and_prepare_eqtlgen_manifest.R finished successfully.\n")
cat("Generated files:\n")
cat("- 00_registry/35_support_download_status.tsv\n")
cat("- 00_registry/35_support_download_log.tsv\n")
cat("- 00_registry/35_eqtlgen_batch1_manifest.tsv\n")
cat("- 00_registry/35_status_summary.tsv\n")
cat("- 02_processed_data/batch1_queries/35_bbj_outcome_snp_placeholder.tsv\n")
cat("- 06_logs/35_download_and_manifest_notes.txt\n")