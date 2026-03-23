# ============================================================
# 00_build_survival_ready_dataset.R
# 用途：
# 1. 固定项目根目录
# 2. 读取 strict-PDAC cohort、MCD subtype、表达矩阵
# 3. 统一样本ID
# 4. 构建 survival-ready 数据表
# 5. 导出审计表、survival-ready 表、对齐后的表达矩阵
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# =========================
# 0. 固定项目根目录
# =========================
project_root <- "/Users/wmz_mac/Desktop/胰腺癌"

if (!dir.exists(project_root)) {
  stop("项目根目录不存在: ", project_root)
}

setwd(project_root)
cat("Current working directory:", getwd(), "\n\n")

# =========================
# 1. 文件路径
# =========================
cohort_file  <- "02_processed_data/TCGA_PAAD/tcga_paad_cohort_master_strict_pdac_ordered.tsv"
subtype_file <- "04_bulk_analysis/02_clustering/tcga_paad_mcd_subtype_phenotype.tsv"
expr_file    <- "02_processed_data/TCGA_PAAD/tcga_paad_expr_strict_pdac_symbol_counts_filtered.rds"

out_dir <- "04_bulk_analysis/03_survival_model"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# 2. 检查输入文件
# =========================
files_to_check <- c(cohort_file, subtype_file, expr_file)
check_res <- data.frame(
  file = files_to_check,
  exists = file.exists(files_to_check),
  stringsAsFactors = FALSE
)

print(check_res)
cat("\n")

if (!all(file.exists(files_to_check))) {
  stop("存在输入文件缺失，请先修正路径。")
}

# =========================
# 3. 读取数据
# =========================
cohort <- fread(cohort_file)
subtype <- fread(subtype_file)
expr <- readRDS(expr_file)

cat("Cohort dim :", dim(cohort), "\n")
cat("Subtype dim:", dim(subtype), "\n")
cat("Expr dim   :", dim(expr), "\n\n")

cat("Cohort columns:\n")
print(colnames(cohort))
cat("\nSubtype columns:\n")
print(colnames(subtype))
cat("\n")

# =========================
# 4. 工具函数
# =========================
pick_first_existing <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

to_sample16 <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  substr(x, 1, 16)
}

to_patient12 <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  substr(x, 1, 12)
}

normalize_event <- function(x) {
  if (is.numeric(x)) {
    return(as.numeric(x))
  }
  
  x_chr <- as.character(x)
  x_chr <- trimws(tolower(x_chr))
  
  out <- rep(NA_real_, length(x_chr))
  out[x_chr %in% c("1", "dead", "deceased", "event", "yes", "true")] <- 1
  out[x_chr %in% c("0", "alive", "living", "censored", "no", "false")] <- 0
  
  suppressWarnings({
    x_num <- as.numeric(x_chr)
  })
  out[is.na(out) & !is.na(x_num)] <- x_num[is.na(out) & !is.na(x_num)]
  
  as.numeric(out)
}

# =========================
# 5. cohort 标准化
# =========================
sample_col_cohort <- pick_first_existing(cohort, c(
  "sample_barcode", "sample_id", "SampleID", "sample", "submitter_id", "barcode", "bcr_patient_barcode"
))

patient_col_cohort <- pick_first_existing(cohort, c(
  "patient_id", "case_id", "patient_barcode", "bcr_patient_barcode"
))

os_time_col <- pick_first_existing(cohort, c(
  "OS.time", "OS_time", "os_time", "overall_survival_time", "days_to_death_or_last_followup",
  "days_to_last_follow_up_or_death", "days_to_event", "time"
))

os_event_col <- pick_first_existing(cohort, c(
  "OS.status", "OS_event", "os_event", "overall_survival_event", "vital_status_binary",
  "event", "status", "vital_status_final"
))

if (is.na(sample_col_cohort)) stop("cohort 中未找到样本ID列")
if (is.na(os_time_col)) stop("cohort 中未找到 OS 时间列")
if (is.na(os_event_col)) stop("cohort 中未找到 OS 结局列")

cohort2 <- cohort %>%
  dplyr::mutate(
    sample_id_raw = .data[[sample_col_cohort]],
    sample_id_16  = to_sample16(.data[[sample_col_cohort]]),
    patient_id_12 = if (!is.na(patient_col_cohort)) to_patient12(.data[[patient_col_cohort]]) else to_patient12(.data[[sample_col_cohort]]),
    OS_time       = suppressWarnings(as.numeric(.data[[os_time_col]])),
    OS_event_raw  = .data[[os_event_col]],
    OS_event      = normalize_event(.data[[os_event_col]])
  )

cat("Cohort standardized preview:\n")
print(head(as.data.frame(cohort2[, c("sample_id_raw", "sample_id_16", "patient_id_12", "OS_time", "OS_event_raw", "OS_event")])))
cat("\n")

# =========================
# 6. subtype 标准化
# =========================
sample_col_subtype <- pick_first_existing(subtype, c(
  "sample_id", "SampleID", "sample", "barcode", "submitter_id", "sample_barcode"
))

patient_col_subtype <- pick_first_existing(subtype, c(
  "patient_id", "case_id", "patient_barcode", "bcr_patient_barcode"
))

subtype_col <- pick_first_existing(subtype, c(
  "mcd_subtype", "MCD_subtype", "subtype", "cluster_label", "Subtype", "cluster", "mcd_cluster"
))

if (is.na(sample_col_subtype) && is.na(patient_col_subtype)) {
  stop("subtype 中未找到样本或患者ID列")
}
if (is.na(subtype_col)) {
  stop("subtype 中未找到亚型列，可用列: ", paste(colnames(subtype), collapse = ", "))
}

subtype2 <- subtype %>%
  dplyr::mutate(
    subtype_id_raw = if (!is.na(sample_col_subtype)) .data[[sample_col_subtype]] else .data[[patient_col_subtype]],
    sample_id_16   = if (!is.na(sample_col_subtype)) to_sample16(.data[[sample_col_subtype]]) else NA_character_,
    patient_id_12  = if (!is.na(patient_col_subtype)) to_patient12(.data[[patient_col_subtype]]) else to_patient12(subtype_id_raw),
    mcd_subtype    = .data[[subtype_col]]
  ) %>%
  dplyr::select(subtype_id_raw, sample_id_16, patient_id_12, mcd_subtype)

cat("Subtype standardized preview:\n")
print(head(as.data.frame(subtype2)))
cat("\n")

# =========================
# 7. 表达矩阵样本ID识别
# =========================
expr_sample_ids_raw <- colnames(expr)
expr_sample_ids_16  <- to_sample16(expr_sample_ids_raw)
expr_patient_ids_12 <- to_patient12(expr_sample_ids_raw)

expr_map <- data.frame(
  expr_colname = expr_sample_ids_raw,
  sample_id_16 = expr_sample_ids_16,
  patient_id_12 = expr_patient_ids_12,
  stringsAsFactors = FALSE
)

cat("Expression sample preview:\n")
print(head(expr_map))
cat("\n")

# =========================
# 8. 判断合并层级
# =========================
n_match_sample_cohort_expr   <- sum(cohort2$sample_id_16 %in% expr_map$sample_id_16)
n_match_sample_subtype_expr  <- sum(subtype2$sample_id_16 %in% expr_map$sample_id_16, na.rm = TRUE)

n_match_patient_cohort_expr  <- sum(cohort2$patient_id_12 %in% expr_map$patient_id_12)
n_match_patient_subtype_expr <- sum(subtype2$patient_id_12 %in% expr_map$patient_id_12)

cat("Match summary:\n")
cat("sample-level  cohort~expr  :", n_match_sample_cohort_expr, "\n")
cat("sample-level  subtype~expr :", n_match_sample_subtype_expr, "\n")
cat("patient-level cohort~expr  :", n_match_patient_cohort_expr, "\n")
cat("patient-level subtype~expr :", n_match_patient_subtype_expr, "\n\n")

use_sample_level <- (n_match_sample_cohort_expr > 0) && (n_match_sample_subtype_expr > 0)

if (use_sample_level) {
  cat("Using sample-level matching (16-char TCGA barcode).\n\n")
  
  dat <- cohort2 %>%
    dplyr::inner_join(
      subtype2 %>% dplyr::select(sample_id_16, mcd_subtype),
      by = "sample_id_16"
    ) %>%
    dplyr::filter(sample_id_16 %in% expr_map$sample_id_16) %>%
    dplyr::left_join(
      expr_map %>% dplyr::select(sample_id_16, expr_colname),
      by = "sample_id_16"
    ) %>%
    dplyr::mutate(sample_id = sample_id_16)
  
} else {
  cat("Using patient-level matching (12-char TCGA patient barcode).\n\n")
  
  dat <- cohort2 %>%
    dplyr::inner_join(
      subtype2 %>% dplyr::select(patient_id_12, mcd_subtype),
      by = "patient_id_12"
    ) %>%
    dplyr::filter(patient_id_12 %in% expr_map$patient_id_12) %>%
    dplyr::left_join(
      expr_map %>% dplyr::select(patient_id_12, expr_colname),
      by = "patient_id_12"
    ) %>%
    dplyr::mutate(sample_id = patient_id_12)
}

dat <- dat %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)

# =========================
# 9. 生存审计
# =========================
audit <- data.frame(
  metric = c(
    "n_cohort_raw",
    "n_subtype_raw",
    "n_expr_samples_raw",
    "n_merged_before_filter",
    "missing_OS_time",
    "missing_OS_event",
    "nonpositive_OS_time",
    "invalid_OS_event",
    "duplicated_sample_id",
    "n_subtype_missing"
  ),
  value = c(
    nrow(cohort),
    nrow(subtype),
    ncol(expr),
    nrow(dat),
    sum(is.na(dat$OS_time)),
    sum(is.na(dat$OS_event)),
    sum(dat$OS_time <= 0, na.rm = TRUE),
    sum(!dat$OS_event %in% c(0, 1), na.rm = TRUE),
    sum(duplicated(dat$sample_id)),
    sum(is.na(dat$mcd_subtype))
  ),
  stringsAsFactors = FALSE
)

fwrite(audit, file.path(out_dir, "00_survival_endpoint_audit.tsv"), sep = "\t")

cat("Audit table:\n")
print(audit)
cat("\n")

# =========================
# 10. 生成 survival-ready 数据
# =========================
surv_ready <- dat %>%
  dplyr::filter(!is.na(OS_time), !is.na(OS_event), !is.na(mcd_subtype)) %>%
  dplyr::filter(OS_time > 0, OS_event %in% c(0, 1)) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)

expr_keep <- surv_ready$expr_colname
expr_keep <- expr_keep[!is.na(expr_keep)]
expr_keep <- intersect(expr_keep, colnames(expr))

surv_ready <- surv_ready %>%
  dplyr::filter(expr_colname %in% expr_keep) %>%
  dplyr::arrange(match(expr_colname, expr_keep))

expr_aligned <- expr[, expr_keep, drop = FALSE]

stopifnot(nrow(surv_ready) == ncol(expr_aligned))
stopifnot(all(surv_ready$expr_colname == colnames(expr_aligned)))

# =========================
# 11. 输出结果
# =========================
fwrite(surv_ready, file.path(out_dir, "01_tcga_strict_pdac_survival_ready.tsv"), sep = "\t")
saveRDS(expr_aligned, file.path(out_dir, "01_tcga_strict_pdac_expr_aligned.rds"))

cat("========================================\n")
cat("Done.\n")
cat("Survival-ready n =", nrow(surv_ready), "\n\n")

cat("Subtype distribution:\n")
print(table(surv_ready$mcd_subtype))
cat("\n")

cat("OS event distribution:\n")
print(table(surv_ready$OS_event))
cat("\n")

cat("Aligned expression dim:\n")
print(dim(expr_aligned))
cat("\n")

cat("Output files:\n")
cat(file.path(out_dir, "00_survival_endpoint_audit.tsv"), "\n")
cat(file.path(out_dir, "01_tcga_strict_pdac_survival_ready.tsv"), "\n")
cat(file.path(out_dir, "01_tcga_strict_pdac_expr_aligned.rds"), "\n")
cat("========================================\n")