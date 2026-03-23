# =========================
# 26_gse57495_refinement_and_cross_external_summary.R
# 目的：
# 1. 汇总 GSE21501 与 GSE57495 的外部验证结果
# 2. 生成跨外部队列对比表
# 3. 固定 bulk 主线最终验证口径
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})
setwd("/Users/wmz_mac/Desktop/胰腺癌")
out_dir <- "04_bulk_analysis/05_publication_summary"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================
# 辅助函数
# ============================

#' 安全读取文件：检查存在性后读入，否则停止并报错
safe_fread <- function(path, label = NULL) {
  tag <- if (!is.null(label)) label else basename(path)
  if (!file.exists(path)) {
    stop(sprintf("[%s] 文件不存在: %s", tag, path))
  }
  dt <- fread(path)
  if (nrow(dt) == 0) {
    warning(sprintf("[%s] 文件为空 (0 行): %s", tag, path))
  }
  message(sprintf("[%s] 读取成功: %d 行 x %d 列", tag, nrow(dt), ncol(dt)))
  return(dt)
}

#' 安全提取标量值：返回第一个匹配值，无匹配时返回 NA 并发出警告
safe_scalar <- function(x, label = "unknown") {
  if (length(x) == 0 || all(is.na(x))) {
    warning(sprintf("safe_scalar [%s]: 取值为空，返回 NA", label))
    return(NA_real_)
  }
  return(x[1])
}

#' 从 metric/value 格式的表中安全提取某一 metric 的值
safe_metric <- function(dt, metric_name, label = "unknown") {
  idx <- which(dt$metric == metric_name)
  if (length(idx) == 0) {
    warning(sprintf("safe_metric [%s]: 未找到 metric='%s'，返回 NA", label, metric_name))
    return(NA_character_)
  }
  return(dt$value[idx[1]])
}

#' 检查数据表是否包含预期列名；缺失则停止
assert_cols <- function(dt, expected_cols, label = "unknown") {
  missing <- setdiff(expected_cols, names(dt))
  if (length(missing) > 0) {
    stop(sprintf("[%s] 缺少必需列: %s\n  现有列: %s",
                 label,
                 paste(missing, collapse = ", "),
                 paste(names(dt), collapse = ", ")))
  }
}

#' 从 timeROC 表中安全取出指定月份的 AUC
safe_auc <- function(roc_dt, months, label = "unknown") {
  assert_cols(roc_dt, c("time_months", "AUC"), label)
  idx <- which(roc_dt$time_months == months)
  if (length(idx) == 0) {
    warning(sprintf("safe_auc [%s]: 未找到 time_months=%d，返回 NA", label, months))
    return(NA_real_)
  }
  return(roc_dt$AUC[idx[1]])
}

# ============================
# 输入路径定义
# ============================
gse21501_score_file <- "04_bulk_analysis/04_external_validation/04_signature_validation/GSE21501/20_GSE21501_univcox_risk_score_confirm.tsv"
gse21501_group_file <- "04_bulk_analysis/04_external_validation/04_signature_validation/GSE21501/20_GSE21501_univcox_external_median_group.tsv"
gse21501_roc_file   <- "04_bulk_analysis/04_external_validation/04_signature_validation/GSE21501/20_GSE21501_timeROC_months.tsv"
gse21501_sum_file   <- "04_bulk_analysis/04_external_validation/04_signature_validation/GSE21501/20_GSE21501_refined_validation_summary.tsv"

gse57495_score_file  <- "04_bulk_analysis/04_external_validation/04_signature_validation/GSE57495/25_GSE57495_univcox_risk_score.tsv"
gse57495_group_file  <- "04_bulk_analysis/04_external_validation/04_signature_validation/GSE57495/25_GSE57495_univcox_risk_group_training_cutoff.tsv"
gse57495_group2_file <- "04_bulk_analysis/04_external_validation/04_signature_validation/GSE57495/25_GSE57495_univcox_risk_group_external_median.tsv"
gse57495_roc_file    <- "04_bulk_analysis/04_external_validation/04_signature_validation/GSE57495/25_GSE57495_timeROC_months.tsv"
gse57495_audit_file  <- "04_bulk_analysis/04_external_validation/04_signature_validation/GSE57495/25_GSE57495_validation_audit.tsv"

# ============================
# 读取所有输入（带安全检查）
# ============================
a_score  <- safe_fread(gse21501_score_file, "GSE21501_score")
a_group  <- safe_fread(gse21501_group_file, "GSE21501_group")
a_roc    <- safe_fread(gse21501_roc_file,   "GSE21501_roc")
a_sum    <- safe_fread(gse21501_sum_file,    "GSE21501_summary")

b_score  <- safe_fread(gse57495_score_file,  "GSE57495_score")
b_group  <- safe_fread(gse57495_group_file,  "GSE57495_group")
b_group2 <- safe_fread(gse57495_group2_file, "GSE57495_group2")
b_roc    <- safe_fread(gse57495_roc_file,    "GSE57495_roc")
b_audit  <- safe_fread(gse57495_audit_file,  "GSE57495_audit")

# ============================
# 列名验证
# ============================
assert_cols(a_score, c("estimate", "p.value"), "GSE21501_score")
assert_cols(a_group, c("estimate", "p.value"), "GSE21501_group")
assert_cols(a_roc,   c("time_months", "AUC"),  "GSE21501_roc")
assert_cols(a_sum,   c("metric", "value"),      "GSE21501_summary")

assert_cols(b_score,  c("estimate", "p.value"), "GSE57495_score")
assert_cols(b_group,  c("estimate", "p.value"), "GSE57495_group")
assert_cols(b_group2, c("estimate", "p.value"), "GSE57495_group2")
assert_cols(b_roc,    c("time_months", "AUC"),  "GSE57495_roc")
assert_cols(b_audit,  c("metric", "value"),     "GSE57495_audit")

# ============================
# 1. 外部队列汇总对比表
# ============================
summary_tab <- data.frame(
  cohort = c("GSE21501", "GSE57495"),
  
  n_samples = c(
    as.numeric(safe_metric(a_sum,   "n_external_samples",  "GSE21501_sum")),
    as.numeric(safe_metric(b_audit, "n_matched_samples",   "GSE57495_audit"))
  ),
  n_events = c(
    as.numeric(safe_metric(a_sum,   "n_external_events",   "GSE21501_sum")),
    as.numeric(safe_metric(b_audit, "n_events",            "GSE57495_audit"))
  ),
  
  continuous_HR = c(
    safe_scalar(a_score$estimate, "GSE21501_cont_HR"),
    safe_scalar(b_score$estimate, "GSE57495_cont_HR")
  ),
  continuous_p = c(
    safe_scalar(a_score$p.value,  "GSE21501_cont_p"),
    safe_scalar(b_score$p.value,  "GSE57495_cont_p")
  ),
  
  grouping_HR = c(
    safe_scalar(a_group$estimate, "GSE21501_grp_HR"),
    safe_scalar(b_group$estimate, "GSE57495_grp_HR")
  ),
  grouping_p = c(
    safe_scalar(a_group$p.value,  "GSE21501_grp_p"),
    safe_scalar(b_group$p.value,  "GSE57495_grp_p")
  ),
  
  auc_12m = c(
    safe_auc(a_roc, 12, "GSE21501"),
    safe_auc(b_roc, 12, "GSE57495")
  ),
  auc_24m = c(
    safe_auc(a_roc, 24, "GSE21501"),
    safe_auc(b_roc, 24, "GSE57495")
  ),
  auc_36m = c(
    safe_auc(a_roc, 36, "GSE21501"),
    safe_auc(b_roc, 36, "GSE57495")
  ),
  
  interpretation = c(
    "directionally consistent but weak external support",
    "externally supportive with significant stratification"
  ),
  
  stringsAsFactors = FALSE
)

fwrite(summary_tab,
       file.path(out_dir, "26_cross_external_validation_summary.tsv"),
       sep = "\t")
message("已输出: 26_cross_external_validation_summary.tsv")

# ============================
# 2. GSE57495 敏感性分析表
# ============================
gse57495_sens <- data.frame(
  analysis = c("training_cutoff", "external_median_cutoff", "continuous_score"),
  HR = c(
    safe_scalar(b_group$estimate,  "GSE57495_sens_train"),
    safe_scalar(b_group2$estimate, "GSE57495_sens_median"),
    safe_scalar(b_score$estimate,  "GSE57495_sens_cont")
  ),
  p_value = c(
    safe_scalar(b_group$p.value,  "GSE57495_sens_train_p"),
    safe_scalar(b_group2$p.value, "GSE57495_sens_median_p"),
    safe_scalar(b_score$p.value,  "GSE57495_sens_cont_p")
  ),
  stringsAsFactors = FALSE
)

fwrite(gse57495_sens,
       file.path(out_dir, "26_GSE57495_sensitivity_summary.tsv"),
       sep = "\t")
message("已输出: 26_GSE57495_sensitivity_summary.tsv")

# ============================
# 3. 最终口径声明
# ============================
final_statement <- data.frame(
  item = c(
    "bulk_mainline_status",
    "external_validation_overall_judgement",
    "recommended_next_priority"
  ),
  value = c(
    "completed",
    "one weak-support external cohort plus one supportive external cohort",
    "start MR analysis"
  ),
  stringsAsFactors = FALSE
)

fwrite(final_statement,
       file.path(out_dir, "26_bulk_stage_final_judgement.tsv"),
       sep = "\t")
message("已输出: 26_bulk_stage_final_judgement.tsv")

# ============================
# 4. 运行环境与版本记录
# ============================
run_info <- data.frame(
  item = c("script", "run_time", "R_version", "data.table_version", "dplyr_version"),
  value = c(
    "26_gse57495_refinement_and_cross_external_summary.R",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("data.table")),
    as.character(packageVersion("dplyr"))
  ),
  stringsAsFactors = FALSE
)

fwrite(run_info,
       file.path(out_dir, "26_run_environment.tsv"),
       sep = "\t")
message("已输出: 26_run_environment.tsv")

cat("\n============================\n")
cat("26_gse57495_refinement_and_cross_external_summary.R 全部完成。\n")
cat("============================\n")