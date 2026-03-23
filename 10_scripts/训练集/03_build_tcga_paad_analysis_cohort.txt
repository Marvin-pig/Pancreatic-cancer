# =========================================================
# 03_build_tcga_paad_analysis_cohort.R
# 作用：
# 1. 读取 TCGA-PAAD coldata 和 sample-clinical map
# 2. 构建 analysis-ready cohort
# 3. 标准化 OS 终点
# 4. 生成严格 PDAC 队列
# =========================================================

project_root <- "/Users/wmz_mac/Desktop/胰腺癌"
options(stringsAsFactors = FALSE)
library(data.table)

raw_root  <- file.path(project_root, "00_raw_data", "TCGA_PAAD")
clin_dir  <- file.path(raw_root, "clinical")
proc_dir  <- file.path(project_root, "02_processed_data", "TCGA_PAAD")
log_dir   <- file.path(project_root, "11_logs")

dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir,  recursive = TRUE, showWarnings = FALSE)

# 同时写入日志文件和控制台
log_file <- file.path(log_dir, paste0("03_cohort_build_", Sys.Date(), ".log"))
con <- file(log_file, open = "wt")
sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

cat("========================================\n")
cat("TCGA-PAAD Cohort Build\n")
cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# -------------------------
# 辅助函数
# -------------------------

# 优先取 .y 列，其次 .x 列；两者均缺时警告并返回 NA 向量
pick_col <- function(df, y_col, x_col = NULL) {
  if (!is.null(y_col) && y_col %in% colnames(df)) return(df[[y_col]])
  if (!is.null(x_col) && x_col %in% colnames(df)) return(df[[x_col]])
  warning(sprintf("pick_col: 未找到列 '%s' 或 '%s'，该列将全部为 NA",
                  y_col, if (is.null(x_col)) "<none>" else x_col))
  return(rep(NA_character_, nrow(df)))
}

# -------------------------
# 读取数据
# -------------------------
coldata_file <- file.path(clin_dir, "tcga_paad_coldata_from_se.tsv")
map_file     <- file.path(clin_dir, "tcga_paad_sample_clinical_map.tsv")

coldata <- fread(coldata_file, sep = "\t", header = TRUE, data.table = FALSE)
map_df  <- fread(map_file,     sep = "\t", header = TRUE, data.table = FALSE)

cat(sprintf("读取 coldata：%d 行 × %d 列\n", nrow(coldata), ncol(coldata)))
cat(sprintf("读取 map_df ：%d 行 × %d 列\n", nrow(map_df),  ncol(map_df)))

# -------------------------
# 1. 只保留 TP 样本
# -------------------------
if ("shortLetterCode" %in% colnames(map_df)) {
  cohort     <- subset(map_df, shortLetterCode == "TP")
  n_tp_raw   <- sum(map_df$shortLetterCode == "TP", na.rm = TRUE)
} else if ("sample_type" %in% colnames(map_df)) {
  cohort     <- subset(map_df, sample_type %in% c("Primary Tumor", "TP"))
  n_tp_raw   <- sum(map_df$sample_type %in% c("Primary Tumor", "TP"), na.rm = TRUE)
} else {
  stop("未找到样本类型列 shortLetterCode / sample_type，请检查输入文件")
}
cat(sprintf("\n[1] TP 样本筛选：原始 %d 条 → 保留 %d 条\n", nrow(map_df), nrow(cohort)))

# -------------------------
# 2. 统一关键临床列
#    优先用 .y（clinical_raw）列，其次 .x
# -------------------------
cohort$primary_diagnosis_final <- pick_col(cohort, "primary_diagnosis.y",        "primary_diagnosis.x")
cohort$ajcc_stage_final        <- pick_col(cohort, "ajcc_pathologic_stage.y",    "ajcc_pathologic_stage.x")
cohort$tumor_grade_final       <- pick_col(cohort, "tumor_grade.y",              "tumor_grade.x")
cohort$vital_status_final      <- pick_col(cohort, "vital_status.y",             "vital_status.x")
cohort$days_to_death_final     <- suppressWarnings(
  as.numeric(pick_col(cohort, "days_to_death.y",           "days_to_death.x")))
cohort$days_to_followup_final  <- suppressWarnings(
  as.numeric(pick_col(cohort, "days_to_last_follow_up.y",  "days_to_last_follow_up.x")))
cohort$age_at_diagnosis_days   <- suppressWarnings(
  as.numeric(pick_col(cohort, "age_at_diagnosis.y",        "age_at_diagnosis.x")))
cohort$gender_final            <- pick_col(cohort, "gender.y",                   "gender.x")

# TCGA age_at_diagnosis 单位为天（出生至诊断），换算为年
cohort$age_at_diagnosis_years <- cohort$age_at_diagnosis_days / 365.25

# -------------------------
# 3. 构建 OS
#    先确定 OS.status，再用状态决定 OS.time，避免 NA 边界歧义
# -------------------------
cohort$OS.status <- ifelse(cohort$vital_status_final == "Dead",  1L,
                    ifelse(cohort$vital_status_final == "Alive", 0L, NA_integer_))

# 死亡者取 days_to_death；存活者取 days_to_last_follow_up
cohort$OS.time <- ifelse(cohort$OS.status == 1L,
                         cohort$days_to_death_final,
                         cohort$days_to_followup_final)

# -------------------------
# 4. 基本清洗
# -------------------------
# 去除无 sample_barcode 及重复 barcode
cohort <- cohort[!is.na(cohort$sample_barcode), ]
cohort <- cohort[!duplicated(cohort$sample_barcode), ]

# patient_id 去重：保留 OS 信息最完整的记录
if (any(duplicated(cohort$patient_id))) {
  n_before <- nrow(cohort)
  cohort$os_info_score <- (!is.na(cohort$OS.time)) + (!is.na(cohort$OS.status))
  cohort <- cohort[order(cohort$patient_id, -cohort$os_info_score), ]
  cohort <- cohort[!duplicated(cohort$patient_id), ]
  cohort$os_info_score <- NULL
  cat(sprintf("[4] patient_id 去重：%d → %d（移除 %d 条）\n",
              n_before, nrow(cohort), n_before - nrow(cohort)))
}

# 去除 OS.status 缺失的样本（vital_status 不明确）
n_before <- nrow(cohort)
cohort <- cohort[!is.na(cohort$OS.status), ]
cat(sprintf("[4] 去除 OS.status 缺失：%d → %d\n", n_before, nrow(cohort)))

# 去除 OS.time 缺失或 ≤ 0 的样本
n_before <- nrow(cohort)
cohort <- cohort[!is.na(cohort$OS.time) & cohort$OS.time > 0, ]
cat(sprintf("[4] 去除 OS.time 无效（NA 或 ≤0）：%d → %d\n", n_before, nrow(cohort)))

# OS.time 值域提示（辅助发现数据异常，如录入错误导致超长随访）
cat(sprintf("[4] OS.time 分布（天）：min=%g, median=%g, max=%g\n",
            min(cohort$OS.time), median(cohort$OS.time), max(cohort$OS.time)))
if (max(cohort$OS.time) > 10950) {  # 30 年
  warning("存在 OS.time > 10950 天（30年）的记录，请核查是否为数据录入错误")
}

# -------------------------
# 5. 诊断分布审计
# -------------------------
diag_table <- as.data.frame(
  sort(table(cohort$primary_diagnosis_final), decreasing = TRUE),
  stringsAsFactors = FALSE
)
colnames(diag_table) <- c("primary_diagnosis_final", "n")
write.table(
  diag_table,
  file      = file.path(proc_dir, "tcga_paad_primary_diagnosis_distribution.tsv"),
  sep       = "\t", quote = FALSE, row.names = FALSE
)
cat("\n[5] 诊断分布（Top 10）：\n")
print(head(diag_table, 10), row.names = FALSE)

# -------------------------
# 6. 严格 PDAC 队列
# -------------------------
pdac_keep_patterns <- c(
  "Infiltrating duct carcinoma",
  "Adenocarcinoma, NOS",
  "Mucinous adenocarcinoma",
  "Adenocarcinoma with mixed subtypes"
)

# 合并为正则，统一大小写不敏感匹配
pdac_pattern <- paste(pdac_keep_patterns, collapse = "|")

strict_pdac <- cohort[
  !is.na(cohort$primary_diagnosis_final) &
    grepl(pdac_pattern, cohort$primary_diagnosis_final, ignore.case = TRUE),
]

# 排除明确不纳入的病理（与纳入保持相同策略：ignore.case = TRUE）
strict_pdac <- strict_pdac[
  !grepl("Neuroendocrine",  strict_pdac$primary_diagnosis_final, ignore.case = TRUE) &
  !grepl("undifferentiated",strict_pdac$primary_diagnosis_final, ignore.case = TRUE),
]

cat(sprintf("\n[6] 严格 PDAC 筛选：all-TP %d → strict_pdac %d\n",
            nrow(cohort), nrow(strict_pdac)))

# -------------------------
# 7. 导出
# -------------------------
keep_cols <- c(
  "sample_barcode", "patient_id",
  "primary_diagnosis_final",
  "ajcc_stage_final",
  "tumor_grade_final",
  "age_at_diagnosis_days",    # 原始单位（天）
  "age_at_diagnosis_years",   # 换算后（年）
  "gender_final",
  "vital_status_final",
  "days_to_death_final",
  "days_to_followup_final",
  "OS.time", "OS.status"
)
keep_cols <- keep_cols[keep_cols %in% colnames(cohort)]

write.table(
  cohort[, keep_cols, drop = FALSE],
  file  = file.path(proc_dir, "tcga_paad_cohort_master.tsv"),
  sep   = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  strict_pdac[, keep_cols[keep_cols %in% colnames(strict_pdac)], drop = FALSE],
  file  = file.path(proc_dir, "tcga_paad_cohort_master_strict_pdac.tsv"),
  sep   = "\t", quote = FALSE, row.names = FALSE
)

# -------------------------
# 8. 审计输出
# -------------------------
audit_df <- data.frame(
  metric = c(
    "n_tp_raw",
    "n_after_barcode_dedup",
    "n_final_cohort_all_tp",
    "n_final_cohort_strict_pdac",
    "n_os_event_all_tp",
    "n_os_event_strict_pdac",
    "os_time_min_days",
    "os_time_median_days",
    "os_time_max_days"
  ),
  value = c(
    n_tp_raw,
    length(unique(cohort$sample_barcode)),
    nrow(cohort),
    nrow(strict_pdac),
    sum(cohort$OS.status == 1L, na.rm = TRUE),
    sum(strict_pdac$OS.status == 1L, na.rm = TRUE),
    min(cohort$OS.time),
    median(cohort$OS.time),
    max(cohort$OS.time)
  )
)
write.table(
  audit_df,
  file  = file.path(proc_dir, "tcga_paad_cohort_build_audit.tsv"),
  sep   = "\t", quote = FALSE, row.names = FALSE
)

cat("\n========================================\n")
cat("Cohort building done.\n")
cat(sprintf("All-TP cohort      n = %d\n", nrow(cohort)))
cat(sprintf("Strict PDAC cohort n = %d\n", nrow(strict_pdac)))
cat(sprintf("OS events (all-TP)      = %d\n", sum(cohort$OS.status == 1L, na.rm = TRUE)))
cat(sprintf("OS events (strict PDAC) = %d\n", sum(strict_pdac$OS.status == 1L, na.rm = TRUE)))
cat("========================================\n")
