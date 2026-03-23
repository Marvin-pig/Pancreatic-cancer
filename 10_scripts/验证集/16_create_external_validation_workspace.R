# =========================
# 16_create_external_validation_workspace.R
# 目的：
# 1. 为外部验证建立标准目录
# 2. 生成队列登记表、处理日志、任务清单
# 3. 为 GEO / ICGC / 后续结果输出预留统一结构
# =========================

suppressPackageStartupMessages({
  library(data.table)
})

# ---------- 根目录 ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "04_bulk_analysis/04_external_validation"

dirs <- c(
  root_dir,
  file.path(root_dir, "00_registry"),
  file.path(root_dir, "01_raw_data"),
  file.path(root_dir, "01_raw_data/GEO"),
  file.path(root_dir, "01_raw_data/ICGC"),
  file.path(root_dir, "02_processed_data"),
  file.path(root_dir, "02_processed_data/GSE21501"),
  file.path(root_dir, "03_qc"),
  file.path(root_dir, "04_signature_validation"),
  file.path(root_dir, "05_figures"),
  file.path(root_dir, "06_tables"),
  file.path(root_dir, "07_logs")
)

for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---------- 队列登记表 ----------
cohort_registry <- data.frame(
  dataset_id = c("GSE21501", "GSE28735", "ICGC_PACA"),
  source = c("GEO", "GEO", "ICGC"),
  priority = c("HIGH", "LOW", "MEDIUM"),
  intended_use = c(
    "primary external prognostic validation",
    "biology/supporting expression comparison",
    "secondary external validation if access available"
  ),
  download_status = c("TODO", "TODO", "TODO"),
  expression_ready = c("NO", "NO", "NO"),
  clinical_ready = c("NO", "NO", "NO"),
  signature_applicable = c("UNKNOWN", "UNKNOWN", "UNKNOWN"),
  notes = c(
    "start here first",
    "not first-choice survival validation cohort",
    "requires access planning"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  cohort_registry,
  file.path(root_dir, "00_registry/16_external_cohort_registry.tsv"),
  sep = "\t"
)

# ---------- 任务清单 ----------
task_list <- data.frame(
  step_id = c("16", "17", "18", "19", "20", "21"),
  task = c(
    "create workspace",
    "download GSE21501",
    "extract expression matrix",
    "build clinical survival table",
    "map signature genes",
    "run external signature validation"
  ),
  status = c("DONE", "TODO", "TODO", "TODO", "TODO", "TODO"),
  stringsAsFactors = FALSE
)

fwrite(
  task_list,
  file.path(root_dir, "00_registry/16_external_task_list.tsv"),
  sep = "\t"
)

# ---------- 日志模板 ----------
log_template <- c(
  "External validation log",
  paste0("Created at: ", Sys.time()),
  "",
  "[Rules]",
  "1. Do not change the frozen training signature coefficients.",
  "2. Do not re-fit model in external cohort.",
  "3. Use training median cutoff unless explicitly justified otherwise.",
  "4. Keep expression and clinical sample IDs strictly matched.",
  "5. Record every manual filtering decision."
)

writeLines(
  log_template,
  file.path(root_dir, "07_logs/16_external_validation_log.txt")
)

cat("External validation workspace created successfully.\n")
cat("Root directory:", root_dir, "\n")