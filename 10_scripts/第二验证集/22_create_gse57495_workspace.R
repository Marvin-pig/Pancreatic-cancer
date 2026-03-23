# =========================
# 22_create_gse57495_workspace.R
# 目的：
# 1. 为 GSE57495 建立标准目录
# 2. 更新外部验证队列登记表
# 3. 建立处理日志
# =========================

suppressPackageStartupMessages({
  library(data.table)
})
setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "04_bulk_analysis/04_external_validation"

dirs <- c(
  file.path(root_dir, "01_raw_data/GEO/GSE57495"),
  file.path(root_dir, "02_processed_data/GSE57495"),
  file.path(root_dir, "03_qc/GSE57495"),
  file.path(root_dir, "04_signature_validation/GSE57495"),
  file.path(root_dir, "07_logs/GSE57495")
)

for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

registry_file <- file.path(root_dir, "00_registry/16_external_cohort_registry.tsv")

if (file.exists(registry_file)) {
  reg <- fread(registry_file)
  if ("GSE57495" %in% reg$dataset_id) {
    reg$download_status[reg$dataset_id == "GSE57495"] <- "IN_PROGRESS"
    fwrite(reg, registry_file, sep = "\t")
  }
}

log_lines <- c(
  "GSE57495 external validation log",
  paste0("Created at: ", Sys.time()),
  "",
  "[Planned usage]",
  "Second external validation cohort for 12-gene MCD-axis-derived signature",
  "",
  "[Rules]",
  "1. Do not refit the model in GSE57495",
  "2. Use frozen coefficients from TCGA training model",
  "3. First build clean expression and survival tables",
  "4. Record all manual filtering decisions"
)

writeLines(
  log_lines,
  file.path(root_dir, "07_logs/GSE57495/22_GSE57495_processing_log.txt")
)

cat("GSE57495 workspace created successfully.\n")