# =========================
# 32_outcome_first_manual_attachment_plan.R
# 目的：
# 1. 固定 PanScan/PanC4 主 outcome 的本地文件位置
# 2. 生成 outcome 文件命名规范
# 3. 生成 outcome 可用性检查表
# 4. 为后续 MR 主分析准备最小 outcome 接口
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

root_dir <- "05_MR"
outcome_dir <- file.path(root_dir, "01_raw_data/outcome/PanScan_PanC4")

# ---- 确保所有输出目录存在 ----
dir.create(outcome_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "00_registry"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "02_processed_data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "06_logs"), recursive = TRUE, showWarnings = FALSE)

# ---------- 1. 固定 outcome 文件规范 ----------
outcome_registry <- data.frame(
  item = c(
    "preferred_file_name",
    "preferred_format",
    "build",
    "trait",
    "ancestry",
    "role"
  ),
  value = c(
    "pancreatic_cancer_sumstats.tsv",
    "tab-delimited text",
    "GRCh37_or_liftable",
    "pancreatic cancer risk",
    "European",
    "primary_outcome"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  outcome_registry,
  file.path(root_dir, "00_registry/32_primary_outcome_registry.tsv"),
  sep = "\t"
)

# ---------- 2. outcome 占位文件说明 ----------
placeholder_lines <- c(
  "Place the primary pancreatic cancer GWAS summary statistics here:",
  file.path(outcome_dir, "pancreatic_cancer_sumstats.tsv"),
  "",
  "Required columns:",
  "SNP",
  "effect_allele",
  "other_allele",
  "beta",
  "se",
  "pval",
  "eaf",
  "outcome",
  "",
  "Optional but recommended:",
  "chr",
  "pos",
  "ncase",
  "ncontrol",
  "samplesize"
)

writeLines(
  placeholder_lines,
  file.path(outcome_dir, "32_PLACE_PRIMARY_OUTCOME_FILE_HERE.txt")
)

# ---------- 3. outcome 模板复制一份到目标目录 ----------
template_src <- file.path(root_dir, "02_processed_data/30_outcome_summary_template.tsv")
template_dst <- file.path(outcome_dir, "32_outcome_summary_template_copy.tsv")

if (file.exists(template_src)) {
  file.copy(template_src, template_dst, overwrite = TRUE)
  cat("Template copied:", template_src, "->", template_dst, "\n")
} else {
  cat("Note: Template source not found, skipping copy:", template_src, "\n")
}

# ---------- 4. outcome 预检查表 ----------
expected_outcome_file <- file.path(outcome_dir, "pancreatic_cancer_sumstats.tsv")

check_tab <- data.frame(
  check_item = c(
    "outcome_file_exists",
    "outcome_template_copy_exists"
  ),
  status = c(
    file.exists(expected_outcome_file),
    file.exists(template_dst)
  ),
  stringsAsFactors = FALSE
)

fwrite(
  check_tab,
  file.path(root_dir, "00_registry/32_primary_outcome_precheck.tsv"),
  sep = "\t"
)

# ---------- 5. outcome 读取检查脚本模板 ----------
check_script_lines <- c(
  "library(data.table)",
  sprintf("f <- '%s'", gsub("\\\\", "/", expected_outcome_file)),
  "if (!file.exists(f)) stop('Primary outcome file not found: ', f)",
  "x <- fread(f)",
  "required <- c('SNP','effect_allele','other_allele','beta','se','pval','eaf','outcome')",
  "missing <- setdiff(required, colnames(x))",
  "if (length(missing) > 0) stop('Missing columns: ', paste(missing, collapse=', '))",
  "cat('Outcome file looks structurally valid.\\n')",
  "cat('Rows:', nrow(x), '\\n')",
  "cat('Columns:', ncol(x), '\\n')",
  "print(head(x))"
)

writeLines(
  check_script_lines,
  file.path(root_dir, "02_processed_data/32_check_primary_outcome_file.R")
)

# ---------- 6. 状态总结 ----------
status_lines <- c(
  "Primary outcome attachment plan",
  paste0("Created at: ", Sys.time()),
  "",
  "[Current status]",
  paste0("Outcome file exists: ", file.exists(expected_outcome_file)),
  "",
  "[Next action for user]",
  "Place PanScan/PanC4 pancreatic cancer GWAS summary statistics file at:",
  expected_outcome_file,
  "",
  "[Then run]",
  "05_MR/02_processed_data/32_check_primary_outcome_file.R"
)

writeLines(
  status_lines,
  file.path(root_dir, "06_logs/32_primary_outcome_attachment_notes.txt")
)

cat("32_outcome_first_manual_attachment_plan.R finished successfully.\n")
cat("Expected primary outcome path:\n")
cat(expected_outcome_file, "\n")