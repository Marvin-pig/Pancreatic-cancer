# =========================
# 30_define_mr_input_file_format_and_manual_data_acquisition_plan.R
# 目的：
# 1. 固定 MR 输入文件格式
# 2. 生成 exposure / outcome 模板文件
# 3. 生成手动数据获取清单
# 4. 为后续 31_run_primary_mr_batch1.R 做准备
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ---------- 设置工作目录（必须最先执行） ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")

root_dir <- "05_MR"

# ---------- 创建所有必要的子目录 ----------
dirs_needed <- c(
  file.path(root_dir, "00_registry"),
  file.path(root_dir, "01_raw_data/exposure/eQTLGen"),
  file.path(root_dir, "01_raw_data/exposure/GTEx_pancreas"),
  file.path(root_dir, "01_raw_data/outcome/PanScan_PanC4"),
  file.path(root_dir, "02_processed_data/batch1_queries"),
  file.path(root_dir, "06_logs")
)
for (d in dirs_needed) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---------- 读取已有配置 ----------
config_file <- file.path(root_dir, "00_registry/29_mr_access_config.tsv")
batch1_file <- file.path(root_dir, "02_processed_data/batch1_queries/29_batch1_priority_genes.tsv")

if (!file.exists(config_file)) {
  stop("配置文件不存在: ", config_file,
       "\n请先运行 29 号脚本或手动创建该文件。")
}
if (!file.exists(batch1_file)) {
  stop("Batch1 基因列表不存在: ", batch1_file,
       "\n请先运行 29 号脚本或手动创建该文件。")
}

config_tab <- fread(config_file)
batch1     <- fread(batch1_file)

cat("已读取配置文件：", nrow(config_tab), "行\n")
cat("已读取 batch1 基因列表：", nrow(batch1), "行\n")

# ---------- 固定输入格式：exposure 标准模板 ----------
exposure_template <- data.frame(
  SNP           = c("rs12345", "rs67890"),
  chr           = c(1, 2),
  pos           = c(1234567, 7654321),
  effect_allele = c("A", "G"),
  other_allele  = c("G", "A"),
  beta          = c(0.10, -0.08),
  se            = c(0.02, 0.03),
  pval          = c(1e-8, 5e-7),
  eaf           = c(0.42, 0.31),
  gene          = c("MET", "MET"),
  exposure      = c("eQTLGen_MET", "eQTLGen_MET"),
  stringsAsFactors = FALSE
)

fwrite(
  exposure_template,
  file.path(root_dir, "02_processed_data/30_exposure_summary_template.tsv"),
  sep = "\t"
)

# ---------- 固定输入格式：outcome 标准模板 ----------
outcome_template <- data.frame(
  SNP           = c("rs12345", "rs67890"),
  effect_allele = c("A", "G"),
  other_allele  = c("G", "A"),
  beta          = c(0.05, -0.03),
  se            = c(0.01, 0.015),
  pval          = c(1e-4, 2e-3),
  eaf           = c(0.41, 0.30),
  outcome       = c("Pancreatic_cancer_PanScanPanC4", "Pancreatic_cancer_PanScanPanC4"),
  stringsAsFactors = FALSE
)

fwrite(
  outcome_template,
  file.path(root_dir, "02_processed_data/30_outcome_summary_template.tsv"),
  sep = "\t"
)

# ---------- 生成 batch1 基因级数据模板 ----------
# 确认 batch1 含有必要列
required_cols <- c("gene", "target_class", "priority", "batch", "run_priority")
missing_cols  <- setdiff(required_cols, colnames(batch1))
if (length(missing_cols) > 0) {
  stop("batch1 文件缺少以下列: ", paste(missing_cols, collapse = ", "))
}

batch1_templates <- batch1 %>%
  select(all_of(required_cols)) %>%
  mutate(
    exposure_file_expected = paste0("01_raw_data/exposure/eQTLGen/", gene, "_cis_eqtl.tsv"),
    gtex_file_expected     = paste0("01_raw_data/exposure/GTEx_pancreas/", gene, "_cis_eqtl.tsv"),
    outcome_file_expected  = "01_raw_data/outcome/PanScan_PanC4/pancreatic_cancer_sumstats.tsv",
    status = "TODO"
  )

fwrite(
  batch1_templates,
  file.path(root_dir, "02_processed_data/30_batch1_expected_input_files.tsv"),
  sep = "\t"
)

# ---------- 手动数据获取清单 ----------
manual_plan <- data.frame(
  item_id = c("30A", "30B", "30C", "30D", "30E"),
  task = c(
    "Download / prepare eQTLGen cis-eQTL data for batch1 genes",
    "Download / prepare GTEx pancreas cis-eQTL data for batch1 genes",
    "Obtain PanScan/PanC4 pancreatic cancer GWAS summary statistics",
    "Prepare OpenGWAS JWT token if BBJ sensitivity outcome will be used",
    "Place files into expected local directories and update config"
  ),
  destination = c(
    "05_MR/01_raw_data/exposure/eQTLGen/",
    "05_MR/01_raw_data/exposure/GTEx_pancreas/",
    "05_MR/01_raw_data/outcome/PanScan_PanC4/",
    "05_MR/00_registry/29_mr_access_config.tsv",
    "05_MR/01_raw_data/..."
  ),
  status = c("TODO", "TODO", "TODO", "TODO", "TODO"),
  stringsAsFactors = FALSE
)

fwrite(
  manual_plan,
  file.path(root_dir, "00_registry/30_manual_data_acquisition_plan.tsv"),
  sep = "\t"
)

# ---------- 自动检查本地文件是否已放置 ----------
expected_files <- c(
  batch1_templates$exposure_file_expected,
  batch1_templates$gtex_file_expected,
  unique(batch1_templates$outcome_file_expected)
)

precheck <- data.frame(
  expected_file = expected_files,
  full_path     = file.path(root_dir, expected_files),
  exists        = file.exists(file.path(root_dir, expected_files)),
  stringsAsFactors = FALSE
)

fwrite(
  precheck,
  file.path(root_dir, "00_registry/30_local_file_presence_check.tsv"),
  sep = "\t"
)

n_found   <- sum(precheck$exists)
n_total   <- nrow(precheck)
cat(sprintf("文件就绪检查：%d / %d 已存在\n", n_found, n_total))

# ---------- 操作说明 ----------
notes <- c(
  "MR input format and manual acquisition notes",
  paste0("Created at: ", Sys.time()),
  "",
  "[Exposure summary statistics required columns]",
  "SNP, chr, pos, effect_allele, other_allele, beta, se, pval, eaf, gene, exposure",
  "",
  "[Outcome summary statistics required columns]",
  "SNP, effect_allele, other_allele, beta, se, pval, eaf, outcome",
  "",
  "[Current status]",
  sprintf("Templates are ready. %d / %d raw data files present.", n_found, n_total),
  "",
  "[Next rule]",
  "Do not run MR until at least one exposure file and one outcome file are present.",
  "",
  "[Files generated by this script]",
  "- 02_processed_data/30_exposure_summary_template.tsv",
  "- 02_processed_data/30_outcome_summary_template.tsv",
  "- 02_processed_data/30_batch1_expected_input_files.tsv",
  "- 00_registry/30_manual_data_acquisition_plan.tsv",
  "- 00_registry/30_local_file_presence_check.tsv",
  "- 06_logs/30_mr_input_format_notes.txt"
)

writeLines(
  notes,
  file.path(root_dir, "06_logs/30_mr_input_format_notes.txt")
)

# ---------- 完成 ----------
cat("\n30_define_mr_input_file_format_and_manual_data_acquisition_plan.R finished successfully.\n")
cat("Generated:\n")
cat("- 30_exposure_summary_template.tsv\n")
cat("- 30_outcome_summary_template.tsv\n")
cat("- 30_batch1_expected_input_files.tsv\n")
cat("- 30_manual_data_acquisition_plan.tsv\n")
cat("- 30_local_file_presence_check.tsv\n")
cat("- 30_mr_input_format_notes.txt\n")