# =========================
# 29_prepare_openGWAS_and_eqtlgen_access.R
# 目的：
# 1. 建立 MR 数据访问与本地文件接口
# 2. 生成 batch1 基因查询模板
# 3. 生成 OpenGWAS / 本地文件配置表
# 4. 生成 MR 运行前检查表
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "05_MR"
dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 固定目录 ----------
dirs <- c(
  file.path(root_dir, "00_registry"),
  file.path(root_dir, "01_raw_data/exposure/eQTLGen"),
  file.path(root_dir, "01_raw_data/exposure/GTEx_pancreas"),
  file.path(root_dir, "01_raw_data/outcome/PanScan_PanC4"),
  file.path(root_dir, "01_raw_data/outcome/OpenGWAS_BBJ"),
  file.path(root_dir, "02_processed_data/batch1_queries"),
  file.path(root_dir, "06_logs")
)

for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---------- 读取 gene execution plan ----------
gene_plan_file <- file.path(root_dir, "00_registry/28_mr_gene_execution_plan.tsv")

if (!file.exists(gene_plan_file)) {
  stop(
    "[ERROR] Gene execution plan not found: ", gene_plan_file, "\n",
    "  Please run script 28 first to generate this file."
  )
}

gene_plan <- fread(gene_plan_file)

# 基本校验：必需列是否存在
required_cols <- c("gene", "batch")
missing_cols <- setdiff(required_cols, colnames(gene_plan))
if (length(missing_cols) > 0) {
  stop(
    "[ERROR] Missing required columns in gene execution plan: ",
    paste(missing_cols, collapse = ", ")
  )
}

batch1 <- gene_plan %>%
  filter(batch == "batch1_priority") %>%
  arrange(gene)

if (nrow(batch1) == 0) {
  stop(
    "[ERROR] No genes found with batch == 'batch1_priority' in:\n  ",
    gene_plan_file, "\n",
    "  Available batch values: ",
    paste(unique(gene_plan$batch), collapse = ", ")
  )
}

cat(sprintf("[INFO] batch1_priority genes: %d\n", nrow(batch1)))

fwrite(
  batch1,
  file.path(root_dir, "02_processed_data/batch1_queries/29_batch1_priority_genes.tsv"),
  sep = "\t"
)

# ---------- OpenGWAS / 本地文件配置 ----------
config_tab <- data.frame(
  config_key = c(
    "opengwas_jwt",
    "bbj_outcome_id",
    "panscan_panc4_summary_stats_file",
    "eqtlgen_source_mode",
    "eqtlgen_local_file",
    "gtex_pancreas_source_mode",
    "gtex_pancreas_local_file"
  ),
  config_value = c(
    "",
    "bbj-a-140",
    "",
    "local_or_manual_download",
    "",
    "local_or_manual_download",
    ""
  ),
  description = c(
    "OpenGWAS JWT token if API access is used",
    "Sensitivity outcome dataset ID in OpenGWAS",
    "Local file path to PanScan/PanC4 summary statistics",
    "Set to 'local' when you have downloaded eQTLGen files",
    "Local eQTLGen cis-eQTL summary stats file or folder",
    "Set to 'local' when GTEx pancreas cis-eQTL files are ready",
    "Local GTEx pancreas cis-eQTL file or folder"
  ),
  stringsAsFactors = FALSE
)

config_path <- file.path(root_dir, "00_registry/29_mr_access_config.tsv")

# 如果配置文件已存在，提示而非覆盖（避免丢失手动填入的值）
if (file.exists(config_path)) {
  cat("[WARN] Config file already exists, overwriting: ", config_path, "\n")
  cat("  If you have manually filled values, back up the file first.\n")
}

fwrite(config_tab, config_path, sep = "\t")

# ---------- batch1 基因查询模板 ----------
query_template <- batch1 %>%
  transmute(
    gene            = gene,
    exposure_primary    = "eQTLGen_whole_blood_cis_eQTL",
    exposure_support    = "GTEx_v8_pancreas_cis_eQTL",
    outcome_primary     = "PanScan_PanC4_European_pancreatic_cancer_GWAS",
    outcome_sensitivity = "BBJ_pancreatic_cancer_bbj-a-140",
    cis_window          = "1Mb",
    p_threshold         = "5e-8",
    ld_r2               = "0.001",
    clump_kb            = "10000",
    status              = "TODO"
  )

fwrite(
  query_template,
  file.path(root_dir, "02_processed_data/batch1_queries/29_batch1_mr_query_template.tsv"),
  sep = "\t"
)

# ---------- MR 运行前检查表 ----------
check_items <- data.frame(
  item = c(
    "batch1_gene_plan_ready",
    "access_config_ready",
    "eqtlgen_data_ready",
    "gtex_pancreas_data_ready",
    "panscan_panc4_outcome_ready",
    "bbj_outcome_ready",
    "jwt_ready_if_needed"
  ),
  status = c(
    "DONE",
    "DONE",
    "TODO",
    "TODO",
    "TODO",
    "TODO",
    "TODO"
  ),
  note = c(
    "29_batch1_priority_genes.tsv generated",
    "29_mr_access_config.tsv generated",
    "Need local eQTLGen cis-eQTL file",
    "Need local GTEx v8 pancreas file",
    "Need local PanScan/PanC4 summary statistics",
    "Can use OpenGWAS bbj-a-140 if access is available",
    "Only needed for OpenGWAS API-based access"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  check_items,
  file.path(root_dir, "00_registry/29_mr_preflight_checklist.tsv"),
  sep = "\t"
)

# ---------- 操作说明 ----------
notes <- c(
  "=== MR access preparation notes ===",
  paste0("Created at: ", Sys.time()),
  paste0("Script: 29_prepare_openGWAS_and_eqtlgen_access.R"),
  paste0("batch1 gene count: ", nrow(batch1)),
  "",
  "[Current status]",
  "MR execution framework is ready.",
  "Raw exposure/outcome summary statistics are NOT yet attached.",
  "",
  "[Immediate next actions]",
  "1. Open and fill 29_mr_access_config.tsv with actual file paths / tokens",
  "2. Place local eQTLGen / GTEx / PanScan-PanC4 files in 01_raw_data/",
  "3. If using OpenGWAS API, prepare JWT token and confirm bbj-a-140 access",
  "4. Re-run preflight checklist to verify all data sources before MR",
  "",
  "[Do NOT run main MR yet]",
  "Run MR only after at least one exposure source AND one outcome source are confirmed."
)

writeLines(
  notes,
  file.path(root_dir, "06_logs/29_mr_access_notes.txt")
)

# ---------- 完成汇报 ----------
cat("\n========================================\n")
cat("29_prepare_openGWAS_and_eqtlgen_access.R finished successfully.\n")
cat("========================================\n")
cat("Files generated:\n")
cat("  - 00_registry/29_mr_access_config.tsv\n")
cat("  - 00_registry/29_mr_preflight_checklist.tsv\n")
cat("  - 02_processed_data/batch1_queries/29_batch1_priority_genes.tsv\n")
cat("  - 02_processed_data/batch1_queries/29_batch1_mr_query_template.tsv\n")
cat("  - 06_logs/29_mr_access_notes.txt\n")