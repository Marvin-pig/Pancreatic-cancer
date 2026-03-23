# =========================
# 66_initialize_drug_sensitivity_and_repositioning_module.R
# 目的：
# 1. 初始化胰腺癌药物敏感性 / 药物重定位模块
# 2. 固定输入来源与分析策略
# 3. 为 67 正式药敏预测分析做准备
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "07_drug"

registry_dir <- file.path(root_dir, "00_registry")
raw_dir      <- file.path(root_dir, "01_raw_data")
proc_dir     <- file.path(root_dir, "02_processed_data")
result_dir   <- file.path(root_dir, "03_results")
figure_dir   <- file.path(root_dir, "04_figures")
table_dir    <- file.path(root_dir, "05_tables")
log_dir      <- file.path(root_dir, "06_logs")

dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("10_scripts/药物", recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Step 1. 冻结药物模块总体方案
# =========================================================
drug_module_plan <- data.frame(
  item = c(
    "module_name",
    "current_stage",
    "primary_goal",
    "primary_analysis_line",
    "secondary_analysis_line",
    "main_dataset",
    "risk_group_source",
    "next_step"
  ),
  value = c(
    "drug_sensitivity_and_repositioning",
    "initialization",
    "complete translational extension of bulk risk model",
    "predicted_drug_sensitivity_between_high_low_risk_groups",
    "drug_repositioning_based_on_risk_associated_gene_programs",
    "TCGA-PAAD training cohort",
    "frozen 12-gene signature with training median cutoff",
    "67_run_predicted_drug_sensitivity_analysis"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  drug_module_plan,
  file.path(registry_dir, "66_drug_module_plan.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 2. 输入来源登记
# 注意：此处先登记，不强行假设你后续文件名完全一致
# =========================================================
drug_input_registry <- data.frame(
  input_type = c(
    "training_expression_matrix",
    "training_clinical_table",
    "frozen_signature_definition",
    "risk_group_assignment",
    "deg_or_program_for_repositioning"
  ),
  expected_source = c(
    "TCGA-PAAD normalized expression matrix",
    "TCGA-PAAD clinical survival table",
    "12-gene frozen signature",
    "training median cutoff based high/low risk grouping",
    "high vs low risk DEGs or malignant program genes"
  ),
  status = c("expected", "expected", "available", "available", "expected"),
  note = c(
    "used for drug sensitivity prediction",
    "used for sample alignment if needed",
    "already frozen in project",
    "already frozen in project",
    "can be imported later from bulk/single-cell downstream results"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  drug_input_registry,
  file.path(registry_dir, "66_drug_input_registry.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 3. 分析策略冻结
# =========================================================
drug_strategy <- data.frame(
  strategy_block = c(
    "A1",
    "A2",
    "A3",
    "B1",
    "B2",
    "B3"
  ),
  strategy_name = c(
    "predicted_drug_sensitivity",
    "group_comparison",
    "main_output",
    "drug_repositioning",
    "gene_program_source",
    "supplementary_output"
  ),
  detail = c(
    "Use transcriptome-based pharmacogenomic prediction on TCGA-PAAD samples",
    "Compare predicted response between high-risk and low-risk groups",
    "Generate ranked candidate drugs with effect direction and visualization",
    "Use risk-associated genes/programs for candidate drug reversal matching",
    "Use bulk DEGs and optionally malignant-pool program genes",
    "Generate repositioning candidate table and biological interpretation notes"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  drug_strategy,
  file.path(registry_dir, "66_drug_analysis_strategy.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 4. 阶段摘要
# =========================================================
stage_summary <- data.frame(
  item = c(
    "stage66_completed",
    "module_plan_written",
    "input_registry_written",
    "strategy_written",
    "next_stage"
  ),
  value = c(
    TRUE,
    file.exists(file.path(registry_dir, "66_drug_module_plan.tsv")),
    file.exists(file.path(registry_dir, "66_drug_input_registry.tsv")),
    file.exists(file.path(registry_dir, "66_drug_analysis_strategy.tsv")),
    "67_run_predicted_drug_sensitivity_analysis"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  stage_summary,
  file.path(registry_dir, "66_stage_summary.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 5. 日志
# =========================================================
writeLines(
  c(
    "66 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Initialize the drug sensitivity and drug repositioning module.",
    "",
    "[Current decision]",
    "Do not start manuscript writing yet.",
    "Finish Stage A by prioritizing drug sensitivity/repositioning results first.",
    "",
    "[Next step]",
    "Run predicted drug sensitivity analysis using frozen TCGA-PAAD risk groups."
  ),
  file.path(log_dir, "66_drug_module_initialization_notes.txt")
)

cat("66_initialize_drug_sensitivity_and_repositioning_module.R finished successfully.\n")