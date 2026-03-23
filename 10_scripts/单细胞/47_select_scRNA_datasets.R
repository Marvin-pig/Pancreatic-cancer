# =========================
# 47_select_scRNA_datasets.R
# 目的：
# 1. 初始化胰腺癌单细胞分析工作区
# 2. 建立 scRNA 数据集筛选 registry
# 3. 固定纳入/排除标准
# 4. 为后续主分析队列和验证队列选择做准备
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "06_scRNA"

registry_dir  <- file.path(root_dir, "00_registry")
raw_dir       <- file.path(root_dir, "01_raw_data")
proc_dir      <- file.path(root_dir, "02_processed_data")
result_dir    <- file.path(root_dir, "03_results")
figure_dir    <- file.path(root_dir, "04_figures")
table_dir     <- file.path(root_dir, "05_tables")
log_dir       <- file.path(root_dir, "06_logs")

dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Step 1. 建立 scRNA 数据集 registry 模板
# =========================================================
candidate_registry <- data.frame(
  dataset_id = c("PDAC_scRNA_01", "PDAC_scRNA_02", "PDAC_scRNA_03"),
  accession = c(NA, NA, NA),
  source = c(NA, NA, NA),
  title = c(NA, NA, NA),
  species = c("Human", "Human", "Human"),
  disease = c("PDAC", "PDAC", "PDAC"),
  sample_type = c(NA, NA, NA),
  tissue_origin = c(NA, NA, NA),
  treatment_status = c(NA, NA, NA),
  platform = c(NA, NA, NA),
  n_samples = c(NA, NA, NA),
  raw_count_available = c(NA, NA, NA),
  processed_object_available = c(NA, NA, NA),
  cell_annotation_available = c(NA, NA, NA),
  malignant_compartment = c(NA, NA, NA),
  immune_compartment = c(NA, NA, NA),
  stromal_compartment = c(NA, NA, NA),
  suitable_for_cellchat = c(NA, NA, NA),
  suitable_for_signature_mapping = c(NA, NA, NA),
  include_as_main = c(NA, NA, NA),
  include_as_validation = c(NA, NA, NA),
  exclusion_reason = c(NA, NA, NA),
  note = c(NA, NA, NA),
  stringsAsFactors = FALSE
)

fwrite(
  candidate_registry,
  file.path(registry_dir, "47_scRNA_dataset_registry.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 2. 固定纳入/排除标准
# =========================================================
inclusion_criteria <- data.frame(
  item = c(
    "disease_scope",
    "species",
    "preferred_sample_type",
    "raw_data_requirement",
    "cell_compartment_requirement",
    "treatment_annotation_requirement",
    "communication_analysis_requirement",
    "signature_mapping_requirement"
  ),
  value = c(
    "Pancreatic ductal adenocarcinoma / pancreatic cancer",
    "Human",
    "Primary tumor preferred",
    "Raw count matrix or standard processed object required",
    "Prefer malignant + immune + stromal compartments",
    "Untreated or clearly annotated treatment status",
    "Should support cell-cell communication analysis",
    "Must allow mapping of 12-gene signature and MCD-related genes"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  inclusion_criteria,
  file.path(registry_dir, "47_scRNA_inclusion_criteria.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 3. 生成初始选择决策表
# =========================================================
selection_decision <- data.frame(
  decision_item = c(
    "main_dataset_selected",
    "validation_dataset_selected",
    "minimum_candidate_dataset_n",
    "current_stage",
    "recommended_next_action"
  ),
  value = c(
    FALSE,
    FALSE,
    3,
    "dataset_selection_initialization",
    "Fill 47_scRNA_dataset_registry.tsv and decide 1 main + 1 validation dataset before Seurat analysis"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  selection_decision,
  file.path(registry_dir, "47_scRNA_selection_decision.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 4. 输出一个空的候选数据集备注文件
# =========================================================
writeLines(
  c(
    "47 scRNA dataset selection notes",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Prepare pancreatic cancer scRNA-seq dataset selection framework",
    "for downstream malignant-cell localisation, TME mapping,",
    "12-gene signature positioning, and cell-cell communication analysis.",
    "",
    "[Expected minimum deliverables]",
    "1. At least 3 candidate scRNA datasets",
    "2. 1 main dataset selected",
    "3. 1 validation dataset selected",
    "",
    "[Do not start Seurat workflow until dataset selection is fixed.]"
  ),
  file.path(log_dir, "47_scRNA_dataset_selection_notes.txt")
)

# =========================================================
# Step 5. 输出总状态
# =========================================================
status_summary <- data.frame(
  item = c(
    "registry_written",
    "criteria_written",
    "selection_decision_written",
    "log_written",
    "ready_for_dataset_filling"
  ),
  value = c(
    file.exists(file.path(registry_dir, "47_scRNA_dataset_registry.tsv")),
    file.exists(file.path(registry_dir, "47_scRNA_inclusion_criteria.tsv")),
    file.exists(file.path(registry_dir, "47_scRNA_selection_decision.tsv")),
    file.exists(file.path(log_dir, "47_scRNA_dataset_selection_notes.txt")),
    TRUE
  ),
  stringsAsFactors = FALSE
)

fwrite(
  status_summary,
  file.path(registry_dir, "47_scRNA_status_summary.tsv"),
  sep = "\t",
  na = "NA"
)

cat("47_select_scRNA_datasets.R finished successfully.\n")
cat("Generated files:\n")
cat("- 00_registry/47_scRNA_dataset_registry.tsv\n")
cat("- 00_registry/47_scRNA_inclusion_criteria.tsv\n")
cat("- 00_registry/47_scRNA_selection_decision.tsv\n")
cat("- 00_registry/47_scRNA_status_summary.tsv\n")
cat("- 06_logs/47_scRNA_dataset_selection_notes.txt\n")
cat("Next step:\n")
cat("- Fill candidate datasets into 47_scRNA_dataset_registry.tsv\n")
cat("- Select 1 main dataset and 1 validation dataset before downstream analysis\n")