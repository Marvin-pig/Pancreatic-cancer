# =========================
# 48_fill_scRNA_dataset_registry.R
# 目的：
# 1. 将候选 PDAC 单细胞数据集正式写入 registry
# 2. 固定主分析队列、验证队列和补充扩展队列
# 3. 为后续下载、Seurat 预处理和细胞通讯分析做准备
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "06_scRNA"

registry_dir <- file.path(root_dir, "00_registry")
log_dir      <- file.path(root_dir, "06_logs")

dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Step 1. 写入候选数据集 registry
# =========================================================
candidate_registry <- data.frame(
  dataset_id = c("PDAC_scRNA_01", "PDAC_scRNA_02", "PDAC_scRNA_03"),
  accession = c("GSE155698", "CRA001160", "GSE205013"),
  source = c("GEO", "CNCB-NGDC/GSA", "GEO"),
  title = c(
    "Multimodal Mapping of the Tumor and Peripheral Blood Immune Landscape in Human Pancreatic Cancer",
    "GSA-PDAC",
    "Effects of Chemotherapy on Human Pancreatic Adenocarcinoma and its Tumor Microenvironment"
  ),
  species = c("Human", "Human", "Human"),
  disease = c("PDAC", "PDAC", "PDAC"),
  sample_type = c("Primary tumor + adjacent normal", "Primary tumor + normal pancreas", "Primary tumor + liver metastasis"),
  tissue_origin = c("Pancreas", "Pancreas", "Pancreas and liver metastasis"),
  treatment_status = c("Prefer treatment-naive / tumor tissue cohort", "Not clearly heterogeneous in public reuse literature", "Mixed pre-treatment and post-chemotherapy"),
  platform = c("scRNA-seq", "scRNA-seq", "scRNA-seq"),
  n_samples = c("16 PDAC + 3 adjacent normal", "24 primary PDAC + 11 normal pancreas", "Primary/metastatic mixed; check series matrix / supplement"),
  raw_count_available = c("YES", "YES", "YES"),
  processed_object_available = c("UNKNOWN", "UNKNOWN", "YES/filtered matrices reported"),
  cell_annotation_available = c("NO", "NO", "NO"),
  malignant_compartment = c("YES", "YES", "YES"),
  immune_compartment = c("YES", "YES", "YES"),
  stromal_compartment = c("YES", "YES", "YES"),
  suitable_for_cellchat = c("YES", "YES", "PARTLY"),
  suitable_for_signature_mapping = c("YES", "YES", "YES"),
  include_as_main = c("YES", "NO", "NO"),
  include_as_validation = c("NO", "YES", "NO"),
  exclusion_reason = c(
    NA,
    NA,
    "Not preferred as main dataset because treatment/metastasis heterogeneity may confound baseline TME interpretation"
  ),
  note = c(
    "Recommended main dataset for malignant/immune/stromal localisation and module scoring",
    "Recommended validation dataset for confirming cell-type localisation and communication findings",
    "Recommended supplementary extension dataset for chemotherapy/metastasis heterogeneity analysis"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  candidate_registry,
  file.path(registry_dir, "47_scRNA_dataset_registry.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 2. 更新选择决策
# =========================================================
selection_decision <- data.frame(
  decision_item = c(
    "main_dataset_selected",
    "main_dataset_accession",
    "validation_dataset_selected",
    "validation_dataset_accession",
    "supplementary_dataset_accession",
    "minimum_candidate_dataset_n",
    "current_stage",
    "recommended_next_action"
  ),
  value = c(
    "TRUE",
    "GSE155698",
    "TRUE",
    "CRA001160",
    "GSE205013",
    "3",
    "dataset_selection_fixed",
    "Start 49: download preparation and dataset-level metadata audit before Seurat preprocessing"
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
# Step 3. 更新状态摘要
# =========================================================
status_summary <- data.frame(
  item = c(
    "registry_written",
    "candidate_dataset_n",
    "main_dataset_selected",
    "validation_dataset_selected",
    "ready_for_download_preparation"
  ),
  value = c(
    TRUE,
    nrow(candidate_registry),
    TRUE,
    TRUE,
    TRUE
  ),
  stringsAsFactors = FALSE
)

fwrite(
  status_summary,
  file.path(registry_dir, "48_scRNA_status_summary.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 4. 记录说明日志
# =========================================================
writeLines(
  c(
    "48 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Fill candidate PDAC scRNA-seq datasets and fix main/validation cohort selection.",
    "",
    "[Selected datasets]",
    "Main dataset: GSE155698",
    "Validation dataset: CRA001160",
    "Supplementary extension dataset: GSE205013",
    "",
    "[Next step]",
    "Start 49: download preparation and dataset-level metadata audit before Seurat preprocessing."
  ),
  file.path(log_dir, "48_fill_scRNA_dataset_registry_notes.txt")
)

cat("48_fill_scRNA_dataset_registry.R finished successfully.\n")
cat("Generated files:\n")
cat("- 00_registry/47_scRNA_dataset_registry.tsv\n")
cat("- 00_registry/47_scRNA_selection_decision.tsv\n")
cat("- 00_registry/48_scRNA_status_summary.tsv\n")
cat("- 06_logs/48_fill_scRNA_dataset_registry_notes.txt\n")
cat("Next step:\n")
cat("- Start 49: download preparation and dataset-level metadata audit\n")
