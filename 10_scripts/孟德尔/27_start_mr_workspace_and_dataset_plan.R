# =========================
# 27_start_mr_workspace_and_dataset_plan.R
# 目的：
# 1. 建立 MR 分析工作区
# 2. 冻结候选暴露对象
# 3. 生成 MR 数据需求清单
# 4. 生成后续 MR 执行路线图
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})
setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "05_MR"
dirs <- c(
  root_dir,
  file.path(root_dir, "00_registry"),
  file.path(root_dir, "01_raw_data"),
  file.path(root_dir, "01_raw_data/exposure"),
  file.path(root_dir, "01_raw_data/outcome"),
  file.path(root_dir, "02_processed_data"),
  file.path(root_dir, "03_results"),
  file.path(root_dir, "04_figures"),
  file.path(root_dir, "05_tables"),
  file.path(root_dir, "06_logs")
)

for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---------- 读取 bulk 已完成结果 ----------
coef_file <- "04_bulk_analysis/03_survival_model/15_final_signature_coefficients.tsv"
candidate_file <- "04_bulk_analysis/03_survival_model/12_prognostic_gene_shortlist.tsv"

coef_df <- if (file.exists(coef_file)) fread(coef_file) else data.frame()
cand_df <- if (file.exists(candidate_file)) fread(candidate_file) else data.frame()

# ---------- 定义 MR 候选对象 ----------
# 第一层：final signature genes
signature_genes <- if ("gene" %in% colnames(coef_df)) unique(coef_df$gene) else character(0)

# 第二层：可考虑纳入的主线机制基因（手动补充，可后续调整）
core_mechanism_genes <- c(
  "GPX4", "SLC31A1", "ABI1", "FLNA",
  "CXCL9", "CXCL10", "CXCL11", "FGL2",
  "IL7R", "CMKLR1", "NOTCH2", "INHBA", "TGFBR3",
  "MET", "TGFBI", "IDO1"
)

mr_targets <- unique(c(signature_genes, core_mechanism_genes))

mr_target_df <- data.frame(
  gene = mr_targets,
  target_class = ifelse(mr_targets %in% signature_genes, "signature_gene", "mechanism_gene"),
  priority = ifelse(mr_targets %in% signature_genes, "HIGH", "MEDIUM"),
  stringsAsFactors = FALSE
)

fwrite(
  mr_target_df,
  file.path(root_dir, "00_registry/27_mr_target_gene_registry.tsv"),
  sep = "\t"
)

# ---------- MR 设计说明 ----------
design_lines <- c(
  "MR project design notes",
  paste0("Created at: ", Sys.time()),
  "",
  "[Primary aim]",
  "Evaluate potential causal effects of genetically predicted gene expression / protein abundance on pancreatic cancer risk or prognosis-related traits.",
  "",
  "[Current candidate targets]",
  paste(mr_targets, collapse = ", "),
  "",
  "[Suggested primary exposure types]",
  "1. cis-eQTL for signature genes / mechanism genes",
  "2. pQTL if available for key proteins",
  "",
  "[Suggested primary outcomes]",
  "1. Pancreatic cancer risk GWAS",
  "2. If feasible, PDAC-specific outcome datasets",
  "",
  "[Core MR quality rules]",
  "1. Instrument SNPs should satisfy genome-wide or predefined significance threshold",
  "2. F-statistic should be checked",
  "3. LD clumping is required",
  "4. Harmonization and palindromic SNP handling must be explicit",
  "5. Sensitivity analyses should include heterogeneity / pleiotropy tests when instrument count permits",
  "6. Colocalization should be considered for high-priority positive findings"
)

writeLines(
  design_lines,
  file.path(root_dir, "06_logs/27_mr_design_notes.txt")
)

# ---------- 数据需求清单 ----------
data_need <- data.frame(
  item = c(
    "exposure_summary_statistics",
    "outcome_summary_statistics",
    "gene_to_exposure_mapping",
    "instrument_selection_rules",
    "harmonization_rules",
    "sensitivity_analysis_plan",
    "colocalization_candidate_list"
  ),
  status = c("TODO", "TODO", "TODO", "TODO", "TODO", "TODO", "TODO"),
  note = c(
    "cis-eQTL or pQTL summary data for selected genes",
    "pancreatic cancer GWAS summary statistics",
    "match genes to available eQTL/pQTL datasets",
    "p threshold, LD clumping, F-statistics",
    "effect allele alignment and palindromic handling",
    "MR-Egger / weighted median / heterogeneity where applicable",
    "apply to top positive MR signals"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  data_need,
  file.path(root_dir, "00_registry/27_mr_data_requirements.tsv"),
  sep = "\t"
)

# ---------- 执行路线图 ----------
roadmap <- data.frame(
  step_id = c("27", "28", "29", "30", "31", "32"),
  step_name = c(
    "build MR workspace and target registry",
    "select exposure and outcome datasets",
    "download / prepare MR summary statistics",
    "run primary MR",
    "run sensitivity analyses",
    "prioritize findings for colocalization / interpretation"
  ),
  status = c("DONE", "TODO", "TODO", "TODO", "TODO", "TODO"),
  stringsAsFactors = FALSE
)

fwrite(
  roadmap,
  file.path(root_dir, "00_registry/27_mr_roadmap.tsv"),
  sep = "\t"
)

cat("27_start_mr_workspace_and_dataset_plan.R finished successfully.\n")
cat("MR workspace created at:", root_dir, "\n")