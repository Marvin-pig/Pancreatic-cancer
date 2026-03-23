# =========================
# 68_drug_repositioning_and_biological_interpretation.R
# 目的：
# 1. 基于 67 的药敏结果，固定 low-risk / high-risk favored drugs
# 2. 结合 bulk 12-gene signature 与 single-cell malignant program
# 3. 生成药物结果的生物学解释表与写作桥接表
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "07_drug"

registry_dir <- file.path(root_dir, "00_registry")
proc_dir     <- file.path(root_dir, "02_processed_data")
result_dir   <- file.path(root_dir, "03_results")
table_dir    <- file.path(root_dir, "05_tables")
log_dir      <- file.path(root_dir, "06_logs")

dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Step 1. 读取 67 结果
# 优先用你当前项目目录里的文件；如果不存在，再尝试你上传到对话的文件名
# =========================================================
cmp_file_candidates <- c(
  file.path(result_dir, "67_candidate_drug_group_comparison.tsv"),
  "/mnt/data/67_candidate_drug_group_comparison.tsv"
)

priority_file_candidates <- c(
  file.path(table_dir, "67_candidate_drug_priority_table.tsv"),
  "/mnt/data/67_candidate_drug_priority_table.tsv"
)

pick_first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

cmp_file <- pick_first_existing(cmp_file_candidates)
priority_file <- pick_first_existing(priority_file_candidates)

if (is.na(cmp_file)) {
  stop("未找到 67_candidate_drug_group_comparison.tsv")
}
if (is.na(priority_file)) {
  stop("未找到 67_candidate_drug_priority_table.tsv")
}

drug_cmp <- fread(cmp_file)
priority_tab <- fread(priority_file)

# =========================================================
# Step 2. 固定药物方向
# 默认 predicted_response 越低 = 预测越敏感
# delta_high_minus_low > 0 代表 low-risk 更敏感
# delta_high_minus_low < 0 代表 high-risk 更敏感
# =========================================================
drug_direction <- drug_cmp %>%
  dplyr::mutate(
    response_direction = dplyr::case_when(
      preferred_in_group == "low_risk_more_sensitive" ~ "low-risk favored",
      preferred_in_group == "high_risk_more_sensitive" ~ "high-risk favored",
      TRUE ~ "no clear preference"
    ),
    significance_tier = dplyr::case_when(
      p_adj < 0.001 ~ "strong",
      p_adj < 0.05 ~ "moderate",
      TRUE ~ "weak_or_ns"
    )
  ) %>%
  dplyr::arrange(p_adj, p_value)

fwrite(
  drug_direction,
  file.path(result_dir, "68_drug_direction_summary.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 3. 生成药物生物学解释模板
# =========================================================
drug_bio <- drug_direction %>%
  dplyr::mutate(
    interpretation_axis = dplyr::case_when(
      drug %in% c("Cisplatin", "Oxaliplatin", "Gemcitabine", "Irinotecan") ~ "cytotoxic_or_DNA_damage_related",
      drug %in% c("Sorafenib") ~ "multi-kinase_stress_signal_related",
      drug %in% c("Trametinib", "Selumetinib") ~ "MAPK_pathway_related",
      drug %in% c("Dasatinib") ~ "SRC_family_kinase_related",
      TRUE ~ "other_or_unspecified"
    ),
    provisional_biological_interpretation = dplyr::case_when(
      response_direction == "low-risk favored" &
        drug %in% c("Cisplatin", "Oxaliplatin", "Gemcitabine", "Irinotecan") ~
        "Lower-risk tumors may retain relatively greater vulnerability to cytotoxic or replication-stress-associated treatments.",
      response_direction == "low-risk favored" &
        drug %in% c("Sorafenib") ~
        "Lower-risk tumors may remain more susceptible to multi-kinase inhibition or broader stress signaling perturbation.",
      response_direction == "high-risk favored" &
        drug %in% c("Trametinib", "Selumetinib") ~
        "Higher-risk tumors may exhibit a stronger dependency on MAPK-related signaling programs.",
      response_direction == "high-risk favored" &
        drug %in% c("Dasatinib") ~
        "Higher-risk tumors may harbor increased vulnerability to kinase signaling dependencies linked to aggressive transcriptional states.",
      TRUE ~
        "Requires additional biological interpretation and cross-checking with bulk/single-cell programs."
    ),
    recommended_use_in_manuscript = dplyr::case_when(
      significance_tier == "strong" ~ "main_or_high_priority_supplement",
      significance_tier == "moderate" ~ "supplement_or_secondary_main",
      TRUE ~ "supplement_only"
    )
  )

fwrite(
  drug_bio,
  file.path(result_dir, "68_drug_biological_interpretation.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 4. 与 bulk / single-cell 主线做桥接
# =========================================================
signature_12 <- c(
  "LCP1", "GBP4", "MET", "SNAPC4", "TGFBI", "SLC25A11",
  "RFC4", "IDO1", "NMUR1", "NDUFA7", "INHBA", "IL15RA"
)

bulk_singlecell_drug_bridge <- data.frame(
  bridge_item = c(
    "bulk_risk_model",
    "frozen_signature_genes",
    "single_cell_localization",
    "drug_result_core_pattern",
    "translational_interpretation",
    "manuscript_caution"
  ),
  bridge_text = c(
    "The frozen 12-gene signature captures a risk-associated transcriptional program in TCGA-PAAD.",
    paste(signature_12, collapse = ", "),
    "Single-cell analysis showed that this program was primarily localized to a subset of ductal-like malignant epithelial cells.",
    "Predicted drug response analysis suggested that several cytotoxic or multi-kinase agents were relatively favored in the low-risk group, whereas Trametinib showed an opposite pattern favoring the high-risk group.",
    "These findings imply that the prognostic program may not only stratify survival risk, but may also reflect distinct therapeutic vulnerabilities across transcriptionally defined risk states.",
    "All drug sensitivity findings should be described as transcriptome-inferred predicted response rather than confirmed clinical efficacy."
  ),
  stringsAsFactors = FALSE
)

fwrite(
  bulk_singlecell_drug_bridge,
  file.path(result_dir, "68_bulk_singlecell_drug_bridge.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 5. 药物重定位策略登记
# =========================================================
repositioning_strategy <- data.frame(
  item = c(
    "current_stage",
    "primary_goal",
    "first_layer_output",
    "second_layer_output",
    "current_decision",
    "next_stage"
  ),
  value = c(
    "drug_interpretation_and_repositioning_setup",
    "connect predicted drug response with bulk/single-cell biology",
    "ranked drug direction summary",
    "biological interpretation and manuscript-ready bridge",
    "prioritize interpretation of candidate-drug panel before broader reversal screening",
    "69_prepare_drug_result_paragraph_and_figure_layout"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  repositioning_strategy,
  file.path(registry_dir, "68_drug_repositioning_strategy.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 6. 写作 digest（修复 slice 冲突）
# =========================================================
top_low <- drug_bio %>%
  dplyr::filter(response_direction == "low-risk favored") %>%
  dplyr::arrange(p_adj) %>%
  dplyr::slice(1)

top_high <- drug_bio %>%
  dplyr::filter(response_direction == "high-risk favored") %>%
  dplyr::arrange(p_adj) %>%
  dplyr::slice(1)

writing_digest <- data.frame(
  item = c(
    "top_low_risk_favored_drug",
    "top_low_risk_favored_p_adj",
    "top_high_risk_favored_drug",
    "top_high_risk_favored_p_adj",
    "overall_interpretation"
  ),
  value = c(
    if (nrow(top_low) > 0) top_low$drug[1] else NA,
    if (nrow(top_low) > 0) top_low$p_adj[1] else NA,
    if (nrow(top_high) > 0) top_high$drug[1] else NA,
    if (nrow(top_high) > 0) top_high$p_adj[1] else NA,
    "Predicted drug sensitivity differed between transcriptomically defined risk groups, suggesting distinct therapeutic vulnerabilities."
  ),
  stringsAsFactors = FALSE
)

fwrite(
  writing_digest,
  file.path(result_dir, "68_drug_biology_writing_digest.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 7. 阶段摘要
# =========================================================
stage_summary <- data.frame(
  item = c(
    "stage68_completed",
    "drug_direction_summary_written",
    "drug_biological_interpretation_written",
    "bulk_singlecell_drug_bridge_written",
    "repositioning_strategy_written",
    "writing_digest_written",
    "next_stage"
  ),
  value = c(
    TRUE,
    file.exists(file.path(result_dir, "68_drug_direction_summary.tsv")),
    file.exists(file.path(result_dir, "68_drug_biological_interpretation.tsv")),
    file.exists(file.path(result_dir, "68_bulk_singlecell_drug_bridge.tsv")),
    file.exists(file.path(registry_dir, "68_drug_repositioning_strategy.tsv")),
    file.exists(file.path(result_dir, "68_drug_biology_writing_digest.tsv")),
    "69_prepare_drug_result_paragraph_and_figure_layout"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  stage_summary,
  file.path(registry_dir, "68_stage_summary.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 8. 日志
# =========================================================
writeLines(
  c(
    "68 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Interpret predicted drug sensitivity and connect it with bulk / single-cell biology.",
    "",
    "[Current key message]",
    "Drug response differences across transcriptomic risk groups may reflect distinct therapeutic vulnerabilities rather than uniform treatment sensitivity.",
    "",
    "[Next step]",
    "Prepare drug result paragraph and figure layout."
  ),
  file.path(log_dir, "68_drug_repositioning_and_interpretation_notes.txt")
)

cat("68_drug_repositioning_and_biological_interpretation.R finished successfully.\n")