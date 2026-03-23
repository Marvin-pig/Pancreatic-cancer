# =========================
# 45_summarise_mr_results_and_assign_reporting_tiers.R
# 目的：
# 1. 汇总 44 的 MR 主结果、异质性、多效性结果
# 2. 对每个基因进行“主文/补充/当前不支持”分层
# 3. 生成可直接用于结果段写作的汇总表
# 4. 为后续 PanScan/PanC4 补强提供优先级列表
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ---------- 路径 ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "05_MR"

result_dir   <- file.path(root_dir, "03_results")
registry_dir <- file.path(root_dir, "00_registry")
proc_dir     <- file.path(root_dir, "02_processed_data")
log_dir      <- file.path(root_dir, "06_logs")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

mr_file      <- file.path(result_dir, "44_mr_main_results.tsv")
het_file     <- file.path(result_dir, "44_mr_heterogeneity.tsv")
pleio_file   <- file.path(result_dir, "44_mr_pleiotropy.tsv")
report_file  <- file.path(registry_dir, "44_reporting_decision.tsv")
status_file  <- file.path(registry_dir, "44_status_summary.tsv")
input_file   <- file.path(registry_dir, "44_gene_input_summary.tsv")

if (!file.exists(mr_file)) stop("缺少 44_mr_main_results.tsv: ", mr_file)
if (!file.exists(report_file)) stop("缺少 44_reporting_decision.tsv: ", report_file)

mr_res <- fread(mr_file)
report_decision <- fread(report_file)

het_res <- if (file.exists(het_file)) fread(het_file) else data.table()
pleio_res <- if (file.exists(pleio_file)) fread(pleio_file) else data.table()
gene_input <- if (file.exists(input_file)) fread(input_file) else data.table()

# ---------- 基础检查 ----------
if (!"gene" %in% names(mr_res)) stop("44_mr_main_results.tsv 缺少 gene 列。")
if (!"gene" %in% names(report_decision)) stop("44_reporting_decision.tsv 缺少 gene 列。")

mr_res$gene <- as.character(mr_res$gene)
report_decision$gene <- as.character(report_decision$gene)

# ---------- 统一方法命名 ----------
mr_res <- mr_res %>%
  mutate(
    method = as.character(method),
    method_simple = case_when(
      str_detect(method, regex("inverse variance weighted|ivw", ignore_case = TRUE)) ~ "IVW",
      str_detect(method, regex("weighted median", ignore_case = TRUE)) ~ "Weighted_median",
      str_detect(method, regex("egger", ignore_case = TRUE)) ~ "MR_Egger",
      str_detect(method, regex("wald ratio", ignore_case = TRUE)) ~ "Wald_ratio",
      TRUE ~ method
    )
  )

# ---------- 每个基因优先方法 ----------
# 优先级：
# 1. IVW
# 2. Weighted median
# 3. MR-Egger
# 4. Wald ratio
method_priority <- c("IVW", "Weighted_median", "MR_Egger", "Wald_ratio")

mr_res <- mr_res %>%
  mutate(
    method_rank = match(method_simple, method_priority),
    method_rank = ifelse(is.na(method_rank), 999, method_rank)
  )

lead_result <- mr_res %>%
  arrange(gene, method_rank, pval) %>%
  group_by(gene) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    gene,
    lead_method = method_simple,
    lead_nsnp = if ("nsnp" %in% names(.)) nsnp else NA_integer_,
    lead_b = if ("b" %in% names(.)) b else NA_real_,
    lead_se = if ("se" %in% names(.)) se else NA_real_,
    lead_pval = if ("pval" %in% names(.)) pval else NA_real_,
    lead_OR = if ("OR" %in% names(.)) OR else if ("b" %in% names(.)) exp(b) else NA_real_,
    lead_OR_lci95 = if ("OR_lci95" %in% names(.)) OR_lci95 else if (all(c("b","se") %in% names(.))) exp(b - 1.96*se) else NA_real_,
    lead_OR_uci95 = if ("OR_uci95" %in% names(.)) OR_uci95 else if (all(c("b","se") %in% names(.))) exp(b + 1.96*se) else NA_real_
  )

# ---------- 多 SNP 基因主结果摘要 ----------
multi_snp_summary <- mr_res %>%
  filter(method_simple %in% c("IVW", "Weighted_median", "MR_Egger")) %>%
  group_by(gene) %>%
  summarise(
    ivw_p = ifelse(any(method_simple == "IVW"), pval[match("IVW", method_simple)], NA_real_),
    wm_p  = ifelse(any(method_simple == "Weighted_median"), pval[match("Weighted_median", method_simple)], NA_real_),
    egger_p = ifelse(any(method_simple == "MR_Egger"), pval[match("MR_Egger", method_simple)], NA_real_),
    ivw_b = ifelse(any(method_simple == "IVW"), b[match("IVW", method_simple)], NA_real_),
    wm_b  = ifelse(any(method_simple == "Weighted_median"), b[match("Weighted_median", method_simple)], NA_real_),
    egger_b = ifelse(any(method_simple == "MR_Egger"), b[match("MR_Egger", method_simple)], NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    sign_consistent_main_methods = case_when(
      !is.na(ivw_b) & !is.na(wm_b) ~ sign(ivw_b) == sign(wm_b),
      TRUE ~ NA
    ),
    any_main_method_sig = case_when(
      (!is.na(ivw_p) & ivw_p < 0.05) |
        (!is.na(wm_p) & wm_p < 0.05) |
        (!is.na(egger_p) & egger_p < 0.05) ~ TRUE,
      TRUE ~ FALSE
    )
  )

# ---------- 异质性摘要 ----------
if (nrow(het_res) > 0) {
  het_res <- het_res %>%
    mutate(
      method = as.character(method),
      method_simple = case_when(
        str_detect(method, regex("inverse variance weighted|ivw", ignore_case = TRUE)) ~ "IVW",
        str_detect(method, regex("egger", ignore_case = TRUE)) ~ "MR_Egger",
        TRUE ~ method
      )
    )
  
  het_summary <- het_res %>%
    group_by(gene) %>%
    summarise(
      ivw_heterogeneity_p = ifelse(any(method_simple == "IVW"), Q_pval[match("IVW", method_simple)], NA_real_),
      egger_heterogeneity_p = ifelse(any(method_simple == "MR_Egger"), Q_pval[match("MR_Egger", method_simple)], NA_real_),
      .groups = "drop"
    )
} else {
  het_summary <- data.frame(
    gene = unique(mr_res$gene),
    ivw_heterogeneity_p = NA_real_,
    egger_heterogeneity_p = NA_real_,
    stringsAsFactors = FALSE
  )
}

# ---------- 多效性摘要 ----------
if (nrow(pleio_res) > 0) {
  pleio_summary <- pleio_res %>%
    transmute(
      gene = as.character(gene),
      egger_intercept = if ("egger_intercept" %in% names(.)) egger_intercept else NA_real_,
      egger_intercept_se = if ("se" %in% names(.)) se else NA_real_,
      egger_intercept_p = if ("pval" %in% names(.)) pval else NA_real_
    )
} else {
  pleio_summary <- data.frame(
    gene = unique(mr_res$gene),
    egger_intercept = NA_real_,
    egger_intercept_se = NA_real_,
    egger_intercept_p = NA_real_,
    stringsAsFactors = FALSE
  )
}

# ---------- 合并总表 ----------
summary_tab <- report_decision %>%
  select(gene, n_mr_keep, reporting_tier, recommended_section) %>%
  left_join(lead_result, by = "gene") %>%
  left_join(multi_snp_summary, by = "gene") %>%
  left_join(het_summary, by = "gene") %>%
  left_join(pleio_summary, by = "gene")

if (nrow(gene_input) > 0 && "gene" %in% names(gene_input)) {
  summary_tab <- summary_tab %>%
    left_join(
      gene_input %>% select(gene, outcome, outcome_id, exposure),
      by = "gene"
    )
}

# ---------- 证据分层规则 ----------
summary_tab <- summary_tab %>%
  mutate(
    evidence_class = case_when(
      n_mr_keep >= 3 &
        !is.na(ivw_p) & ivw_p < 0.05 &
        (is.na(wm_p) | wm_p < 0.05 | is.na(sign_consistent_main_methods) | sign_consistent_main_methods) &
        (is.na(egger_intercept_p) | egger_intercept_p >= 0.05) ~ "supported_multi_snp_signal",
      
      n_mr_keep >= 3 &
        !is.na(ivw_p) & ivw_p >= 0.05 &
        (is.na(wm_p) | wm_p >= 0.05) ~ "multi_snp_but_no_support",
      
      n_mr_keep == 2 &
        !is.na(lead_pval) & lead_pval < 0.05 ~ "few_snp_exploratory_positive",
      
      n_mr_keep == 2 &
        (is.na(lead_pval) | lead_pval >= 0.05) ~ "few_snp_exploratory_negative",
      
      n_mr_keep == 1 &
        !is.na(lead_pval) & lead_pval < 0.05 ~ "single_snp_positive",
      
      n_mr_keep == 1 &
        (is.na(lead_pval) | lead_pval >= 0.05) ~ "single_snp_negative",
      
      TRUE ~ "insufficient_or_unclassified"
    ),
    writing_priority = case_when(
      evidence_class == "supported_multi_snp_signal" ~ "main_text_priority",
      evidence_class == "multi_snp_but_no_support" ~ "main_text_negative_or_neutral",
      evidence_class %in% c("few_snp_exploratory_positive", "few_snp_exploratory_negative") ~ "supplementary_exploratory",
      evidence_class %in% c("single_snp_positive", "single_snp_negative") ~ "supplementary_wald_ratio",
      TRUE ~ "not_for_emphasis"
    )
  )

# ---------- 需要补强的优先级 ----------
summary_tab <- summary_tab %>%
  mutate(
    next_action_priority = case_when(
      writing_priority %in% c("main_text_priority", "main_text_negative_or_neutral") ~ "prioritise_primary_outcome_replication",
      writing_priority == "supplementary_exploratory" ~ "consider_primary_outcome_replication_if_biologically_important",
      writing_priority == "supplementary_wald_ratio" ~ "do_not_overinterpret_need_more_iv_or_other_evidence",
      TRUE ~ "low_priority"
    )
  )

# ---------- 结果段写作用摘要 ----------
writing_digest <- summary_tab %>%
  transmute(
    gene,
    n_mr_keep,
    lead_method,
    lead_OR = round(lead_OR, 3),
    lead_OR_lci95 = round(lead_OR_lci95, 3),
    lead_OR_uci95 = round(lead_OR_uci95, 3),
    lead_pval = signif(lead_pval, 3),
    ivw_p = signif(ivw_p, 3),
    wm_p = signif(wm_p, 3),
    egger_p = signif(egger_p, 3),
    ivw_heterogeneity_p = signif(ivw_heterogeneity_p, 3),
    egger_intercept_p = signif(egger_intercept_p, 3),
    evidence_class,
    writing_priority,
    next_action_priority
  )

# ---------- 主文候选 / 补充候选 ----------
main_text_candidates <- summary_tab %>%
  filter(writing_priority %in% c("main_text_priority", "main_text_negative_or_neutral"))

supplementary_candidates <- summary_tab %>%
  filter(writing_priority %in% c("supplementary_exploratory", "supplementary_wald_ratio"))

# ---------- 生成一句话解释模板 ----------
interpret_template <- function(n_mr_keep, evidence_class, gene, lead_method, lead_OR, lci, uci, pval) {
  if (is.na(lead_method)) {
    return(paste0(gene, ": 当前缺少可解释的 MR 结果。"))
  }
  
  if (evidence_class == "supported_multi_snp_signal") {
    return(sprintf(
      "%s 在 BBJ pancreatic cancer 路线中显示多 SNP MR 支持的因果信号（%s, OR=%.3f, 95%%CI %.3f-%.3f, P=%s）。",
      gene, lead_method, lead_OR, lci, uci, format(signif(pval, 3), scientific = FALSE)
    ))
  }
  
  if (evidence_class == "multi_snp_but_no_support") {
    return(sprintf(
      "%s 具备多 SNP MR 分析条件，但在 BBJ 路线中未观察到稳定的因果支持（%s, OR=%.3f, 95%%CI %.3f-%.3f, P=%s）。",
      gene, lead_method, lead_OR, lci, uci, format(signif(pval, 3), scientific = FALSE)
    ))
  }
  
  if (evidence_class == "few_snp_exploratory_positive") {
    return(sprintf(
      "%s 仅保留少量工具变量，当前结果提示探索性阳性信号，但证据层级有限（%s, OR=%.3f, 95%%CI %.3f-%.3f, P=%s）。",
      gene, lead_method, lead_OR, lci, uci, format(signif(pval, 3), scientific = FALSE)
    ))
  }
  
  if (evidence_class == "few_snp_exploratory_negative") {
    return(sprintf(
      "%s 仅保留少量工具变量，当前探索性分析未见显著因果支持（%s, OR=%.3f, 95%%CI %.3f-%.3f, P=%s）。",
      gene, lead_method, lead_OR, lci, uci, format(signif(pval, 3), scientific = FALSE)
    ))
  }
  
  if (evidence_class == "single_snp_positive") {
    return(sprintf(
      "%s 当前仅能进行单 SNP Wald ratio 分析，结果提示阳性信号，但需谨慎解释（OR=%.3f, 95%%CI %.3f-%.3f, P=%s）。",
      gene, lead_OR, lci, uci, format(signif(pval, 3), scientific = FALSE)
    ))
  }
  
  if (evidence_class == "single_snp_negative") {
    return(sprintf(
      "%s 当前仅能进行单 SNP Wald ratio 分析，未见显著因果支持（OR=%.3f, 95%%CI %.3f-%.3f, P=%s）。",
      gene, lead_OR, lci, uci, format(signif(pval, 3), scientific = FALSE)
    ))
  }
  
  paste0(gene, ": 当前证据不足，暂不建议强调。")
}

text_digest <- summary_tab %>%
  rowwise() %>%
  mutate(
    interpretation_sentence = interpret_template(
      n_mr_keep, evidence_class, gene, lead_method,
      lead_OR, lead_OR_lci95, lead_OR_uci95, lead_pval
    )
  ) %>%
  ungroup() %>%
  select(gene, interpretation_sentence)

# ---------- 输出 ----------
fwrite(
  summary_tab,
  file.path(result_dir, "45_mr_summary_master_table.tsv"),
  sep = "\t", na = "NA"
)

fwrite(
  writing_digest,
  file.path(result_dir, "45_mr_writing_digest.tsv"),
  sep = "\t", na = "NA"
)

fwrite(
  main_text_candidates,
  file.path(result_dir, "45_mr_main_text_candidates.tsv"),
  sep = "\t", na = "NA"
)

fwrite(
  supplementary_candidates,
  file.path(result_dir, "45_mr_supplementary_candidates.tsv"),
  sep = "\t", na = "NA"
)

fwrite(
  text_digest,
  file.path(result_dir, "45_mr_interpretation_sentences.tsv"),
  sep = "\t", na = "NA"
)

# ---------- 总状态 ----------
status_summary <- data.frame(
  item = c(
    "n_genes_with_mr_results",
    "n_main_text_priority",
    "n_main_text_negative_or_neutral",
    "n_supplementary_exploratory",
    "n_supplementary_wald_ratio",
    "recommended_next_action"
  ),
  value = c(
    nrow(summary_tab),
    sum(summary_tab$writing_priority == "main_text_priority", na.rm = TRUE),
    sum(summary_tab$writing_priority == "main_text_negative_or_neutral", na.rm = TRUE),
    sum(summary_tab$writing_priority == "supplementary_exploratory", na.rm = TRUE),
    sum(summary_tab$writing_priority == "supplementary_wald_ratio", na.rm = TRUE),
    ifelse(
      sum(summary_tab$writing_priority %in% c("main_text_priority", "main_text_negative_or_neutral"), na.rm = TRUE) > 0,
      "Proceed to primary outcome replication (PanScan/PanC4) and draft MR result text",
      "Current BBJ MR is mainly supplementary; prioritise primary outcome replication"
    )
  ),
  stringsAsFactors = FALSE
)

fwrite(
  status_summary,
  file.path(registry_dir, "45_status_summary.tsv"),
  sep = "\t", na = "NA"
)

# ---------- 日志 ----------
log_lines <- c(
  "=== Code 45 Summary ===",
  paste0("Created at: ", Sys.time()),
  "",
  "[Goal]",
  "Summarise MR results into manuscript-oriented evidence tiers and reporting priorities.",
  "",
  "[Outputs]",
  "- 45_mr_summary_master_table.tsv",
  "- 45_mr_writing_digest.tsv",
  "- 45_mr_main_text_candidates.tsv",
  "- 45_mr_supplementary_candidates.tsv",
  "- 45_mr_interpretation_sentences.tsv",
  "- 45_status_summary.tsv"
)

writeLines(
  log_lines,
  file.path(log_dir, "45_summarise_mr_results_notes.txt")
)

cat("45_summarise_mr_results_and_assign_reporting_tiers.R finished successfully.\n")
cat("Generated files:\n")
cat("- 03_results/45_mr_summary_master_table.tsv\n")
cat("- 03_results/45_mr_writing_digest.tsv\n")
cat("- 03_results/45_mr_main_text_candidates.tsv\n")
cat("- 03_results/45_mr_supplementary_candidates.tsv\n")
cat("- 03_results/45_mr_interpretation_sentences.tsv\n")
cat("- 00_registry/45_status_summary.tsv\n")