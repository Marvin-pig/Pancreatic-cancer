# =========================
# 21_build_publication_ready_result_summary_tables.R
# 目的：
# 1. 汇总训练集与 GSE21501 外部验证的核心结果
# 2. 生成可直接用于论文 Results / Supplementary Table 的汇总表
# 3. 为后续结果写作与作图准备统一表格
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})
setwd("/Users/wmz_mac/Desktop/胰腺癌")

out_dir <- "04_bulk_analysis/05_publication_summary"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 输入 ----------
train_cindex_file   <- "04_bulk_analysis/03_survival_model/14_cindex_training.tsv"
train_timeROC_file  <- "04_bulk_analysis/03_survival_model/14_timeROC_training.tsv"
train_multi_file    <- "04_bulk_analysis/03_survival_model/14_multivcox_clinical_and_risk.tsv"
train_formula_file  <- "04_bulk_analysis/03_survival_model/15_final_signature_model.tsv"

ext_summary_file    <- "04_bulk_analysis/04_external_validation/04_signature_validation/GSE21501/20_GSE21501_refined_validation_summary.tsv"
ext_timeROC_file    <- "04_bulk_analysis/04_external_validation/04_signature_validation/GSE21501/20_GSE21501_timeROC_months.tsv"
ext_group_file      <- "04_bulk_analysis/04_external_validation/04_signature_validation/GSE21501/20_GSE21501_univcox_external_median_group.tsv"
ext_score_file      <- "04_bulk_analysis/04_external_validation/04_signature_validation/GSE21501/20_GSE21501_univcox_risk_score_confirm.tsv"

coef_file           <- "04_bulk_analysis/03_survival_model/15_final_signature_coefficients.tsv"

# ---------- 读取 ----------
train_cindex  <- fread(train_cindex_file)
train_roc     <- fread(train_timeROC_file)
train_multi   <- fread(train_multi_file)
train_formula <- fread(train_formula_file)

ext_summary   <- fread(ext_summary_file)
ext_roc       <- fread(ext_timeROC_file)
ext_group     <- fread(ext_group_file)
ext_score     <- fread(ext_score_file)

coef_df       <- fread(coef_file)

# ---------- 表1：signature 系数表 ----------
coef_tab <- coef_df %>%
  mutate(
    direction = ifelse(coefficient > 0, "Risk", "Protective")
  ) %>%
  arrange(desc(abs(coefficient)))

fwrite(coef_tab,
       file.path(out_dir, "21_signature_coefficients_publication.tsv"),
       sep = "\t")

# ---------- 表2：训练集 vs 外部验证结果总表 ----------
risk_multi <- train_multi %>% filter(term == "risk_score")

summary_tab <- data.frame(
  cohort = c("TCGA_train", "GSE21501_external"),
  n_samples = c(168, ext_summary$value[ext_summary$metric == "n_external_samples"]),
  n_events = c(93, ext_summary$value[ext_summary$metric == "n_external_events"]),
  c_index = c(train_cindex$c_index[1], NA),
  auc_1 = c(train_roc$AUC[train_roc$time_days == 365], ext_roc$AUC[ext_roc$time_months == 12]),
  auc_2 = c(train_roc$AUC[train_roc$time_days == 730], ext_roc$AUC[ext_roc$time_months == 24]),
  auc_3 = c(train_roc$AUC[train_roc$time_days == 1095], ext_roc$AUC[ext_roc$time_months == 36]),
  continuous_HR = c(risk_multi$estimate[1], ext_score$estimate[1]),
  continuous_p = c(risk_multi$p.value[1], ext_score$p.value[1]),
  stringsAsFactors = FALSE
)

fwrite(summary_tab,
       file.path(out_dir, "21_train_vs_external_validation_summary.tsv"),
       sep = "\t")

# ---------- 表3：外部验证敏感性分析 ----------
external_sensitivity_tab <- data.frame(
  analysis = c(
    "training_cutoff_grouping",
    "external_median_grouping",
    "continuous_risk_score"
  ),
  HR = c(
    NA,
    ext_group$estimate[1],
    ext_score$estimate[1]
  ),
  p_value = c(
    NA,
    ext_group$p.value[1],
    ext_score$p.value[1]
  ),
  note = c(
    "distribution shift caused severe imbalance",
    "balanced grouping but non-significant",
    "direction consistent but non-significant"
  ),
  stringsAsFactors = FALSE
)

fwrite(external_sensitivity_tab,
       file.path(out_dir, "21_external_validation_sensitivity_summary.tsv"),
       sep = "\t")

# ---------- 表4：模型信息 ----------
model_tab <- train_formula
fwrite(model_tab,
       file.path(out_dir, "21_final_model_information.tsv"),
       sep = "\t")

cat("21_build_publication_ready_result_summary_tables.R finished successfully.\n")