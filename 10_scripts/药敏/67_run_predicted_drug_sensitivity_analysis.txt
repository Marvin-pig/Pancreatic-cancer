# =========================
# 67_run_predicted_drug_sensitivity_analysis.R
# 目的：
# 1. 基于 TCGA-PAAD 冻结风险分组，使用 oncoPredict 进行 predicted drug sensitivity
# 2. 比较 high-risk vs low-risk 在候选药物上的预测反应差异
# 3. 输出候选药物优先级表与箱线图
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
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

# =========================================================
# Step 1. 固定你的输入文件路径
# =========================================================
expr_file <- "/Users/wmz_mac/Desktop/胰腺癌/07_drug/02_processed_data/67b_tcga_paad_expr_for_drug.tsv"
risk_file <- "/Users/wmz_mac/Desktop/胰腺癌/07_drug/02_processed_data/67b_tcga_paad_risk_group_for_drug.tsv"

if (!file.exists(expr_file)) {
  stop("未找到表达矩阵文件: ", expr_file)
}
if (!file.exists(risk_file)) {
  stop("未找到风险分组文件: ", risk_file)
}

# =========================================================
# Step 2. 读取并重建 test_expr
# 关键：先合并重复 gene，再设置 rownames
# =========================================================
expr <- data.table::fread(
  expr_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

gene_col <- names(expr)[1]

gene_vec <- as.character(expr[[gene_col]])
gene_vec <- trimws(gene_vec)
gene_vec <- sub("\\..*$", "", gene_vec)     # 去版本号
gene_vec <- toupper(gene_vec)               # 统一大写 symbol

expr_only <- expr[, setdiff(names(expr), gene_col), drop = FALSE]

for (j in seq_len(ncol(expr_only))) {
  expr_only[[j]] <- suppressWarnings(as.numeric(expr_only[[j]]))
}

expr_clean <- data.frame(
  gene_name = gene_vec,
  expr_only,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expr_clean <- expr_clean %>%
  filter(!is.na(gene_name), gene_name != "")

# 合并重复 gene
expr_clean <- expr_clean %>%
  group_by(gene_name) %>%
  summarise(
    across(everything(), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

expr_mat_df <- as.data.frame(expr_clean, stringsAsFactors = FALSE)
rownames(expr_mat_df) <- expr_mat_df$gene_name
expr_mat_df$gene_name <- NULL

test_expr <- as.matrix(expr_mat_df)
mode(test_expr) <- "numeric"

# 去掉全 NA / 零方差基因
keep_gene <- apply(test_expr, 1, function(x) {
  sum(!is.na(x)) > 0 && stats::var(x, na.rm = TRUE) > 0
})
test_expr <- test_expr[keep_gene, , drop = FALSE]

cat("==== rebuilt test_expr rownames preview ====\n")
print(head(rownames(test_expr), 10))
cat("==== n test genes ====\n")
print(nrow(test_expr))

# =========================================================
# Step 3. 读取风险分组
# =========================================================
risk <- data.table::fread(
  risk_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

if (!all(c("sample_id", "risk_group") %in% colnames(risk))) {
  stop("risk_file 必须包含列：sample_id 和 risk_group")
}

risk <- risk %>%
  mutate(
    sample_id = as.character(sample_id),
    risk_group = tolower(as.character(risk_group))
  ) %>%
  filter(risk_group %in% c("high", "low")) %>%
  distinct(sample_id, .keep_all = TRUE)

common_samples <- intersect(colnames(test_expr), risk$sample_id)

if (length(common_samples) < 20) {
  stop("test_expr 与风险分组重叠样本过少，请检查 sample_id 一致性。")
}

test_expr <- test_expr[, common_samples, drop = FALSE]
risk2 <- risk[match(common_samples, risk$sample_id), ]

alignment_registry <- data.frame(
  item = c(
    "n_expr_samples",
    "n_risk_samples",
    "n_common_samples",
    "risk_high_n",
    "risk_low_n"
  ),
  value = c(
    ncol(expr_mat_df),
    nrow(risk),
    ncol(test_expr),
    sum(risk2$risk_group == "high"),
    sum(risk2$risk_group == "low")
  ),
  stringsAsFactors = FALSE
)

fwrite(
  alignment_registry,
  file.path(registry_dir, "67_sample_alignment_registry.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 4. 读取 oncoPredict 完整训练数据
# 关键：禁止 short/demo 版本
# =========================================================
if (!requireNamespace("oncoPredict", quietly = TRUE)) {
  stop("未安装 oncoPredict，请先安装。")
}

gdsc_expr_candidates <- c(
  "/Users/wmz_mac/Desktop/胰腺癌/07_drug/01_raw_data/GDSC2_Expr.rds",
  "/Users/wmz_mac/Desktop/胰腺癌/07_drug/01_raw_data/GDSC2_Expr_full.rds",
  "/Users/wmz_mac/Desktop/胰腺癌/oncoPredict-main/vignettes/GDSC2_Expr.rds"
)

gdsc_res_candidates <- c(
  "/Users/wmz_mac/Desktop/胰腺癌/07_drug/01_raw_data/GDSC2_Res.rds",
  "/Users/wmz_mac/Desktop/胰腺癌/07_drug/01_raw_data/GDSC2_Res_full.rds",
  "/Users/wmz_mac/Desktop/胰腺癌/oncoPredict-main/vignettes/GDSC2_Res.rds"
)

pick_first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

training_expr_file <- pick_first_existing(gdsc_expr_candidates)
training_res_file  <- pick_first_existing(gdsc_res_candidates)

if (is.na(training_expr_file) || is.na(training_res_file)) {
  stop(
    paste0(
      "未找到完整 GDSC 训练数据。\n",
      "请确保以下文件至少存在一组：\n",
      "- /Users/wmz_mac/Desktop/胰腺癌/07_drug/01_raw_data/GDSC2_Expr.rds\n",
      "- /Users/wmz_mac/Desktop/胰腺癌/07_drug/01_raw_data/GDSC2_Res.rds\n",
      "或放在 oncoPredict-main/vignettes 目录下。"
    )
  )
}

training_expr <- readRDS(training_expr_file)
training_res  <- readRDS(training_res_file)

if (!is.matrix(training_expr)) training_expr <- as.matrix(training_expr)
if (!is.matrix(training_res))  training_res  <- as.matrix(training_res)

mode(training_expr) <- "numeric"
mode(training_res)  <- "numeric"

cat("==== training_expr_file ====\n")
print(training_expr_file)
cat("==== training_res_file ====\n")
print(training_res_file)
cat("==== dim(training_expr) ====\n")
print(dim(training_expr))
cat("==== dim(training_res) ====\n")
print(dim(training_res))

# 强制阻止 short/demo 版本
if (nrow(training_expr) < 5000) {
  stop(
    paste0(
      "当前 training_expr 只有 ", nrow(training_expr),
      " 个基因，明显不是完整 GDSC 训练表达矩阵。请改用完整 GDSC2_Expr 文件。"
    )
  )
}

# =========================================================
# Step 5. 清洗 training_expr rownames，并统一 gene symbol 大写
# =========================================================
tr_genes <- rownames(training_expr)
tr_genes <- sub("\\..*$", "", tr_genes)
tr_genes <- toupper(tr_genes)
rownames(training_expr) <- tr_genes

training_expr <- training_expr[
  !is.na(rownames(training_expr)) &
    rownames(training_expr) != "" &
    !duplicated(rownames(training_expr)),
  , drop = FALSE
]

test_expr <- test_expr[
  !is.na(rownames(test_expr)) &
    rownames(test_expr) != "" &
    !duplicated(rownames(test_expr)),
  , drop = FALSE
]

cat("==== training rownames preview ====\n")
print(head(rownames(training_expr), 10))
cat("==== test rownames preview ====\n")
print(head(rownames(test_expr), 10))

common_genes <- intersect(rownames(training_expr), rownames(test_expr))
cat("==== number of common genes ====\n")
print(length(common_genes))

if (length(common_genes) < 3000) {
  stop(
    paste0(
      "training/test 共用基因太少（", length(common_genes),
      "）。当前仍不适合跑 oncoPredict，请检查训练矩阵是否为完整 GDSC2_Expr。"
    )
  )
}

training_expr2 <- training_expr[common_genes, , drop = FALSE]
test_expr2     <- test_expr[common_genes, , drop = FALSE]

# =========================================================
# Step 6. 对齐 training_expr 与 training_res 样本
# =========================================================
common_train_samples <- intersect(colnames(training_expr2), rownames(training_res))

if (length(common_train_samples) >= 20) {
  training_expr2 <- training_expr2[, common_train_samples, drop = FALSE]
  training_res   <- training_res[common_train_samples, , drop = FALSE]
} else {
  common_train_samples <- intersect(colnames(training_expr2), colnames(training_res))
  if (length(common_train_samples) < 20) {
    stop("training_expr 与 training_res 的样本无法对齐，请检查 GDSC2_Expr / GDSC2_Res 的格式。")
  }
  training_expr2 <- training_expr2[, common_train_samples, drop = FALSE]
  training_res   <- training_res[, common_train_samples, drop = FALSE]
}

cat("==== aligned training samples ====\n")
print(length(common_train_samples))

# =========================================================
# Step 7. 跑 oncoPredict
# =========================================================
predicted_all <- oncoPredict::calcPhenotype(
  trainingExprData = training_expr2,
  trainingPtype = training_res,
  testExprData = test_expr2,
  batchCorrect = "eb",
  powerTransformPhenotype = TRUE,
  removeLowVaryingGenes = 0.2,
  minNumSamples = 10,
  selection = 1,
  printOutput = TRUE,
  removeLowVaringGenesFrom = "rawData",
  folder = FALSE
)

if (is.null(predicted_all)) {
  stop("oncoPredict 未返回有效结果。")
}

predicted_all <- as.data.frame(predicted_all)

# =========================================================
# =========================================================
# =========================================================
# =========================================================
# Step 8. 统一输出方向
# 常见情况：行=药物，列=样本；若相反则转置
# =========================================================
predicted_all <- as.data.frame(predicted_all, stringsAsFactors = FALSE)

# 清掉可能残留的旧列，避免 pivot_longer 报错
drop_cols <- intersect(c("drug", "drug_full", "drug_base"), colnames(predicted_all))
if (length(drop_cols) > 0) {
  predicted_all <- predicted_all[, setdiff(colnames(predicted_all), drop_cols), drop = FALSE]
}

sample_overlap_in_col <- sum(colnames(predicted_all) %in% risk2$sample_id)
sample_overlap_in_row <- sum(rownames(predicted_all) %in% risk2$sample_id)

if (sample_overlap_in_row > sample_overlap_in_col) {
  predicted_all <- as.data.frame(t(as.matrix(predicted_all)), stringsAsFactors = FALSE)
}

predicted_all$drug_full <- rownames(predicted_all)

# 去掉末尾的 _数字，得到基础药物名
predicted_all$drug_base <- sub("_[0-9]+$", "", predicted_all$drug_full)

fwrite(
  predicted_all,
  file.path(proc_dir, "67_predicted_drug_response_raw.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 9. 候选药物集合
# =========================================================
candidate_drugs <- c(
  "Erlotinib",
  "Gemcitabine",
  "Paclitaxel",
  "Cisplatin",
  "Oxaliplatin",
  "5-Fluorouracil",
  "Irinotecan",
  "Trametinib",
  "Selumetinib",
  "Bortezomib",
  "Sorafenib",
  "Dasatinib"
)

fwrite(
  data.frame(drug = candidate_drugs, stringsAsFactors = FALSE),
  file.path(registry_dir, "67_candidate_drug_list.tsv"),
  sep = "\t", na = "NA"
)

# 检查候选药物是否在输出里出现
drug_match_registry <- data.frame(
  candidate_drug = candidate_drugs,
  matched = candidate_drugs %in% unique(predicted_all$drug_base),
  stringsAsFactors = FALSE
)

fwrite(
  drug_match_registry,
  file.path(registry_dir, "67_candidate_drug_match_registry.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 10. 转成长表并与风险分组合并
# 注意：按 drug_base 匹配，而不是 drug_full
# =========================================================
pred_long <- predicted_all %>%
  filter(drug_base %in% candidate_drugs) %>%
  pivot_longer(
    cols = -c(drug_full, drug_base),
    names_to = "sample_id",
    values_to = "predicted_response"
  ) %>%
  left_join(risk2 %>% select(sample_id, risk_group), by = "sample_id") %>%
  filter(!is.na(risk_group))

fwrite(
  pred_long,
  file.path(proc_dir, "67_candidate_drug_response_long.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 11. high vs low 比较
# 先按 drug_full 分别比较，再按 drug_base 汇总
# =========================================================
drug_compare_full <- bind_rows(lapply(unique(pred_long$drug_full), function(d) {
  sub <- pred_long %>% filter(drug_full == d)
  
  x_high <- sub$predicted_response[sub$risk_group == "high"]
  x_low  <- sub$predicted_response[sub$risk_group == "low"]
  
  if (length(x_high) < 3 || length(x_low) < 3) return(NULL)
  
  wt <- wilcox.test(x_high, x_low)
  
  data.frame(
    drug_full = d,
    drug_base = unique(sub$drug_base)[1],
    n_high = length(x_high),
    n_low = length(x_low),
    median_high = median(x_high, na.rm = TRUE),
    median_low = median(x_low, na.rm = TRUE),
    delta_high_minus_low = median(x_high, na.rm = TRUE) - median(x_low, na.rm = TRUE),
    p_value = wt$p.value,
    stringsAsFactors = FALSE
  )
}))

if (nrow(drug_compare_full) == 0) {
  stop("候选药物比较结果为空，请检查候选药物名与 oncoPredict 输出。")
}

# 强制转成普通 data.frame，避免 Rle/list 兼容问题
drug_compare_full <- as.data.frame(drug_compare_full, stringsAsFactors = FALSE)

drug_compare_full$drug_full <- as.character(drug_compare_full$drug_full)
drug_compare_full$drug_base <- as.character(drug_compare_full$drug_base)
drug_compare_full$n_high <- as.numeric(drug_compare_full$n_high)
drug_compare_full$n_low <- as.numeric(drug_compare_full$n_low)
drug_compare_full$median_high <- as.numeric(drug_compare_full$median_high)
drug_compare_full$median_low <- as.numeric(drug_compare_full$median_low)
drug_compare_full$delta_high_minus_low <- as.numeric(drug_compare_full$delta_high_minus_low)
drug_compare_full$p_value <- as.numeric(drug_compare_full$p_value)

drug_compare_full$p_adj <- p.adjust(drug_compare_full$p_value, method = "BH")
drug_compare_full$preferred_in_group <- ifelse(
  drug_compare_full$delta_high_minus_low < 0,
  "high_risk_more_sensitive",
  ifelse(
    drug_compare_full$delta_high_minus_low > 0,
    "low_risk_more_sensitive",
    "no_difference"
  )
)

drug_compare_full <- drug_compare_full[order(drug_compare_full$p_adj, drug_compare_full$p_value), ]

fwrite(
  drug_compare_full,
  file.path(result_dir, "67_candidate_drug_group_comparison_full.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 12. 对同一种基础药物汇总
# 例如 Oxaliplatin 可能有多个版本，取 p 最小的一条作为代表
# =========================================================
split_list <- split(drug_compare_full, drug_compare_full$drug_base)

drug_compare <- do.call(
  rbind,
  lapply(split_list, function(df) {
    df <- df[order(df$p_adj, df$p_value), , drop = FALSE]
    df[1, , drop = FALSE]
  })
)

drug_compare <- as.data.frame(drug_compare, stringsAsFactors = FALSE)
drug_compare$drug <- drug_compare$drug_base
drug_compare <- drug_compare[order(drug_compare$p_adj, drug_compare$p_value), ]

fwrite(
  drug_compare,
  file.path(result_dir, "67_candidate_drug_group_comparison.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 13. 候选药物优先级表
# =========================================================
priority_table <- drug_compare[order(drug_compare$p_adj, -abs(drug_compare$delta_high_minus_low)), , drop = FALSE]
priority_table$priority_rank <- seq_len(nrow(priority_table))

fwrite(
  priority_table,
  file.path(table_dir, "67_candidate_drug_priority_table.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 14. Top 药物箱线图
# =========================================================
top_drugs <- head(priority_table$drug, 6)

plot_df <- pred_long %>%
  filter(drug_base %in% top_drugs)

if (nrow(plot_df) > 0) {
  p <- ggplot(plot_df, aes(x = risk_group, y = predicted_response, fill = risk_group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
    facet_wrap(~ drug_base, scales = "free_y") +
    theme_bw(base_size = 12) +
    labs(
      title = "Predicted drug response between high- and low-risk groups",
      x = "Risk group",
      y = "Predicted response"
    )
  
  ggsave(
    filename = file.path(figure_dir, "67_top_candidate_drug_boxplots.pdf"),
    plot = p,
    width = 12,
    height = 8
  )
}

# =========================================================
# Step 15. 写作 digest
# =========================================================
top_hit <- priority_table[1, , drop = FALSE]

writing_digest <- data.frame(
  item = c(
    "n_candidate_drugs_tested",
    "best_ranked_drug",
    "best_ranked_drug_full_id",
    "best_ranked_drug_p_adj",
    "best_ranked_drug_delta_high_minus_low",
    "best_ranked_direction"
  ),
  value = c(
    nrow(priority_table),
    top_hit$drug[1],
    top_hit$drug_full[1],
    top_hit$p_adj[1],
    top_hit$delta_high_minus_low[1],
    top_hit$preferred_in_group[1]
  ),
  stringsAsFactors = FALSE
)

fwrite(
  writing_digest,
  file.path(result_dir, "67_drug_sensitivity_writing_digest.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 16. 阶段摘要
# =========================================================
stage_summary <- data.frame(
  item = c(
    "stage67_completed",
    "sample_alignment_written",
    "raw_prediction_written",
    "match_registry_written",
    "group_comparison_full_written",
    "group_comparison_written",
    "priority_table_written",
    "boxplot_written",
    "writing_digest_written",
    "next_stage"
  ),
  value = c(
    TRUE,
    file.exists(file.path(registry_dir, "67_sample_alignment_registry.tsv")),
    file.exists(file.path(proc_dir, "67_predicted_drug_response_raw.tsv")),
    file.exists(file.path(registry_dir, "67_candidate_drug_match_registry.tsv")),
    file.exists(file.path(result_dir, "67_candidate_drug_group_comparison_full.tsv")),
    file.exists(file.path(result_dir, "67_candidate_drug_group_comparison.tsv")),
    file.exists(file.path(table_dir, "67_candidate_drug_priority_table.tsv")),
    file.exists(file.path(figure_dir, "67_top_candidate_drug_boxplots.pdf")),
    file.exists(file.path(result_dir, "67_drug_sensitivity_writing_digest.tsv")),
    "68_drug_repositioning_and_biological_interpretation"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  stage_summary,
  file.path(registry_dir, "67_stage_summary.tsv"),
  sep = "\t", na = "NA"
)

writeLines(
  c(
    "67 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Run predicted drug sensitivity analysis using oncoPredict and frozen TCGA-PAAD risk groups.",
    "",
    "[Next step]",
    "Interpret candidate drugs and connect them with bulk / single-cell biology."
  ),
  file.path(log_dir, "67_run_predicted_drug_sensitivity_notes.txt")
)

cat("67_run_predicted_drug_sensitivity_analysis.R finished successfully.\n")