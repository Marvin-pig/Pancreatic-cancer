# =========================
# 16_external_validation_prepare_and_apply_signature.R
# 目的：
# 1. 读取外部队列表达矩阵与临床表
# 2. 检查 signature genes 是否齐全
# 3. 按冻结公式计算外部 risk score（可选冻结 z-score 统计量）
# 4. 使用训练集 median cutoff 分组
# 5. 输出外部验证基础结果
# =========================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(survival)
  library(survminer)
  library(broom)
})

# ---------- 路径配置 ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")
ext_expr_file   <- "03_external_validation/YOUR_EXTERNAL_EXPR_MATRIX.tsv"
ext_cli_file    <- "03_external_validation/YOUR_EXTERNAL_CLINICAL.tsv"
coef_file       <- "04_bulk_analysis/03_survival_model/15_final_signature_coefficients.tsv"
train_risk_file <- "04_bulk_analysis/03_survival_model/13_risk_score_training.tsv"

# [可选] 若训练阶段保存了各基因 z-score 参数，填写路径；否则置 NULL 使用外部队列内部归一化
# 文件须含列：gene, mean_expr, sd_expr
train_zscore_stat_file <- NULL  # 例："04_bulk_analysis/03_survival_model/13_zscore_stats.tsv"

out_dir <- "4_bulk_analysis/04_external_validation/01_signature_validation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 读取数据 ----------
expr_ext   <- fread(ext_expr_file)
cli_ext    <- fread(ext_cli_file)
coef_df    <- fread(coef_file)
train_risk <- fread(train_risk_file)

# ---------- 列名校验 ----------
stopifnot(all(c("sample_id", "OS_time", "OS_event") %in% colnames(cli_ext)))
stopifnot(all(c("gene", "coefficient") %in% colnames(coef_df)))
stopifnot("risk_score" %in% colnames(train_risk))

# ---------- 临床数据值域验证 ----------
# 修复：仅检查列名存在是不够的，需验证值的合法性
invalid_time  <- cli_ext$OS_time <= 0 | is.na(cli_ext$OS_time)
invalid_event <- !cli_ext$OS_event %in% c(0, 1) | is.na(cli_ext$OS_event)

if (any(invalid_time)) {
  n_bad <- sum(invalid_time)
  warning(sprintf("发现 %d 个 OS_time <= 0 或 NA 的样本，已自动剔除。", n_bad))
  cli_ext <- cli_ext[!invalid_time, ]
}
if (any(invalid_event)) {
  n_bad <- sum(invalid_event)
  warning(sprintf("发现 %d 个 OS_event 非 0/1 或 NA 的样本，已自动剔除。", n_bad))
  cli_ext <- cli_ext[!invalid_event, ]
}
stopifnot(nrow(cli_ext) > 0)

# ---------- 表达矩阵格式化 ----------
if (!("gene" %in% colnames(expr_ext))) {
  colnames(expr_ext)[1] <- "gene"
}
expr_mat <- as.data.frame(expr_ext)
rownames(expr_mat) <- expr_mat$gene
expr_mat$gene <- NULL
expr_mat <- as.matrix(expr_mat)

# ---------- Signature genes 匹配 ----------
sig_genes     <- coef_df$gene
genes_found   <- intersect(sig_genes, rownames(expr_mat))
genes_missing <- setdiff(sig_genes, rownames(expr_mat))

cat(sprintf("Signature genes: %d total | %d found | %d missing\n",
            length(sig_genes), length(genes_found), length(genes_missing)))

if (length(genes_missing) > 0) {
  cat("Missing genes:", paste(genes_missing, collapse = ", "), "\n")
  # 修复：明确警告 cutoff 的可比性问题
  pct_missing <- length(genes_missing) / length(sig_genes) * 100
  warning(sprintf(
    "%.1f%% 的 signature genes 在外部队列中缺失。\n",
    pct_missing
  ), "缺失基因将以 0（z-score 均值）代替，risk score 分布可能与训练集有偏移，\n",
  "建议谨慎解读分层结果。")
}

fwrite(data.table(gene = genes_found),
       file.path(out_dir, "16_signature_genes_found.tsv"), sep = "\t")
fwrite(data.table(gene = genes_missing),
       file.path(out_dir, "16_signature_genes_missing.tsv"), sep = "\t")

if (length(genes_found) < 8) {
  stop(sprintf("可匹配 signature genes 仅 %d 个，过少，终止验证。", length(genes_found)))
}

# ---------- 样本对齐 ----------
common_samples <- intersect(colnames(expr_mat), cli_ext$sample_id)
cat("Matched external samples:", length(common_samples), "\n")
if (length(common_samples) < 30) {
  stop("外部队列表达矩阵与临床交集样本过少（< 30）。")
}

expr_mat <- expr_mat[, common_samples, drop = FALSE]
cli_ext  <- cli_ext %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

# 关键对齐校验
stopifnot(identical(colnames(expr_mat), cli_ext$sample_id))

# ---------- 提取 & 归一化 ----------
# 先对已找到的基因建子矩阵（genes × samples）
expr_sub <- expr_mat[genes_found, , drop = FALSE]

if (!is.null(train_zscore_stat_file)) {
  # === 冻结归一化（推荐）：使用训练集的 mean/SD ===
  # 修复：外部验证应使用训练集统计量，避免外部队列内部分布影响
  train_stat <- fread(train_zscore_stat_file)
  stopifnot(all(c("gene", "mean_expr", "sd_expr") %in% colnames(train_stat)))
  
  stat_use <- train_stat %>% filter(gene %in% genes_found)
  row_order <- match(genes_found, stat_use$gene)
  stat_use  <- stat_use[row_order, ]
  stopifnot(all(stat_use$gene == genes_found))
  
  expr_sub_z_genes <- (expr_sub - stat_use$mean_expr) /
    ifelse(stat_use$sd_expr == 0, 1, stat_use$sd_expr)
  # 结果仍为 genes × samples；转置为 samples × genes
  expr_sub_z <- t(expr_sub_z_genes)
  
} else {
  # === 队列内归一化（退而求其次）===
  # 修复：用 scale(t(...)) 替代双重 t(apply(...))，逻辑等价但更易读
  # scale() 默认对列（=基因）做 z-score，输入 t(expr_sub) = samples × genes
  expr_sub_z <- scale(t(expr_sub))  # samples × genes
  # 处理 SD=0 的基因（常数列）
  zero_sd <- apply(expr_sub_z, 2, function(x) all(is.nan(x) | is.na(x)))
  if (any(zero_sd)) {
    warning(sprintf("%d 个基因表达值为常数，z-score 设为 0。",
                    sum(zero_sd)))
    expr_sub_z[, zero_sd] <- 0
  }
}

# ---------- 处理缺失基因（补 0）----------
# 若有缺失基因，在矩阵中补充全零列（z-score 均值），保持与训练时系数维度一致
if (length(genes_missing) > 0) {
  zero_cols <- matrix(0,
                      nrow = nrow(expr_sub_z),
                      ncol = length(genes_missing),
                      dimnames = list(rownames(expr_sub_z), genes_missing))
  expr_sub_z <- cbind(expr_sub_z, zero_cols)
}

# ---------- Risk score 计算 ----------
# 修复：显式按 sig_genes 顺序对齐列，确保与 coef_df 顺序一致，%*% 按位置而非名称匹配
coef_ordered <- coef_df$coefficient
names(coef_ordered) <- coef_df$gene
# 保证 expr_sub_z 的列顺序与 coef_ordered 完全一致
stopifnot(all(sig_genes %in% colnames(expr_sub_z)))
expr_aligned  <- expr_sub_z[, sig_genes, drop = FALSE]
coef_vec_aln  <- coef_ordered[sig_genes]
stopifnot(identical(colnames(expr_aligned), names(coef_vec_aln)))  # 最终对齐校验

risk_score <- as.numeric(expr_aligned %*% coef_vec_aln)

# 导出实际使用的基因与系数（包含缺失基因，标记为补零）
genes_used_df <- data.frame(
  gene        = sig_genes,
  coefficient = coef_vec_aln,
  imputed     = sig_genes %in% genes_missing
)
fwrite(genes_used_df,
       file.path(out_dir, "16_external_genes_used_coefficients.tsv"), sep = "\t")

# ---------- 固定 cutoff 分组 ----------
train_cutoff <- median(train_risk$risk_score, na.rm = TRUE)
cat(sprintf("Training median cutoff: %.4f\n", train_cutoff))

# 修复：显式设置 factor levels，确保 KM 图曲线、风险表顺序可预期
risk_ext <- cli_ext %>%
  mutate(
    risk_score = risk_score,
    risk_group = factor(
      ifelse(risk_score >= train_cutoff, "High", "Low"),
      levels = c("Low", "High")   # Low 在前：KM 图中 Low 曲线先绘制，视觉上更清晰
    )
  )

cat(table(risk_ext$risk_group), "\n")

fwrite(risk_ext,
       file.path(out_dir, "16_external_risk_score.tsv"), sep = "\t")

# ---------- KM 曲线 ----------
fit_km <- survfit(Surv(OS_time, OS_event) ~ risk_group, data = risk_ext)

pdf(file.path(out_dir, "16_external_km_risk_group.pdf"), width = 7.2, height = 6.5)
p <- ggsurvplot(
  fit_km,
  data          = risk_ext,
  pval          = TRUE,
  pval.method   = TRUE,
  risk.table    = TRUE,
  conf.int      = FALSE,
  palette       = c("#2E86AB", "#E84855"),  # 修复：固定颜色；Low=蓝，High=红
  xlab          = "Time (days)",
  ylab          = "Overall survival probability",
  legend.title  = "Risk group",
  legend.labs   = c("Low", "High"),         # 与 factor levels 顺序一致
  risk.table.height = 0.25,
  ggtheme       = theme_classic()
)
print(p)
dev.off()

# ---------- 单因素 Cox ----------
cox_risk <- coxph(Surv(OS_time, OS_event) ~ risk_score, data = risk_ext)
cox_res  <- broom::tidy(cox_risk, exponentiate = TRUE, conf.int = TRUE)

fwrite(cox_res,
       file.path(out_dir, "16_external_univcox_risk_score.tsv"), sep = "\t")
print(cox_res)

cat("\n[DONE] 16_external_validation_prepare_and_apply_signature.R finished.\n")
cat(sprintf("Output directory: %s\n", out_dir))
