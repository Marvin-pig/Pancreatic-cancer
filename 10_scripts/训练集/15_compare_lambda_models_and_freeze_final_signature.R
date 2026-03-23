# =========================
# 15_compare_lambda_models_and_freeze_final_signature.R
# 目的：
# 1. 重新基于 13_lasso_input_matrix.tsv 比较 lambda.min 与 lambda.1se
# 2. 输出两个版本的非零基因与训练集表现
# 3. 冻结最终模型版本
# 4. 输出供外部验证直接调用的最终系数表
# =========================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(survival)
  library(glmnet)
  library(survminer)
  library(broom)
})

# ── 阈值常量（集中定义，便于维护）──────────────────────────────────────────
COEF_ZERO_THRESH  <- 1e-10   # LASSO 系数判零阈值（避免浮点精度误差）
P_THRESHOLD       <- 0.01    # lambda.1se 显著性门槛
HR_RATIO_MIN      <- 0.85    # lambda.1se HR 不得低于 lambda.min HR 的此比例
# ─────────────────────────────────────────────────────────────────────────────
setwd("/Users/wmz_mac/Desktop/胰腺癌")
infile  <- "04_bulk_analysis/03_survival_model/13_lasso_input_matrix.tsv"
out_dir <- "04_bulk_analysis/03_survival_model"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── 读取数据 ─────────────────────────────────────────────────────────────────
dat <- fread(infile)

required_cols <- c("sample_id", "OS_time", "OS_event")
missing_cols  <- setdiff(required_cols, colnames(dat))
if (length(missing_cols) > 0) {
  stop("输入文件缺少必要列：", paste(missing_cols, collapse = ", "))
}

# 去除非表达量列（sample_id_16 列若不存在也不报错）
meta_cols  <- intersect(c("sample_id", "sample_id_16", "OS_time", "OS_event"),
                        colnames(dat))
expr_cols  <- setdiff(colnames(dat), meta_cols)
if (length(expr_cols) == 0) stop("未找到任何基因表达列，请检查输入文件。")

x_mat  <- as.matrix(dat[, ..expr_cols])
y_surv <- Surv(dat$OS_time, dat$OS_event)

# ── LASSO Cox 交叉验证 ───────────────────────────────────────────────────────
set.seed(20260318)
cvfit <- cv.glmnet(
  x         = x_mat,
  y         = y_surv,
  family    = "cox",
  alpha     = 1,
  nfolds    = 10,
  standardize = FALSE
)

# ── 提取非零系数 ─────────────────────────────────────────────────────────────
# 修复：使用阈值而非 != 0，避免浮点精度导致的假非零
extract_nonzero <- function(cvfit, s_value, thresh = COEF_ZERO_THRESH) {
  coef_obj <- coef(cvfit, s = s_value)               # dgCMatrix (p × 1)
  coef_vec <- drop(as.matrix(coef_obj))              # 转为命名数值向量
  idx      <- which(abs(coef_vec) > thresh)
  if (length(idx) == 0) return(data.frame())
  data.frame(
    gene        = names(coef_vec)[idx],
    coefficient = coef_vec[idx],
    model       = s_value,
    stringsAsFactors = FALSE
  )
}

# ── 计算风险评分 ─────────────────────────────────────────────────────────────
# 修复：增加基因名存在性校验，明确矩阵乘法语义
score_samples <- function(x_mat, coef_df) {
  missing_genes <- setdiff(coef_df$gene, colnames(x_mat))
  if (length(missing_genes) > 0) {
    stop("x_mat 中缺少以下基因列：", paste(missing_genes, collapse = ", "))
  }
  coef_vec <- setNames(coef_df$coefficient, coef_df$gene)
  # drop() 将 n×1 矩阵结果显式降为向量
  drop(x_mat[, coef_df$gene, drop = FALSE] %*% coef_vec)
}

# ── 评估模型（训练集 CoxPH）─────────────────────────────────────────────────
# 修复：
#   1. 移除从未使用的 risk_group 列
#   2. 将 n_genes 改为参数传入，逻辑集中在调用处
eval_model <- function(dat, score, label, n_genes) {
  tmp <- dat %>%
    select(OS_time, OS_event) %>%
    mutate(risk_score = score)
  
  fit <- coxph(Surv(OS_time, OS_event) ~ risk_score, data = tmp)
  res <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  res$model  <- label
  res$n_genes <- n_genes
  res
}

# ── 提取两个模型的系数 ───────────────────────────────────────────────────────
nz_min <- extract_nonzero(cvfit, "lambda.min")
nz_1se <- extract_nonzero(cvfit, "lambda.1se")

fwrite(
  bind_rows(nz_min, nz_1se),
  file.path(out_dir, "15_lasso_nonzero_genes_compare.tsv"),
  sep = "\t"
)

if (nrow(nz_min) == 0) {
  stop("lambda.min 下无非零基因，LASSO 拟合可能异常，请检查输入数据。")
}
if (nrow(nz_1se) == 0) {
  warning("lambda.1se 下无非零基因，将直接使用 lambda.min 模型。")
}

# ── 计算训练集评分与 CoxPH 性能 ─────────────────────────────────────────────
score_min <- score_samples(x_mat, nz_min)
model_perf <- eval_model(dat, score_min, "lambda.min", nrow(nz_min))

if (nrow(nz_1se) > 0) {
  score_1se <- score_samples(x_mat, nz_1se)
  res_1se   <- eval_model(dat, score_1se, "lambda.1se", nrow(nz_1se))
  model_perf <- bind_rows(model_perf, res_1se)
}

fwrite(model_perf,
       file.path(out_dir, "15_lambda_model_performance.tsv"),
       sep = "\t")
print(model_perf)

# ── 选择最终模型 ─────────────────────────────────────────────────────────────
# 规则：
#   若 lambda.1se 模型存在，且满足：
#     (1) p < P_THRESHOLD（显著）
#     (2) 基因数更少（更稀疏）
#     (3) HR 不低于 lambda.min HR 的 HR_RATIO_MIN 倍（效果相近）
#   则优先选 lambda.1se；否则保留 lambda.min
#
# 修复：
#   1. 改用 "lambda.1se" %in% model_perf$model 替代脆弱的 nrow == 2 判断
#   2. 对 hr_min、hr_1se 均加 NA 保护

final_model <- "lambda.min"

if ("lambda.1se" %in% model_perf$model) {
  hr_min <- model_perf$estimate[model_perf$model == "lambda.min"]
  hr_1se <- model_perf$estimate[model_perf$model == "lambda.1se"]
  p_1se  <- model_perf$p.value[model_perf$model == "lambda.1se"]
  n_min  <- model_perf$n_genes[model_perf$model == "lambda.min"]
  n_1se  <- model_perf$n_genes[model_perf$model == "lambda.1se"]
  
  select_1se <- (
    !is.na(p_1se)  && p_1se  < P_THRESHOLD   &&
      !is.na(hr_min) && !is.na(hr_1se)          &&
      n_1se < n_min                             &&
      hr_1se >= HR_RATIO_MIN * hr_min
  )
  
  if (isTRUE(select_1se)) {
    final_model <- "lambda.1se"
    cat(sprintf(
      "\n[选择 lambda.1se] p=%.4f < %.2f，基因数 %d < %d，HR ratio=%.3f >= %.2f\n",
      p_1se, P_THRESHOLD, n_1se, n_min, hr_1se / hr_min, HR_RATIO_MIN
    ))
  } else {
    cat(sprintf(
      "\n[保留 lambda.min] lambda.1se 未满足选择条件（p=%.4f, genes: %d vs %d, HR ratio=%.3f）\n",
      ifelse(is.na(p_1se), NaN, p_1se),
      ifelse(is.na(n_1se), -1L, n_1se), n_min,
      ifelse(is.na(hr_1se) || is.na(hr_min), NaN, hr_1se / hr_min)
    ))
  }
}

# ── 输出最终签名系数 ─────────────────────────────────────────────────────────
final_coef <- bind_rows(nz_min, nz_1se) %>%
  filter(model == final_model) %>%
  arrange(desc(abs(coefficient)))

fwrite(final_coef,
       file.path(out_dir, "15_final_signature_coefficients.tsv"),
       sep = "\t")

# 记录模型元信息
formula_text <- paste0(
  sprintf("%.6f", final_coef$coefficient), " * ", final_coef$gene,
  collapse = " + "
)
final_info <- data.frame(
  final_model      = final_model,
  n_genes          = nrow(final_coef),
  cutoff_rule      = "training median",
  lambda_value     = cvfit[[final_model]],   # 记录实际 lambda 数值，便于复现
  risk_score_formula = paste0("risk_score = ", formula_text),
  stringsAsFactors = FALSE
)
fwrite(final_info,
       file.path(out_dir, "15_final_signature_model.tsv"),
       sep = "\t")

cat("\nFinal model selected:", final_model, "\n")
cat("Number of genes:", nrow(final_coef), "\n")
print(final_coef)
