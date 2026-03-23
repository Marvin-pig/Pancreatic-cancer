# =========================
# 13_lasso_cox_model_train.R
# =========================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(survival)
  library(glmnet)
  library(survminer)
  library(ggplot2)
  library(broom)
})

# ---------- 路径 ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")
surv_file     <- "04_bulk_analysis/03_survival_model/01_tcga_strict_pdac_survival_ready.tsv"
expr_file     <- "02_processed_data/TCGA_PAAD/tcga_paad_expr_strict_pdac_symbol_counts_filtered.rds"
shortlist_file <- "04_bulk_analysis/03_survival_model/12_prognostic_gene_shortlist.tsv"
out_dir       <- "04_bulk_analysis/03_survival_model"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 工具函数 ----------
tcga16 <- function(x) substr(as.character(x), 1, 16)

zscore_vec <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

# 【修复5】稀疏矩阵系数提取，加 rownames 防御
extract_nonzero <- function(coef_obj, model_label = "") {
  coef_dense <- as.numeric(coef_obj)
  idx <- which(coef_dense != 0)
  if (length(idx) == 0) {
    message("No non-zero coefficients at ", model_label)
    return(data.frame(gene = character(), coefficient = numeric(),
                      stringsAsFactors = FALSE))
  }
  gene_names <- rownames(coef_obj)
  if (is.null(gene_names) || length(gene_names) == 0) {
    stop("coef() returned object with no rownames – check glmnet version.")
  }
  data.frame(
    gene        = gene_names[idx],
    coefficient = coef_dense[idx],
    stringsAsFactors = FALSE
  )
}

# ---------- 读取 ----------
surv_dat  <- fread(surv_file)
expr      <- readRDS(expr_file)
shortlist <- fread(shortlist_file)

if (!is.matrix(expr)) expr <- as.matrix(expr)

stopifnot(all(c("sample_id", "OS_time", "OS_event") %in% colnames(surv_dat)))
stopifnot("gene" %in% colnames(shortlist))

# ---------- 样本对齐 ----------
surv_dat$sample_id_16 <- tcga16(surv_dat$sample_id)
expr_sample_16 <- tcga16(colnames(expr))
common_samples <- intersect(surv_dat$sample_id_16, expr_sample_16)

cat("Matched samples:", length(common_samples), "\n")
if (length(common_samples) < 50) {
  stop("Too few matched samples after ID harmonization.")
}

# 对齐 surv_dat
surv_dat <- surv_dat %>%
  filter(sample_id_16 %in% common_samples) %>%
  distinct(sample_id_16, .keep_all = TRUE)

# 【修复2】match() NA 安全检查
idx_surv <- match(common_samples, surv_dat$sample_id_16)
if (any(is.na(idx_surv))) {
  stop("match() returned NA for surv_dat – ID mismatch after deduplication: ",
       paste(common_samples[is.na(idx_surv)], collapse = ", "))
}
surv_dat <- surv_dat[idx_surv, ]

# 对齐 expr
keep_expr <- !duplicated(expr_sample_16) & expr_sample_16 %in% common_samples
expr <- expr[, keep_expr, drop = FALSE]
expr_sample_16 <- tcga16(colnames(expr))

# 【修复2】match() NA 安全检查
idx_expr <- match(common_samples, expr_sample_16)
if (any(is.na(idx_expr))) {
  stop("match() returned NA for expr – ID mismatch after deduplication: ",
       paste(common_samples[is.na(idx_expr)], collapse = ", "))
}
expr <- expr[, idx_expr, drop = FALSE]

# 最终对齐验证
stopifnot(all(tcga16(colnames(expr)) == surv_dat$sample_id_16))

# ---------- 提取 shortlist 基因 ----------
genes_use <- intersect(unique(shortlist$gene), rownames(expr))
cat("Shortlist genes found in expression matrix:", length(genes_use), "\n")
print(genes_use)
if (length(genes_use) < 5) {
  stop("Too few shortlist genes found in expression matrix.")
}

expr_sub <- expr[genes_use, , drop = FALSE]

# ---------- 转为样本×基因矩阵 ----------
# 【修复1】check.names=FALSE 保留基因原始名称（含连字符等特殊字符）
x_df <- as.data.frame(t(expr_sub), check.names = FALSE)
x_df$sample_id_16 <- tcga16(rownames(x_df))

# ---------- 与生存数据合并 ----------
dat <- surv_dat %>%
  select(sample_id, sample_id_16, OS_time, OS_event) %>%
  left_join(x_df, by = "sample_id_16")

if (any(is.na(dat$OS_time)) || any(is.na(dat$OS_event))) {
  stop("Missing OS_time or OS_event after merge.")
}

expr_cols <- setdiff(colnames(dat), c("sample_id", "sample_id_16", "OS_time", "OS_event"))

# ---------- 标准化表达 ----------
# 【修复3】用 mutate(across(...)) 替代 tibble 上的 [<- 赋值
dat <- dat %>%
  mutate(across(all_of(expr_cols), zscore_vec))

# 保存 LASSO 输入矩阵
fwrite(as.data.frame(dat), file.path(out_dir, "13_lasso_input_matrix.tsv"), sep = "\t")

# ---------- glmnet 输入 ----------
x_mat <- as.matrix(dat[, ..expr_cols])
y_surv <- Surv(dat$OS_time, dat$OS_event)

# ---------- LASSO-Cox ----------
set.seed(20260318)
cvfit <- cv.glmnet(
  x         = x_mat,
  y         = y_surv,
  family    = "cox",
  alpha     = 1,
  nfolds    = 10,
  standardize = FALSE
)

pdf(file.path(out_dir, "13_lasso_cv_plot.pdf"), width = 7, height = 5.5)
plot(cvfit)
dev.off()

# ---------- 提取非零基因 ----------
coef_min <- coef(cvfit, s = "lambda.min")
coef_1se <- coef(cvfit, s = "lambda.1se")

extract_nonzero <- function(coef_obj, model_name) {
  idx <- which(as.numeric(coef_obj) != 0)
  
  if (length(idx) == 0) {
    message("No non-zero coefficients at ", model_name)
    return(data.frame(
      gene = character(0),
      coefficient = numeric(0),
      model = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    gene = rownames(coef_obj)[idx],
    coefficient = as.numeric(coef_obj)[idx],
    model = model_name,
    stringsAsFactors = FALSE
  )
}

nz_min <- extract_nonzero(coef_min, "lambda.min")
nz_1se <- extract_nonzero(coef_1se, "lambda.1se")

nz_all <- bind_rows(nz_min, nz_1se)
fwrite(nz_all, file.path(out_dir, "13_lasso_nonzero_genes.tsv"), sep = "\t")

cat("Non-zero genes at lambda.min:", nrow(nz_min), "\n")
cat("Non-zero genes at lambda.1se:", nrow(nz_1se), "\n")

if (nrow(nz_min) == 0) {
  stop("No non-zero genes at lambda.min. Consider using a broader candidate set.")
}

# ---------- Risk score ----------
# 验证基因名一致性（修复1的下游保障）
missing_in_xmat <- setdiff(nz_min$gene, colnames(x_mat))
if (length(missing_in_xmat) > 0) {
  stop("Gene name mismatch between LASSO coefficients and x_mat columns: ",
       paste(missing_in_xmat, collapse = ", "),
       "\nThis usually means gene names were sanitized. Check check.names=FALSE.")
}

coef_vec <- setNames(nz_min$coefficient, nz_min$gene)
risk_score <- as.numeric(x_mat[, names(coef_vec), drop = FALSE] %*% coef_vec)

# ---------- 【修复4】显式设定因子水平，保证图例顺序正确 ----------
risk_dat <- dat %>%
  select(sample_id, sample_id_16, OS_time, OS_event) %>%
  mutate(
    risk_score = risk_score,
    risk_group = factor(
      ifelse(risk_score >= median(risk_score, na.rm = TRUE), "High", "Low"),
      levels = c("Low", "High")   # Low 为参照，High 为高风险
    )
  )

fwrite(as.data.frame(risk_dat),
       file.path(out_dir, "13_risk_score_training.tsv"), sep = "\t")

# ---------- 训练集 KM ----------
fit_km <- survfit(Surv(OS_time, OS_event) ~ risk_group, data = risk_dat)

pdf(file.path(out_dir, "13_km_training_risk_group.pdf"), width = 7.2, height = 6.2)
p <- ggsurvplot(
  fit_km,
  data              = risk_dat,
  pval              = TRUE,
  risk.table        = TRUE,
  conf.int          = FALSE,
  xlab              = "Time (days)",
  ylab              = "Overall survival probability",
  legend.title      = "Risk group",
  legend.labs       = c("Low", "High"),    # 与 factor levels 顺序一致
  palette           = c("#2E9FDF", "#E7B800"),
  risk.table.height = 0.25,
  ggtheme           = theme_bw(base_size = 12)
)
print(p)
dev.off()

# ---------- 单因素 Cox：risk score ----------
cox_risk <- coxph(Surv(OS_time, OS_event) ~ risk_score, data = risk_dat)
cox_risk_res <- tidy(cox_risk, exponentiate = TRUE, conf.int = TRUE)
fwrite(cox_risk_res,
       file.path(out_dir, "13_univcox_risk_score.tsv"), sep = "\t")
print(cox_risk_res)

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "13_sessionInfo.txt"))

cat("\n13_lasso_cox_model_train.R finished successfully.\n")

