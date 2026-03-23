# =========================
# 19a_build_train_gene_norm_params.R
# 目的：
# 1. 基于训练集 LASSO 输入矩阵
# 2. 提取最终 12 个 signature genes 的训练集 mean / sd
# 3. 生成外部验证所需的标准化参数文件
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ---------- 路径 ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")
lasso_input_file <- "04_bulk_analysis/03_survival_model/13_lasso_input_matrix.tsv"
coef_file        <- "04_bulk_analysis/03_survival_model/15_final_signature_coefficients.tsv"
out_file         <- "04_bulk_analysis/03_survival_model/15_train_gene_norm_params.tsv"

# ---------- 读取 ----------
dat <- fread(lasso_input_file)
coef_df <- fread(coef_file)

stopifnot("gene" %in% colnames(coef_df))

sig_genes <- unique(coef_df$gene)

# 训练集 LASSO 输入矩阵里，gene 列应该已经是 z-score 前的表达还是 z-score 后？
# 你当前 13_lasso_input_matrix.tsv 里保存的是用于 LASSO 的标准化矩阵
# 为了让外部验证严格对齐当前模型实现，我们就按这个文件中的列统计 mean/sd。
#
# 如果某基因在该文件中不存在，直接报错
missing_genes <- setdiff(sig_genes, colnames(dat))
if (length(missing_genes) > 0) {
  stop(
    "以下 signature genes 不在 13_lasso_input_matrix.tsv 中：",
    paste(missing_genes, collapse = ", ")
  )
}

norm_df <- data.frame(
  gene = sig_genes,
  mean_expr = sapply(sig_genes, function(g) mean(dat[[g]], na.rm = TRUE)),
  sd_expr   = sapply(sig_genes, function(g) sd(dat[[g]], na.rm = TRUE)),
  stringsAsFactors = FALSE
)

# 防止 sd=0
norm_df$sd_expr[is.na(norm_df$sd_expr) | norm_df$sd_expr == 0] <- 1

fwrite(norm_df, out_file, sep = "\t")

cat("Training normalization parameter file generated:\n")
cat(out_file, "\n\n")
print(norm_df)