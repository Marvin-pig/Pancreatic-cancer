# =========================================================
# 04_prepare_tcga_paad_expression_for_analysis.R
# 作用：
# 1. 读取原始 counts matrix
# 2. 读取 cohort_master / strict_pdac
# 3. 裁剪表达矩阵到最终队列
# 4. 检查样本一一对应
# 5. 导出 analysis-ready counts matrix
# =========================================================

# ---------------------------------------------------------
# 0. 环境配置
# ---------------------------------------------------------
project_root <- Sys.getenv("PROJECT_ROOT",
                            unset = "/Users/wmz_mac/Desktop/胰腺癌")

library(data.table)

expr_dir <- file.path(project_root, "00_raw_data",      "TCGA_PAAD", "expression")
proc_dir <- file.path(project_root, "02_processed_data", "TCGA_PAAD")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

expr_file       <- file.path(expr_dir, "tcga_paad_counts_matrix_selected_assay.rds")
gene_anno_f     <- file.path(expr_dir, "tcga_paad_gene_annotation.tsv")
cohort_all_f    <- file.path(proc_dir, "tcga_paad_cohort_master.tsv")
cohort_strict_f <- file.path(proc_dir, "tcga_paad_cohort_master_strict_pdac.tsv")

# ---------------------------------------------------------
# 1. 文件存在性检查
# ---------------------------------------------------------
for (f in c(expr_file, gene_anno_f, cohort_all_f, cohort_strict_f)) {
  if (!file.exists(f)) stop("找不到文件: ", f)
}

# ---------------------------------------------------------
# 2. 读取数据
# ---------------------------------------------------------
expr_mat   <- as.matrix(readRDS(expr_file))
gene_anno  <- fread(gene_anno_f,     data.table = FALSE)
cohort_all <- fread(cohort_all_f,    data.table = FALSE)
cohort_str <- fread(cohort_strict_f, data.table = FALSE)

cat("原始表达矩阵维度（基因 x 样本）:", dim(expr_mat), "\n")

# ---------------------------------------------------------
# 3. 验证表达矩阵列名
# ---------------------------------------------------------
expr_barcodes <- colnames(expr_mat)
if (is.null(expr_barcodes) || length(expr_barcodes) == 0) {
  stop("表达矩阵没有列名，无法匹配样本。")
}

# ---------------------------------------------------------
# 4. 验证队列文件无重复 barcode
# ---------------------------------------------------------
if (anyDuplicated(cohort_all$sample_barcode)) {
  stop("cohort_all 中存在重复的 sample_barcode: ",
       paste(cohort_all$sample_barcode[duplicated(cohort_all$sample_barcode)],
             collapse = ", "))
}
if (anyDuplicated(cohort_str$sample_barcode)) {
  stop("cohort_str 中存在重复的 sample_barcode: ",
       paste(cohort_str$sample_barcode[duplicated(cohort_str$sample_barcode)],
             collapse = ", "))
}

# ---------------------------------------------------------
# 5. 验证基因注释与表达矩阵行名一致
# ---------------------------------------------------------
if (!is.null(gene_anno$gene_id) && !is.null(rownames(expr_mat))) {
  missing_genes <- setdiff(rownames(expr_mat), gene_anno$gene_id)
  if (length(missing_genes) > 0) {
    warning(sprintf("表达矩阵中有 %d 个基因在注释文件中缺失。", length(missing_genes)))
  }
}

# ---------------------------------------------------------
# 6. all-TP 队列匹配
# ---------------------------------------------------------
missing_all <- setdiff(cohort_all$sample_barcode, expr_barcodes)
if (length(missing_all) > 0) {
  warning(sprintf(
    "all-TP cohort：%d 个样本在表达矩阵中缺失，已跳过:\n  %s",
    length(missing_all), paste(missing_all, collapse = "\n  ")
  ))
}

keep_all  <- intersect(cohort_all$sample_barcode, expr_barcodes)
expr_all  <- expr_mat[, keep_all, drop = FALSE]
cohort_all2 <- cohort_all[match(colnames(expr_all), cohort_all$sample_barcode), ]
stopifnot(identical(colnames(expr_all), cohort_all2$sample_barcode))

# ---------------------------------------------------------
# 7. strict-PDAC 队列匹配
# ---------------------------------------------------------
missing_str <- setdiff(cohort_str$sample_barcode, expr_barcodes)
if (length(missing_str) > 0) {
  warning(sprintf(
    "strict-PDAC cohort：%d 个样本在表达矩阵中缺失，已跳过:\n  %s",
    length(missing_str), paste(missing_str, collapse = "\n  ")
  ))
}

keep_str  <- intersect(cohort_str$sample_barcode, expr_barcodes)
expr_str  <- expr_mat[, keep_str, drop = FALSE]
cohort_str2 <- cohort_str[match(colnames(expr_str), cohort_str$sample_barcode), ]
stopifnot(identical(colnames(expr_str), cohort_str2$sample_barcode))

# ---------------------------------------------------------
# 8. 导出表达矩阵与队列文件
# ---------------------------------------------------------
saveRDS(expr_all, file = file.path(proc_dir, "tcga_paad_expr_all_tp_counts.rds"))
saveRDS(expr_str, file = file.path(proc_dir, "tcga_paad_expr_strict_pdac_counts.rds"))

write.table(
  cohort_all2,
  file  = file.path(proc_dir, "tcga_paad_cohort_master_ordered.tsv"),
  sep   = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  cohort_str2,
  file  = file.path(proc_dir, "tcga_paad_cohort_master_strict_pdac_ordered.tsv"),
  sep   = "\t", quote = FALSE, row.names = FALSE
)

# ---------------------------------------------------------
# 9. 导出审计日志
# ---------------------------------------------------------
audit_df <- data.frame(
  metric = c(
    "raw_genes",
    "raw_samples",
    "cohort_all_tp_expected",
    "cohort_strict_pdac_expected",
    "matched_all_tp_samples",
    "matched_strict_pdac_samples",
    "dropped_all_tp_samples",
    "dropped_strict_pdac_samples"
  ),
  value = c(
    nrow(expr_mat),
    ncol(expr_mat),
    nrow(cohort_all),
    nrow(cohort_str),
    ncol(expr_all),
    ncol(expr_str),
    nrow(cohort_all) - ncol(expr_all),
    nrow(cohort_str)  - ncol(expr_str)
  )
)
write.table(
  audit_df,
  file  = file.path(proc_dir, "tcga_paad_expression_prepare_audit.tsv"),
  sep   = "\t", quote = FALSE, row.names = FALSE
)

# ---------------------------------------------------------
# 10. 汇总输出
# ---------------------------------------------------------
cat("\n========== 表达矩阵准备完成 ==========\n")
cat(sprintf("原始基因数:             %d\n", nrow(expr_mat)))
cat(sprintf("原始样本数:             %d\n", ncol(expr_mat)))
cat(sprintf("all-TP   队列期望/实际: %d / %d（丢失 %d）\n",
            nrow(cohort_all), ncol(expr_all), nrow(cohort_all) - ncol(expr_all)))
cat(sprintf("strict-PDAC 期望/实际: %d / %d（丢失 %d）\n",
            nrow(cohort_str), ncol(expr_str), nrow(cohort_str) - ncol(expr_str)))
cat("=======================================\n\n")

sessionInfo()
