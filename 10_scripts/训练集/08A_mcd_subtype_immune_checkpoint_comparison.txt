# =========================================================
# 08A_mcd_subtype_immune_checkpoint_comparison.R
# 作用：
# 1. 读取表达矩阵与 MCD subtype phenotype
# 2. 提取 immune checkpoint genes
# 3. 比较 subtype 间表达差异（KW + 事后 Wilcoxon）
# 4. 输出结果文件与箱线图
# =========================================================

# ---- 0. 路径配置 ----
project_root <- Sys.getenv("PROJECT_ROOT", unset = "/Users/wmz_mac/Desktop/胰腺癌")
proc_dir     <- file.path(project_root, "02_processed_data", "TCGA_PAAD")
clust_dir    <- file.path(project_root, "04_bulk_analysis", "02_clustering")
immune_dir   <- file.path(project_root, "04_bulk_analysis", "05_immune")
dir.create(immune_dir, recursive = TRUE, showWarnings = FALSE)

expr_file  <- file.path(proc_dir,   "tcga_paad_expr_all_tp_symbol_counts_filtered.rds")
pheno_file <- file.path(clust_dir,  "tcga_paad_mcd_subtype_phenotype.tsv")

# ---- 1. 加载包 ----
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- 2. 读取数据 ----
expr  <- readRDS(expr_file)
pheno <- fread(pheno_file, data.table = FALSE)

# Bug修复①：确保 expr 为 matrix，data.frame 的行索引行为不同
if (!is.matrix(expr)) expr <- as.matrix(expr)

stopifnot(!is.null(rownames(expr)))
stopifnot(!anyDuplicated(rownames(expr)))
stopifnot(!anyDuplicated(colnames(expr)))

# ---- 3. 识别样本列 ----
sample_col <- intersect(c("sample_id", "sample_barcode"), colnames(pheno))[1]
if (is.na(sample_col)) {
  stop("phenotype 中未找到样本列（需含 sample_id 或 sample_barcode）。")
}
if (!"subtype_label" %in% colnames(pheno)) {
  stop("phenotype 中缺少 subtype_label 列。")
}

# ---- 4. 对齐样本顺序（含缺失样本处理） ----
# Bug修复②：match 可能返回 NA，需显式处理而非让 stopifnot 失败
matched_idx <- match(colnames(expr), pheno[[sample_col]])
n_missing   <- sum(is.na(matched_idx))

if (n_missing > 0) {
  warning(sprintf(
    "%d 个表达矩阵样本在 phenotype 中未找到，已排除：%s",
    n_missing,
    paste(colnames(expr)[is.na(matched_idx)], collapse = ", ")
  ))
  keep        <- !is.na(matched_idx)
  expr        <- expr[, keep, drop = FALSE]
  matched_idx <- matched_idx[keep]
}

pheno <- pheno[matched_idx, , drop = FALSE]
stopifnot(identical(colnames(expr), pheno[[sample_col]]))

# ---- 5. 提取 immune checkpoint genes ----
icp_genes <- c("CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3",
               "TIGIT", "HAVCR2", "IDO1", "SIGLEC15")
icp_genes <- intersect(icp_genes, rownames(expr))

if (length(icp_genes) == 0) stop("表达矩阵中未找到任何 immune checkpoint genes。")
message(sprintf("检测到 %d 个 ICP genes: %s",
                length(icp_genes), paste(icp_genes, collapse = ", ")))

expr_icp <- expr[icp_genes, , drop = FALSE]

# ---- 6. 构建长表 ----
long_df <- do.call(rbind, lapply(icp_genes, function(g) {
  data.frame(
    sample_id   = colnames(expr_icp),
    gene_symbol = g,
    expression  = as.numeric(expr_icp[g, ]),
    stringsAsFactors = FALSE
  )
}))

long_df <- merge(
  long_df,
  pheno[, c(sample_col, "subtype_label"), drop = FALSE],
  by.x  = "sample_id",
  by.y  = sample_col,
  all.x = TRUE,
  sort  = FALSE
)

# Bug修复③：过滤 merge 后 subtype_label 为 NA 的行
n_na_sub <- sum(is.na(long_df$subtype_label))
if (n_na_sub > 0) {
  warning(sprintf("过滤 %d 行无 subtype_label 的记录。", n_na_sub))
  long_df <- long_df[!is.na(long_df$subtype_label), ]
}

# 检查 subtype 数量（KW 需要至少 2 组）
n_groups <- length(unique(long_df$subtype_label))
if (n_groups < 2) stop(sprintf("subtype_label 仅有 %d 个水平，无法做组间比较。", n_groups))

# 添加 log2 变换列（用于可视化）
long_df$log2_expr <- log2(long_df$expression + 1)

fwrite(
  long_df,
  file.path(immune_dir, "tcga_paad_mcd_subtype_icp_expression.tsv"),
  sep = "\t", quote = FALSE
)

# ---- 7. Kruskal-Wallis 检验 ----
# 修复：用 icp_genes 保证顺序稳定
kw_res <- do.call(rbind, lapply(icp_genes, function(g) {
  tmp <- long_df[long_df$gene_symbol == g, , drop = FALSE]
  fit <- kruskal.test(expression ~ subtype_label, data = tmp)
  data.frame(
    gene_symbol = g,
    statistic   = as.numeric(fit$statistic),
    df          = as.numeric(fit$parameter),
    p_value     = fit$p.value,
    stringsAsFactors = FALSE
  )
}))

kw_res$p_adj_BH <- p.adjust(kw_res$p_value, method = "BH")
kw_res <- kw_res[order(kw_res$p_value), ]

fwrite(
  kw_res,
  file.path(immune_dir, "tcga_paad_mcd_subtype_icp_comparison.tsv"),
  sep = "\t", quote = FALSE
)

# ---- 8. Subtype 均值表达 ----
mean_res <- do.call(rbind, lapply(icp_genes, function(g) {
  tmp <- long_df[long_df$gene_symbol == g, , drop = FALSE]
  agg <- aggregate(expression ~ subtype_label, data = tmp, FUN = mean)
  data.frame(
    gene_symbol     = g,
    subtype_label   = agg$subtype_label,
    mean_expression = agg$expression,
    stringsAsFactors = FALSE
  )
}))
mean_res <- mean_res[order(mean_res$gene_symbol, mean_res$subtype_label), ]
rownames(mean_res) <- NULL

fwrite(
  mean_res,
  file.path(immune_dir, "tcga_paad_mcd_subtype_icp_cluster_means.tsv"),
  sep = "\t", quote = FALSE
)

# ---- 9. 两两 Wilcoxon 事后检验（仅对 KW BH<0.05 的基因） ----
sig_genes <- kw_res$gene_symbol[kw_res$p_adj_BH < 0.05]

if (length(sig_genes) > 0) {
  pw_res <- do.call(rbind, lapply(sig_genes, function(g) {
    tmp    <- long_df[long_df$gene_symbol == g, ]
    groups <- sort(unique(tmp$subtype_label))
    pairs  <- combn(groups, 2, simplify = FALSE)
    do.call(rbind, lapply(pairs, function(pr) {
      x  <- tmp$expression[tmp$subtype_label == pr[1]]
      y  <- tmp$expression[tmp$subtype_label == pr[2]]
      wt <- wilcox.test(x, y, exact = FALSE)
      data.frame(
        gene_symbol = g,
        group1      = pr[1],
        group2      = pr[2],
        W_statistic = as.numeric(wt$statistic),
        p_value     = wt$p.value,
        stringsAsFactors = FALSE
      )
    }))
  }))
  pw_res$p_adj_BH <- p.adjust(pw_res$p_value, method = "BH")
  fwrite(
    pw_res,
    file.path(immune_dir, "tcga_paad_mcd_subtype_icp_posthoc_wilcoxon.tsv"),
    sep = "\t", quote = FALSE
  )
  message(sprintf("事后检验完成，BH<0.05 的显著基因：%s",
                  paste(sig_genes, collapse = ", ")))
} else {
  message("无 KW BH<0.05 的显著基因，跳过事后检验。")
}

# ---- 10. 可视化：分面箱线图 ----
# 将 KW 结果标签拼入 long_df 用于 facet 标题
kw_label_df <- data.frame(
  gene_symbol = kw_res$gene_symbol,
  facet_label = sprintf("%s\n(KW p_adj=%.3g%s)",
                        kw_res$gene_symbol,
                        kw_res$p_adj_BH,
                        ifelse(kw_res$p_adj_BH < 0.05, " *", "")),
  stringsAsFactors = FALSE
)

plot_df <- merge(long_df, kw_label_df, by = "gene_symbol")

# 固定分面顺序（按 p 值排序）
ordered_labels <- kw_label_df$facet_label[match(kw_res$gene_symbol, kw_label_df$gene_symbol)]
plot_df$facet_label <- factor(plot_df$facet_label, levels = ordered_labels)

ncol_plot <- 3
nrow_plot  <- ceiling(length(icp_genes) / ncol_plot)

p <- ggplot(plot_df, aes(x = subtype_label, y = log2_expr, fill = subtype_label)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.35, color = "grey30") +
  facet_wrap(~ facet_label, scales = "free_y", ncol = ncol_plot) +
  labs(
    title = "Immune Checkpoint Gene Expression by MCD Subtype (TCGA-PAAD)",
    x     = "MCD Subtype",
    y     = expression(log[2](Count + 1)),
    fill  = "Subtype"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1, size = 9),
    strip.text      = element_text(size = 8),
    legend.position = "bottom"
  )

fig_h <- 4 * nrow_plot

ggsave(
  file.path(immune_dir, "tcga_paad_mcd_subtype_icp_boxplot.pdf"),
  plot = p, width = 12, height = fig_h, useDingbats = FALSE
)
ggsave(
  file.path(immune_dir, "tcga_paad_mcd_subtype_icp_boxplot.png"),
  plot = p, width = 12, height = fig_h, dpi = 150
)
message(sprintf("箱线图已保存至 %s", immune_dir))

# ---- 11. 汇总打印 ----
print(kw_res)
message("Immune checkpoint comparison finished.")
