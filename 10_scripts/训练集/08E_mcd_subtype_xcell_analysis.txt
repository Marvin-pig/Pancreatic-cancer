# =========================================================
# 08E_mcd_subtype_xcell_analysis.R
# 作用：
# 1. 读取 all-TP 表达矩阵与 MCD phenotype
# 2. logCPM 标准化
# 3. 运行 xCell
# 4. 比较 MCD subtype 间 xCell score 差异
# 5. 输出结果表、热图和候选箱线图
# =========================================================
project_root <- Sys.getenv("PROJECT_ROOT", unset = "/Users/wmz_mac/Desktop/胰腺癌")
proc_dir     <- file.path(project_root, "02_processed_data", "TCGA_PAAD")
clust_dir    <- file.path(project_root, "04_bulk_analysis", "02_clustering")
immune_dir   <- file.path(project_root, "04_bulk_analysis", "05_immune")
dir.create(immune_dir, recursive = TRUE, showWarnings = FALSE)
expr_file  <- file.path(proc_dir,  "tcga_paad_expr_all_tp_symbol_counts_filtered.rds")
pheno_file <- file.path(clust_dir, "tcga_paad_mcd_subtype_phenotype.tsv")
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(pheatmap)
  library(rstatix)
  library(xCell)
})
# -------------------------------
# 1. 读取数据
# -------------------------------
expr  <- readRDS(expr_file)
pheno <- fread(pheno_file, data.table = FALSE)
if (!is.matrix(expr)) expr <- as.matrix(expr)
stopifnot(!is.null(rownames(expr)))
stopifnot(!anyDuplicated(rownames(expr)))
stopifnot(!anyDuplicated(colnames(expr)))
sample_col <- intersect(c("sample_id", "sample_barcode"), colnames(pheno))[1]
if (is.na(sample_col)) stop("phenotype 中未找到 sample_id / sample_barcode")
if (!"subtype_label" %in% colnames(pheno)) stop("phenotype 中缺少 subtype_label")
# -------------------------------
# 2. 样本对齐
# -------------------------------
matched_idx <- match(colnames(expr), pheno[[sample_col]])
if (any(is.na(matched_idx))) {
  keep        <- !is.na(matched_idx)
  expr        <- expr[, keep, drop = FALSE]
  matched_idx <- matched_idx[keep]
}
pheno <- pheno[matched_idx, , drop = FALSE]
stopifnot(identical(colnames(expr), pheno[[sample_col]]))
# -------------------------------
# 3. 表达标准化
# 说明：
# [修复 Bug 1] xCellAnalysis 在 rnaseq=TRUE 时内部执行 log2(expr+1)。
# 若传入已 log 变换的 logCPM，将造成双重 log 变换，严重扭曲评分。
# 正确做法：传入线性尺度 CPM，由 xCell 自行 log 变换（rnaseq=TRUE）。
# logCPM 仍保留，供统计分析和可视化使用。
# -------------------------------
# 线性 CPM，供 xCell 使用
cpm_mat <- edgeR::cpm(expr, log = FALSE)

keep_gene <- !is.na(rownames(cpm_mat)) & rownames(cpm_mat) != ""
cpm_mat   <- cpm_mat[keep_gene, , drop = FALSE]

# logCPM，供后续可视化/统计（若需要）
logcpm <- log2(cpm_mat + 1)

# 如果仍有重复 symbol，保留第一条
if (anyDuplicated(rownames(cpm_mat))) {
  keep_rows <- !duplicated(rownames(cpm_mat))
  cpm_mat   <- cpm_mat[keep_rows, , drop = FALSE]
  logcpm    <- logcpm[keep_rows, , drop = FALSE]
}
# -------------------------------
# 4. 运行 xCell
# 返回通常是 cell type x sample
# rnaseq = TRUE：xCell 内部执行 log2(cpm+1)，与原设计一致
# -------------------------------
xcell_mat <- xCell::xCellAnalysis(
  expr   = cpm_mat,
  rnaseq = TRUE
)
# 转为 sample x feature
# 注：as.data.frame 默认 check.names=TRUE，
# 会将 "CD8+ T-cells" 等转为 "CD8..T.cells"（以点替换非法字符）
xcell_df <- as.data.frame(t(xcell_mat), stringsAsFactors = FALSE)
xcell_df$sample_id <- rownames(xcell_df)
rownames(xcell_df) <- NULL
# 合并 phenotype
xcell_df <- merge(
  xcell_df,
  pheno[, c(sample_col, "subtype_label"), drop = FALSE],
  by.x   = "sample_id",
  by.y   = sample_col,
  all.x  = TRUE,
  sort   = FALSE
)
xcell_df$subtype_label <- factor(
  xcell_df$subtype_label,
  levels = c("MCD-low", "MCD-intermediate", "MCD-high")
)
# 保存完整分数表
fwrite(
  xcell_df,
  file.path(immune_dir, "tcga_paad_mcd_subtype_xcell_scores.tsv"),
  sep   = "\t",
  quote = FALSE
)
# -------------------------------
# 5. 统计比较：全局 Kruskal-Wallis
# -------------------------------
features <- setdiff(colnames(xcell_df), c("sample_id", "subtype_label"))
kw_res <- do.call(rbind, lapply(features, function(f) {
  tmp <- xcell_df[, c("subtype_label", f), drop = FALSE]
  colnames(tmp)[2] <- "value"
  tmp <- tmp[!is.na(tmp$value), , drop = FALSE]
  if (nrow(tmp) == 0 || length(unique(tmp$value)) < 2) {
    return(data.frame(
      feature   = f,
      statistic = NA_real_,
      df        = NA_real_,
      p_value   = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  fit <- kruskal.test(value ~ subtype_label, data = tmp)
  data.frame(
    feature   = f,
    statistic = as.numeric(fit$statistic),
    df        = as.numeric(fit$parameter),
    p_value   = fit$p.value,
    stringsAsFactors = FALSE
  )
}))
kw_res$p_adj_BH <- p.adjust(kw_res$p_value, method = "BH")
kw_res <- kw_res[order(kw_res$p_value), ]
fwrite(
  kw_res,
  file.path(immune_dir, "tcga_paad_mcd_subtype_xcell_kw.tsv"),
  sep   = "\t",
  quote = FALSE
)
# -------------------------------
# 6. 事后比较：Dunn test
# [修复 Bug 4] 使用 rbindlist 处理 tibble 列表，更稳健
# -------------------------------
pw_list <- lapply(features, function(f) {
  tmp <- xcell_df[, c("subtype_label", f), drop = FALSE]
  colnames(tmp)[2] <- "value"
  tmp <- tmp[!is.na(tmp$value), , drop = FALSE]
  if (nrow(tmp) == 0 || length(unique(tmp$value)) < 2) return(NULL)
  out <- rstatix::dunn_test(
    data            = tmp,
    formula         = value ~ subtype_label,
    p.adjust.method = "BH"
  )
  out$feature <- f
  as.data.frame(out)          # 转为普通 data.frame，避免 tibble 属性冲突
})
pw_list <- Filter(Negate(is.null), pw_list)   # 移除 NULL 元素
if (length(pw_list) > 0) {
  pw_res <- data.table::rbindlist(pw_list, fill = TRUE)
  fwrite(
    pw_res,
    file.path(immune_dir, "tcga_paad_mcd_subtype_xcell_posthoc.tsv"),
    sep   = "\t",
    quote = FALSE
  )
} else {
  warning("posthoc 结果为空，未输出 tcga_paad_mcd_subtype_xcell_posthoc.tsv")
}
# -------------------------------
# 7. 选择重点细胞群
# 逻辑：
# 优先保留和当前主线最相关的 immune/stromal cells
# 若匹配不足，则用最显著的前若干特征补足
#
# [修复 Bug 3] as.data.frame(check.names=TRUE) 会将 "B-cells" 变为 "B.cells"，
# 连字符和空格均转为 "."，因此正则中以 "[. -]" 或直接省略分隔符来兼容。
# -------------------------------
focus_patterns <- c(
  "CD8", "CD4", "Treg",
  "Th1", "Th2",
  "NK", "NKT",
  "Mono", "Macrophage", "M1", "M2",
  "DC",                                   # 匹配 aDC/iDC/pDC/cDC/DC
  "Fibro", "Endothel", "Stroma",
  "ImmuneScore", "MicroenvironmentScore",
  "B[.-]cell", "Plasma", "Neutrophil"    # "B-cell"→"B.cell"，用 [.-] 兼容
)
focus_regex <- paste(focus_patterns, collapse = "|")
focus_kw <- kw_res[
  !is.na(kw_res$p_adj_BH) &
    kw_res$p_adj_BH < 0.05 &
    grepl(focus_regex, kw_res$feature, ignore.case = TRUE, perl = TRUE),
  ,
  drop = FALSE
]
selected_features <- unique(head(focus_kw$feature, 12))
# 如果重点特征不足，则用整体最显著结果补足
if (length(selected_features) < 8) {
  fill_features <- kw_res$feature[
    !is.na(kw_res$p_adj_BH) &
      kw_res$p_adj_BH < 0.05
  ]
  selected_features <- unique(c(selected_features, head(fill_features, 12)))
  selected_features <- selected_features[seq_len(min(length(selected_features), 12))]
}
fwrite(
  data.frame(feature = selected_features),
  file.path(immune_dir, "tcga_paad_mcd_subtype_xcell_selected_features.tsv"),
  sep   = "\t",
  quote = FALSE
)
# -------------------------------
# 8. 热图：subtype 中位数 + row z-score
# [修复 Bug 2] match() 在亚型缺失时返回 NA，导致全-NA 行进入热图矩阵。
# 修复：过滤掉 NA 行后再转置/缩放。
# -------------------------------
if (length(selected_features) >= 2) {
  heat_df <- aggregate(
    . ~ subtype_label,
    data = xcell_df[, c("subtype_label", selected_features), drop = FALSE],
    FUN  = median
  )
  # 固定列顺序（仅保留实际存在的亚型）
  level_order <- c("MCD-low", "MCD-intermediate", "MCD-high")
  row_idx     <- match(level_order, heat_df$subtype_label)
  row_idx     <- row_idx[!is.na(row_idx)]          # [修复] 跳过缺失亚型
  heat_df     <- heat_df[row_idx, , drop = FALSE]
  rownames(heat_df) <- heat_df$subtype_label
  heat_mat <- as.matrix(heat_df[, -1, drop = FALSE])
  # 转置后：feature x subtype
  heat_mat <- t(heat_mat)
  # row z-score（仅在列数 >= 2 时有意义）
  if (ncol(heat_mat) >= 2) {
    heat_scaled <- t(scale(t(heat_mat)))
    heat_scaled[!is.finite(heat_scaled)] <- 0      # 处理 NaN / Inf
  } else {
    heat_scaled <- heat_mat
    warning("亚型数量不足 2，跳过 row z-score 标准化。")
  }
  pdf(
    file.path(immune_dir, "tcga_paad_mcd_subtype_xcell_heatmap.pdf"),
    width  = 8.5,
    height = max(5, 0.35 * nrow(heat_scaled) + 2.5)
  )
  pheatmap(
    heat_scaled,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    main         = "xCell cell enrichment by MCD subtype\n(row z-score across subtypes)"
  )
  dev.off()
} else {
  warning("可用于热图的 selected_features 不足，未输出热图。")
}
# -------------------------------
# 9. 候选箱线图
# 优先展示最显著的 6 个重点特征
# -------------------------------
plot_features <- head(selected_features, 6)
if (length(plot_features) > 0) {
  long_df <- do.call(rbind, lapply(plot_features, function(f) {
    data.frame(
      sample_id     = xcell_df$sample_id,
      subtype_label = xcell_df$subtype_label,
      feature       = f,
      score         = xcell_df[[f]],
      stringsAsFactors = FALSE
    )
  }))
  long_df$feature <- factor(long_df$feature, levels = plot_features)
  p <- ggplot(long_df, aes(x = subtype_label, y = score, fill = subtype_label)) +
    geom_boxplot(outlier.size = 0.4, width = 0.62, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.35, alpha = 0.30, color = "grey30") +
    facet_wrap(~ feature, scales = "free_y", ncol = 3) +
    labs(
      title = "Selected xCell scores by MCD subtype",
      x     = "MCD subtype",
      y     = "xCell score"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x    = element_text(angle = 30, hjust = 1),
      legend.position = "none",
      strip.text     = element_text(size = 9)
    )
  ggsave(
    file.path(immune_dir, "tcga_paad_mcd_subtype_xcell_selected_boxplot.pdf"),
    p,
    width       = 10,
    height      = 6,
    useDingbats = FALSE
  )
} else {
  warning("没有可绘图的 plot_features，未输出 selected_boxplot。")
}
# -------------------------------
# 10. 控制台输出摘要
# -------------------------------
message("08E xCell analysis finished.")
message("Top significant xCell features:")
print(head(kw_res, 20))
message("Selected features for heatmap/boxplot:")
print(selected_features)
