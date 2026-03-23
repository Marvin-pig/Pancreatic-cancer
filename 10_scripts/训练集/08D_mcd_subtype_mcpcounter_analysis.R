# =========================================================
# 08D_mcd_subtype_mcpcounter_analysis.R
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
  library(MCPcounter)
  library(ggplot2)
  library(pheatmap)
  library(rstatix)
})

# ---- 1. 读取数据 ----
expr  <- readRDS(expr_file)
pheno <- fread(pheno_file, data.table = FALSE)

if (!is.matrix(expr)) expr <- as.matrix(expr)
stopifnot(!is.null(rownames(expr)))
stopifnot(!anyDuplicated(rownames(expr)))
stopifnot(!anyDuplicated(colnames(expr)))

sample_col <- intersect(c("sample_id", "sample_barcode"), colnames(pheno))[1]
if (is.na(sample_col)) stop("phenotype 中未找到 sample_id / sample_barcode")
if (!"subtype_label" %in% colnames(pheno)) stop("phenotype 中缺少 subtype_label")

# ---- 2. 样本对齐 ----
matched_idx <- match(colnames(expr), pheno[[sample_col]])
if (any(is.na(matched_idx))) {
  keep        <- !is.na(matched_idx)
  expr        <- expr[, keep, drop = FALSE]
  matched_idx <- matched_idx[keep]
}
pheno <- pheno[matched_idx, , drop = FALSE]
stopifnot(identical(colnames(expr), pheno[[sample_col]]))

# ---- 3. logCPM ----
logcpm    <- edgeR::cpm(expr, log = TRUE, prior.count = 1)
keep_gene <- !is.na(rownames(logcpm)) & rownames(logcpm) != ""
logcpm    <- logcpm[keep_gene, , drop = FALSE]

# ---- 4. MCPcounter ----
mcpc_mat <- MCPcounter::MCPcounter.estimate(
  expression   = logcpm,
  featuresType = "HUGO_symbols"
)

# [Fix 7] 显式声明 check.names = TRUE，空格转为点号，避免后续列名解析问题
mcpc_df <- as.data.frame(t(mcpc_mat), check.names = TRUE, stringsAsFactors = FALSE)
mcpc_df$sample_id <- rownames(mcpc_df)
rownames(mcpc_df) <- NULL

mcpc_df <- merge(
  mcpc_df,
  pheno[, c(sample_col, "subtype_label"), drop = FALSE],
  by.x  = "sample_id",
  by.y  = sample_col,
  all.x = TRUE,
  sort  = FALSE
)

subtype_levels        <- c("MCD-low", "MCD-intermediate", "MCD-high")
mcpc_df$subtype_label <- factor(mcpc_df$subtype_label, levels = subtype_levels)

fwrite(
  mcpc_df,
  file.path(immune_dir, "tcga_paad_mcd_subtype_mcpcounter_scores.tsv"),
  sep = "\t", quote = FALSE
)

# ---- 5. 统计比较 ----
features <- setdiff(colnames(mcpc_df), c("sample_id", "subtype_label"))

kw_res <- do.call(rbind, lapply(features, function(f) {
  tmp <- mcpc_df[, c("subtype_label", f), drop = FALSE]
  colnames(tmp)[2] <- "value"
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
  file.path(immune_dir, "tcga_paad_mcd_subtype_mcpcounter_kw.tsv"),
  sep = "\t", quote = FALSE
)

# [Fix 5] 将 tibble 转为 data.frame，避免 rbind 时列类型冲突
pw_res <- do.call(rbind, lapply(features, function(f) {
  tmp <- mcpc_df[, c("subtype_label", f), drop = FALSE]
  colnames(tmp)[2] <- "value"
  out         <- rstatix::dunn_test(tmp, value ~ subtype_label, p.adjust.method = "BH")
  out$feature <- f
  as.data.frame(out)           # tibble → data.frame
}))

fwrite(
  pw_res,
  file.path(immune_dir, "tcga_paad_mcd_subtype_mcpcounter_posthoc.tsv"),
  sep = "\t", quote = FALSE
)

# ---- 6. 热图（subtype 中位数，row z-score）----
# [Fix 1] median 加 na.rm = TRUE，防止 NA 传播
heat_df <- aggregate(
  . ~ subtype_label,
  data = mcpc_df[, c("subtype_label", features)],
  FUN  = function(x) median(x, na.rm = TRUE)
)

# [Fix 3] 显式按 factor levels 排序，保证列顺序为 low → mid → high
heat_df           <- heat_df[match(subtype_levels, as.character(heat_df$subtype_label)), ]
rownames(heat_df) <- as.character(heat_df$subtype_label)
heat_mat          <- as.matrix(heat_df[, -1, drop = FALSE])

# row z-score：scale() 对列（即每个 feature 跨 3 个亚型）标准化，t() 后 feature 为行
mat_scaled <- t(scale(heat_mat))
# [Fix 2] 统一替换 NaN / NA / Inf / -Inf
mat_scaled[!is.finite(mat_scaled)] <- 0

# [Fix 6] 添加列注释，标明亚型
ann_col <- data.frame(
  Subtype   = factor(colnames(mat_scaled), levels = subtype_levels),
  row.names = colnames(mat_scaled)
)

pdf(file.path(immune_dir, "tcga_paad_mcd_subtype_mcpcounter_heatmap.pdf"),
    width = 8.5, height = 5.2)
pheatmap(
  mat_scaled,
  cluster_rows   = TRUE,
  cluster_cols   = FALSE,
  annotation_col = ann_col,           # 亚型颜色条
  main = "MCPcounter cell abundance by MCD subtype\n(row z-score across subtypes)"
)
dev.off()

# ---- 7. 候选箱线图 ----
sig_features <- head(kw_res$feature[kw_res$p_adj_BH < 0.05], 6)

if (length(sig_features) > 0) {
  long_df <- do.call(rbind, lapply(sig_features, function(f) {
    data.frame(
      sample_id     = mcpc_df$sample_id,
      subtype_label = mcpc_df$subtype_label,
      feature       = f,
      score         = mcpc_df[[f]],
      stringsAsFactors = FALSE
    )
  }))
  long_df$feature <- factor(long_df$feature, levels = sig_features)
  
  p <- ggplot(long_df, aes(x = subtype_label, y = score, fill = subtype_label)) +
    geom_boxplot(outlier.size = 0.4, width = 0.62, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.35, alpha = 0.30, color = "grey30") +
    facet_wrap(~ feature, scales = "free_y", ncol = 3) +
    labs(
      title = "Selected MCPcounter cell abundance scores by MCD subtype",
      x     = "MCD subtype",
      y     = "MCPcounter score"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x     = element_text(angle = 30, hjust = 1),
      legend.position = "none",
      strip.text      = element_text(size = 9)
    )
  
  # [Fix 4] 移除已废弃的 useDingbats 参数
  ggsave(
    file.path(immune_dir, "tcga_paad_mcd_subtype_mcpcounter_selected_boxplot.pdf"),
    p, width = 10, height = 6
  )
} else {
  message("没有 BH<0.05 的显著细胞群，跳过候选箱线图。")
}

message("08D MCPcounter analysis finished.")
print(kw_res)
