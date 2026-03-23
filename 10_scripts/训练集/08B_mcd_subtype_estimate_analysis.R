# =========================================================
# 08B_mcd_subtype_estimate_analysis.R
# =========================================================

# ---- 0. 路径配置 ----
project_root <- Sys.getenv("PROJECT_ROOT", unset = "/Users/wmz_mac/Desktop/胰腺癌")
proc_dir     <- file.path(project_root, "02_processed_data", "TCGA_PAAD")
clust_dir    <- file.path(project_root, "04_bulk_analysis", "02_clustering")
immune_dir   <- file.path(project_root, "04_bulk_analysis", "05_immune")
dir.create(immune_dir, recursive = TRUE, showWarnings = FALSE)

expr_file  <- file.path(proc_dir,  "tcga_paad_expr_all_tp_symbol_counts_filtered.rds")
pheno_file <- file.path(clust_dir, "tcga_paad_mcd_subtype_phenotype.tsv")

# ---- 1. 加载包 ----
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(estimate)
  library(ggplot2)
  library(ggpubr)    # stat_pvalue_manual
  library(rstatix)
})

# ---- 2. 读取数据 ----
expr  <- readRDS(expr_file)
pheno <- fread(pheno_file, data.table = FALSE)

if (!is.matrix(expr)) expr <- as.matrix(expr)
stopifnot(!is.null(rownames(expr)))
stopifnot(!anyDuplicated(rownames(expr)))
stopifnot(!anyDuplicated(colnames(expr)))

sample_col <- intersect(c("sample_id", "sample_barcode"), colnames(pheno))[1]
if (is.na(sample_col)) {
  stop("phenotype 中未找到样本列（需含 sample_id 或 sample_barcode）。")
}
if (!"subtype_label" %in% colnames(pheno)) {
  stop("phenotype 中缺少 subtype_label 列。")
}

# ---- 3. 样本对齐 ----
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
message(sprintf("样本对齐完成：%d 个样本进入分析。", ncol(expr)))

# ---- 4. 表达标准化：logCPM ----
logcpm <- edgeR::cpm(expr, log = TRUE, prior.count = 1)

keep_gene <- !is.na(rownames(logcpm)) & rownames(logcpm) != ""
logcpm    <- logcpm[keep_gene, , drop = FALSE]

# 保存原始样本名，用于后续修复 ESTIMATE 的名称转换
# ESTIMATE 内部调用 make.names()，会将 TCGA 条形码中的 '-' 转为 '.'
orig_sample_ids <- colnames(logcpm)

# ---- 5. 写入 ESTIMATE 输入 ----
estimate_input_txt  <- file.path(immune_dir, "tcga_paad_estimate_input_logcpm.txt")
estimate_input_gct  <- file.path(immune_dir, "tcga_paad_estimate_input_common_genes.gct")
estimate_output_gct <- file.path(immune_dir, "tcga_paad_estimate_scores.gct")

write.table(
  logcpm,
  file      = estimate_input_txt,
  sep       = "\t",
  quote     = FALSE,
  col.names = NA   # ESTIMATE 要求行名列的列头为空
)

estimate::filterCommonGenes(
  input.f  = estimate_input_txt,
  output.f = estimate_input_gct,
  id       = "GeneSymbol"
)

# FIX #1: RNA-seq 数据应使用 "illumina" 平台，而非 "affymetrix"
estimate::estimateScore(
  input.ds  = estimate_input_gct,
  output.ds = estimate_output_gct,
  platform  = "illumina"   # 修正：TCGA 为 Illumina RNA-seq，不是微阵列
)

# ---- 6. 读取 ESTIMATE 输出 ----
score_raw <- fread(estimate_output_gct, skip = 2, data.table = FALSE)

# FIX #5: 按名称重命名，不依赖列位置
names(score_raw)[names(score_raw) == "NAME"]        <- "feature"
names(score_raw)[names(score_raw) == "Description"] <- "Description"

score_mat <- as.matrix(score_raw[, !colnames(score_raw) %in% c("feature", "Description"),
                                 drop = FALSE])
rownames(score_mat) <- score_raw$feature

score_df <- as.data.frame(t(score_mat), stringsAsFactors = FALSE)

# FIX #2: 将 ESTIMATE 输出的 make.names() 名称映射回原始 TCGA 条形码
# make.names() 将 '-' 转为 '.'，需要反向查找原始 ID
score_df$sample_id <- orig_sample_ids[
  match(rownames(score_df), make.names(orig_sample_ids))
]

unmatched_n <- sum(is.na(score_df$sample_id))
if (unmatched_n > 0) {
  warning(sprintf(
    "%d 个 ESTIMATE 输出样本无法映射回原始 ID，将被排除。",
    unmatched_n
  ))
  score_df <- score_df[!is.na(score_df$sample_id), , drop = FALSE]
}
rownames(score_df) <- NULL

# 转数值
score_cols <- setdiff(colnames(score_df), "sample_id")
for (cc in score_cols) {
  score_df[[cc]] <- as.numeric(score_df[[cc]])
}

# FIX #3: 补算 TumorPurity 并夹紧到 [0, 1]
if (!"TumorPurity" %in% colnames(score_df) && "ESTIMATEScore" %in% colnames(score_df)) {
  raw_purity <- cos(0.6049872018 + 0.0001467884 * score_df$ESTIMATEScore)
  score_df$TumorPurity <- pmax(0, pmin(1, raw_purity))  # 夹紧至 [0, 1]
}

# ---- 7. 合并 phenotype ----
plot_df <- merge(
  score_df,
  pheno[, c(sample_col, "subtype_label"), drop = FALSE],
  by.x  = "sample_id",
  by.y  = sample_col,
  all.x = TRUE,
  sort  = FALSE
)

# 检查合并质量
n_no_subtype <- sum(is.na(plot_df$subtype_label))
if (n_no_subtype > 0) {
  warning(sprintf("%d 个样本合并后 subtype_label 为 NA，请检查样本 ID 匹配。", n_no_subtype))
}
message(sprintf("合并后数据：%d 行，%d 列。", nrow(plot_df), ncol(plot_df)))

plot_df$subtype_label <- factor(
  plot_df$subtype_label,
  levels = c("MCD-low", "MCD-intermediate", "MCD-high")
)

fwrite(
  plot_df,
  file.path(immune_dir, "tcga_paad_mcd_subtype_estimate_scores.tsv"),
  sep   = "\t",
  quote = FALSE
)

# ---- 8. 统计比较 ----
target_scores <- intersect(
  c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity"),
  colnames(plot_df)
)

kw_res <- do.call(rbind, lapply(target_scores, function(s) {
  tmp <- plot_df[!is.na(plot_df$subtype_label), c("sample_id", "subtype_label", s),
                 drop = FALSE]
  colnames(tmp)[3] <- "value"
  fit <- kruskal.test(value ~ subtype_label, data = tmp)
  data.frame(
    feature   = s,
    statistic = as.numeric(fit$statistic),
    df        = as.numeric(fit$parameter),
    p_value   = fit$p.value,
    stringsAsFactors = FALSE
  )
}))
kw_res$p_adj_BH <- p.adjust(kw_res$p_value, method = "BH")

fwrite(
  kw_res,
  file.path(immune_dir, "tcga_paad_mcd_subtype_estimate_kw.tsv"),
  sep   = "\t",
  quote = FALSE
)

# Dunn post hoc
pw_res <- do.call(rbind, lapply(target_scores, function(s) {
  tmp <- plot_df[!is.na(plot_df$subtype_label), c("sample_id", "subtype_label", s),
                 drop = FALSE]
  colnames(tmp)[3] <- "value"
  out <- rstatix::dunn_test(tmp, value ~ subtype_label, p.adjust.method = "BH")
  out$feature <- s
  out
}))

fwrite(
  pw_res,
  file.path(immune_dir, "tcga_paad_mcd_subtype_estimate_posthoc.tsv"),
  sep   = "\t",
  quote = FALSE
)

# ---- 9. 作图 ----
# FIX #4: 在箱线图上叠加显著性标注（仅显示 p.adj < 0.05 的配对）
plot_one_score <- function(df, score_name, pw_df, out_file) {
  # 过滤出当前指标的显著配对，计算标注 y 位置
  pw_sub <- pw_df[pw_df$feature == score_name & pw_df$p.adj < 0.05, , drop = FALSE]
  
  y_vals   <- df[[score_name]][!is.na(df[[score_name]]) & !is.na(df$subtype_label)]
  y_max    <- max(y_vals, na.rm = TRUE)
  y_range  <- diff(range(y_vals, na.rm = TRUE))
  
  p <- ggplot(df[!is.na(df$subtype_label), ],
              aes(x = subtype_label, y = .data[[score_name]], fill = subtype_label)) +
    geom_boxplot(outlier.size = 0.6, width = 0.65, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.7, alpha = 0.35, color = "grey30") +
    labs(
      title = paste0(score_name, " by MCD subtype"),
      x     = "MCD subtype",
      y     = score_name
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title    = element_text(hjust = 0.5),
      axis.text.x   = element_text(angle = 30, hjust = 1),
      legend.position = "none"
    )
  
  # 若有显著配对，添加显著性括号
  if (nrow(pw_sub) > 0) {
    # 依次为每对分配递增的 y 位置，避免重叠
    pw_sub <- pw_sub[order(pw_sub$p.adj), ]
    pw_sub$y.position <- y_max + y_range * seq(0.08, by = 0.10,
                                               length.out = nrow(pw_sub))
    p <- p + ggpubr::stat_pvalue_manual(
      pw_sub,
      label         = "p.adj.signif",
      tip.length    = 0.01,
      hide.ns       = FALSE
    )
  }
  
  ggsave(out_file, p, width = 5.2, height = 4.8, useDingbats = FALSE)
}

for (s in target_scores) {
  plot_one_score(
    plot_df,
    s,
    pw_res,
    file.path(immune_dir, paste0("tcga_paad_", s, "_by_mcd_subtype.pdf"))
  )
}

message("08B ESTIMATE analysis finished.")
print(kw_res)
