# =========================================================
# 08C_mcd_subtype_immune_function_ssgsea.R
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
  library(GSVA)
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

# ---- 4. 定义核心 immune/stromal function gene sets ----
gene_sets <- list(
  IFN_Gamma_Response    = c("STAT1","IFNG","CXCL9","CXCL10","CXCL11","IDO1","IRF1","HLA-DRA"),
  Cytolytic_Activity    = c("GZMA","GZMB","PRF1","NKG7","GNLY"),
  Antigen_Presentation  = c("B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"),
  Tcell_Activation      = c("CD3D","CD3E","CD2","LCK","TRBC1","IL7R"),
  Tcell_Exhaustion      = c("PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","TOX","ENTPD1"),
  Treg_Signature        = c("FOXP3","IL2RA","TIGIT","CTLA4","IKZF2"),
  Chemokine_Recruitment = c("CCL5","CXCL9","CXCL10","CXCL11","CCL2","CCL19"),
  CAF_Stromal_Activation= c("COL1A1","COL1A2","COL3A1","TAGLN","ACTA2","FAP","PDGFRB","POSTN")
)

gene_sets <- lapply(gene_sets, function(gs) intersect(gs, rownames(logcpm)))
gene_sets <- gene_sets[sapply(gene_sets, length) >= 3]
if (length(gene_sets) == 0) stop("所有 gene set 均未能映射到表达矩阵。")

# ---- 5. ssGSEA（兼容 GSVA 新旧 API） ----
# FIX 1: GSVA >= 1.50.0 改用 param 对象方式调用
#         abs.ranking 默认 FALSE（标准有符号排名），不建议改为 TRUE
gsva_new_api <- tryCatch(
  utils::packageVersion("GSVA") >= "1.50.0",
  error = function(e) FALSE
)

expr_mat <- as.matrix(logcpm)

if (gsva_new_api) {
  message("检测到 GSVA >= 1.50.0，使用新式 ssgseaParam() API")
  param        <- GSVA::ssgseaParam(
    exprData = expr_mat,
    geneSets = gene_sets,
    normalize = TRUE        # 等价于旧版默认行为
  )
  ssgsea_score <- GSVA::gsva(param, verbose = FALSE)
} else {
  message("检测到 GSVA < 1.50.0，使用旧式 gsva() API")
  ssgsea_score <- GSVA::gsva(
    expr          = expr_mat,
    gset.idx.list = gene_sets,
    method        = "ssgsea",
    kcdf          = "Gaussian",
    abs.ranking   = FALSE,   # FIX 3: 标准 ssGSEA 用有符号排名
    verbose       = FALSE
  )
}

# ---- 整理分数表 ----
ssgsea_df            <- as.data.frame(t(ssgsea_score), stringsAsFactors = FALSE)
ssgsea_df$sample_id  <- rownames(ssgsea_df)
rownames(ssgsea_df)  <- NULL

ssgsea_df <- merge(
  ssgsea_df,
  pheno[, c(sample_col, "subtype_label"), drop = FALSE],
  by.x  = "sample_id",
  by.y  = sample_col,
  all.x = TRUE,
  sort  = FALSE
)
ssgsea_df$subtype_label <- factor(
  ssgsea_df$subtype_label,
  levels = c("MCD-low", "MCD-intermediate", "MCD-high")
)

fwrite(
  ssgsea_df,
  file.path(immune_dir, "tcga_paad_mcd_subtype_ssgsea_scores.tsv"),
  sep = "\t", quote = FALSE
)

# ---- 6. 统计比较 ----
features <- setdiff(colnames(ssgsea_df), c("sample_id", "subtype_label"))

kw_res <- do.call(rbind, lapply(features, function(f) {
  tmp <- ssgsea_df[, c("subtype_label", f), drop = FALSE]
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
  file.path(immune_dir, "tcga_paad_mcd_subtype_ssgsea_kw.tsv"),
  sep = "\t", quote = FALSE
)

pw_res <- do.call(rbind, lapply(features, function(f) {
  tmp <- ssgsea_df[, c("subtype_label", f), drop = FALSE]
  colnames(tmp)[2] <- "value"
  out <- rstatix::dunn_test(tmp, value ~ subtype_label, p.adjust.method = "BH")
  out$feature <- f
  out
}))

fwrite(
  pw_res,
  file.path(immune_dir, "tcga_paad_mcd_subtype_ssgsea_posthoc.tsv"),
  sep = "\t", quote = FALSE
)

# ---- 7. 热图 ----
# heat_mat: subtype（行）× gene_set（列）
heat_df <- aggregate(
  . ~ subtype_label,
  data = ssgsea_df[, c("subtype_label", features)],
  FUN  = median
)
rownames(heat_df) <- heat_df$subtype_label
heat_mat          <- as.matrix(heat_df[, -1, drop = FALSE])

# FIX 2: 正确的行标准化方向
# heat_mat 列 = gene_set → scale() 对每列（gene_set）跨 subtype 标准化
# t() → gene_set（行）× subtype（列），即 pheatmap 需要的方向
mat_scaled <- t(scale(heat_mat))

# 极端情况：某 gene_set 在三个 subtype 中位数完全相同 → sd=0 → NaN → 置 0
mat_scaled[is.nan(mat_scaled)] <- 0

pdf(file.path(immune_dir, "tcga_paad_mcd_subtype_ssgsea_heatmap.pdf"),
    width = 8, height = 4.8)
pheatmap(
  mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,   # subtype 保持 low / intermediate / high 顺序
  main         = "Immune / stromal ssGSEA scores by MCD subtype\n(row z-score across subtypes)"
)
dev.off()

# ---- 8. 显著特征箱线图 ----
sig_features <- head(kw_res$feature[kw_res$p_adj_BH < 0.05], 6)

if (length(sig_features) > 0) {
  long_df <- do.call(rbind, lapply(sig_features, function(f) {
    data.frame(
      sample_id     = ssgsea_df$sample_id,
      subtype_label = ssgsea_df$subtype_label,
      feature       = f,
      score         = ssgsea_df[[f]],
      stringsAsFactors = FALSE
    )
  }))
  long_df$feature <- factor(long_df$feature, levels = sig_features)
  
  p <- ggplot(long_df, aes(x = subtype_label, y = score, fill = subtype_label)) +
    geom_boxplot(outlier.size = 0.4, width = 0.62, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.35, alpha = 0.30, color = "grey30") +
    facet_wrap(~ feature, scales = "free_y", ncol = 3) +
    labs(
      title = "Selected immune/stromal ssGSEA scores by MCD subtype",
      x     = "MCD subtype",
      y     = "ssGSEA score"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x    = element_text(angle = 30, hjust = 1),
      legend.position = "none",
      strip.text     = element_text(size = 9)
    )
  
  ggsave(
    file.path(immune_dir, "tcga_paad_mcd_subtype_ssgsea_selected_boxplot.pdf"),
    p, width = 10, height = 6, useDingbats = FALSE
  )
} else {
  message("没有 BH 校正后显著的 gene set（p_adj_BH < 0.05），跳过箱线图。")
}

message("08C ssGSEA analysis finished.")
print(kw_res)
