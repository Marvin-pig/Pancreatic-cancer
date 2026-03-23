# =========================================================
# 09A_mcd_subtype_deg_analysis.R
# 作用：
#   1. 读取 all-TP filtered count matrix 与 MCD phenotype
#   2. 进行 limma-voom 差异分析
#   3. 输出三组 subtype 对比的 DEG 结果
#   4. 绘制主对比（MCD-high vs MCD-low）的 volcano 和 heatmap
# =========================================================

project_root <- Sys.getenv("PROJECT_ROOT", unset = "/Users/wmz_mac/Desktop/胰腺癌")
proc_dir  <- file.path(project_root, "02_processed_data", "TCGA_PAAD")
clust_dir <- file.path(project_root, "04_bulk_analysis", "02_clustering")
deg_dir   <- file.path(project_root, "04_bulk_analysis", "06_deg")
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)

expr_file  <- file.path(proc_dir,  "tcga_paad_expr_all_tp_symbol_counts_filtered.rds")
pheno_file <- file.path(clust_dir, "tcga_paad_mcd_subtype_phenotype.tsv")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(pheatmap)
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

# 修复 2：显式检查长度，避免依赖 [1] 隐式返回 NA
sample_col_candidates <- intersect(c("sample_id", "sample_barcode"), colnames(pheno))
if (length(sample_col_candidates) == 0) {
  stop("phenotype 中未找到 sample_id / sample_barcode")
}
sample_col <- sample_col_candidates[1]

if (!"subtype_label" %in% colnames(pheno)) {
  stop("phenotype 中缺少 subtype_label")
}

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

pheno$subtype_label <- factor(
  pheno$subtype_label,
  levels = c("MCD-low", "MCD-intermediate", "MCD-high")
)

# -------------------------------
# 3. edgeR + voom
# -------------------------------
dge    <- DGEList(counts = expr)
dge    <- calcNormFactors(dge)
design <- model.matrix(~ 0 + subtype_label, data = pheno)

# 修复 1：将设计矩阵列名中的连字符替换为下划线，避免 makeContrasts 将其解析为减法
colnames(design) <- gsub("-", "_",
                         gsub("^subtype_label", "", colnames(design)))
# 列名现在为：MCD_low / MCD_intermediate / MCD_high

v    <- voom(dge, design = design, plot = FALSE)
fit  <- lmFit(v, design)

contrast_mat <- makeContrasts(
  high_vs_low              = MCD_high - MCD_low,
  high_vs_intermediate     = MCD_high - MCD_intermediate,
  intermediate_vs_low      = MCD_intermediate - MCD_low,
  levels                   = design
)

fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

# -------------------------------
# 4. 导出 DEG 结果
# -------------------------------
export_deg <- function(coef_name, out_file) {
  res <- topTable(
    fit2,
    coef    = coef_name,
    number  = Inf,
    sort.by = "P"
  )
  res$gene_symbol <- rownames(res)
  rownames(res)   <- NULL
  res <- res[, c("gene_symbol", setdiff(colnames(res), "gene_symbol"))]
  fwrite(res, out_file, sep = "\t", quote = FALSE)
  return(res)
}

deg_high_low <- export_deg(
  "high_vs_low",
  file.path(deg_dir, "tcga_paad_mcd_high_vs_low_deg.tsv")
)

deg_high_intermediate <- export_deg(
  "high_vs_intermediate",
  file.path(deg_dir, "tcga_paad_mcd_high_vs_intermediate_deg.tsv")
)

deg_intermediate_low <- export_deg(
  "intermediate_vs_low",
  file.path(deg_dir, "tcga_paad_mcd_intermediate_vs_low_deg.tsv")
)

# 汇总表
deg_summary <- data.frame(
  contrast = c(
    "MCD-high vs MCD-low",
    "MCD-high vs MCD-intermediate",
    "MCD-intermediate vs MCD-low"
  ),
  n_FDR_0.05 = c(
    sum(deg_high_low$adj.P.Val          < 0.05, na.rm = TRUE),
    sum(deg_high_intermediate$adj.P.Val < 0.05, na.rm = TRUE),
    sum(deg_intermediate_low$adj.P.Val  < 0.05, na.rm = TRUE)
  ),
  n_FDR_0.05_logFC_1 = c(
    sum(deg_high_low$adj.P.Val          < 0.05 & abs(deg_high_low$logFC)          >= 1, na.rm = TRUE),
    sum(deg_high_intermediate$adj.P.Val < 0.05 & abs(deg_high_intermediate$logFC) >= 1, na.rm = TRUE),
    sum(deg_intermediate_low$adj.P.Val  < 0.05 & abs(deg_intermediate_low$logFC)  >= 1, na.rm = TRUE)
  )
)

fwrite(
  deg_summary,
  file.path(deg_dir, "tcga_paad_mcd_deg_summary.tsv"),
  sep = "\t", quote = FALSE
)

# -------------------------------
# 5. Volcano: MCD-high vs MCD-low
# -------------------------------
vol_df <- deg_high_low
vol_df$significance <- "Not significant"
vol_df$significance[vol_df$adj.P.Val < 0.05 & vol_df$logFC >=  1] <- "Up in MCD-high"
vol_df$significance[vol_df$adj.P.Val < 0.05 & vol_df$logFC <= -1] <- "Up in MCD-low"

# 修复 3：显式指定颜色，避免 ggplot2 自动配色不稳定
vol_colors <- c(
  "Up in MCD-high"  = "#d73027",
  "Up in MCD-low"   = "#4575b4",
  "Not significant" = "grey60"
)

p_vol <- ggplot(vol_df, aes(x = logFC, y = -log10(adj.P.Val), color = significance)) +
  geom_point(size = 0.8, alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  scale_color_manual(values = vol_colors) +
  theme_bw(base_size = 11) +
  labs(
    title = "DEG volcano: MCD-high vs MCD-low",
    x     = "log2 Fold Change",
    y     = "-log10(FDR)",
    color = NULL
  )

ggsave(
  file.path(deg_dir, "tcga_paad_mcd_high_vs_low_deg_volcano.pdf"),
  p_vol, width = 6.2, height = 5.2
)

# -------------------------------
# 6. Heatmap: top DEG
# -------------------------------
sig_df   <- deg_high_low[deg_high_low$adj.P.Val < 0.05 & abs(deg_high_low$logFC) >= 1, ]
sig_df   <- sig_df[order(sig_df$adj.P.Val), ]
top_up   <- head(sig_df[sig_df$logFC > 0, "gene_symbol"], 25)
top_down <- head(sig_df[sig_df$logFC < 0, "gene_symbol"], 25)
heat_genes <- unique(c(top_up, top_down))
heat_genes <- heat_genes[heat_genes %in% rownames(expr)]

if (length(heat_genes) >= 2) {
  logcpm   <- cpm(dge, log = TRUE, prior.count = 1)
  heat_mat <- logcpm[heat_genes, , drop = FALSE]

  ann_col <- data.frame(
    subtype_label = pheno$subtype_label,
    row.names     = colnames(heat_mat)
  )

  # 修复 5：显式传入注释色板，与 subtype_label 因子水平对应
  ann_colors <- list(
    subtype_label = c(
      "MCD-low"          = "#4575b4",
      "MCD-intermediate" = "#fee090",
      "MCD-high"         = "#d73027"
    )
  )

  pdf(file.path(deg_dir, "tcga_paad_mcd_high_vs_low_deg_heatmap.pdf"),
      width = 9, height = 8)
  pheatmap(
    heat_mat,
    annotation_col  = ann_col,
    annotation_colors = ann_colors,
    scale           = "row",
    cluster_cols    = TRUE,
    cluster_rows    = TRUE,
    show_rownames   = TRUE,
    show_colnames   = FALSE,
    main            = "Top DEG heatmap: MCD-high vs MCD-low"
  )
  invisible(dev.off())   # 修复 4：抑制 "null device 1" 控制台输出
}

message("09A DEG analysis finished.")
print(deg_summary)
