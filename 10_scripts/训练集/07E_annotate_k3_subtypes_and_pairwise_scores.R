# =========================================================
# 07E_annotate_k3_subtypes_and_pairwise_scores.R
# 作用：
# 1. 读取 K3 signature score
# 2. 根据综合 MCD 活性注释 C1/C2/C3
# 3. 输出 subtype annotation
# 4. 对三类 signature 做两两比较
#
# 修复记录：
# - 修复 pairwise 分组变量应使用 subtype_label 而非原始 cluster 编号
# - 增加 K3 簇数验证，避免硬编码 3 个标签时静默出错
# - merge 后按 sample_id 恢复原始行顺序
# - 统一使用 as.data.frame(pw$p.value) 替代冗余 as.table 转换
# =========================================================

project_root <- Sys.getenv("PROJECT_ROOT", unset = "/Users/wmz_mac/Desktop/胰腺癌")
clust_dir    <- file.path(project_root, "04_bulk_analysis", "02_clustering")
score_file   <- file.path(clust_dir, "tcga_paad_k3_signature_scores.tsv")

library(data.table)

dat <- fread(score_file, data.table = FALSE)

# ---- 列完整性检查 ----
required_cols <- c("sample_id", "cluster", "ferroptosis", "cuproptosis", "disulfidptosis")
miss <- setdiff(required_cols, colnames(dat))
if (length(miss) > 0) {
  stop("缺少必要列：", paste(miss, collapse = ", "))
}

# ---- 确认恰好 3 个簇（K3 前提） ----
n_clusters <- length(unique(dat$cluster))
if (n_clusters != 3L) {
  stop(
    "期望 K=3 聚类，但检测到 ", n_clusters, " 个簇。",
    " 请确认输入文件来自 K3 聚类结果。"
  )
}

# ---- 定义综合 MCD score ----
dat$MCD_overall_score <- rowMeans(
  dat[, c("ferroptosis", "cuproptosis", "disulfidptosis")],
  na.rm = TRUE
)

# ---- 记录原始行顺序，merge 后恢复 ----
dat$.row_order <- seq_len(nrow(dat))

# ---- 按 cluster 计算均值并决定注释 ----
cluster_mean <- aggregate(
  MCD_overall_score ~ cluster,
  data = dat,
  FUN  = mean
)
colnames(cluster_mean)[2] <- "mean_MCD_overall_score"

# 升序排列 → 依次标记 low / intermediate / high
cluster_mean <- cluster_mean[order(cluster_mean$mean_MCD_overall_score), , drop = FALSE]
cluster_mean$subtype_label <- c("MCD-low", "MCD-intermediate", "MCD-high")

# ---- 映射回样本，恢复原始顺序 ----
dat <- merge(dat, cluster_mean[, c("cluster", "subtype_label", "mean_MCD_overall_score")],
             by = "cluster", all.x = TRUE, sort = FALSE)
dat <- dat[order(dat$.row_order), ]
dat$.row_order <- NULL

# ---- 输出样本级 annotation ----
annotation_df <- dat[, c(
  "sample_id", "cluster", "subtype_label",
  "MCD_overall_score", "ferroptosis", "cuproptosis", "disulfidptosis"
)]
fwrite(annotation_df,
       file.path(clust_dir, "tcga_paad_k3_cluster_annotation.tsv"),
       sep = "\t", quote = FALSE)

# ---- 输出 cluster 汇总表 ----
fwrite(cluster_mean,
       file.path(clust_dir, "tcga_paad_k3_cluster_annotation_summary.tsv"),
       sep = "\t", quote = FALSE)

# ---- pairwise Wilcoxon（分组使用 subtype_label，结果语义清晰） ----
# 确保因子顺序与 MCD 活性一致
dat$subtype_label <- factor(
  dat$subtype_label,
  levels = c("MCD-low", "MCD-intermediate", "MCD-high")
)

score_cols <- c("ferroptosis", "cuproptosis", "disulfidptosis", "MCD_overall_score")

pairwise_res <- lapply(score_cols, function(sc) {
  pw  <- pairwise.wilcox.test(
    x              = dat[[sc]],
    g              = dat$subtype_label,   # 修复：使用语义标签而非原始编号
    p.adjust.method = "BH"
  )
  # pw$p.value 本身即 matrix，无需 as.table 中转
  mat <- as.data.frame(pw$p.value, stringsAsFactors = FALSE)
  mat$group1 <- rownames(mat)
  mat <- reshape(mat,
                 varying   = setdiff(colnames(mat), "group1"),
                 v.names   = "p_adj",
                 timevar   = "group2",
                 times     = setdiff(colnames(mat), "group1"),
                 direction = "long")
  rownames(mat) <- NULL
  mat <- mat[, c("group1", "group2", "p_adj")]
  mat <- mat[!is.na(mat$p_adj), , drop = FALSE]
  mat$signature <- sc
  mat
})

pairwise_res <- do.call(rbind, pairwise_res)
pairwise_res <- pairwise_res[, c("signature", "group1", "group2", "p_adj")]

fwrite(pairwise_res,
       file.path(clust_dir, "tcga_paad_k3_signature_score_pairwise.tsv"),
       sep = "\t", quote = FALSE)

print(cluster_mean)
print(pairwise_res)
message("K3 subtype annotation and pairwise comparison finished.")
