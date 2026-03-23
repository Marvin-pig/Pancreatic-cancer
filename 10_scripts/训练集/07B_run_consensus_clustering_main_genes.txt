# =========================================================
# 07B_run_consensus_clustering_main_genes.R
# =========================================================

# ---------- 包安装（仅首次运行时安装缺失的包）----------

# CRAN 包
cran_pkgs <- c("data.table")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Bioconductor 包
bioc_pkgs <- c("ConsensusClusterPlus")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# ---------- 加载包 ----------
library(data.table)
library(ConsensusClusterPlus)

# ========= 以下为原有脚本逻辑，保持不变 =========
project_root <- Sys.getenv("PROJECT_ROOT", unset = "/Users/wmz_mac/Desktop/胰腺癌")
# ... （其余代码不变）

# =========================================================
# 07B_run_consensus_clustering_main_genes.R
# 作用：
# 1. 读取 main gene clustering matrix
# 2. 过滤零方差基因
# 3. 对每个基因做 Z-score 标准化
# 4. 运行 ConsensusClusterPlus (K=2~4)
# 5. 计算 ICL 并基于 PAC 推荐最优 K
# 6. 导出 cluster assignment（所有 K + 最优 K）
# 7. 导出审计表 + sessionInfo
# =========================================================

project_root <- Sys.getenv("PROJECT_ROOT", unset = "/Users/wmz_mac/Desktop/胰腺癌")
clust_dir    <- file.path(project_root, "04_bulk_analysis", "02_clustering")
out_dir      <- file.path(clust_dir, "ConsensusClusterPlus_main")
dir.create(clust_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir,   recursive = TRUE, showWarnings = FALSE)

expr_file <- file.path(clust_dir, "tcga_paad_main_gene_expr_for_clustering.rds")

library(data.table)
library(ConsensusClusterPlus)


expr <- readRDS(expr_file)

# ---------- 基础检查 ----------
stopifnot(is.matrix(expr) || is.data.frame(expr))
expr <- as.matrix(expr)
stopifnot(!is.null(rownames(expr)))
stopifnot(!is.null(colnames(expr)))
stopifnot(!anyDuplicated(rownames(expr)))
stopifnot(!anyDuplicated(colnames(expr)))

n_genes_raw    <- nrow(expr)
n_samples_raw  <- ncol(expr)
message(sprintf("输入矩阵：%d 基因 × %d 样本", n_genes_raw, n_samples_raw))

# ---------- 过滤零方差基因（z-score 前必须做）----------
gene_var <- apply(expr, 1, var, na.rm = TRUE)
zero_var_genes <- sum(gene_var == 0 | is.na(gene_var))
if (zero_var_genes > 0) {
  message(sprintf("过滤零方差基因：%d 个", zero_var_genes))
  expr <- expr[gene_var > 0 & !is.na(gene_var), , drop = FALSE]
}
n_genes_filtered <- nrow(expr)

# ---------- row-wise z-score ----------
expr_scaled <- t(scale(t(expr)))

# 记录并处理残余 NA（理论上过滤后不应再出现）
na_count <- sum(is.na(expr_scaled))
if (na_count > 0) {
  warning(sprintf("z-score 后仍有 %d 个 NA，已替换为 0，请检查输入数据", na_count))
  expr_scaled[is.na(expr_scaled)] <- 0
}

saveRDS(
  expr_scaled,
  file.path(clust_dir, "consensuscluster_input_scaled_main_expr.rds")
)

# ---------- 运行 ConsensusClusterPlus ----------
maxK <- 4
ccp_res <- ConsensusClusterPlus(
  d           = expr_scaled,
  maxK        = maxK,
  reps        = 1000,
  pItem       = 0.8,
  pFeature    = 1,
  clusterAlg  = "hc",
  innerLinkage = "ward.D2",   # 修正：ward.D2 比 average 更适合 expression 聚类
  finalLinkage = "ward.D2",
  distance    = "pearson",
  seed        = 20260317,
  plot        = "png",
  title       = out_dir,
  verbose     = TRUE
)

saveRDS(
  ccp_res,
  file.path(clust_dir, "tcga_paad_consensusclusterplus_main_results.rds")
)

# ---------- 计算 ICL ----------
icl <- calcICL(ccp_res, title = out_dir, plot = "png")
saveRDS(
  icl,
  file.path(clust_dir, "tcga_paad_consensusclusterplus_main_icl.rds")
)

# ---------- 基于 PAC 推荐最优 K ----------
# PAC = 共识矩阵中 [0.1, 0.9] 区间内的比例（越小越好）
calc_pac <- function(ccp_res, k, lower = 0.1, upper = 0.9) {
  cm <- ccp_res[[k]]$consensusMatrix
  # 仅取上三角（避免重复计数对角线）
  vals <- cm[upper.tri(cm)]
  mean(vals > lower & vals < upper)
}

pac_values <- sapply(2:maxK, function(k) calc_pac(ccp_res, k))
names(pac_values) <- paste0("K", 2:maxK)
best_k <- as.integer(sub("K", "", names(which.min(pac_values))))
message(sprintf("PAC 推荐最优 K = %d", best_k))
message("各 K 的 PAC 值：")
print(round(pac_values, 4))

# ---------- 导出 cluster assignment ----------
for (k in 2:maxK) {
  cls <- data.frame(
    sample_id = colnames(expr_scaled),
    cluster   = paste0("C", ccp_res[[k]]$consensusClass),
    stringsAsFactors = FALSE
  )
  fwrite(
    cls,
    file.path(clust_dir, paste0("tcga_paad_cluster_assignment_k", k, ".tsv")),
    sep = "\t", quote = FALSE
  )
}

# 单独保存最优 K 结果，方便下游直接引用
best_cls <- data.frame(
  sample_id = colnames(expr_scaled),
  cluster   = paste0("C", ccp_res[[best_k]]$consensusClass),
  stringsAsFactors = FALSE
)
fwrite(
  best_cls,
  file.path(clust_dir, "tcga_paad_cluster_assignment_best_k.tsv"),
  sep = "\t", quote = FALSE
)

# ---------- 审计表 ----------
audit <- data.frame(
  metric = c(
    "n_genes_raw",
    "n_genes_after_filter",
    "n_zero_var_genes_removed",
    "n_na_replaced_post_zscore",
    "n_samples",
    "maxK",
    "best_k_by_PAC",
    paste0("PAC_K", 2:maxK),
    "reps",
    "pItem",
    "pFeature",
    "clusterAlg",
    "innerLinkage",
    "finalLinkage",
    "distance",
    "seed"
  ),
  value = c(
    n_genes_raw,
    n_genes_filtered,
    zero_var_genes,
    na_count,
    n_samples_raw,
    maxK,
    best_k,
    round(pac_values, 6),
    1000,
    0.8,
    1,
    "hc",
    "ward.D2",
    "ward.D2",
    "pearson",
    20260317
  ),
  stringsAsFactors = FALSE
)
fwrite(
  audit,
  file.path(clust_dir, "tcga_paad_consensus_clustering_audit.tsv"),
  sep = "\t", quote = FALSE
)

# ---------- 保存 sessionInfo ----------
writeLines(
  capture.output(sessionInfo()),
  file.path(clust_dir, "tcga_paad_consensus_clustering_sessionInfo.txt")
)

message("Consensus clustering completed successfully.")
message(sprintf("最优 K = %d，结果已写入 %s", best_k, clust_dir))
print(audit)
