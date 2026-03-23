# =========================================================
# 07D_k3_signature_score_comparison.R
# 作用：
# 1. 读取 all-TP 表达矩阵
# 2. 读取 analysis-ready gene set
# 3. 读取 K3 cluster assignment
# 4. 计算 ferroptosis / cuproptosis / disulfidptosis signature score
# 5. 输出每个样本的 score 表
# 6. 比较 K3 各 cluster 的 score 差异（Kruskal-Wallis + Dunn 事后检验）
# =========================================================

project_root <- Sys.getenv("PROJECT_ROOT", unset = "/Users/wmz_mac/Desktop/胰腺癌")
proc_dir     <- file.path(project_root, "02_processed_data",  "TCGA_PAAD")
clust_dir    <- file.path(project_root, "04_bulk_analysis",   "02_clustering")

expr_file <- file.path(proc_dir,  "tcga_paad_expr_all_tp_symbol_counts_filtered.rds")
gene_file <- file.path(proc_dir,  "tcga_paad_metabolic_cell_death_gene_set_analysis_ready.tsv")
k3_file   <- file.path(clust_dir, "tcga_paad_cluster_assignment_k3.tsv")

# ---- 依赖包 --------------------------------------------------------
library(data.table)
# dunn.test 用于 Kruskal-Wallis 事后两两比较（需要安装）
if (!requireNamespace("dunn.test", quietly = TRUE)) {
  install.packages("dunn.test")
}
library(dunn.test)

# ---- 输入文件存在性检查 -------------------------------------------
for (f in c(expr_file, gene_file, k3_file)) {
  if (!file.exists(f)) stop("输入文件不存在: ", f)
}

# ---- 确保输出目录存在 ---------------------------------------------
dir.create(clust_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 读取数据 -----------------------------------------------------
expr    <- readRDS(expr_file)
gene_df <- fread(gene_file, data.table = FALSE)
k3      <- fread(k3_file,   data.table = FALSE)

# ---- 表达矩阵基本校验 ---------------------------------------------
stopifnot(!is.null(rownames(expr)))
stopifnot(!anyDuplicated(rownames(expr)))
stopifnot(!anyDuplicated(colnames(expr)))

# ---- k3 列校验 ---------------------------------------------------
stopifnot("sample_id" %in% colnames(k3))
stopifnot("cluster"   %in% colnames(k3))

# ---- 仅用 main tier gene set --------------------------------------
gene_main <- gene_df[gene_df$tier == "main", ]

# ---- 构建每类唯一 gene symbol，并过滤到表达矩阵中存在的基因 --------
split_genes <- split(gene_main$gene_symbol, gene_main$death_type)
split_genes <- lapply(split_genes, unique)

# 报告过滤前后基因数量（P1 修复：不再静默丢弃基因）
split_genes <- lapply(names(split_genes), function(tp) {
  g_all      <- split_genes[[tp]]
  g_filtered <- intersect(g_all, rownames(expr))
  message(sprintf("[基因过滤] %s: %d 个 → %d 个（矩阵中存在）",
                  tp, length(g_all), length(g_filtered)))
  if (length(g_filtered) == 0) {
    stop("death_type 无可用基因（全部不在表达矩阵中）: ", tp)
  }
  g_filtered
})
names(split_genes) <- unique(gene_main$death_type)

# ---- 表达矩阵预处理 ------------------------------------------------
# 修复 P2：若为原始计数（文件名含 counts），先 log1p 再 z-score
# 如已是 log-normalized，注释掉下一行
expr <- log1p(expr)

expr_z <- t(scale(t(expr)))         # 行（基因）标准化
# 修复：仅将全零方差（恒定表达）基因的 NA 赋 0，并给出警告
n_na_genes <- sum(is.na(expr_z[, 1]))
if (n_na_genes > 0) {
  warning(sprintf("%d 个基因方差为 0（恒定表达），z-score 置为 0", n_na_genes))
}
expr_z[is.na(expr_z)] <- 0

# ---- 计算 signature score（每样本各类型均值）----------------------
score_list <- lapply(names(split_genes), function(tp) {
  genes <- split_genes[[tp]]
  score <- colMeans(expr_z[genes, , drop = FALSE])
  data.frame(
    sample_id  = names(score),
    death_type = tp,
    score      = as.numeric(score)
  )
})
score_long <- do.call(rbind, score_list)

# ---- 长表转宽表（P3 修复：替代脆弱的 reshape()）------------------
score_wide <- data.frame(
  sample_id = unique(score_long$sample_id),
  stringsAsFactors = FALSE
)
for (tp in names(split_genes)) {
  sub_df <- score_long[score_long$death_type == tp, c("sample_id", "score")]
  colnames(sub_df)[2] <- tp
  score_wide <- merge(score_wide, sub_df, by = "sample_id", all.x = TRUE)
}

# ---- 验证宽表列名（P0 修复）--------------------------------------
score_cols <- names(split_genes)   # 动态获取，不再硬编码
missing_cols <- setdiff(score_cols, colnames(score_wide))
if (length(missing_cols) > 0) {
  stop("score_wide 中缺少以下列: ", paste(missing_cols, collapse = ", "))
}

# ---- 合并 K3 cluster assignment ----------------------------------
dat <- merge(k3, score_wide, by = "sample_id", all.x = TRUE, sort = FALSE)

# 报告合并后 NA 数量
n_na <- sum(is.na(dat[[score_cols[1]]]))
if (n_na > 0) {
  warning(sprintf("合并后有 %d 个样本 score 为 NA（在 score_wide 中无匹配）", n_na))
}

# ---- 保存每样本 score 宽表 ----------------------------------------
fwrite(
  dat,
  file.path(clust_dir, "tcga_paad_k3_signature_scores.tsv"),
  sep   = "\t",
  quote = FALSE
)

# ---- Kruskal-Wallis 检验（P1 修复：去除 NA 后再检验）-------------
kw_res <- lapply(score_cols, function(sc) {
  sub <- dat[!is.na(dat[[sc]]) & !is.na(dat$cluster), ]
  
  # 校验每组样本量 ≥ 2
  grp_n <- table(sub$cluster)
  if (any(grp_n < 2)) {
    warning(sprintf("[%s] 某些 cluster 样本量 < 2，Kruskal-Wallis 结果不可信", sc))
  }
  
  fit <- kruskal.test(sub[[sc]] ~ factor(sub$cluster))
  data.frame(
    signature  = sc,
    statistic  = fit$statistic,
    df         = fit$parameter,
    p_value    = fit$p.value,
    p_adj_BH   = NA_real_   # 后续多重校正
  )
})
kw_res <- do.call(rbind, kw_res)
kw_res$p_adj_BH <- p.adjust(kw_res$p_value, method = "BH")

fwrite(
  kw_res,
  file.path(clust_dir, "tcga_paad_k3_signature_score_kruskal.tsv"),
  sep   = "\t",
  quote = FALSE
)

# ---- Dunn 事后两两比较（P2 修复：新增）---------------------------
dunn_res_list <- lapply(score_cols, function(sc) {
  sub <- dat[!is.na(dat[[sc]]) & !is.na(dat$cluster), ]
  res <- dunn.test(
    x      = sub[[sc]],
    g      = factor(sub$cluster),
    method = "bh",   # Benjamini-Hochberg 校正
    altp   = TRUE    # 双侧 p 值
  )
  data.frame(
    signature    = sc,
    comparison   = res$comparisons,
    statistic    = res$Z,
    p_value      = res$altP,
    p_adj_BH     = res$altP.adjusted
  )
})
dunn_res <- do.call(rbind, dunn_res_list)

fwrite(
  dunn_res,
  file.path(clust_dir, "tcga_paad_k3_signature_score_dunn_posthoc.tsv"),
  sep   = "\t",
  quote = FALSE
)

# ---- cluster 均值汇总（P1 修复：na.rm = TRUE）--------------------
mean_res <- do.call(rbind, lapply(score_cols, function(sc) {
  sub <- dat[!is.na(dat[[sc]]) & !is.na(dat$cluster), ]
  agg <- aggregate(
    sub[[sc]],
    by  = list(cluster = sub$cluster),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  data.frame(
    signature  = sc,
    cluster    = agg$cluster,
    mean_score = agg$x
  )
}))

fwrite(
  mean_res,
  file.path(clust_dir, "tcga_paad_k3_signature_score_cluster_means.tsv"),
  sep   = "\t",
  quote = FALSE
)

# ---- 打印结果 -----------------------------------------------------
message("\n===== Kruskal-Wallis 检验结果 =====")
print(kw_res)

message("\n===== Dunn 事后两两比较 =====")
print(dunn_res)

message("\n===== 各 cluster 均值 =====")
print(mean_res)

message("\nK3 signature score comparison finished.")
