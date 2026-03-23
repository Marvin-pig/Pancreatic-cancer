# =========================================================
# 07A_prepare_main_gene_matrix_for_clustering.R
# 作用：
# 1. 读取 all-TP 表达矩阵
# 2. 读取 analysis-ready gene set
# 3. 提取 main gene set
# 4. 按 gene_symbol 去重，生成 clustering 输入矩阵
# 5. 输出 audit
# =========================================================

project_root <- Sys.getenv("PROJECT_ROOT", unset = "/Users/wmz_mac/Desktop/胰腺癌")
proc_dir  <- file.path(project_root, "02_processed_data", "TCGA_PAAD")
gene_dir  <- file.path(project_root, "03_gene_sets")
clust_dir <- file.path(project_root, "04_bulk_analysis", "02_clustering")
dir.create(clust_dir, recursive = TRUE, showWarnings = FALSE)

expr_file <- file.path(proc_dir, "tcga_paad_expr_all_tp_symbol_counts_filtered.rds")
gene_file <- file.path(proc_dir, "tcga_paad_metabolic_cell_death_gene_set_analysis_ready.tsv")
# BUG FIX #2: dup_file 若暂不参与主流程，先移除加载，避免逻辑歧义
# dup_file  <- file.path(gene_dir, "metabolic_cell_death_gene_set_duplicate_resolution.tsv")

library(data.table)

expr    <- readRDS(expr_file)
gene_df <- fread(gene_file, data.table = FALSE)

# ---- basic check ----
# BUG FIX #4: 明确检查 expr 类型与 rownames
if (!is.matrix(expr) && !is.data.frame(expr)) {
  stop("expr 对象类型异常，应为 matrix 或 data.frame，实际为：", class(expr))
}
stopifnot(!is.null(rownames(expr)))
stopifnot(!anyDuplicated(rownames(expr)))
stopifnot(!anyDuplicated(colnames(expr)))

# BUG FIX #3: 验证必要列存在
required_cols <- c("gene_symbol", "tier", "death_type", "subclass", "evidence_tier")
missing_cols  <- setdiff(required_cols, colnames(gene_df))
if (length(missing_cols) > 0) {
  stop("gene_df 缺少以下必要列：", paste(missing_cols, collapse = ", "))
}

# ---- select main genes ----
main_df <- gene_df[gene_df$tier == "main", , drop = FALSE]
if (nrow(main_df) == 0) {
  stop("未找到 main gene set，请检查 gene_file 中 tier 列是否含 'main'。")
}

# ---- duplicate handling for clustering ----
# 对 clustering 输入，gene_symbol 必须唯一
# 若同一 gene_symbol 对应多条记录（如 SLC7A11），只保留首条
main_unique_symbols <- unique(main_df$gene_symbol)

missing_main <- setdiff(main_unique_symbols, rownames(expr))
if (length(missing_main) > 0) {
  stop(
    "以下 main genes 不在表达矩阵 rownames 中（共 ", length(missing_main), " 个）：\n",
    paste(missing_main, collapse = ", ")
  )
}

expr_main <- expr[main_unique_symbols, , drop = FALSE]

# ---- reorder symbols by death_type then gene_symbol for display convenience ----
main_unique_map <- main_df[
  !duplicated(main_df$gene_symbol),
  c("gene_symbol", "death_type", "subclass", "evidence_tier"),
  drop = FALSE
]
main_unique_map <- main_unique_map[
  order(main_unique_map$death_type, main_unique_map$gene_symbol), ,
  drop = FALSE
]
expr_main <- expr_main[main_unique_map$gene_symbol, , drop = FALSE]

# ---- audit ----
audit <- data.frame(
  metric = c(
    "n_main_records",
    "n_main_unique_symbols",
    "n_expr_rows",
    "n_expr_cols",
    "n_duplicate_main_records_removed"
  ),
  value = c(
    nrow(main_df),
    length(main_unique_symbols),
    nrow(expr_main),
    ncol(expr_main),
    nrow(main_df) - length(main_unique_symbols)
  )
)

# ---- write outputs ----
fwrite(
  main_unique_map,
  file.path(clust_dir, "tcga_paad_main_gene_set_unique_symbols.tsv"),
  sep = "\t", quote = FALSE
)

saveRDS(
  expr_main,
  file.path(clust_dir, "tcga_paad_main_gene_expr_for_clustering.rds")
)

# BUG FIX #1: fwrite 不支持 row.names，需手动将 rownames 转为列再写出
expr_main_tsv <- as.data.frame(expr_main, check.names = FALSE)
expr_main_tsv <- cbind(gene_symbol = rownames(expr_main_tsv), expr_main_tsv)
fwrite(
  expr_main_tsv,
  file.path(clust_dir, "tcga_paad_main_gene_expr_for_clustering.tsv"),
  sep = "\t", quote = FALSE
)

fwrite(
  audit,
  file.path(clust_dir, "tcga_paad_clustering_input_audit.tsv"),
  sep = "\t", quote = FALSE
)

message("Main gene matrix for clustering prepared successfully.")
print(audit)
