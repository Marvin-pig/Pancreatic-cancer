# =========================================================
# 06B_map_starter_gene_sets_to_tcga.R
# 作用：
# 1. 读取 starter gene set master
# 2. 与 TCGA filtered expression matrix 做交集
# 3. 输出 main / extended 两套可分析 gene set
# =========================================================

# ---------- 路径配置（优先读取环境变量，方便跨机器运行）----------
project_root <- Sys.getenv("PROJECT_ROOT", unset = "/Users/wmz_mac/Desktop/胰腺癌")
gene_dir     <- file.path(project_root, "03_gene_sets")
proc_dir     <- file.path(project_root, "02_processed_data", "TCGA_PAAD")

master_file    <- file.path(gene_dir,  "metabolic_cell_death_gene_sets_master_starter.tsv")
expr_all_file  <- file.path(proc_dir,  "tcga_paad_expr_all_tp_symbol_counts_filtered.rds")
expr_str_file  <- file.path(proc_dir,  "tcga_paad_expr_strict_pdac_symbol_counts_filtered.rds")

# ---------- 依赖 ----------
library(data.table)

# ---------- 输入文件校验 ----------
input_files <- c(master_file, expr_all_file, expr_str_file)
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs) > 0) {
  stop("以下输入文件不存在，请检查路径：\n  ",
       paste(missing_inputs, collapse = "\n  "))
}

# 确保输出目录存在
dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 读入数据 ----------
message("[1/4] 读取 master gene set ...")
gene_df <- fread(master_file, data.table = FALSE)

message("[2/4] 读取表达矩阵 ...")
expr_all <- readRDS(expr_all_file)
expr_str <- readRDS(expr_str_file)

# ---------- 列名校验 ----------
required_cols <- c("gene_symbol", "death_type", "include_main", "include_extended")
missing_cols  <- setdiff(required_cols, colnames(gene_df))
if (length(missing_cols) > 0) {
  stop("master 文件缺少必需列：", paste(missing_cols, collapse = ", "))
}

# ---------- 基因符号规范化 ----------
gene_df$gene_symbol <- toupper(trimws(gene_df$gene_symbol))
expr_all_genes       <- toupper(rownames(expr_all))
expr_str_genes       <- toupper(rownames(expr_str))

# ---------- 重复基因符号检查 ----------
dup_genes <- gene_df$gene_symbol[duplicated(gene_df$gene_symbol)]
if (length(dup_genes) > 0) {
  warning("发现重复 gene_symbol，请确认是否有意为之：\n  ",
          paste(unique(dup_genes), collapse = ", "))
}

# ---------- include_main / include_extended 类型统一 ----------
gene_df$include_main     <- as.integer(gene_df$include_main)
gene_df$include_extended <- as.integer(gene_df$include_extended)

# ---------- 映射标记 ----------
gene_df$in_all_tp_filtered    <- gene_df$gene_symbol %in% expr_all_genes
gene_df$in_strict_pdac_filtered <- gene_df$gene_symbol %in% expr_str_genes

# ---------- 子集过滤 ----------
gene_df_detected <- gene_df[gene_df$in_all_tp_filtered, , drop = FALSE]
main_set         <- gene_df_detected[gene_df_detected$include_main     == 1L, , drop = FALSE]
extended_set     <- gene_df_detected[gene_df_detected$include_extended == 1L, , drop = FALSE]
missing_df       <- gene_df[!gene_df$in_all_tp_filtered, , drop = FALSE]

# ---------- 写出结果 ----------
message("[3/4] 写出结果文件 ...")

out_files <- list(
  mapped   = list(df = gene_df,      name = "metabolic_cell_death_gene_sets_master_starter_mapped.tsv"),
  main     = list(df = main_set,     name = "metabolic_cell_death_gene_sets_main_detected.tsv"),
  extended = list(df = extended_set, name = "metabolic_cell_death_gene_sets_extended_detected.tsv"),
  missing  = list(df = missing_df,   name = "metabolic_cell_death_gene_sets_missing_in_tcga.tsv")
)

for (item in out_files) {
  fwrite(item$df,
         file  = file.path(gene_dir, item$name),
         sep   = "\t",
         quote = FALSE)
}

# ---------- 分类型审计 ----------
message("[4/4] 生成审计表 ...")

audit_by_type <- do.call(rbind, lapply(split(gene_df, gene_df$death_type), function(x) {
  data.frame(
    death_type              = unique(x$death_type),
    n_input                 = nrow(x),
    n_detected_all_tp       = sum(x$in_all_tp_filtered),
    n_detected_strict_pdac  = sum(x$in_strict_pdac_filtered),
    n_main_detected         = sum(x$in_all_tp_filtered & (x$include_main     == 1L)),
    n_extended_detected     = sum(x$in_all_tp_filtered & (x$include_extended == 1L))
  )
}))

fwrite(audit_by_type,
       file  = file.path(gene_dir, "metabolic_cell_death_gene_sets_mapping_audit.tsv"),
       sep   = "\t",
       quote = FALSE)

# ---------- 汇总输出 ----------
message("\n========== Gene Set Mapping Summary ==========")
message(sprintf("  Input genes       : %d", nrow(gene_df)))
message(sprintf("  Detected (all-TP) : %d (%.1f%%)",
                sum(gene_df$in_all_tp_filtered),
                100 * mean(gene_df$in_all_tp_filtered)))
message(sprintf("  Main set          : %d", nrow(main_set)))
message(sprintf("  Extended set      : %d", nrow(extended_set)))
message(sprintf("  Missing           : %d", nrow(missing_df)))
message("==============================================\n")

print(audit_by_type)
message("Gene set mapping finished.")
