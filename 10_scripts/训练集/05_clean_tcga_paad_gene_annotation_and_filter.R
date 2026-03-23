# =========================================================
# 05_clean_tcga_paad_gene_annotation_and_filter.R
# =========================================================

project_root <- "/Users/wmz_mac/Desktop/胰腺癌"
options(stringsAsFactors = FALSE)
library(data.table)

expr_dir <- file.path(project_root, "00_raw_data",       "TCGA_PAAD", "expression")
proc_dir <- file.path(project_root, "02_processed_data", "TCGA_PAAD")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

gene_anno_file <- file.path(expr_dir, "tcga_paad_gene_annotation.tsv")
expr_all_file  <- file.path(proc_dir, "tcga_paad_expr_all_tp_counts.rds")
expr_str_file  <- file.path(proc_dir, "tcga_paad_expr_strict_pdac_counts.rds")

expr_all  <- as.matrix(readRDS(expr_all_file))
expr_str  <- as.matrix(readRDS(expr_str_file))
gene_anno <- fread(gene_anno_file, data.table = FALSE)

# [Fix 3] 确保两矩阵行名完全一致，否则后续 valid_idx 会错位
if (!identical(rownames(expr_all), rownames(expr_str))) {
  stop("rownames of expr_all and expr_str differ – check upstream pipeline.")
}

# -------------------------
# 1. 清洗 Ensembl ID
# -------------------------
raw_ids       <- rownames(expr_all)
ensembl_clean <- sub("\\..*$", "", raw_ids)

candidate_gene_id_cols <- c("gene_id", "ensembl_gene_id", "Geneid")
candidate_symbol_cols  <- c("gene_name", "external_gene_name", "gene_symbol", "Symbol")
candidate_type_cols    <- c("gene_type", "gene_biotype", "gene_type.x")

gene_id_col <- candidate_gene_id_cols[candidate_gene_id_cols %in% colnames(gene_anno)][1]
symbol_col  <- candidate_symbol_cols [candidate_symbol_cols  %in% colnames(gene_anno)][1]
type_col    <- candidate_type_cols   [candidate_type_cols    %in% colnames(gene_anno)][1]

if (is.na(gene_id_col)) stop("gene annotation 中未找到 gene_id 列")
if (is.na(symbol_col))  stop("gene annotation 中未找到 gene symbol 列")

gene_anno$ensembl_clean <- sub("\\..*$", "", gene_anno[[gene_id_col]])

# [Fix 1] type_col 缺失时优雅降级，填充 NA 而非报错
if (!is.na(type_col)) {
  anno_use <- gene_anno[, c("ensembl_clean", symbol_col, type_col), drop = FALSE]
  colnames(anno_use) <- c("ensembl_clean", "gene_symbol", "gene_type")
} else {
  message("Warning: gene_type column not found – filling with NA")
  anno_use <- gene_anno[, c("ensembl_clean", symbol_col), drop = FALSE]
  colnames(anno_use) <- c("ensembl_clean", "gene_symbol")
  anno_use$gene_type <- NA_character_
}

# [Fix 2] 按 ensembl_clean 去重，防止 merge 产生一对多重复行
anno_use <- anno_use[!duplicated(anno_use$ensembl_clean), ]

# -------------------------
# 2. 构建 annotation table
# -------------------------
annot_df <- data.frame(
  raw_ensembl   = raw_ids,
  ensembl_clean = ensembl_clean,
  stringsAsFactors = FALSE
)
annot_df <- merge(annot_df, anno_use, by = "ensembl_clean",
                  all.x = TRUE, sort = FALSE)
# 恢复原始行顺序
annot_df <- annot_df[match(raw_ids, annot_df$raw_ensembl), ]

# -------------------------
# 3. 去除无效 symbol
# -------------------------
valid_idx <- !is.na(annot_df$gene_symbol) &
  nzchar(annot_df$gene_symbol) &        # 替代 != ""，更严格
  annot_df$gene_symbol != "NA"

expr_all_valid <- expr_all[valid_idx, , drop = FALSE]
expr_str_valid <- expr_str[valid_idx, , drop = FALSE]
annot_valid    <- annot_df[valid_idx, ]

# -------------------------
# 4. 去重：保留每个 symbol 均值表达最高的行
# -------------------------
annot_valid$mean_expr_all <- rowMeans(expr_all_valid)
annot_valid$row_index     <- seq_len(nrow(annot_valid))
annot_valid <- annot_valid[order(annot_valid$gene_symbol,
                                 -annot_valid$mean_expr_all), ]
keep_idx    <- !duplicated(annot_valid$gene_symbol)
annot_dedup <- annot_valid[keep_idx, ]

expr_all_dedup <- expr_all_valid[annot_dedup$row_index, , drop = FALSE]
expr_str_dedup <- expr_str_valid[annot_dedup$row_index, , drop = FALSE]
rownames(expr_all_dedup) <- annot_dedup$gene_symbol
rownames(expr_str_dedup) <- annot_dedup$gene_symbol

# [Fix 4] 清理辅助列，避免污染导出的 annotation 文件
annot_dedup <- annot_dedup[, setdiff(colnames(annot_dedup),
                                     c("mean_expr_all", "row_index")),
                           drop = FALSE]
rownames(annot_dedup) <- NULL

# -------------------------
# 5. 低表达过滤
# 规则：count >= 10 的样本数至少占 20%
# -------------------------
min_samples_all <- ceiling(ncol(expr_all_dedup) * 0.20)
min_samples_str <- ceiling(ncol(expr_str_dedup) * 0.20)

keep_gene_all <- rowSums(expr_all_dedup >= 10) >= min_samples_all
keep_gene_str <- rowSums(expr_str_dedup >= 10) >= min_samples_str

expr_all_filt <- expr_all_dedup[keep_gene_all, , drop = FALSE]
expr_str_filt <- expr_str_dedup[keep_gene_str, , drop = FALSE]

annot_all_filt <- annot_dedup[match(rownames(expr_all_filt), annot_dedup$gene_symbol), ]
annot_str_filt <- annot_dedup[match(rownames(expr_str_filt), annot_dedup$gene_symbol), ]
rownames(annot_all_filt) <- NULL
rownames(annot_str_filt) <- NULL

# -------------------------
# 6. 导出
# -------------------------
saveRDS(expr_all_dedup, file.path(proc_dir, "tcga_paad_expr_all_tp_symbol_counts_unfiltered.rds"))
saveRDS(expr_str_dedup, file.path(proc_dir, "tcga_paad_expr_strict_pdac_symbol_counts_unfiltered.rds"))
saveRDS(expr_all_filt,  file.path(proc_dir, "tcga_paad_expr_all_tp_symbol_counts_filtered.rds"))
saveRDS(expr_str_filt,  file.path(proc_dir, "tcga_paad_expr_strict_pdac_symbol_counts_filtered.rds"))

write.table(annot_dedup,
            file.path(proc_dir, "tcga_paad_gene_annotation_deduplicated.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(annot_all_filt,
            file.path(proc_dir, "tcga_paad_gene_annotation_all_tp_filtered.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(annot_str_filt,
            file.path(proc_dir, "tcga_paad_gene_annotation_strict_pdac_filtered.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# [Fix 5] audit 用 sum(valid_idx) 而非已 subset 的中间变量
audit_df <- data.frame(
  metric = c(
    "raw_gene_rows",
    "valid_symbol_rows",
    "deduplicated_symbol_rows",
    "filtered_all_tp_rows",
    "filtered_strict_pdac_rows",
    "all_tp_samples",
    "strict_pdac_samples"
  ),
  value = c(
    nrow(expr_all),
    sum(valid_idx),
    nrow(expr_all_dedup),
    nrow(expr_all_filt),
    nrow(expr_str_filt),
    ncol(expr_all_filt),
    ncol(expr_str_filt)
  )
)
write.table(audit_df,
            file.path(proc_dir, "tcga_paad_gene_cleaning_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Gene cleaning finished.\n")
cat("All-TP filtered dim:",      dim(expr_all_filt), "\n")
cat("Strict-PDAC filtered dim:", dim(expr_str_filt), "\n")
