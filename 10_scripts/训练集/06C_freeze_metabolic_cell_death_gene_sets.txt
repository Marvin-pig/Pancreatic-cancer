# =========================================================
# 06C_freeze_metabolic_cell_death_gene_sets.R
# 作用：
# 1. 读取 starter gene set master
# 2. 标准化字段并严格校验
# 3. 映射到 TCGA all-TP / strict-PDAC filtered expression
# 4. 输出 analysis-ready frozen gene set
# 5. 输出 duplicate resolution / audit / log
# =========================================================

project_root <- Sys.getenv("PROJECT_ROOT", unset = "/Users/wmz_mac/Desktop/胰腺癌")
gene_dir <- file.path(project_root, "03_gene_sets")
proc_dir <- file.path(project_root, "02_processed_data", "TCGA_PAAD")
log_dir  <- file.path(project_root, "11_logs")
dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir,  recursive = TRUE, showWarnings = FALSE)

master_file   <- file.path(gene_dir, "metabolic_cell_death_gene_sets_master_starter.tsv")
expr_all_file <- file.path(proc_dir, "tcga_paad_expr_all_tp_symbol_counts_filtered.rds")
expr_str_file <- file.path(proc_dir, "tcga_paad_expr_strict_pdac_symbol_counts_filtered.rds")

library(data.table)

# ---------- input check ----------
input_files <- c(master_file, expr_all_file, expr_str_file)
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs) > 0) {
  stop("以下输入文件不存在：\n  ", paste(missing_inputs, collapse = "\n  "))
}

# ---------- read ----------
gene_df  <- fread(master_file, data.table = FALSE)
expr_all <- readRDS(expr_all_file)
expr_str <- readRDS(expr_str_file)

# ---------- required columns ----------
required_cols <- c("gene_symbol", "death_type", "include_main", "include_extended")
missing_cols <- setdiff(required_cols, colnames(gene_df))
if (length(missing_cols) > 0) {
  stop("master 文件缺少必需列：", paste(missing_cols, collapse = ", "))
}

# ---------- validate expression rownames ----------
if (is.null(rownames(expr_all)) || is.null(rownames(expr_str))) {
  stop("表达矩阵缺少 rownames(gene symbols)。")
}
if (anyDuplicated(rownames(expr_all))) stop("expr_all rownames 存在重复 gene symbol。")
if (anyDuplicated(rownames(expr_str))) stop("expr_str rownames 存在重复 gene symbol。")
if (any(is.na(rownames(expr_all)) | trimws(rownames(expr_all)) == "")) {
  stop("expr_all rownames 存在空 gene symbol。")
}
if (any(is.na(rownames(expr_str)) | trimws(rownames(expr_str)) == "")) {
  stop("expr_str rownames 存在空 gene symbol。")
}

# ---------- standardize ----------
gene_df$gene_symbol <- toupper(trimws(gene_df$gene_symbol))
gene_df$death_type  <- tolower(trimws(gene_df$death_type))

valid_death_types <- c("ferroptosis", "cuproptosis", "disulfidptosis")
if (!all(gene_df$death_type %in% valid_death_types)) {
  bad_types <- unique(gene_df$death_type[!gene_df$death_type %in% valid_death_types])
  stop("death_type 存在非法值：", paste(bad_types, collapse = ", "))
}

gene_df$include_main     <- as.integer(trimws(as.character(gene_df$include_main)))
gene_df$include_extended <- as.integer(trimws(as.character(gene_df$include_extended)))

if (any(is.na(gene_df$include_main)))     stop("include_main 存在非法值或缺失。")
if (any(is.na(gene_df$include_extended))) stop("include_extended 存在非法值或缺失。")
if (!all(gene_df$include_main     %in% c(0L, 1L))) stop("include_main 必须为 0/1。")
if (!all(gene_df$include_extended %in% c(0L, 1L))) stop("include_extended 必须为 0/1。")
if (any(is.na(gene_df$gene_symbol) | gene_df$gene_symbol == "")) {
  stop("gene_symbol 存在空值。")
}

# ---------- 校验 include_main 与 include_extended 不应同时为 1 ----------
conflict_rows <- gene_df$include_main == 1L & gene_df$include_extended == 1L
if (any(conflict_rows)) {
  warning("以下基因同时标记 include_main=1 和 include_extended=1，将归为 main：",
          paste(unique(gene_df$gene_symbol[conflict_rows]), collapse = ", "))
}

# ---------- map ----------
expr_all_genes <- toupper(rownames(expr_all))
expr_str_genes <- toupper(rownames(expr_str))

gene_df$in_all_tp_filtered      <- gene_df$gene_symbol %in% expr_all_genes
gene_df$in_strict_pdac_filtered <- gene_df$gene_symbol %in% expr_str_genes
gene_df$in_both_cohorts         <- gene_df$in_all_tp_filtered & gene_df$in_strict_pdac_filtered

# ---------- tier ----------
gene_df$tier <- ifelse(gene_df$include_main == 1L, "main",
                       ifelse(gene_df$include_extended == 1L, "extended", "exclude"))

# ---------- duplicate resolution ----------
dup_tab <- as.data.table(gene_df)[
  , .(
    n_records      = .N,
    death_type_list = paste(sort(unique(death_type)), collapse = ";"),
    tier_list      = paste(sort(unique(tier)), collapse = ";")
  ),
  by = gene_symbol
]
dup_tab[, duplicate_flag   := n_records > 1]
dup_tab[, duplicate_policy := ifelse(
  duplicate_flag,
  "keep_multi_label_for_gene_set; deduplicate_before_model_if_needed",
  "unique"
)]

gene_df <- merge(gene_df, dup_tab[, .(gene_symbol, duplicate_flag, duplicate_policy)],
                 by = "gene_symbol", all.x = TRUE, sort = FALSE)

# ---------- analysis-ready outputs ----------
analysis_ready <- gene_df[gene_df$tier %in% c("main", "extended"), , drop = FALSE]
analysis_ready <- analysis_ready[
  order(analysis_ready$death_type, analysis_ready$tier, analysis_ready$gene_symbol), ]

main_detected_alltp <- analysis_ready[
  analysis_ready$tier == "main"     & analysis_ready$in_all_tp_filtered, , drop = FALSE]
ext_detected_alltp  <- analysis_ready[
  analysis_ready$tier == "extended" & analysis_ready$in_all_tp_filtered, , drop = FALSE]
main_detected_both  <- analysis_ready[
  analysis_ready$tier == "main"     & analysis_ready$in_both_cohorts, , drop = FALSE]
ext_detected_both   <- analysis_ready[
  analysis_ready$tier == "extended" & analysis_ready$in_both_cohorts, , drop = FALSE]

# 修复：仅报告 analysis-ready 集合中未在 all-TP 中检测到的基因（排除 exclude 层）
missing_alltp <- analysis_ready[!analysis_ready$in_all_tp_filtered, , drop = FALSE]

# ---------- write ----------
fwrite(
  gene_df,
  file.path(gene_dir, "metabolic_cell_death_gene_sets_master_starter_checked.tsv"),
  sep = "\t", quote = FALSE
)
fwrite(
  dup_tab[order(-as.integer(duplicate_flag), gene_symbol), ],
  file.path(gene_dir, "metabolic_cell_death_gene_set_duplicate_resolution.tsv"),
  sep = "\t", quote = FALSE
)
fwrite(
  analysis_ready,
  file.path(proc_dir, "tcga_paad_metabolic_cell_death_gene_set_analysis_ready.tsv"),
  sep = "\t", quote = FALSE
)
fwrite(
  main_detected_alltp,
  file.path(proc_dir, "tcga_paad_metabolic_cell_death_main_detected_all_tp.tsv"),
  sep = "\t", quote = FALSE
)
fwrite(
  ext_detected_alltp,
  file.path(proc_dir, "tcga_paad_metabolic_cell_death_extended_detected_all_tp.tsv"),
  sep = "\t", quote = FALSE
)
fwrite(
  main_detected_both,
  file.path(proc_dir, "tcga_paad_metabolic_cell_death_main_detected_both_cohorts.tsv"),
  sep = "\t", quote = FALSE
)
fwrite(
  ext_detected_both,
  file.path(proc_dir, "tcga_paad_metabolic_cell_death_extended_detected_both_cohorts.tsv"),
  sep = "\t", quote = FALSE
)
fwrite(
  missing_alltp,
  file.path(proc_dir, "tcga_paad_metabolic_cell_death_missing_in_all_tp.tsv"),
  sep = "\t", quote = FALSE
)

# ---------- audit ----------
audit_by_type <- as.data.table(gene_df)[
  , .(
    n_input                    = .N,
    n_detected_all_tp          = sum(in_all_tp_filtered),
    n_detected_strict_pdac     = sum(in_strict_pdac_filtered),
    n_detected_both            = sum(in_both_cohorts),
    n_main                     = sum(tier == "main"),
    n_extended                 = sum(tier == "extended"),
    n_main_detected_all_tp     = sum(tier == "main"     & in_all_tp_filtered),
    n_extended_detected_all_tp = sum(tier == "extended" & in_all_tp_filtered),
    # 修复：计算含重复记录的唯一基因数，而非重复记录总行数
    n_duplicate_genes          = uniqueN(gene_symbol[duplicate_flag])
  ),
  by = death_type
]

fwrite(
  audit_by_type,
  file.path(proc_dir, "tcga_paad_metabolic_cell_death_gene_set_mapping_audit_frozen.tsv"),
  sep = "\t", quote = FALSE
)

# ---------- log ----------
log_lines <- c(
  paste0("date\t",                    as.character(Sys.time())),
  paste0("input_master\t",            master_file),
  paste0("input_expr_all\t",          expr_all_file),
  paste0("input_expr_strict\t",       expr_str_file),
  paste0("n_input\t",                 nrow(gene_df)),
  paste0("n_analysis_ready\t",        nrow(analysis_ready)),
  paste0("n_main_detected_all_tp\t",  nrow(main_detected_alltp)),
  paste0("n_extended_detected_all_tp\t", nrow(ext_detected_alltp)),
  paste0("n_main_detected_both\t",    nrow(main_detected_both)),
  paste0("n_extended_detected_both\t", nrow(ext_detected_both)),
  paste0("n_missing_all_tp\t",        nrow(missing_alltp)),
  # 修复：与 audit 保持一致，记录唯一重复基因数
  paste0("n_duplicate_genes\t",       sum(dup_tab$duplicate_flag)),
  "",
  "--- sessionInfo ---",
  capture.output(sessionInfo())
)

writeLines(
  log_lines,
  con = file.path(log_dir, "phase3_gene_set_freeze.log")
)

message("Gene set freeze completed successfully.")
print(audit_by_type)
