# =========================================================
# 02_download_tcga_paad_full.R
# 作用：
# 1. 下载 TCGA-PAAD RNA-seq STAR-counts
# 2. 下载 TCGA-PAAD clinical
# 3. 下载 TCGA-PAAD mutation MAF
# 4. 生成 expression-clinical mapping
# 5. 生成 sample audit
# 6. 保存日志和 sessionInfo
# =========================================================

project_root <- "/Users/wmz_mac/Desktop/胰腺癌"

# =========================
# 0. 基础设置
# =========================
options(stringsAsFactors = FALSE)

if (!dir.exists(project_root)) {
  stop("项目根目录不存在：", project_root)
}

raw_root <- file.path(project_root, "00_raw_data", "TCGA_PAAD")
expr_dir <- file.path(raw_root, "expression")
clin_dir <- file.path(raw_root, "clinical")
mut_dir  <- file.path(raw_root, "mutation")
meta_dir <- file.path(project_root, "01_metadata")
log_dir  <- file.path(project_root, "11_logs")

dir.create(expr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(clin_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mut_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir,  recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "phase2_tcga_download.log")

log_msg <- function(...) {
  txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(txt, "\n")
  write(txt, file = log_file, append = TRUE)
}

# 将 data.frame 中的 list 列转换成可写出的字符列
flatten_df_list_cols <- function(df) {
  df[] <- lapply(df, function(x) {
    if (is.list(x)) {
      vapply(x, function(elem) {
        if (length(elem) == 0 || all(is.na(elem))) {
          NA_character_
        } else {
          paste(as.character(elem), collapse = ";")
        }
      }, character(1))
    } else {
      x
    }
  })
  df
}

log_msg("START phase2_tcga_download")

# =========================
# 1. 安装 / 加载依赖
# =========================
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_pkgs <- c("data.table")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

bioc_pkgs <- c("TCGAbiolinks", "SummarizedExperiment")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(TCGAbiolinks)
library(SummarizedExperiment)
library(data.table)

# =========================
# 2. 下载 expression: STAR - Counts
# =========================
log_msg("Query TCGA-PAAD RNA-seq STAR-counts ...")

query_exp <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts"
)

log_msg("Expression query completed.")

old_wd <- getwd()
setwd(raw_root)
on.exit(setwd(old_wd), add = TRUE)

log_msg("Downloading expression data ...")
GDCdownload(query_exp, method = "api", files.per.chunk = 20)

log_msg("Preparing expression SummarizedExperiment ...")
se_exp <- GDCprepare(query_exp)

saveRDS(se_exp, file = file.path(expr_dir, "tcga_paad_star_counts_se.rds"))

# 查看 assay 名称
assay_names <- assayNames(se_exp)
writeLines(
  assay_names,
  con = file.path(expr_dir, "tcga_paad_assay_names.txt")
)

log_msg("Assay names detected: ", paste(assay_names, collapse = ", "))

# 优先选择 raw counts，避免优先落到 TPM/FPKM
preferred_assays <- c(
  "unstranded",
  "stranded_first",
  "stranded_second",
  "raw_count",
  "count"
)

matched_assays <- preferred_assays[preferred_assays %in% assay_names]
if (length(matched_assays) > 0) {
  selected_assay <- matched_assays[1]
} else {
  selected_assay <- assay_names[1]
  log_msg("WARNING: preferred assay not found, fallback to first assay: ", selected_assay)
}

expr_mat <- assay(se_exp, selected_assay)
expr_mat <- as.matrix(expr_mat)

row_anno <- as.data.frame(rowData(se_exp), stringsAsFactors = FALSE)
col_anno <- as.data.frame(colData(se_exp), stringsAsFactors = FALSE)

# 修复 list 列，避免 write.table 出错
row_anno <- flatten_df_list_cols(row_anno)
col_anno <- flatten_df_list_cols(col_anno)

# 补充样本和患者 ID
if ("barcode" %in% colnames(col_anno)) {
  col_anno$sample_barcode <- as.character(col_anno$barcode)
} else {
  col_anno$sample_barcode <- rownames(col_anno)
}
col_anno$patient_id <- substr(col_anno$sample_barcode, 1, 12)

# 保存 expression 矩阵
saveRDS(expr_mat, file = file.path(expr_dir, "tcga_paad_counts_matrix_selected_assay.rds"))

write.table(
  expr_mat,
  file = gzfile(file.path(expr_dir, "tcga_paad_counts_matrix.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

write.table(
  row_anno,
  file = file.path(expr_dir, "tcga_paad_gene_annotation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  col_anno,
  file = file.path(clin_dir, "tcga_paad_coldata_from_se.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

log_msg("Expression done. selected_assay=", selected_assay,
        "; genes=", nrow(expr_mat),
        "; samples=", ncol(expr_mat))

# =========================
# 3. 下载 clinical
# =========================
log_msg("Downloading indexed clinical ...")

clinical_raw <- GDCquery_clinic(project = "TCGA-PAAD", type = "clinical")
clinical_raw <- as.data.frame(clinical_raw, stringsAsFactors = FALSE)
clinical_raw <- flatten_df_list_cols(clinical_raw)

write.table(
  clinical_raw,
  file = file.path(clin_dir, "tcga_paad_clinical_raw.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(
  colnames(clinical_raw),
  con = file.path(clin_dir, "tcga_paad_clinical_columns.txt")
)

log_msg("Clinical done. rows=", nrow(clinical_raw), "; cols=", ncol(clinical_raw))

# =========================
# 4. expression-clinical mapping
# =========================
log_msg("Building sample-clinical mapping ...")

clinical_id_col <- NULL
candidate_id_cols <- c("submitter_id", "bcr_patient_barcode", "case_submitter_id")

for (cc in candidate_id_cols) {
  if (cc %in% colnames(clinical_raw)) {
    clinical_id_col <- cc
    break
  }
}

if (is.null(clinical_id_col)) {
  stop("clinical_raw 中没有找到可用于 patient 对应的列：submitter_id / bcr_patient_barcode / case_submitter_id")
}

# 增加匹配标记，便于准确统计映射情况
clinical_raw$clinical_matched_flag <- TRUE

sample_clinical_map <- merge(
  col_anno,
  clinical_raw,
  by.x = "patient_id",
  by.y = clinical_id_col,
  all.x = TRUE
)

write.table(
  sample_clinical_map,
  file = file.path(clin_dir, "tcga_paad_sample_clinical_map.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

log_msg("Sample-clinical mapping done. rows=", nrow(sample_clinical_map))

# =========================
# 5. 下载 mutation MAF
# =========================
log_msg("Query mutation MAF ...")

maf_success <- TRUE
maf_error_msg <- NA_character_

tryCatch({
  
  query_maf <- GDCquery(
    project = "TCGA-PAAD",
    data.category = "Simple Nucleotide Variation",
    access = "open",
    data.type = "Masked Somatic Mutation",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  )
  
  log_msg("Downloading mutation MAF ...")
  GDCdownload(query_maf, method = "api", files.per.chunk = 20)
  
  log_msg("Preparing mutation MAF ...")
  maf_raw <- GDCprepare(query_maf)
  maf_df <- as.data.frame(maf_raw, stringsAsFactors = FALSE)
  maf_df <- flatten_df_list_cols(maf_df)
  
  write.table(
    maf_df,
    file = gzfile(file.path(mut_dir, "tcga_paad_maf_raw.tsv.gz")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  log_msg("Mutation done. rows=", nrow(maf_df), "; cols=", ncol(maf_df))
  
}, error = function(e) {
  maf_success <<- FALSE
  maf_error_msg <<- conditionMessage(e)
  log_msg("WARNING: MAF step failed: ", maf_error_msg)
})

# =========================
# 6. 生成 sample audit
# =========================
log_msg("Generating sample audit ...")

sample_type_col <- NULL
candidate_sample_type_cols <- c(
  "shortLetterCode",
  "short_letter_code",
  "sample_type",
  "definition"
)

for (cc in candidate_sample_type_cols) {
  if (cc %in% colnames(col_anno)) {
    sample_type_col <- cc
    break
  }
}

if (is.null(sample_type_col)) {
  col_anno$sample_type_final <- NA_character_
} else {
  col_anno$sample_type_final <- as.character(col_anno[[sample_type_col]])
}

n_expr_samples      <- ncol(expr_mat)
n_unique_patients   <- length(unique(col_anno$patient_id))
n_primary_tumor     <- sum(col_anno$sample_type_final %in% c("TP", "Primary Tumor"), na.rm = TRUE)
n_solid_normal      <- sum(col_anno$sample_type_final %in% c("NT", "Solid Tissue Normal"), na.rm = TRUE)
n_clinical_patients <- length(unique(clinical_raw[[clinical_id_col]]))
n_mapped_samples    <- sum(!is.na(sample_clinical_map$clinical_matched_flag))
n_unmapped_samples  <- sum(is.na(sample_clinical_map$clinical_matched_flag))

sample_type_table <- as.data.frame(table(col_anno$sample_type_final), stringsAsFactors = FALSE)
colnames(sample_type_table) <- c("sample_type_final", "n_samples")

write.table(
  sample_type_table,
  file = file.path(meta_dir, "tcga_paad_sample_type_table.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sample_audit <- data.frame(
  metric = c(
    "selected_expression_assay",
    "n_expression_samples",
    "n_unique_patients_expression",
    "n_primary_tumor",
    "n_solid_tissue_normal",
    "n_clinical_patients",
    "n_rows_in_sample_clinical_map",
    "n_mapped_samples",
    "n_unmapped_samples",
    "maf_download_success"
  ),
  value = c(
    selected_assay,
    n_expr_samples,
    n_unique_patients,
    n_primary_tumor,
    n_solid_normal,
    n_clinical_patients,
    nrow(sample_clinical_map),
    n_mapped_samples,
    n_unmapped_samples,
    maf_success
  ),
  stringsAsFactors = FALSE
)

write.table(
  sample_audit,
  file = file.path(meta_dir, "tcga_paad_sample_audit.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 保存前几行预览，便于快速核查
fwrite(
  as.data.table(head(clinical_raw, 5)),
  file = file.path(clin_dir, "tcga_paad_clinical_head5.tsv"),
  sep = "\t"
)

fwrite(
  as.data.table(head(col_anno, 5)),
  file = file.path(clin_dir, "tcga_paad_coldata_head5.tsv"),
  sep = "\t"
)

log_msg("Sample audit saved.")

# =========================
# 7. 保存 sessionInfo
# =========================
capture.output(
  sessionInfo(),
  file = file.path(log_dir, "phase2_sessionInfo.txt")
)

# =========================
# 8. 完成提示
# =========================
log_msg("END phase2_tcga_download")
cat("\nTCGA-PAAD download and audit finished.\n")
cat("Key files generated:\n")
cat("1) 00_raw_data/TCGA_PAAD/expression/tcga_paad_star_counts_se.rds\n")
cat("2) 00_raw_data/TCGA_PAAD/expression/tcga_paad_counts_matrix.tsv.gz\n")
cat("3) 00_raw_data/TCGA_PAAD/clinical/tcga_paad_clinical_raw.tsv\n")
cat("4) 00_raw_data/TCGA_PAAD/clinical/tcga_paad_sample_clinical_map.tsv\n")
cat("5) 01_metadata/tcga_paad_sample_audit.tsv\n")
cat("6) 11_logs/phase2_tcga_download.log\n")

if (!maf_success) {
  cat("\nNote: MAF step failed, but expression + clinical should still be usable.\n")
  cat("MAF error:\n", maf_error_msg, "\n")
}