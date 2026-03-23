# =========================
# 24_build_gse57495_expression_and_clinical_tables_fix.R
# 修正版：
# 1. 修复 GPL / pheno 表重复列名导致 transmute 报错
# 2. 下载 GPL15048 注释
# 3. 建立 probe -> gene symbol 映射
# 4. 构建 gene-level expression 和 survival-ready clinical table
# =========================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ---------- 路径 ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")
raw_dir <- "04_bulk_analysis/04_external_validation/01_raw_data/GEO/GSE57495"
out_dir <- "04_bulk_analysis/04_external_validation/02_processed_data/GSE57495"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr_file  <- file.path(raw_dir, "23_GSE57495_expression_probe_level.tsv")
pheno_file <- file.path(raw_dir, "23_GSE57495_pheno_raw.tsv")

# ---------- 读取 ----------
expr_df  <- fread(expr_file)
pheno_df <- fread(pheno_file)

# FIX 1: 防止重复列名导致 dplyr 失败
colnames(pheno_df) <- make.unique(colnames(pheno_df))

stopifnot("probe_id" %in% colnames(expr_df))
stopifnot("geo_accession" %in% colnames(pheno_df))

cat("expr_df dim:", dim(expr_df), "\n")
cat("pheno_df dim:", dim(pheno_df), "\n")
cat("Any duplicated pheno colnames:", anyDuplicated(colnames(pheno_df)), "\n")

# ---------- 下载 GPL15048 ----------
cat("Downloading GPL15048 annotation...\n")
gpl <- getGEO("GPL15048", AnnotGPL = TRUE, destdir = raw_dir)
saveRDS(gpl, file.path(raw_dir, "24_GPL15048_gpl.rds"))

gpl_tab <- Table(gpl)
gpl_tab <- as.data.frame(gpl_tab, check.names = FALSE, stringsAsFactors = FALSE)

# FIX 2: 防止 GPL 注释表重复列名导致 dplyr 失败
colnames(gpl_tab) <- make.unique(colnames(gpl_tab))

fwrite(
  gpl_tab,
  file.path(out_dir, "24_GSE57495_GPL15048_annotation.tsv"),
  sep = "\t"
)

cat("GPL15048 annotation dim:", dim(gpl_tab), "\n")
cat("Any duplicated GPL colnames:", anyDuplicated(colnames(gpl_tab)), "\n")

# ---------- 自动识别 probe 与 gene symbol 列 ----------
probe_col_candidates <- c("ID", "ID_REF", "ProbeName", "probe_id", "SPOT_ID")
gene_col_candidates  <- c(
  "Gene Symbol", "GENE_SYMBOL", "Gene symbol", "Symbol",
  "GENE", "GeneSymbol", "gene_assignment", "Gene Assignment"
)

probe_col <- intersect(probe_col_candidates, colnames(gpl_tab))
gene_col  <- intersect(gene_col_candidates, colnames(gpl_tab))

if (length(probe_col) == 0) {
  stop("GPL annotation: no probe ID column found.")
}
if (length(gene_col) == 0) {
  stop(
    paste0(
      "GPL annotation: no gene symbol column found.\nAvailable columns:\n",
      paste(colnames(gpl_tab), collapse = ", ")
    )
  )
}

probe_col <- probe_col[1]
gene_col  <- gene_col[1]

cat("Using probe column:", probe_col, "\n")
cat("Using gene column :", gene_col, "\n")

# ---------- 建立 probe -> gene 映射 ----------
map_df <- data.frame(
  probe_id = as.character(gpl_tab[[probe_col]]),
  gene_raw = as.character(gpl_tab[[gene_col]]),
  stringsAsFactors = FALSE
)

map_df$gene_raw <- str_trim(map_df$gene_raw)
map_df$gene_raw[map_df$gene_raw %in% c("", "NA", "NULL", "---")] <- NA

split_gene_symbol <- function(x) {
  x <- strsplit(x, "///|//|;|,")[[1]]
  x <- str_trim(x)
  x <- x[nchar(x) > 0]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

map_df$gene_symbol <- sapply(map_df$gene_raw, function(x) {
  if (is.na(x)) return(NA_character_)
  split_gene_symbol(x)
})

map_df <- map_df %>%
  filter(!is.na(probe_id), !is.na(gene_symbol), gene_symbol != "")

fwrite(
  map_df,
  file.path(out_dir, "24_GSE57495_probe2gene_mapping.tsv"),
  sep = "\t"
)

cat("Mapped probes with gene symbols:", nrow(map_df), "\n")

# ---------- 表达矩阵转 gene-level ----------
expr_mat <- as.data.frame(expr_df, check.names = FALSE, stringsAsFactors = FALSE)
rownames(expr_mat) <- expr_mat$probe_id
expr_mat$probe_id <- NULL
expr_mat <- as.matrix(expr_mat)
mode(expr_mat) <- "numeric"

common_probes <- intersect(rownames(expr_mat), map_df$probe_id)
cat("Common probes between expression and annotation:", length(common_probes), "\n")

if (length(common_probes) < 1000) {
  stop("Too few common probes between expression matrix and GPL annotation.")
}

expr_mat <- expr_mat[common_probes, , drop = FALSE]
map_sub  <- map_df %>%
  filter(probe_id %in% common_probes) %>%
  distinct(probe_id, .keep_all = TRUE)

map_sub <- map_sub[match(rownames(expr_mat), map_sub$probe_id), ]
stopifnot(all(rownames(expr_mat) == map_sub$probe_id))

# 每个 gene 保留平均表达最高 probe
probe_mean <- rowMeans(expr_mat, na.rm = TRUE)

collapse_df <- data.frame(
  probe_id = rownames(expr_mat),
  gene_symbol = map_sub$gene_symbol,
  probe_mean = probe_mean,
  stringsAsFactors = FALSE
) %>%
  group_by(gene_symbol) %>%
  slice_max(order_by = probe_mean, n = 1, with_ties = FALSE) %>%
  ungroup()

expr_gene <- expr_mat[collapse_df$probe_id, , drop = FALSE]
rownames(expr_gene) <- collapse_df$gene_symbol

expr_gene_df <- data.frame(
  gene = rownames(expr_gene),
  expr_gene,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

fwrite(
  expr_gene_df,
  file.path(out_dir, "24_GSE57495_expression_gene_level.tsv"),
  sep = "\t"
)

cat("Gene-level expression dim:", dim(expr_gene_df), "\n")

# ---------- 构建 survival-ready 临床表 ----------
required_cols <- c("geo_accession", "overall survival (month):ch1", "vital.status:ch1")
missing_req <- setdiff(required_cols, colnames(pheno_df))

if (length(missing_req) > 0) {
  stop("Missing required pheno columns: ", paste(missing_req, collapse = ", "))
}

cli <- data.frame(
  sample_id = pheno_df$geo_accession,
  OS_time = suppressWarnings(as.numeric(pheno_df[["overall survival (month):ch1"]])),
  OS_event = case_when(
    toupper(as.character(pheno_df[["vital.status:ch1"]])) == "DEAD"  ~ 1,
    toupper(as.character(pheno_df[["vital.status:ch1"]])) == "ALIVE" ~ 0,
    TRUE ~ NA_real_
  ),
  Stage = if ("Stage:ch1" %in% colnames(pheno_df)) as.character(pheno_df[["Stage:ch1"]]) else NA_character_,
  source_name = if ("source_name_ch1" %in% colnames(pheno_df)) as.character(pheno_df[["source_name_ch1"]]) else NA_character_,
  stringsAsFactors = FALSE
)

cli <- cli %>%
  filter(!is.na(OS_time), !is.na(OS_event)) %>%
  filter(OS_time > 0, OS_event %in% c(0, 1))

fwrite(
  cli,
  file.path(out_dir, "24_GSE57495_clinical_survival_ready.tsv"),
  sep = "\t"
)

cat("Survival-ready clinical samples:", nrow(cli), "\n")

# ---------- 审计 ----------
audit_df <- data.frame(
  metric = c(
    "n_expr_samples_raw",
    "n_pheno_samples_raw",
    "n_annotation_rows",
    "n_mapped_probes",
    "n_common_probes_expr_annotation",
    "n_gene_level_features",
    "n_survival_ready_samples",
    "n_events"
  ),
  value = c(
    ncol(expr_mat),
    nrow(pheno_df),
    nrow(gpl_tab),
    nrow(map_df),
    length(common_probes),
    nrow(expr_gene_df),
    nrow(cli),
    sum(cli$OS_event, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

fwrite(
  audit_df,
  file.path(out_dir, "24_GSE57495_processing_audit.tsv"),
  sep = "\t"
)

cat("\n24_build_gse57495_expression_and_clinical_tables_fix.R finished successfully.\n")