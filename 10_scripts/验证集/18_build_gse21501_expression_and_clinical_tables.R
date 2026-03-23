# =========================
# 18_build_gse21501_expression_and_clinical_tables.R
# 目的：
# 1. 下载 GPL4133 平台注释（有缓存则复用）
# 2. 建立 GSE21501 probe -> gene symbol 映射
# 3. 将 probe-level 表达矩阵转换为 gene-level
# 4. 从 pheno_raw 中构建 survival-ready 临床表
# 5. 输出审计文件
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

raw_dir <- "04_bulk_analysis/04_external_validation/01_raw_data/GEO/GSE21501"
out_dir <- "04_bulk_analysis/04_external_validation/02_processed_data/GSE21501"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr_file  <- file.path(raw_dir, "17_GSE21501_expression_probe_level.tsv")
pheno_file <- file.path(raw_dir, "17_GSE21501_pheno_raw.tsv")

# ---------- 检查输入文件 ----------
for (f in c(expr_file, pheno_file)) {
  if (!file.exists(f)) stop("Input file not found: ", f)
}

# ---------- 读取数据 ----------
expr_df  <- fread(expr_file, check.names = FALSE)
pheno_df <- fread(pheno_file, check.names = FALSE)

# FIX 1: pheno 重复列名会导致 dplyr::transmute 报错
colnames(pheno_df) <- make.unique(colnames(pheno_df), sep = "__dup")

stopifnot("probe_id" %in% colnames(expr_df))
stopifnot("geo_accession" %in% colnames(pheno_df))

cat("expr_df dim:", dim(expr_df), "\n")
cat("pheno_df dim:", dim(pheno_df), "\n")

# 原始表达样本数（减去 probe_id）
n_expr_samples_raw <- ncol(expr_df) - 1L

# ---------- 下载 GPL4133 注释（缓存优先） ----------
gpl_rds <- file.path(raw_dir, "18_GPL4133_gpl.rds")

if (file.exists(gpl_rds)) {
  cat("Loading cached GPL4133 from:", gpl_rds, "\n")
  gpl <- readRDS(gpl_rds)
} else {
  cat("Downloading GPL4133 annotation...\n")
  gpl <- tryCatch(
    getGEO("GPL4133", AnnotGPL = TRUE, destdir = raw_dir),
    error = function(e) stop("getGEO('GPL4133') failed: ", conditionMessage(e))
  )
  saveRDS(gpl, gpl_rds)
  cat("GPL4133 saved to:", gpl_rds, "\n")
}

gpl_tab <- as.data.frame(Table(gpl), check.names = FALSE, stringsAsFactors = FALSE)

fwrite(
  gpl_tab,
  file.path(out_dir, "18_GSE21501_GPL4133_annotation.tsv"),
  sep = "\t"
)

cat("GPL4133 annotation dim:", dim(gpl_tab), "\n")

# ---------- 自动识别 probe / gene symbol 列 ----------
probe_col_candidates <- c("ID", "ID_REF", "ProbeName", "probe_id", "SPOT_ID")
gene_col_candidates  <- c("Gene Symbol", "GENE_SYMBOL", "Gene symbol", "Symbol", "GENE", "GeneSymbol")

probe_col_matches <- intersect(probe_col_candidates, colnames(gpl_tab))
gene_col_matches  <- intersect(gene_col_candidates,  colnames(gpl_tab))

if (length(probe_col_matches) == 0) {
  stop(
    "GPL annotation: no probe ID column found.\nAvailable columns:\n",
    paste(colnames(gpl_tab), collapse = ", ")
  )
}
if (length(gene_col_matches) == 0) {
  stop(
    "GPL annotation: no gene symbol column found.\nAvailable columns:\n",
    paste(colnames(gpl_tab), collapse = ", ")
  )
}

probe_col <- probe_col_matches[1]
gene_col  <- gene_col_matches[1]

if (length(probe_col_matches) > 1) {
  warning("Multiple probe columns found; using: ", probe_col)
}
if (length(gene_col_matches) > 1) {
  warning("Multiple gene columns found; using: ", gene_col)
}

cat("Using probe column:", probe_col, "\n")
cat("Using gene column :", gene_col, "\n")

# ---------- probe -> gene symbol 映射 ----------
first_gene_symbol <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x %in% c("", "NA", "NULL", "---", "N/A", "n/a", "na")] <- NA_character_
  
  # 一个 probe 对多个 symbol 时，仅取第一个
  part1 <- str_split_fixed(x, "///|//|;|,", n = 2)[, 1]
  part1 <- str_trim(part1)
  part1[part1 == ""] <- NA_character_
  part1
}

map_df <- gpl_tab %>%
  transmute(
    probe_id    = as.character(.data[[probe_col]]),
    gene_symbol = first_gene_symbol(.data[[gene_col]])
  ) %>%
  filter(!is.na(probe_id), probe_id != "") %>%
  filter(!is.na(gene_symbol), gene_symbol != "")

fwrite(
  map_df,
  file.path(out_dir, "18_GSE21501_probe2gene_mapping.tsv"),
  sep = "\t"
)

cat("Mapped probes with gene symbols:", nrow(map_df), "\n")
cat("Unique genes in annotation:", dplyr::n_distinct(map_df$gene_symbol), "\n")

# ---------- probe-level 表达矩阵转 matrix ----------
expr_mat <- as.data.frame(expr_df, check.names = FALSE, stringsAsFactors = FALSE)
rownames(expr_mat) <- expr_mat$probe_id
expr_mat$probe_id <- NULL
expr_mat <- as.matrix(expr_mat)

storage.mode(expr_mat) <- "double"

na_frac <- mean(is.na(expr_mat))
cat(sprintf("NA fraction in expression matrix: %.4f\n", na_frac))

if (na_frac > 0.2) {
  warning("More than 20% NA in expression matrix after numeric conversion.")
}

# ---------- 表达矩阵与注释表匹配 ----------
common_probes <- intersect(rownames(expr_mat), map_df$probe_id)
cat("Common probes between expression and annotation:", length(common_probes), "\n")

if (length(common_probes) == 0) {
  stop("Zero common probes between expression matrix and GPL annotation.")
}
if (length(common_probes) < 1000) {
  stop("Too few common probes: ", length(common_probes))
}

expr_mat <- expr_mat[common_probes, , drop = FALSE]

map_sub <- map_df %>%
  filter(probe_id %in% common_probes) %>%
  distinct(probe_id, .keep_all = TRUE)

map_sub <- map_sub[match(rownames(expr_mat), map_sub$probe_id), ]

if (!identical(rownames(expr_mat), map_sub$probe_id)) {
  stop("Probe alignment failed between expression matrix and annotation.")
}

# ---------- 每个 gene 保留平均表达最高的 probe ----------
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

if (anyDuplicated(collapse_df$gene_symbol)) {
  stop("Duplicate gene symbols remained after collapsing probes.")
}

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
  file.path(out_dir, "18_GSE21501_expression_gene_level.tsv"),
  sep = "\t"
)

cat("Gene-level expression dim:", dim(expr_gene_df), "\n")

# ---------- 构建 survival-ready 临床表 ----------
required_cols <- c("geo_accession", "os time:ch2", "os event:ch2")
missing_req <- setdiff(required_cols, colnames(pheno_df))

if (length(missing_req) > 0) {
  stop(
    "Missing required pheno columns: ",
    paste(missing_req, collapse = ", "),
    "\nAvailable columns:\n",
    paste(colnames(pheno_df), collapse = ", ")
  )
}

col_or_na <- function(df, col) {
  if (col %in% colnames(df)) as.character(df[[col]]) else rep(NA_character_, nrow(df))
}

os_time_raw  <- pheno_df[["os time:ch2"]]
os_event_raw <- pheno_df[["os event:ch2"]]

os_time_na  <- sum(is.na(suppressWarnings(as.numeric(os_time_raw))))
os_event_na <- sum(is.na(suppressWarnings(as.numeric(os_event_raw))))

if (os_time_na > 0)  cat("[WARN] os time:ch2  -> NA count:", os_time_na, "\n")
if (os_event_na > 0) cat("[WARN] os event:ch2 -> NA count:", os_event_na, "\n")

# FIX 2: 不再用 transmute 直接处理可能异常的 pheno，而是先显式抽列
cli_raw <- data.frame(
  sample_id      = pheno_df[["geo_accession"]],
  OS_time        = suppressWarnings(as.numeric(pheno_df[["os time:ch2"]])),
  OS_event       = suppressWarnings(as.numeric(pheno_df[["os event:ch2"]])),
  GEO_risk_group = col_or_na(pheno_df, "risk group:ch2"),
  T_stage        = col_or_na(pheno_df, "t stage:ch2"),
  N_stage        = col_or_na(pheno_df, "n stage:ch2"),
  tissue_ch2     = col_or_na(pheno_df, "tissue:ch2"),
  stringsAsFactors = FALSE
)

cat("Pheno samples before survival filter:", nrow(cli_raw), "\n")

cli <- cli_raw %>%
  filter(!is.na(OS_time), !is.na(OS_event)) %>%
  filter(OS_time > 0, OS_event %in% c(0, 1))

cat("Survival-ready clinical samples:", nrow(cli),
    "(dropped:", nrow(cli_raw) - nrow(cli), ")\n")

fwrite(
  cli,
  file.path(out_dir, "18_GSE21501_clinical_survival_ready.tsv"),
  sep = "\t"
)

# ---------- 样本 ID 交叉检查 ----------
expr_samples <- colnames(expr_gene_df)[-1]
cli_samples  <- cli$sample_id
shared_samples <- intersect(expr_samples, cli_samples)

cat("Samples in expression:", length(expr_samples), "\n")
cat("Samples in clinical:  ", length(cli_samples), "\n")
cat("Overlapping samples:  ", length(shared_samples), "\n")

if (length(shared_samples) == 0) {
  warning("No overlapping sample IDs between expression and clinical data.")
}

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
    "n_overlapping_expr_cli_samples"
  ),
  value = c(
    n_expr_samples_raw,
    nrow(pheno_df),
    nrow(gpl_tab),
    nrow(map_df),
    length(common_probes),
    nrow(expr_gene_df),
    nrow(cli),
    length(shared_samples)
  ),
  stringsAsFactors = FALSE
)

fwrite(
  audit_df,
  file.path(out_dir, "18_GSE21501_processing_audit.tsv"),
  sep = "\t"
)

cat("\n18_build_gse21501_expression_and_clinical_tables.R finished successfully.\n")