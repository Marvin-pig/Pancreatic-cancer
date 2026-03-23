# =========================
# 38_reconstruct_eqtlgen_beta_se_eaf_and_write_exposure_files.R
# 目的：
# 1. 读取 37 拆出的基因级原始候选表
# 2. 尝试读取/下载 eQTLGen allele frequency 文件
# 3. 合并 EAF
# 4. 在 zscore + samplesize + eaf 可用时近似重建 beta / se
# 5. 仅对结构完整的基因写出正式 exposure 文件
# 6. 生成完整审计表，区分：
#    - 可进入 MR
#    - 仅部分可用
#    - 暂不可用
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ---------- 路径 ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "05_MR"

batch1_file     <- file.path(root_dir, "02_processed_data/batch1_queries/29_batch1_priority_genes.tsv")
schema_map_file <- file.path(root_dir, "00_registry/37_eqtlgen_schema_map.tsv")
row_count_file  <- file.path(root_dir, "00_registry/37_eqtlgen_gene_row_counts.tsv")
cand_dir        <- file.path(root_dir, "02_processed_data/37_candidate_gene_tables")
eqtlgen_dir     <- file.path(root_dir, "01_raw_data/exposure/eQTLGen")
out_registry    <- file.path(root_dir, "00_registry")
out_logs        <- file.path(root_dir, "06_logs")

dir.create(out_registry, recursive = TRUE, showWarnings = FALSE)
dir.create(out_logs, recursive = TRUE, showWarnings = FALSE)
dir.create(eqtlgen_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(batch1_file)) stop("缺少 batch1_file: ", batch1_file)
if (!file.exists(schema_map_file)) stop("缺少 schema_map_file: ", schema_map_file)
if (!file.exists(row_count_file)) stop("缺少 row_count_file: ", row_count_file)

batch1    <- fread(batch1_file)
schema_map<- fread(schema_map_file)
row_count <- fread(row_count_file)

if (!"gene" %in% names(batch1)) stop("29_batch1_priority_genes.tsv 缺少 gene 列。")
genes <- unique(as.character(batch1$gene))

# ---------- 读取 schema mapping ----------
get_detected_col <- function(target_name) {
  x <- schema_map %>% filter(target_field == target_name) %>% pull(detected_column)
  if (length(x) == 0) return(NA_character_)
  x[1]
}

col_gene_symbol <- get_detected_col("gene_symbol")
col_gene_id     <- get_detected_col("gene_id")
col_snp         <- get_detected_col("SNP")
col_chr         <- get_detected_col("chr")
col_pos         <- get_detected_col("pos")
col_ea          <- get_detected_col("effect_allele")
col_oa          <- get_detected_col("other_allele")
col_p           <- get_detected_col("pval")
col_z           <- get_detected_col("zscore")
col_n           <- get_detected_col("samplesize")

required_raw_cols <- c(col_snp, col_chr, col_pos, col_ea, col_oa, col_p, col_z, col_n)
required_raw_cols <- required_raw_cols[!is.na(required_raw_cols)]

# ---------- Step 1: allele frequency 文件准备 ----------
af_file <- file.path(eqtlgen_dir, "36_eqtlgen_allele_frequencies.txt.gz")
af_url  <- "https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/eQTLGen_AF.txt.gz"

safe_download <- function(url, dest) {
  ok <- FALSE
  msg <- ""
  if (!file.exists(dest)) {
    tryCatch({
      download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
      ok <- file.exists(dest)
    }, error = function(e) {
      msg <<- conditionMessage(e)
    })
  } else {
    ok <- TRUE
    msg <- "already_exists"
  }
  list(ok = ok, msg = msg)
}

af_download_status <- data.frame(
  file = af_file,
  exists_before = file.exists(af_file),
  attempted_download = FALSE,
  exists_after = file.exists(af_file),
  message = "",
  stringsAsFactors = FALSE
)

if (!file.exists(af_file)) {
  af_download_status$attempted_download <- TRUE
  res <- safe_download(af_url, af_file)
  af_download_status$exists_after <- file.exists(af_file)
  af_download_status$message <- res$msg
} else {
  af_download_status$message <- "already_exists"
}

fwrite(
  af_download_status,
  file.path(out_registry, "38_af_download_status.tsv"),
  sep = "\t"
)

# ---------- Step 2: 读取 allele frequency ----------
af_dt <- NULL
af_snp_col <- NA_character_
af_eaf_col <- NA_character_
af_a1_col  <- NA_character_
af_a2_col  <- NA_character_

if (file.exists(af_file)) {
  af_header <- tryCatch(names(fread(af_file, nrows = 0)), error = function(e) character(0))
  
  detect_one <- function(cands, header_vec) {
    hit <- cands[cands %in% header_vec]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }
  
  af_snp_col <- detect_one(c("SNP", "rsid", "variant_id"), af_header)
  af_eaf_col <- detect_one(c("AF", "EAF", "eaf", "freq", "frequency", "effect_allele_freq"), af_header)
  af_a1_col  <- detect_one(c("A1", "AssessedAllele", "effect_allele", "EA"), af_header)
  af_a2_col  <- detect_one(c("A2", "OtherAllele", "other_allele", "OA"), af_header)
  
  af_select <- unique(na.omit(c(af_snp_col, af_eaf_col, af_a1_col, af_a2_col)))
  
  if (length(af_select) >= 2) {
    af_dt <- tryCatch(
      fread(af_file, select = af_select),
      error = function(e) NULL
    )
  }
}

af_schema <- data.frame(
  target_field = c("SNP", "EAF", "A1", "A2"),
  detected_column = c(af_snp_col, af_eaf_col, af_a1_col, af_a2_col),
  stringsAsFactors = FALSE
)

fwrite(
  af_schema,
  file.path(out_registry, "38_af_schema_map.tsv"),
  sep = "\t"
)

# ---------- Step 3: 核心函数 ----------
# 近似重建公式（仅作当前项目推进用途）：
# se ≈ 1 / sqrt(2 * eaf * (1 - eaf) * (n + z^2))
# beta = z * se
#
# 使用前提：
# 1) zscore 非缺失
# 2) samplesize 非缺失
# 3) eaf 在 (0,1)
#
# 这是近似重建，不是官方原始 beta/se

reconstruct_beta_se <- function(z, n, eaf) {
  valid <- !is.na(z) & !is.na(n) & !is.na(eaf) & eaf > 0 & eaf < 1 & n > 0
  se <- rep(NA_real_, length(z))
  beta <- rep(NA_real_, length(z))
  
  se[valid]   <- 1 / sqrt(2 * eaf[valid] * (1 - eaf[valid]) * (n[valid] + z[valid]^2))
  beta[valid] <- z[valid] * se[valid]
  
  list(beta = beta, se = se, valid = valid)
}

# 对齐 EAF：如果 AF 文件提供的是 A1 frequency，而 A1 不一定等于当前表的 effect_allele，
# 则需要根据等位基因方向进行转换。
align_eaf <- function(raw_dt, af_dt, raw_snp, raw_ea, raw_oa, af_snp, af_eaf, af_a1 = NA, af_a2 = NA) {
  x <- copy(raw_dt)
  x[, .__rowid__ := .I]
  
  af2 <- copy(af_dt)
  setnames(af2, af_snp, ".__join_snp__")
  setnames(x, raw_snp, ".__join_snp__")
  
  x <- merge(x, af2, by = ".__join_snp__", all.x = TRUE)
  
  eaf_vec <- rep(NA_real_, nrow(x))
  
  if (!is.na(af_a1) && af_a1 %in% names(x) && !is.na(af_eaf) && af_eaf %in% names(x)) {
    raw_ea_vec <- as.character(x[[raw_ea]])
    raw_oa_vec <- as.character(x[[raw_oa]])
    af_a1_vec  <- as.character(x[[af_a1]])
    
    af_freq <- suppressWarnings(as.numeric(x[[af_eaf]]))
    
    # 若 AF 的 A1 与当前 effect_allele 一致，则直接用
    same_as_ea <- !is.na(raw_ea_vec) & !is.na(af_a1_vec) & raw_ea_vec == af_a1_vec
    # 若 AF 的 A1 与当前 other_allele 一致，则取 1 - AF
    same_as_oa <- !is.na(raw_oa_vec) & !is.na(af_a1_vec) & raw_oa_vec == af_a1_vec
    
    eaf_vec[same_as_ea] <- af_freq[same_as_ea]
    eaf_vec[same_as_oa] <- 1 - af_freq[same_as_oa]
  } else if (!is.na(af_eaf) && af_eaf %in% names(x)) {
    eaf_vec <- suppressWarnings(as.numeric(x[[af_eaf]]))
  }
  
  x$eaf_reconstructed <- eaf_vec
  
  # 还原 SNP 列名
  setnames(x, ".__join_snp__", raw_snp)
  x[, .__rowid__ := NULL]
  x
}

# ---------- Step 4: 循环处理每个基因 ----------
gene_audit_list <- list()
row_audit_list  <- list()

for (g in genes) {
  raw_file <- file.path(cand_dir, paste0(g, "_eqtlgen_raw.tsv"))
  
  gene_audit <- data.frame(
    gene = g,
    raw_file_exists = file.exists(raw_file),
    raw_n = NA_integer_,
    eaf_joined_n = 0L,
    reconstructable_n = 0L,
    complete_n = 0L,
    written_official_exposure = FALSE,
    official_exposure_file = file.path(eqtlgen_dir, paste0(g, "_cis_eqtl.tsv")),
    reason = "",
    stringsAsFactors = FALSE
  )
  
  if (!file.exists(raw_file)) {
    gene_audit$reason <- "raw_gene_table_missing"
    gene_audit_list[[g]] <- gene_audit
    next
  }
  
  x <- fread(raw_file)
  gene_audit$raw_n <- nrow(x)
  
  missing_cols <- setdiff(required_raw_cols, names(x))
  if (length(missing_cols) > 0) {
    gene_audit$reason <- paste0("missing_raw_columns:", paste(missing_cols, collapse = ","))
    gene_audit_list[[g]] <- gene_audit
    next
  }
  
  # 合并 EAF
  if (!is.null(af_dt) && !is.na(af_snp_col) && !is.na(af_eaf_col)) {
    x2 <- align_eaf(
      raw_dt  = x,
      af_dt   = af_dt,
      raw_snp = col_snp,
      raw_ea  = col_ea,
      raw_oa  = col_oa,
      af_snp  = af_snp_col,
      af_eaf  = af_eaf_col,
      af_a1   = af_a1_col,
      af_a2   = af_a2_col
    )
  } else {
    x2 <- copy(x)
    x2$eaf_reconstructed <- NA_real_
  }
  
  gene_audit$eaf_joined_n <- sum(!is.na(x2$eaf_reconstructed))
  
  # 近似重建 beta / se
  z_vec <- suppressWarnings(as.numeric(x2[[col_z]]))
  n_vec <- suppressWarnings(as.numeric(x2[[col_n]]))
  eaf_vec <- suppressWarnings(as.numeric(x2$eaf_reconstructed))
  
  rec <- reconstruct_beta_se(z = z_vec, n = n_vec, eaf = eaf_vec)
  x2$beta_reconstructed <- rec$beta
  x2$se_reconstructed   <- rec$se
  
  gene_audit$reconstructable_n <- sum(rec$valid, na.rm = TRUE)
  
  # 生成标准 exposure 表
  std <- data.table(
    SNP           = as.character(x2[[col_snp]]),
    chr           = suppressWarnings(as.integer(x2[[col_chr]])),
    pos           = suppressWarnings(as.integer(x2[[col_pos]])),
    effect_allele = as.character(x2[[col_ea]]),
    other_allele  = as.character(x2[[col_oa]]),
    beta          = suppressWarnings(as.numeric(x2$beta_reconstructed)),
    se            = suppressWarnings(as.numeric(x2$se_reconstructed)),
    pval          = suppressWarnings(as.numeric(x2[[col_p]])),
    eaf           = suppressWarnings(as.numeric(x2$eaf_reconstructed)),
    gene          = g,
    exposure      = paste0("eQTLGen_", g)
  )
  
  # 基本 QC
  std <- std %>%
    filter(
      !is.na(SNP),
      !is.na(effect_allele),
      !is.na(other_allele),
      effect_allele %in% c("A","C","G","T"),
      other_allele  %in% c("A","C","G","T"),
      effect_allele != other_allele,
      !is.na(pval),
      pval > 0,
      pval <= 1
    ) %>%
    distinct()
  
  # 完整行
  complete_idx <- with(std,
                       !is.na(SNP) &
                         !is.na(chr) &
                         !is.na(pos) &
                         !is.na(effect_allele) &
                         !is.na(other_allele) &
                         !is.na(beta) &
                         !is.na(se) &
                         !is.na(pval) &
                         !is.na(eaf) &
                         !is.na(gene) &
                         !is.na(exposure)
  )
  
  std_complete <- std[complete_idx, ]
  gene_audit$complete_n <- nrow(std_complete)
  
  # 写出中间结果
  fwrite(
    x2,
    file.path(cand_dir, paste0(g, "_with_reconstructed_fields.tsv")),
    sep = "\t",
    na = "NA"
  )
  
  fwrite(
    std,
    file.path(cand_dir, paste0(g, "_candidate_exposure_all_rows.tsv")),
    sep = "\t",
    na = "NA"
  )
  
  if (nrow(std_complete) > 0) {
    fwrite(
      std_complete,
      file.path(eqtlgen_dir, paste0(g, "_cis_eqtl.tsv")),
      sep = "\t",
      na = "NA"
    )
    gene_audit$written_official_exposure <- TRUE
    gene_audit$reason <- "official_exposure_written_from_approx_reconstruction"
  } else {
    gene_audit$reason <- "no_complete_rows_after_eaf_join_and_beta_se_reconstruction"
  }
  
  # 行级审计摘要
  row_audit <- data.frame(
    gene = g,
    raw_rows = nrow(x),
    rows_after_basic_qc = nrow(std),
    complete_rows = nrow(std_complete),
    stringsAsFactors = FALSE
  )
  
  gene_audit_list[[g]] <- gene_audit
  row_audit_list[[g]]  <- row_audit
}

gene_audit_tab <- bind_rows(gene_audit_list)
row_audit_tab  <- bind_rows(row_audit_list)

# ---------- Step 5: 全局总结 ----------
status_summary <- data.frame(
  item = c(
    "af_file_exists",
    "n_genes_total",
    "n_genes_with_raw_table",
    "n_genes_written_official_exposure",
    "recommended_next_action"
  ),
  value = c(
    file.exists(af_file),
    length(genes),
    sum(gene_audit_tab$raw_file_exists, na.rm = TRUE),
    sum(gene_audit_tab$written_official_exposure, na.rm = TRUE),
    ifelse(
      sum(gene_audit_tab$written_official_exposure, na.rm = TRUE) > 0,
      "Proceed to 39: extract instruments, clump, and prepare BBJ outcome retrieval",
      "Check AF file schema / reconstruction assumptions before MR"
    )
  ),
  stringsAsFactors = FALSE
)

# ---------- 输出 ----------
fwrite(
  gene_audit_tab,
  file.path(out_registry, "38_gene_reconstruction_audit.tsv"),
  sep = "\t"
)

fwrite(
  row_audit_tab,
  file.path(out_registry, "38_gene_row_qc_summary.tsv"),
  sep = "\t"
)

fwrite(
  status_summary,
  file.path(out_registry, "38_status_summary.tsv"),
  sep = "\t"
)

# ---------- 日志 ----------
notes <- c(
  "=== Code 38 Summary ===",
  paste0("Created at: ", Sys.time()),
  "",
  "[Goal]",
  "Reconstruct approximate beta / se using zscore + samplesize + eaf.",
  "",
  "[Important warning]",
  "This is approximate reconstruction, not official original beta/se.",
  "Use for project progression and method testing; document clearly in Methods if retained.",
  "",
  "[Output tables]",
  "- 38_af_download_status.tsv",
  "- 38_af_schema_map.tsv",
  "- 38_gene_reconstruction_audit.tsv",
  "- 38_gene_row_qc_summary.tsv",
  "- 38_status_summary.tsv",
  "",
  "[Per-gene intermediate files]",
  "- 02_processed_data/37_candidate_gene_tables/<GENE>_with_reconstructed_fields.tsv",
  "- 02_processed_data/37_candidate_gene_tables/<GENE>_candidate_exposure_all_rows.tsv",
  "",
  "[Official exposure output path if successful]",
  "- 01_raw_data/exposure/eQTLGen/<GENE>_cis_eqtl.tsv"
)

writeLines(
  notes,
  file.path(out_logs, "38_reconstruct_beta_se_eaf_notes.txt")
)

cat("38_reconstruct_eqtlgen_beta_se_eaf_and_write_exposure_files.R finished successfully.\n")
cat("Generated files:\n")
cat("- 00_registry/38_af_download_status.tsv\n")
cat("- 00_registry/38_af_schema_map.tsv\n")
cat("- 00_registry/38_gene_reconstruction_audit.tsv\n")
cat("- 00_registry/38_gene_row_qc_summary.tsv\n")
cat("- 00_registry/38_status_summary.tsv\n")
cat("- 02_processed_data/37_candidate_gene_tables/<GENE>_with_reconstructed_fields.tsv\n")
cat("- 02_processed_data/37_candidate_gene_tables/<GENE>_candidate_exposure_all_rows.tsv\n")
cat("- 01_raw_data/exposure/eQTLGen/<GENE>_cis_eqtl.tsv (if complete rows exist)\n")
cat("- 06_logs/38_reconstruct_beta_se_eaf_notes.txt\n")