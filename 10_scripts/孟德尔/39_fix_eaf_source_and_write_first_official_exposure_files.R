# =========================
# 39_fix_eaf_source_and_write_first_official_exposure_files.R
# 目的：
# 1. 使用正确的 eQTLGen AF 文件来源修复 EAF
# 2. 自动识别 AF 文件字段并与 37 的 raw gene tables 合并
# 3. 重新近似重建 beta / se
# 4. 尽可能写出第一批正式 exposure 文件
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

batch1     <- fread(batch1_file)
schema_map <- fread(schema_map_file)
row_count  <- fread(row_count_file)

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

# ---------- Step 1: 修复 AF 来源 ----------
# 正确的 eQTLGen AF 文件名（当前公开目录可见）
af_primary_path <- file.path(eqtlgen_dir, "2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz")
af_alt_path     <- file.path(eqtlgen_dir, "39_eqtlgen_AF_fixed_source.txt.gz")

af_urls <- c(
  "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz",
  "https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz"
)

safe_download <- function(url, dest) {
  ok <- FALSE
  msg <- ""
  tryCatch({
    download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
    ok <- file.exists(dest)
  }, error = function(e) {
    msg <<- conditionMessage(e)
  })
  list(ok = ok, msg = msg)
}

download_attempts <- list()

if (!file.exists(af_primary_path)) {
  for (i in seq_along(af_urls)) {
    tmp_dest <- if (i == 1) af_primary_path else af_alt_path
    res <- safe_download(af_urls[i], tmp_dest)
    download_attempts[[i]] <- data.frame(
      attempt = i,
      url = af_urls[i],
      dest = tmp_dest,
      ok = isTRUE(res$ok),
      message = res$msg,
      stringsAsFactors = FALSE
    )
    if (isTRUE(res$ok)) break
  }
} else {
  download_attempts[[1]] <- data.frame(
    attempt = 0,
    url = "local_existing",
    dest = af_primary_path,
    ok = TRUE,
    message = "already_exists",
    stringsAsFactors = FALSE
  )
}

af_download_status <- bind_rows(download_attempts)

# 统一最终 AF 路径
af_file <- NA_character_
if (file.exists(af_primary_path)) {
  af_file <- af_primary_path
} else if (file.exists(af_alt_path)) {
  af_file <- af_alt_path
}

fwrite(
  af_download_status,
  file.path(out_registry, "39_af_download_attempts.tsv"),
  sep = "\t"
)

# ---------- Step 2: 自动识别 AF schema ----------
detect_one <- function(cands, header_vec) {
  hit <- cands[cands %in% header_vec]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

af_dt <- NULL
af_snp_col <- NA_character_
af_eaf_col <- NA_character_
af_a1_col  <- NA_character_
af_a2_col  <- NA_character_
af_chr_col <- NA_character_
af_pos_col <- NA_character_

if (!is.na(af_file) && file.exists(af_file)) {
  af_header <- tryCatch(names(fread(af_file, nrows = 0)), error = function(e) character(0))
  
  af_snp_col <- detect_one(c("SNP", "rsid", "ID"), af_header)
  # 这个 AF 文件更可能以 AlleleB / MAF 命名
  af_eaf_col <- detect_one(c("AF", "EAF", "eaf", "AlleleBFrequency", "AlleleB_Frequency", "MAF"), af_header)
  af_a1_col  <- detect_one(c("A1", "AlleleA", "RefAllele", "AssessedAllele"), af_header)
  af_a2_col  <- detect_one(c("A2", "AlleleB", "AltAllele", "OtherAllele"), af_header)
  af_chr_col <- detect_one(c("chr", "CHR", "Chromosome"), af_header)
  af_pos_col <- detect_one(c("pos", "BP", "Position"), af_header)
  
  # 至少把可能有用的列都读进来
  af_select <- unique(na.omit(c(af_snp_col, af_eaf_col, af_a1_col, af_a2_col, af_chr_col, af_pos_col)))
  if (length(af_select) >= 2) {
    af_dt <- tryCatch(fread(af_file, select = af_select), error = function(e) NULL)
  }
}

af_schema <- data.frame(
  target_field = c("SNP", "EAF_or_MAF", "AlleleA_or_A1", "AlleleB_or_A2", "chr", "pos"),
  detected_column = c(af_snp_col, af_eaf_col, af_a1_col, af_a2_col, af_chr_col, af_pos_col),
  stringsAsFactors = FALSE
)

fwrite(
  af_schema,
  file.path(out_registry, "39_af_schema_map.tsv"),
  sep = "\t"
)

# ---------- Step 3: 工具函数 ----------
# 近似重建:
# se ≈ 1 / sqrt(2 * eaf * (1 - eaf) * (n + z^2))
# beta = z * se
reconstruct_beta_se <- function(z, n, eaf) {
  valid <- !is.na(z) & !is.na(n) & !is.na(eaf) & eaf > 0 & eaf < 1 & n > 0
  se <- rep(NA_real_, length(z))
  beta <- rep(NA_real_, length(z))
  
  se[valid]   <- 1 / sqrt(2 * eaf[valid] * (1 - eaf[valid]) * (n[valid] + z[valid]^2))
  beta[valid] <- z[valid] * se[valid]
  
  list(beta = beta, se = se, valid = valid)
}

# 关键：AF 文件可能是 AlleleB 的频率，不一定是当前 effect allele 的频率
# 这版优先策略：
# 1) rsID 精确匹配
# 2) 若 AF 提供 AlleleB / A2，且当前 effect_allele == AlleleB，则 EAF = freq
# 3) 若 current effect_allele == AlleleA，则 EAF = 1 - freq
# 4) 若只拿到 MAF 而拿不到等位基因方向，则只在能和 raw alleles 明确对上时使用
align_eaf_from_af <- function(raw_dt, af_dt,
                              raw_snp, raw_chr, raw_pos, raw_ea, raw_oa,
                              af_snp = NA, af_eaf = NA, af_a1 = NA, af_a2 = NA,
                              af_chr = NA, af_pos = NA) {
  x <- copy(raw_dt)
  
  if (is.null(af_dt) || is.na(af_eaf)) {
    x$eaf_reconstructed <- NA_real_
    x$eaf_match_mode <- "no_af_data"
    return(x)
  }
  
  af2 <- copy(af_dt)
  
  # 优先 rsID
  if (!is.na(af_snp) && af_snp %in% names(af2) && raw_snp %in% names(x)) {
    setnames(af2, af_snp, ".__join_snp__", skip_absent = TRUE)
    setnames(x,   raw_snp, ".__join_snp__", skip_absent = TRUE)
    x <- merge(x, af2, by = ".__join_snp__", all.x = TRUE)
    setnames(x, ".__join_snp__", raw_snp, skip_absent = TRUE)
  } else if (!is.na(af_chr) && !is.na(af_pos) && af_chr %in% names(af2) && af_pos %in% names(af2)) {
    setnames(af2, af_chr, ".__join_chr__", skip_absent = TRUE)
    setnames(af2, af_pos, ".__join_pos__", skip_absent = TRUE)
    setnames(x, raw_chr,  ".__join_chr__", skip_absent = TRUE)
    setnames(x, raw_pos,  ".__join_pos__", skip_absent = TRUE)
    x <- merge(x, af2, by = c(".__join_chr__", ".__join_pos__"), all.x = TRUE)
    setnames(x, ".__join_chr__", raw_chr, skip_absent = TRUE)
    setnames(x, ".__join_pos__", raw_pos, skip_absent = TRUE)
  } else {
    x$eaf_reconstructed <- NA_real_
    x$eaf_match_mode <- "no_join_key"
    return(x)
  }
  
  freq_vec <- suppressWarnings(as.numeric(x[[af_eaf]]))
  eaf_vec  <- rep(NA_real_, nrow(x))
  mode_vec <- rep("unresolved", nrow(x))
  
  raw_ea_vec <- as.character(x[[raw_ea]])
  raw_oa_vec <- as.character(x[[raw_oa]])
  
  # 如果 AF 提供 allele identity
  if (!is.na(af_a2) && af_a2 %in% names(x)) {
    af_b_vec <- as.character(x[[af_a2]])
    ea_is_b <- !is.na(raw_ea_vec) & !is.na(af_b_vec) & raw_ea_vec == af_b_vec
    oa_is_b <- !is.na(raw_oa_vec) & !is.na(af_b_vec) & raw_oa_vec == af_b_vec
    
    eaf_vec[ea_is_b] <- freq_vec[ea_is_b]
    mode_vec[ea_is_b] <- "effect_allele_is_AlleleB"
    
    eaf_vec[oa_is_b] <- 1 - freq_vec[oa_is_b]
    mode_vec[oa_is_b] <- "other_allele_is_AlleleB_invert"
  }
  
  # 若仍未解决，尝试用 af_a1
  if (!is.na(af_a1) && af_a1 %in% names(x)) {
    af_a_vec <- as.character(x[[af_a1]])
    unresolved <- is.na(eaf_vec)
    
    ea_is_a <- unresolved & !is.na(raw_ea_vec) & !is.na(af_a_vec) & raw_ea_vec == af_a_vec
    oa_is_a <- unresolved & !is.na(raw_oa_vec) & !is.na(af_a_vec) & raw_oa_vec == af_a_vec
    
    # 若 af_eaf 是 AlleleB/Alt 频率，则 effect_allele==AlleleA 时应取 1-freq
    eaf_vec[ea_is_a] <- 1 - freq_vec[ea_is_a]
    mode_vec[ea_is_a] <- "effect_allele_is_AlleleA_invert"
    
    eaf_vec[oa_is_a] <- freq_vec[oa_is_a]
    mode_vec[oa_is_a] <- "other_allele_is_AlleleA_keep"
  }
  
  # 若没有 allele identity，但频率列名其实是 MAF，则保守不用
  x$eaf_reconstructed <- eaf_vec
  x$eaf_match_mode <- mode_vec
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
  
  x2 <- align_eaf_from_af(
    raw_dt = x,
    af_dt = af_dt,
    raw_snp = col_snp,
    raw_chr = col_chr,
    raw_pos = col_pos,
    raw_ea = col_ea,
    raw_oa = col_oa,
    af_snp = af_snp_col,
    af_eaf = af_eaf_col,
    af_a1 = af_a1_col,
    af_a2 = af_a2_col,
    af_chr = af_chr_col,
    af_pos = af_pos_col
  )
  
  gene_audit$eaf_joined_n <- sum(!is.na(x2$eaf_reconstructed))
  
  z_vec   <- suppressWarnings(as.numeric(x2[[col_z]]))
  n_vec   <- suppressWarnings(as.numeric(x2[[col_n]]))
  eaf_vec <- suppressWarnings(as.numeric(x2$eaf_reconstructed))
  
  rec <- reconstruct_beta_se(z = z_vec, n = n_vec, eaf = eaf_vec)
  x2$beta_reconstructed <- rec$beta
  x2$se_reconstructed   <- rec$se
  
  gene_audit$reconstructable_n <- sum(rec$valid, na.rm = TRUE)
  
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
  
  fwrite(
    x2,
    file.path(cand_dir, paste0(g, "_39_with_fixed_eaf.tsv")),
    sep = "\t",
    na = "NA"
  )
  
  fwrite(
    std,
    file.path(cand_dir, paste0(g, "_39_candidate_exposure_all_rows.tsv")),
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
    gene_audit$reason <- "official_exposure_written_after_fixing_af_source"
  } else {
    gene_audit$reason <- "no_complete_rows_after_fixed_af_source"
  }
  
  row_audit <- data.frame(
    gene = g,
    raw_rows = nrow(x),
    eaf_joined_rows = sum(!is.na(x2$eaf_reconstructed)),
    reconstructable_rows = sum(rec$valid, na.rm = TRUE),
    rows_after_basic_qc = nrow(std),
    complete_rows = nrow(std_complete),
    stringsAsFactors = FALSE
  )
  
  gene_audit_list[[g]] <- gene_audit
  row_audit_list[[g]]  <- row_audit
}

gene_audit_tab <- bind_rows(gene_audit_list)
row_audit_tab  <- bind_rows(row_audit_list)

status_summary <- data.frame(
  item = c(
    "af_file_exists",
    "n_genes_total",
    "n_genes_with_raw_table",
    "n_genes_written_official_exposure",
    "recommended_next_action"
  ),
  value = c(
    !is.na(af_file) && file.exists(af_file),
    length(genes),
    sum(gene_audit_tab$raw_file_exists, na.rm = TRUE),
    sum(gene_audit_tab$written_official_exposure, na.rm = TRUE),
    ifelse(
      sum(gene_audit_tab$written_official_exposure, na.rm = TRUE) > 0,
      "Proceed to 40: extract instruments and prepare BBJ outcome retrieval",
      "If still zero, inspect 39_af_schema_map.tsv and manually verify AF allele columns"
    )
  ),
  stringsAsFactors = FALSE
)

fwrite(
  gene_audit_tab,
  file.path(out_registry, "39_gene_reconstruction_audit.tsv"),
  sep = "\t"
)

fwrite(
  row_audit_tab,
  file.path(out_registry, "39_gene_row_qc_summary.tsv"),
  sep = "\t"
)

fwrite(
  status_summary,
  file.path(out_registry, "39_status_summary.tsv"),
  sep = "\t"
)

notes <- c(
  "=== Code 39 Summary ===",
  paste0("Created at: ", Sys.time()),
  "",
  "[Goal]",
  "Fix EAF source using the current eQTLGen AF filename and rewrite official exposure files if possible.",
  "",
  "[Important warning]",
  "beta/se are still approximate reconstructions from z, N, and EAF.",
  "",
  "[Outputs]",
  "- 39_af_download_attempts.tsv",
  "- 39_af_schema_map.tsv",
  "- 39_gene_reconstruction_audit.tsv",
  "- 39_gene_row_qc_summary.tsv",
  "- 39_status_summary.tsv",
  "- 02_processed_data/37_candidate_gene_tables/<GENE>_39_with_fixed_eaf.tsv",
  "- 02_processed_data/37_candidate_gene_tables/<GENE>_39_candidate_exposure_all_rows.tsv",
  "- 01_raw_data/exposure/eQTLGen/<GENE>_cis_eqtl.tsv (if successful)"
)

writeLines(
  notes,
  file.path(out_logs, "39_fix_eaf_source_notes.txt")
)

cat("39_fix_eaf_source_and_write_first_official_exposure_files.R finished successfully.\n")
cat("Generated files:\n")
cat("- 00_registry/39_af_download_attempts.tsv\n")
cat("- 00_registry/39_af_schema_map.tsv\n")
cat("- 00_registry/39_gene_reconstruction_audit.tsv\n")
cat("- 00_registry/39_gene_row_qc_summary.tsv\n")
cat("- 00_registry/39_status_summary.tsv\n")
cat("- 02_processed_data/37_candidate_gene_tables/<GENE>_39_with_fixed_eaf.tsv\n")
cat("- 02_processed_data/37_candidate_gene_tables/<GENE>_39_candidate_exposure_all_rows.tsv\n")
cat("- 01_raw_data/exposure/eQTLGen/<GENE>_cis_eqtl.tsv (if successful)\n")
cat("- 06_logs/39_fix_eaf_source_notes.txt\n")