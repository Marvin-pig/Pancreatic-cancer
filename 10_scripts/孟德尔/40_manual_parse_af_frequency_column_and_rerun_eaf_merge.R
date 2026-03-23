# =========================
# 40_manual_parse_af_frequency_column_and_rerun_eaf_merge.R
# 修正版
# 目的：
# 1. 完整预览 eQTLGen AF 文件的表头与前几行
# 2. 手动指定频率列 = AlleleB_all，且该频率对应 AlleleB
# 3. 重新合并 EAF，并近似重建 beta / se
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
cand_dir        <- file.path(root_dir, "02_processed_data/37_candidate_gene_tables")
eqtlgen_dir     <- file.path(root_dir, "01_raw_data/exposure/eQTLGen")
out_registry    <- file.path(root_dir, "00_registry")
out_logs        <- file.path(root_dir, "06_logs")

dir.create(out_registry, recursive = TRUE, showWarnings = FALSE)
dir.create(out_logs, recursive = TRUE, showWarnings = FALSE)
dir.create(eqtlgen_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(batch1_file)) stop("缺少 batch1_file: ", batch1_file)
if (!file.exists(schema_map_file)) stop("缺少 schema_map_file: ", schema_map_file)

batch1     <- fread(batch1_file)
schema_map <- fread(schema_map_file)

if (!"gene" %in% names(batch1)) stop("29_batch1_priority_genes.tsv 缺少 gene 列。")
genes <- unique(as.character(batch1$gene))

# ---------- 读取 37 的 schema mapping ----------
get_detected_col <- function(target_name) {
  x <- schema_map %>% filter(target_field == target_name) %>% pull(detected_column)
  if (length(x) == 0) return(NA_character_)
  x[1]
}

col_snp <- get_detected_col("SNP")
col_chr <- get_detected_col("chr")
col_pos <- get_detected_col("pos")
col_ea  <- get_detected_col("effect_allele")
col_oa  <- get_detected_col("other_allele")
col_p   <- get_detected_col("pval")
col_z   <- get_detected_col("zscore")
col_n   <- get_detected_col("samplesize")

required_raw_cols <- c(col_snp, col_chr, col_pos, col_ea, col_oa, col_p, col_z, col_n)
required_raw_cols <- required_raw_cols[!is.na(required_raw_cols)]

# ---------- AF 文件路径 ----------
af_candidates <- c(
  file.path(eqtlgen_dir, "2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz"),
  file.path(eqtlgen_dir, "39_eqtlgen_AF_fixed_source.txt.gz"),
  file.path(eqtlgen_dir, "36_eqtlgen_allele_frequencies.txt.gz")
)

af_file <- af_candidates[file.exists(af_candidates)][1]
if (length(af_file) == 0 || is.na(af_file)) {
  stop("未找到 AF 文件。请先确认 39 的下载结果。")
}

# =========================================================
# A. 导出 AF 文件表头与前 20 行预览
# =========================================================
af_header <- names(fread(af_file, nrows = 0))
af_preview <- fread(af_file, nrows = 20)

fwrite(
  data.frame(column_name = af_header, stringsAsFactors = FALSE),
  file.path(out_registry, "40_af_full_header.tsv"),
  sep = "\t"
)

fwrite(
  af_preview,
  file.path(out_registry, "40_af_preview_20rows.tsv"),
  sep = "\t"
)

# =========================================================
# B. 手动指定 AF 字段语义（已根据你的预览结果固定）
# =========================================================
freq_col_name <- "AlleleB_all"
freq_is_for   <- "AlleleB"

alleleA_col_name <- if ("AlleleA" %in% af_header) "AlleleA" else NA_character_
alleleB_col_name <- if ("AlleleB" %in% af_header) "AlleleB" else NA_character_

# 修正后的 chr/pos 识别
af_chr_col <- if ("hg19_chr" %in% af_header) {
  "hg19_chr"
} else if ("chr" %in% af_header) {
  "chr"
} else if ("CHR" %in% af_header) {
  "CHR"
} else {
  NA_character_
}

af_pos_col <- if ("hg19_pos" %in% af_header) {
  "hg19_pos"
} else if ("pos" %in% af_header) {
  "pos"
} else if ("BP" %in% af_header) {
  "BP"
} else {
  NA_character_
}

manual_af_config <- data.frame(
  item = c("af_file", "freq_col_name", "freq_is_for", "alleleA_col_name", "alleleB_col_name", "af_chr_col", "af_pos_col"),
  value = c(af_file, freq_col_name, freq_is_for, alleleA_col_name, alleleB_col_name, af_chr_col, af_pos_col),
  stringsAsFactors = FALSE
)

fwrite(
  manual_af_config,
  file.path(out_registry, "40_manual_af_config.tsv"),
  sep = "\t"
)

manual_config_ready <- !is.na(freq_col_name) && freq_col_name %in% af_header &&
  !is.na(freq_is_for) && freq_is_for %in% c("AlleleA", "AlleleB")

if (!manual_config_ready) {
  status_summary <- data.frame(
    item = c(
      "af_file_exists",
      "header_exported",
      "preview_exported",
      "manual_config_ready",
      "recommended_next_action"
    ),
    value = c(
      file.exists(af_file),
      TRUE,
      TRUE,
      FALSE,
      "Check 40_manual_af_config.tsv and ensure freq_col_name / freq_is_for are valid"
    ),
    stringsAsFactors = FALSE
  )
  
  fwrite(
    status_summary,
    file.path(out_registry, "40_status_summary.tsv"),
    sep = "\t"
  )
  
  writeLines(
    c(
      "40 code stopped after AF preview export.",
      paste0("Created at: ", Sys.time()),
      "",
      "Manual AF config is not ready."
    ),
    file.path(out_logs, "40_manual_af_parse_notes.txt")
  )
  
  stop("manual_config_ready = FALSE，请检查 AF 配置。")
}

# =========================================================
# C. 读取 AF 所需列
# =========================================================
af_select <- unique(na.omit(c(
  "SNP",
  freq_col_name,
  alleleA_col_name,
  alleleB_col_name,
  af_chr_col,
  af_pos_col
)))

af_dt <- fread(af_file, select = af_select)

# =========================================================
# D. 工具函数
# =========================================================
reconstruct_beta_se <- function(z, n, eaf) {
  valid <- !is.na(z) & !is.na(n) & !is.na(eaf) & eaf > 0 & eaf < 1 & n > 0
  se <- rep(NA_real_, length(z))
  beta <- rep(NA_real_, length(z))
  
  se[valid]   <- 1 / sqrt(2 * eaf[valid] * (1 - eaf[valid]) * (n[valid] + z[valid]^2))
  beta[valid] <- z[valid] * se[valid]
  
  list(beta = beta, se = se, valid = valid)
}

align_eaf_manual <- function(raw_dt, af_dt,
                             raw_snp, raw_chr, raw_pos, raw_ea, raw_oa,
                             freq_col_name, freq_is_for,
                             alleleA_col_name = NA, alleleB_col_name = NA,
                             af_chr_col = NA, af_pos_col = NA) {
  x <- copy(raw_dt)
  af2 <- copy(af_dt)
  
  # 优先 SNP 合并
  if ("SNP" %in% names(af2) && raw_snp %in% names(x)) {
    setnames(x, raw_snp, ".__join_snp__", skip_absent = TRUE)
    setnames(af2, "SNP", ".__join_snp__", skip_absent = TRUE)
    x <- merge(x, af2, by = ".__join_snp__", all.x = TRUE)
    setnames(x, ".__join_snp__", raw_snp, skip_absent = TRUE)
  } else if (!is.na(af_chr_col) && !is.na(af_pos_col)) {
    setnames(x, raw_chr, ".__join_chr__", skip_absent = TRUE)
    setnames(x, raw_pos, ".__join_pos__", skip_absent = TRUE)
    setnames(af2, af_chr_col, ".__join_chr__", skip_absent = TRUE)
    setnames(af2, af_pos_col, ".__join_pos__", skip_absent = TRUE)
    x <- merge(x, af2, by = c(".__join_chr__", ".__join_pos__"), all.x = TRUE)
    setnames(x, ".__join_chr__", raw_chr, skip_absent = TRUE)
    setnames(x, ".__join_pos__", raw_pos, skip_absent = TRUE)
  } else {
    x$eaf_reconstructed <- NA_real_
    x$eaf_match_mode <- "no_join_key"
    return(x)
  }
  
  freq_vec   <- suppressWarnings(as.numeric(x[[freq_col_name]]))
  raw_ea_vec <- as.character(x[[raw_ea]])
  raw_oa_vec <- as.character(x[[raw_oa]])
  
  eaf_vec  <- rep(NA_real_, nrow(x))
  mode_vec <- rep("unresolved", nrow(x))
  
  if (!is.na(alleleA_col_name) && !is.na(alleleB_col_name) &&
      alleleA_col_name %in% names(x) && alleleB_col_name %in% names(x)) {
    
    a_vec <- as.character(x[[alleleA_col_name]])
    b_vec <- as.character(x[[alleleB_col_name]])
    
    if (freq_is_for == "AlleleB") {
      ea_is_b <- !is.na(raw_ea_vec) & !is.na(b_vec) & raw_ea_vec == b_vec
      oa_is_b <- !is.na(raw_oa_vec) & !is.na(b_vec) & raw_oa_vec == b_vec
      ea_is_a <- !is.na(raw_ea_vec) & !is.na(a_vec) & raw_ea_vec == a_vec
      oa_is_a <- !is.na(raw_oa_vec) & !is.na(a_vec) & raw_oa_vec == a_vec
      
      eaf_vec[ea_is_b] <- freq_vec[ea_is_b]
      mode_vec[ea_is_b] <- "effect_allele_is_AlleleB_keep"
      
      eaf_vec[oa_is_b] <- 1 - freq_vec[oa_is_b]
      mode_vec[oa_is_b] <- "other_allele_is_AlleleB_invert"
      
      unresolved1 <- is.na(eaf_vec) & ea_is_a
      eaf_vec[unresolved1] <- 1 - freq_vec[unresolved1]
      mode_vec[unresolved1] <- "effect_allele_is_AlleleA_invert"
      
      unresolved2 <- is.na(eaf_vec) & oa_is_a
      eaf_vec[unresolved2] <- freq_vec[unresolved2]
      mode_vec[unresolved2] <- "other_allele_is_AlleleA_keep"
    }
    
    if (freq_is_for == "AlleleA") {
      ea_is_a <- !is.na(raw_ea_vec) & !is.na(a_vec) & raw_ea_vec == a_vec
      oa_is_a <- !is.na(raw_oa_vec) & !is.na(a_vec) & raw_oa_vec == a_vec
      ea_is_b <- !is.na(raw_ea_vec) & !is.na(b_vec) & raw_ea_vec == b_vec
      oa_is_b <- !is.na(raw_oa_vec) & !is.na(b_vec) & raw_oa_vec == b_vec
      
      eaf_vec[ea_is_a] <- freq_vec[ea_is_a]
      mode_vec[ea_is_a] <- "effect_allele_is_AlleleA_keep"
      
      eaf_vec[oa_is_a] <- 1 - freq_vec[oa_is_a]
      mode_vec[oa_is_a] <- "other_allele_is_AlleleA_invert"
      
      unresolved1 <- is.na(eaf_vec) & ea_is_b
      eaf_vec[unresolved1] <- 1 - freq_vec[unresolved1]
      mode_vec[unresolved1] <- "effect_allele_is_AlleleB_invert"
      
      unresolved2 <- is.na(eaf_vec) & oa_is_b
      eaf_vec[unresolved2] <- freq_vec[unresolved2]
      mode_vec[unresolved2] <- "other_allele_is_AlleleB_keep"
    }
  }
  
  x$eaf_reconstructed <- eaf_vec
  x$eaf_match_mode <- mode_vec
  x
}

# =========================================================
# E. 重跑 EAF 合并与 beta/se 重建
# =========================================================
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
  
  x2 <- align_eaf_manual(
    raw_dt = x,
    af_dt = af_dt,
    raw_snp = col_snp,
    raw_chr = col_chr,
    raw_pos = col_pos,
    raw_ea = col_ea,
    raw_oa = col_oa,
    freq_col_name = freq_col_name,
    freq_is_for = freq_is_for,
    alleleA_col_name = alleleA_col_name,
    alleleB_col_name = alleleB_col_name,
    af_chr_col = af_chr_col,
    af_pos_col = af_pos_col
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
  ) %>%
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
  
  complete_idx <- with(
    std,
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
    file.path(cand_dir, paste0(g, "_40_with_manual_eaf.tsv")),
    sep = "\t",
    na = "NA"
  )
  
  fwrite(
    std,
    file.path(cand_dir, paste0(g, "_40_candidate_exposure_all_rows.tsv")),
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
    gene_audit$reason <- "official_exposure_written_after_manual_af_parse"
  } else {
    gene_audit$reason <- "no_complete_rows_after_manual_af_parse"
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
    "manual_config_ready",
    "n_genes_total",
    "n_genes_written_official_exposure",
    "recommended_next_action"
  ),
  value = c(
    file.exists(af_file),
    TRUE,
    length(genes),
    sum(gene_audit_tab$written_official_exposure, na.rm = TRUE),
    ifelse(
      sum(gene_audit_tab$written_official_exposure, na.rm = TRUE) > 0,
      "Proceed to 41: audit written exposure files and extract instruments",
      "Still zero: inspect 40_gene_reconstruction_audit.tsv and 40_gene_row_qc_summary.tsv"
    )
  ),
  stringsAsFactors = FALSE
)

fwrite(
  gene_audit_tab,
  file.path(out_registry, "40_gene_reconstruction_audit.tsv"),
  sep = "\t"
)

fwrite(
  row_audit_tab,
  file.path(out_registry, "40_gene_row_qc_summary.tsv"),
  sep = "\t"
)

fwrite(
  status_summary,
  file.path(out_registry, "40_status_summary.tsv"),
  sep = "\t"
)

writeLines(
  c(
    "40 manual AF parsing notes",
    paste0("Created at: ", Sys.time()),
    "",
    "[Manual config used]",
    paste0("freq_col_name = ", freq_col_name),
    paste0("freq_is_for   = ", freq_is_for),
    paste0("af_chr_col    = ", af_chr_col),
    paste0("af_pos_col    = ", af_pos_col),
    "",
    "[Outputs]",
    "- 40_af_full_header.tsv",
    "- 40_af_preview_20rows.tsv",
    "- 40_manual_af_config.tsv",
    "- 40_gene_reconstruction_audit.tsv",
    "- 40_gene_row_qc_summary.tsv",
    "- 40_status_summary.tsv"
  ),
  file.path(out_logs, "40_manual_af_parse_notes.txt")
)

cat("40_manual_parse_af_frequency_column_and_rerun_eaf_merge.R finished successfully.\n")
cat("Generated files:\n")
cat("- 00_registry/40_af_full_header.tsv\n")
cat("- 00_registry/40_af_preview_20rows.tsv\n")
cat("- 00_registry/40_manual_af_config.tsv\n")
cat("- 00_registry/40_gene_reconstruction_audit.tsv\n")
cat("- 00_registry/40_gene_row_qc_summary.tsv\n")
cat("- 00_registry/40_status_summary.tsv\n")
cat("- 02_processed_data/37_candidate_gene_tables/<GENE>_40_with_manual_eaf.tsv\n")
cat("- 02_processed_data/37_candidate_gene_tables/<GENE>_40_candidate_exposure_all_rows.tsv\n")
cat("- 01_raw_data/exposure/eQTLGen/<GENE>_cis_eqtl.tsv (if successful)\n")