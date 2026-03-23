# =========================
# 41_audit_exposure_files_and_prepare_instruments.R
# 目的：
# 1. 读取 40 成功写出的正式 exposure 文件
# 2. 检查列结构是否符合模板
# 3. 进行 instrument 提取前 QC
# 4. 按 p 值阈值筛选候选 SNP
# 5. 生成每个基因的 instrument 审计表与统一 SNP 清单
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ---------- 路径 ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "05_MR"

batch1_file   <- file.path(root_dir, "02_processed_data/batch1_queries/29_batch1_priority_genes.tsv")
audit40_file  <- file.path(root_dir, "00_registry/40_gene_reconstruction_audit.tsv")
eqtlgen_dir   <- file.path(root_dir, "01_raw_data/exposure/eQTLGen")
proc_dir      <- file.path(root_dir, "02_processed_data")
result_dir    <- file.path(root_dir, "03_results")
registry_dir  <- file.path(root_dir, "00_registry")
log_dir       <- file.path(root_dir, "06_logs")

dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(batch1_file)) stop("缺少 batch1_file: ", batch1_file)
if (!file.exists(audit40_file)) stop("缺少 40_gene_reconstruction_audit.tsv: ", audit40_file)

batch1 <- fread(batch1_file)
audit40 <- fread(audit40_file)

if (!"gene" %in% names(batch1)) stop("29_batch1_priority_genes.tsv 缺少 gene 列。")
if (!"gene" %in% names(audit40)) stop("40_gene_reconstruction_audit.tsv 缺少 gene 列。")

genes_all <- unique(as.character(batch1$gene))

# ---------- 仅保留 40 已成功写出 exposure 的基因 ----------
success_genes <- audit40 %>%
  filter(written_official_exposure == TRUE) %>%
  pull(gene) %>%
  as.character() %>%
  unique()

if (length(success_genes) == 0) {
  stop("40 结果显示没有任何成功写出的正式 exposure 文件，不能进入 41。")
}

# ---------- 参数 ----------
# 主阈值与宽松阈值都输出，便于你后续比较
p_threshold_strict <- 5e-8
p_threshold_relax  <- 1e-5

# 经验性 F 统计量
# 对单 SNP 可近似用 (beta/se)^2
calc_F_stat <- function(beta, se) {
  out <- rep(NA_real_, length(beta))
  valid <- !is.na(beta) & !is.na(se) & se > 0
  out[valid] <- (beta[valid] / se[valid])^2
  out
}

# 回文 SNP 判断
is_palindromic <- function(a1, a2) {
  a1 <- toupper(a1)
  a2 <- toupper(a2)
  (a1 == "A" & a2 == "T") |
    (a1 == "T" & a2 == "A") |
    (a1 == "C" & a2 == "G") |
    (a1 == "G" & a2 == "C")
}

# 中间频率区间的回文 SNP 更危险
is_ambiguous_palindromic <- function(a1, a2, eaf) {
  pal <- is_palindromic(a1, a2)
  mid <- !is.na(eaf) & eaf > 0.42 & eaf < 0.58
  pal & mid
}

required_cols <- c(
  "SNP", "chr", "pos", "effect_allele", "other_allele",
  "beta", "se", "pval", "eaf", "gene", "exposure"
)

schema_audit_list <- list()
row_qc_list <- list()
instrument_summary_list <- list()
strict_snp_list <- list()
relax_snp_list <- list()

for (g in success_genes) {
  f <- file.path(eqtlgen_dir, paste0(g, "_cis_eqtl.tsv"))
  
  schema_row <- data.frame(
    gene = g,
    exposure_file = f,
    file_exists = file.exists(f),
    n_rows = NA_integer_,
    all_required_cols_present = FALSE,
    missing_cols = "",
    stringsAsFactors = FALSE
  )
  
  if (!file.exists(f)) {
    schema_row$missing_cols <- "file_missing"
    schema_audit_list[[g]] <- schema_row
    next
  }
  
  x <- fread(f)
  schema_row$n_rows <- nrow(x)
  
  miss <- setdiff(required_cols, names(x))
  schema_row$all_required_cols_present <- length(miss) == 0
  schema_row$missing_cols <- ifelse(length(miss) == 0, "", paste(miss, collapse = ","))
  
  schema_audit_list[[g]] <- schema_row
  
  if (length(miss) > 0) next
  
  # ---------- 基础 QC ----------
  x <- x %>%
    mutate(
      SNP = as.character(SNP),
      chr = suppressWarnings(as.integer(chr)),
      pos = suppressWarnings(as.integer(pos)),
      effect_allele = toupper(as.character(effect_allele)),
      other_allele  = toupper(as.character(other_allele)),
      beta = suppressWarnings(as.numeric(beta)),
      se   = suppressWarnings(as.numeric(se)),
      pval = suppressWarnings(as.numeric(pval)),
      eaf  = suppressWarnings(as.numeric(eaf)),
      gene = as.character(gene),
      exposure = as.character(exposure)
    ) %>%
    filter(
      !is.na(SNP),
      SNP != "",
      !is.na(chr),
      !is.na(pos),
      !is.na(effect_allele),
      !is.na(other_allele),
      effect_allele %in% c("A", "C", "G", "T"),
      other_allele %in% c("A", "C", "G", "T"),
      effect_allele != other_allele,
      !is.na(beta),
      !is.na(se),
      se > 0,
      !is.na(pval),
      pval > 0,
      pval <= 1,
      !is.na(eaf),
      eaf > 0,
      eaf < 1
    ) %>%
    distinct(SNP, .keep_all = TRUE) %>%
    mutate(
      F_stat = calc_F_stat(beta, se),
      palindromic = is_palindromic(effect_allele, other_allele),
      ambiguous_palindromic = is_ambiguous_palindromic(effect_allele, other_allele, eaf)
    )
  
  # 保存 QC 后的 exposure
  qc_file <- file.path(proc_dir, paste0("41_", g, "_exposure_qc.tsv"))
  fwrite(x, qc_file, sep = "\t", na = "NA")
  
  # ---------- 行级 QC 汇总 ----------
  row_qc <- data.frame(
    gene = g,
    n_after_basic_qc = nrow(x),
    n_palindromic = sum(x$palindromic, na.rm = TRUE),
    n_ambiguous_palindromic = sum(x$ambiguous_palindromic, na.rm = TRUE),
    n_F_gt_10 = sum(!is.na(x$F_stat) & x$F_stat > 10),
    median_F = median(x$F_stat, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  row_qc_list[[g]] <- row_qc
  
  # ---------- 候选 instruments ----------
  # 先去掉 ambiguous palindromic，再按阈值筛选
  x_clean <- x %>%
    filter(!ambiguous_palindromic)
  
  inst_strict <- x_clean %>%
    filter(pval < p_threshold_strict)
  
  inst_relax <- x_clean %>%
    filter(pval < p_threshold_relax)
  
  # 保存每个基因的 strict/relax 候选表
  fwrite(
    inst_strict,
    file.path(proc_dir, paste0("41_", g, "_instrument_candidates_strict.tsv")),
    sep = "\t", na = "NA"
  )
  
  fwrite(
    inst_relax,
    file.path(proc_dir, paste0("41_", g, "_instrument_candidates_relax.tsv")),
    sep = "\t", na = "NA"
  )
  
  # 汇总
  inst_summary <- data.frame(
    gene = g,
    n_qc_rows = nrow(x),
    n_after_remove_ambiguous_palindromic = nrow(x_clean),
    n_strict_p_lt_5e8 = nrow(inst_strict),
    n_relax_p_lt_1e5 = nrow(inst_relax),
    n_strict_F_gt_10 = sum(!is.na(inst_strict$F_stat) & inst_strict$F_stat > 10),
    n_relax_F_gt_10 = sum(!is.na(inst_relax$F_stat) & inst_relax$F_stat > 10),
    best_p = ifelse(nrow(x_clean) > 0, min(x_clean$pval, na.rm = TRUE), NA_real_),
    min_F = ifelse(nrow(x_clean) > 0, min(x_clean$F_stat, na.rm = TRUE), NA_real_),
    median_F = ifelse(nrow(x_clean) > 0, median(x_clean$F_stat, na.rm = TRUE), NA_real_),
    stringsAsFactors = FALSE
  )
  instrument_summary_list[[g]] <- inst_summary
  
  if (nrow(inst_strict) > 0) {
    strict_snp_list[[g]] <- inst_strict %>%
      transmute(
        gene, SNP, chr, pos, effect_allele, other_allele, beta, se, pval, eaf, F_stat, exposure
      )
  }
  
  if (nrow(inst_relax) > 0) {
    relax_snp_list[[g]] <- inst_relax %>%
      transmute(
        gene, SNP, chr, pos, effect_allele, other_allele, beta, se, pval, eaf, F_stat, exposure
      )
  }
}

schema_audit_tab <- bind_rows(schema_audit_list)
row_qc_tab <- bind_rows(row_qc_list)
instrument_summary_tab <- bind_rows(instrument_summary_list)
strict_snp_tab <- bind_rows(strict_snp_list)
relax_snp_tab <- bind_rows(relax_snp_list)

# ---------- 输出统一 SNP 清单 ----------
if (nrow(strict_snp_tab) > 0) {
  fwrite(
    strict_snp_tab,
    file.path(result_dir, "41_instrument_candidates_strict_all_genes.tsv"),
    sep = "\t", na = "NA"
  )
  
  fwrite(
    strict_snp_tab %>% distinct(SNP),
    file.path(result_dir, "41_bbj_outcome_snp_list_strict.tsv"),
    sep = "\t", na = "NA"
  )
}

if (nrow(relax_snp_tab) > 0) {
  fwrite(
    relax_snp_tab,
    file.path(result_dir, "41_instrument_candidates_relax_all_genes.tsv"),
    sep = "\t", na = "NA"
  )
  
  fwrite(
    relax_snp_tab %>% distinct(SNP),
    file.path(result_dir, "41_bbj_outcome_snp_list_relax.tsv"),
    sep = "\t", na = "NA"
  )
}

# ---------- 总状态 ----------
status_summary <- data.frame(
  item = c(
    "n_success_exposure_files",
    "n_schema_pass_files",
    "n_genes_with_strict_candidates",
    "n_genes_with_relax_candidates",
    "recommended_next_action"
  ),
  value = c(
    length(success_genes),
    sum(schema_audit_tab$all_required_cols_present, na.rm = TRUE),
    sum(instrument_summary_tab$n_strict_p_lt_5e8 > 0, na.rm = TRUE),
    sum(instrument_summary_tab$n_relax_p_lt_1e5 > 0, na.rm = TRUE),
    ifelse(
      sum(instrument_summary_tab$n_relax_p_lt_1e5 > 0, na.rm = TRUE) > 0,
      "Proceed to 42: LD clumping for candidate instruments and prepare BBJ outcome extraction",
      "No usable candidates after QC: inspect 41_instrument_summary.tsv and reconsider thresholds"
    )
  ),
  stringsAsFactors = FALSE
)

# ---------- 保存表 ----------
fwrite(
  schema_audit_tab,
  file.path(registry_dir, "41_exposure_schema_audit.tsv"),
  sep = "\t", na = "NA"
)

fwrite(
  row_qc_tab,
  file.path(registry_dir, "41_exposure_row_qc_summary.tsv"),
  sep = "\t", na = "NA"
)

fwrite(
  instrument_summary_tab,
  file.path(registry_dir, "41_instrument_summary.tsv"),
  sep = "\t", na = "NA"
)

fwrite(
  status_summary,
  file.path(registry_dir, "41_status_summary.tsv"),
  sep = "\t", na = "NA"
)

# ---------- 日志 ----------
log_lines <- c(
  "=== Code 41 Summary ===",
  paste0("Created at: ", Sys.time()),
  "",
  "[Goal]",
  "Audit official exposure files and prepare candidate instruments before clumping.",
  "",
  "[Thresholds]",
  paste0("Strict threshold: p < ", format(p_threshold_strict, scientific = TRUE)),
  paste0("Relax threshold : p < ", format(p_threshold_relax, scientific = TRUE)),
  "",
  "[Outputs]",
  "- 41_exposure_schema_audit.tsv",
  "- 41_exposure_row_qc_summary.tsv",
  "- 41_instrument_summary.tsv",
  "- 41_status_summary.tsv",
  "- 41_<GENE>_exposure_qc.tsv",
  "- 41_<GENE>_instrument_candidates_strict.tsv",
  "- 41_<GENE>_instrument_candidates_relax.tsv",
  "- 41_bbj_outcome_snp_list_strict.tsv",
  "- 41_bbj_outcome_snp_list_relax.tsv"
)

writeLines(
  log_lines,
  file.path(log_dir, "41_audit_exposure_and_prepare_instruments_notes.txt")
)

cat("41_audit_exposure_files_and_prepare_instruments.R finished successfully.\n")
cat("Generated files:\n")
cat("- 00_registry/41_exposure_schema_audit.tsv\n")
cat("- 00_registry/41_exposure_row_qc_summary.tsv\n")
cat("- 00_registry/41_instrument_summary.tsv\n")
cat("- 00_registry/41_status_summary.tsv\n")
cat("- 02_processed_data/41_<GENE>_exposure_qc.tsv\n")
cat("- 02_processed_data/41_<GENE>_instrument_candidates_strict.tsv\n")
cat("- 02_processed_data/41_<GENE>_instrument_candidates_relax.tsv\n")
cat("- 03_results/41_instrument_candidates_strict_all_genes.tsv\n")
cat("- 03_results/41_instrument_candidates_relax_all_genes.tsv\n")
cat("- 03_results/41_bbj_outcome_snp_list_strict.tsv\n")
cat("- 03_results/41_bbj_outcome_snp_list_relax.tsv\n")