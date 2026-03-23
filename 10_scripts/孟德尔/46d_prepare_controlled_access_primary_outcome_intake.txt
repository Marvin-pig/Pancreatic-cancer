# =========================
# 46d_prepare_controlled_access_primary_outcome_intake.R
# 目的：
# 1. 固定真正的 PanScan / PanC4 primary outcome 数据来源（dbGaP 受控）
# 2. 生成申请与追踪清单
# 3. 为手动获取后的真实 summary stats 文件做 intake / 标准化 / 结构检查
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "05_MR"

registry_dir    <- file.path(root_dir, "00_registry")
outcome_dir     <- file.path(root_dir, "01_raw_data/outcome/PanScan_PanC4")
processed_dir   <- file.path(root_dir, "02_processed_data")
log_dir         <- file.path(root_dir, "06_logs")

dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(outcome_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# A. 固定真正 primary outcome 来源
# =========================================================
primary_source_registry <- data.frame(
  item = c(
    "target_publication",
    "pubmed_id",
    "doi",
    "target_trait",
    "target_ancestry",
    "preferred_local_file",
    "dbgap_panscan_accession",
    "dbgap_panc4_accession",
    "preferred_route",
    "fallback_route"
  ),
  value = c(
    "Genome-wide meta-analysis identifies five new susceptibility loci for pancreatic cancer",
    "29422604",
    "10.1038/s41467-018-02942-5",
    "pancreatic cancer risk",
    "European",
    "05_MR/01_raw_data/outcome/PanScan_PanC4/pancreatic_cancer_sumstats.tsv",
    "phs000206.v5.p3",
    "phs000648.v1.p1",
    "controlled_access_dbGaP",
    "BBJ_bbj-a-140"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  primary_source_registry,
  file.path(registry_dir, "46d_primary_outcome_true_source_registry.tsv"),
  sep = "\t"
)

# =========================================================
# B. 生成受控获取追踪表
# =========================================================
access_tracking <- data.frame(
  task_id = c("46D_A", "46D_B", "46D_C", "46D_D", "46D_E", "46D_F"),
  task = c(
    "Confirm institutional route for dbGaP controlled-access request",
    "Track access request for phs000206.v5.p3 (PanScan)",
    "Track access request for phs000648.v1.p1 (PanC4)",
    "Record whether full summary statistics or per-study result tables were obtained",
    "Place obtained file into local intake path",
    "Run schema normalization and structural validation"
  ),
  status = c("TODO", "TODO", "TODO", "TODO", "TODO", "TODO"),
  note = c(
    "May require PI/institutional data access approval",
    "PanScan controlled data",
    "PanC4 controlled data",
    "Need summary-level GWAS result table usable for MR",
    "Manual file placement required",
    "Use this script after file arrival"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  access_tracking,
  file.path(registry_dir, "46d_controlled_access_tracking.tsv"),
  sep = "\t"
)

# =========================================================
# C. 你获得真实文件后，把路径填这里
# =========================================================
manual_true_outcome_file <- ""
# 示例：
# manual_true_outcome_file <- "/Users/你的用户名/Downloads/PanScan_PanC4_pancreatic_cancer_sumstats.txt.gz"

target_outcome_file <- file.path(outcome_dir, "pancreatic_cancer_sumstats.tsv")

# =========================================================
# D. outcome 结构与列映射工具
# =========================================================
required_cols <- c("SNP", "effect_allele", "other_allele", "beta", "se", "pval", "eaf", "outcome")

detect_one <- function(cands, nms) {
  hit <- cands[cands %in% nms]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

standardize_outcome_columns <- function(x) {
  nms <- names(x)
  
  col_snp   <- detect_one(c("SNP", "rsid", "rsID", "MarkerName", "ID", "variant_id"), nms)
  col_ea    <- detect_one(c("effect_allele", "EA", "A1", "Allele2", "ALT", "EffectAllele"), nms)
  col_oa    <- detect_one(c("other_allele", "NEA", "A2", "Allele1", "REF", "OtherAllele"), nms)
  col_beta  <- detect_one(c("beta", "BETA", "Effect", "effect", "Beta"), nms)
  col_se    <- detect_one(c("se", "SE", "StdErr", "stderr", "Std_Error"), nms)
  col_p     <- detect_one(c("pval", "P", "P-value", "p_value", "PVALUE", "Pvalue"), nms)
  col_eaf   <- detect_one(c("eaf", "EAF", "Freq1", "AF_Allele2", "EffectAlleleFreq", "effect_allele_frequency"), nms)
  col_chr   <- detect_one(c("chr", "CHR", "Chromosome"), nms)
  col_pos   <- detect_one(c("pos", "BP", "Position", "base_pair_location"), nms)
  col_ncase <- detect_one(c("ncase", "N_CASE", "Cases"), nms)
  col_nctrl <- detect_one(c("ncontrol", "N_CONTROL", "Controls"), nms)
  col_n     <- detect_one(c("samplesize", "N", "TotalN"), nms)
  
  map_tab <- data.frame(
    target_field = c("SNP","effect_allele","other_allele","beta","se","pval","eaf","chr","pos","ncase","ncontrol","samplesize"),
    detected_column = c(col_snp,col_ea,col_oa,col_beta,col_se,col_p,col_eaf,col_chr,col_pos,col_ncase,col_nctrl,col_n),
    stringsAsFactors = FALSE
  )
  
  y <- copy(x)
  
  if (!is.na(col_snp))   y[, SNP := as.character(get(col_snp))]
  if (!is.na(col_ea))    y[, effect_allele := toupper(as.character(get(col_ea)))]
  if (!is.na(col_oa))    y[, other_allele  := toupper(as.character(get(col_oa)))]
  if (!is.na(col_beta))  y[, beta := suppressWarnings(as.numeric(get(col_beta)))]
  if (!is.na(col_se))    y[, se := suppressWarnings(as.numeric(get(col_se)))]
  if (!is.na(col_p))     y[, pval := suppressWarnings(as.numeric(get(col_p)))]
  if (!is.na(col_eaf))   y[, eaf := suppressWarnings(as.numeric(get(col_eaf)))]
  if (!is.na(col_chr))   y[, chr := suppressWarnings(as.integer(get(col_chr)))]
  if (!is.na(col_pos))   y[, pos := suppressWarnings(as.integer(get(col_pos)))]
  if (!is.na(col_ncase)) y[, ncase := suppressWarnings(as.numeric(get(col_ncase)))]
  if (!is.na(col_nctrl)) y[, ncontrol := suppressWarnings(as.numeric(get(col_nctrl)))]
  if (!is.na(col_n))     y[, samplesize := suppressWarnings(as.numeric(get(col_n)))]
  
  if (!"outcome" %in% names(y)) {
    y[, outcome := "Pancreatic_cancer_PanScanPanC4"]
  }
  
  list(data = y, map = map_tab)
}

check_outcome_structure <- function(x) {
  missing_cols <- setdiff(required_cols, names(x))
  valid_rows <- rep(TRUE, nrow(x))
  
  if ("SNP" %in% names(x)) valid_rows <- valid_rows & !is.na(x$SNP) & x$SNP != ""
  if ("effect_allele" %in% names(x)) valid_rows <- valid_rows & x$effect_allele %in% c("A","C","G","T")
  if ("other_allele" %in% names(x)) valid_rows <- valid_rows & x$other_allele %in% c("A","C","G","T")
  if ("beta" %in% names(x)) valid_rows <- valid_rows & !is.na(x$beta)
  if ("se" %in% names(x)) valid_rows <- valid_rows & !is.na(x$se) & x$se > 0
  if ("pval" %in% names(x)) valid_rows <- valid_rows & !is.na(x$pval) & x$pval > 0 & x$pval <= 1
  if ("eaf" %in% names(x)) valid_rows <- valid_rows & !is.na(x$eaf) & x$eaf > 0 & x$eaf < 1
  if ("outcome" %in% names(x)) valid_rows <- valid_rows & !is.na(x$outcome) & x$outcome != ""
  
  list(
    required_columns_ok = length(missing_cols) == 0,
    missing_columns = paste(missing_cols, collapse = ","),
    n_valid_rows = sum(valid_rows, na.rm = TRUE)
  )
}

# =========================================================
# E. 若真实文件已获得，执行 intake
# =========================================================
intake_status <- data.frame(
  item = c(
    "manual_true_outcome_file_provided",
    "manual_true_outcome_file_exists",
    "target_written",
    "required_columns_ok",
    "n_rows_raw",
    "n_valid_rows"
  ),
  value = c(FALSE, FALSE, FALSE, FALSE, NA, NA),
  stringsAsFactors = FALSE
)

if (nzchar(manual_true_outcome_file)) {
  intake_status$value[intake_status$item == "manual_true_outcome_file_provided"] <- TRUE
  intake_status$value[intake_status$item == "manual_true_outcome_file_exists"] <- file.exists(manual_true_outcome_file)
  
  if (!file.exists(manual_true_outcome_file)) {
    warning("手动指定的真实 outcome 文件不存在：", manual_true_outcome_file)
  } else {
    raw_dat <- tryCatch(fread(manual_true_outcome_file), error = function(e) NULL)
    
    if (is.null(raw_dat)) {
      stop("无法读取真实 outcome 文件，请检查格式：", manual_true_outcome_file)
    }
    
    fwrite(
      data.frame(original_colnames = names(raw_dat), stringsAsFactors = FALSE),
      file.path(registry_dir, "46d_true_outcome_original_colnames.tsv"),
      sep = "\t"
    )
    
    fwrite(
      raw_dat[1:min(20, nrow(raw_dat)), ],
      file.path(registry_dir, "46d_true_outcome_head_preview.tsv"),
      sep = "\t"
    )
    
    std <- standardize_outcome_columns(raw_dat)
    std_dat <- std$data
    map_tab <- std$map
    
    fwrite(
      map_tab,
      file.path(registry_dir, "46d_true_outcome_column_map.tsv"),
      sep = "\t"
    )
    
    qc <- check_outcome_structure(std_dat)
    
    fwrite(
      data.frame(
        item = c("required_columns_ok", "missing_columns", "n_rows_raw", "n_valid_rows"),
        value = c(qc$required_columns_ok, qc$missing_columns, nrow(std_dat), qc$n_valid_rows),
        stringsAsFactors = FALSE
      ),
      file.path(registry_dir, "46d_true_outcome_attachment_check.tsv"),
      sep = "\t"
    )
    
    intake_status$value[intake_status$item == "n_rows_raw"] <- nrow(std_dat)
    intake_status$value[intake_status$item == "required_columns_ok"] <- qc$required_columns_ok
    intake_status$value[intake_status$item == "n_valid_rows"] <- qc$n_valid_rows
    
    if (isTRUE(qc$required_columns_ok) && qc$n_valid_rows > 0) {
      fwrite(std_dat, target_outcome_file, sep = "\t", na = "NA")
      intake_status$value[intake_status$item == "target_written"] <- file.exists(target_outcome_file)
    }
  }
}

fwrite(
  intake_status,
  file.path(registry_dir, "46d_true_outcome_intake_status.tsv"),
  sep = "\t"
)

# =========================================================
# F. 日志
# =========================================================
notes <- c(
  "46d true primary outcome pursuit notes",
  paste0("Created at: ", Sys.time()),
  "",
  "[True target source]",
  "PanScan = phs000206.v5.p3",
  "PanC4   = phs000648.v1.p1",
  "Publication PMID = 29422604",
  "",
  "[Preferred local target]",
  "05_MR/01_raw_data/outcome/PanScan_PanC4/pancreatic_cancer_sumstats.tsv",
  "",
  "[Required columns]",
  "SNP, effect_allele, other_allele, beta, se, pval, eaf, outcome",
  "",
  "[Rule]",
  "Do not proceed to 47 unless the true outcome file is structurally valid and written to target path."
)

writeLines(
  notes,
  file.path(log_dir, "46d_true_primary_outcome_pursuit_notes.txt")
)

cat("46d_prepare_controlled_access_primary_outcome_intake.R finished successfully.\n")
cat("Generated files:\n")
cat("- 00_registry/46d_primary_outcome_true_source_registry.tsv\n")
cat("- 00_registry/46d_controlled_access_tracking.tsv\n")
cat("- 00_registry/46d_true_outcome_intake_status.tsv\n")
cat("- 06_logs/46d_true_primary_outcome_pursuit_notes.txt\n")
cat("If manual_true_outcome_file is provided and readable, additional files will be generated:\n")
cat("- 00_registry/46d_true_outcome_original_colnames.tsv\n")
cat("- 00_registry/46d_true_outcome_head_preview.tsv\n")
cat("- 00_registry/46d_true_outcome_column_map.tsv\n")
cat("- 00_registry/46d_true_outcome_attachment_check.tsv\n")
cat("- 01_raw_data/outcome/PanScan_PanC4/pancreatic_cancer_sumstats.tsv\n")