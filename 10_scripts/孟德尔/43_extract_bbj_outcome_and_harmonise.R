Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiI2MTk1MDE5OTBAcXEuY29tIiwiaWF0IjoxNzczOTc4NjM5LCJleHAiOjE3NzUxODgyMzl9.exWVTcTBoD5wCJz2mwm0ikv9rmHkzPOjhauRkA5aj_Jbl0pQpavaXK7oRWTOYhQ7f0lQs1v-c8LCWyFwNtlJbvFimMksWeuPjusBig410hR5IePfLXL44gbYpdENhQAwaNVuK6ELFGj-lSH1jPSWiNxpgP6ZwiIVf3Epwm6ubx5f8lggrVHpZw-N7DzXHyg1IfSE0HV3YKGo-CLp24wtF7XlzpQ-v6Wls83pQXemKcMoBcwGFqkxE5uG1Kmo5YCHA0YrLXG_H0tYh34Nn2Z0RUZlt0M8bUhKP5_jMGJk397T3mJV_vtEQqf-4nUnbdOh_QNvw4gcyddEmobKTvJbxQ")

# =========================
# 43_extract_bbj_outcome_and_harmonise.R
# 目的：
# 1. 读取 42b 的 primary clumped instruments
# 2. 抓取 BBJ pancreatic cancer outcome 数据
# 3. 完成 exposure-outcome harmonization
# 4. 输出每个基因的 harmonized 数据集与总审计表
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(TwoSampleMR)
  library(ieugwasr)
  library(httr)
})

# ---------- 路径 ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "05_MR"

result_dir   <- file.path(root_dir, "03_results")
registry_dir <- file.path(root_dir, "00_registry")
proc_dir     <- file.path(root_dir, "02_processed_data")
log_dir      <- file.path(root_dir, "06_logs")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

primary_file <- file.path(result_dir, "42b_clumped_instruments_primary_all_genes.tsv")
choice_file  <- file.path(registry_dir, "42b_primary_instrument_choice.tsv")

if (!file.exists(primary_file)) stop("缺少 primary instruments 文件: ", primary_file)
if (!file.exists(choice_file))  stop("缺少 primary choice 文件: ", choice_file)

primary_dt <- fread(primary_file)
choice_dt  <- fread(choice_file)

# ---------- JWT ----------
jwt <- Sys.getenv("OPENGWAS_JWT", unset = "")
if (jwt == "") {
  stop(
    "未检测到 OPENGWAS_JWT。\n",
    "请先运行：\n",
    'Sys.setenv(OPENGWAS_JWT = "你的JWT")\n',
    "然后重跑 43。"
  )
}

# ---------- API 连通性测试 ----------
api_test <- tryCatch({
  r <- httr::GET(
    url = "https://api.opengwas.io/api/user",
    httr::add_headers(Authorization = paste("Bearer", jwt))
  )
  list(
    status_code = httr::status_code(r),
    ok = httr::status_code(r) == 200
  )
}, error = function(e) {
  list(status_code = NA_integer_, ok = FALSE)
})

fwrite(
  data.frame(
    item = c("api_user_status_code", "api_user_ok"),
    value = c(api_test$status_code, api_test$ok),
    stringsAsFactors = FALSE
  ),
  file.path(registry_dir, "43_api_auth_check.tsv"),
  sep = "\t", na = "NA"
)

if (!isTRUE(api_test$ok)) {
  stop("JWT 认证失败，无法进入 43。")
}

# ---------- outcome 设置 ----------
# 这里沿用你当前项目路线中的 BBJ pancreatic cancer outcome
bbj_outcome_id <- "bbj-a-140"

# extract_outcome_data 参数
use_proxies      <- TRUE
proxy_rsq        <- 0.8
align_alleles    <- 1
allow_palindrome <- 1
maf_threshold    <- 0.3
split_size       <- 10000
proxy_split_size <- 500

# harmonise_data 参数
# action = 2 是常用稳健设置：
# 尝试处理回文 SNP，并在必要时丢弃无法判断方向的位点
harmonise_action <- 2

# ---------- 预检查 ----------
required_exp_cols <- c(
  "gene", "SNP", "beta", "se", "effect_allele", "other_allele", "eaf", "pval", "exposure"
)
miss_exp <- setdiff(required_exp_cols, names(primary_dt))
if (length(miss_exp) > 0) {
  stop("42b primary 文件缺少必要列: ", paste(miss_exp, collapse = ", "))
}

# 仅保留真正有 primary instruments 的基因
genes <- choice_dt %>%
  filter(chosen_n > 0, chosen_source != "none") %>%
  pull(gene) %>%
  as.character() %>%
  unique()

if (length(genes) == 0) {
  stop("42b 结果没有任何可用 primary instruments，不能进入 43。")
}

# ---------- 格式化 exposure ----------
# TwoSampleMR 期望 exposure 列名以 .exposure 结尾
exp_dat <- primary_dt %>%
  filter(gene %in% genes) %>%
  mutate(
    SNP = as.character(SNP),
    beta.exposure = as.numeric(beta),
    se.exposure = as.numeric(se),
    effect_allele.exposure = as.character(effect_allele),
    other_allele.exposure = as.character(other_allele),
    eaf.exposure = as.numeric(eaf),
    pval.exposure = as.numeric(pval),
    exposure = as.character(exposure),
    id.exposure = as.character(gene),
    samplesize.exposure = NA_real_,
    chr.exposure = if ("chr" %in% names(.)) as.integer(chr) else NA_integer_,
    pos.exposure = if ("pos" %in% names(.)) as.integer(pos) else NA_integer_,
    mr_keep.exposure = TRUE,
    data_source.exposure = "eQTLGen_reconstructed"
  ) %>%
  select(
    SNP,
    beta.exposure,
    se.exposure,
    effect_allele.exposure,
    other_allele.exposure,
    eaf.exposure,
    pval.exposure,
    exposure,
    id.exposure,
    samplesize.exposure,
    chr.exposure,
    pos.exposure,
    mr_keep.exposure,
    data_source.exposure
  ) %>%
  distinct(id.exposure, SNP, .keep_all = TRUE)

fwrite(
  exp_dat,
  file.path(result_dir, "43_primary_exposure_for_harmonisation.tsv"),
  sep = "\t", na = "NA"
)

# ---------- 抓取 outcome ----------
all_snps <- unique(exp_dat$SNP)

out_dat <- tryCatch({
  TwoSampleMR::extract_outcome_data(
    snps = all_snps,
    outcomes = bbj_outcome_id,
    proxies = use_proxies,
    rsq = proxy_rsq,
    align_alleles = align_alleles,
    palindromes = allow_palindrome,
    maf_threshold = maf_threshold,
    opengwas_jwt = jwt,
    splitsize = split_size,
    proxy_splitsize = proxy_split_size
  )
}, error = function(e) {
  msg <- conditionMessage(e)
  writeLines(
    c(
      "43 outcome extraction failed",
      paste0("Created at: ", Sys.time()),
      paste0("Error: ", msg)
    ),
    file.path(log_dir, "43_outcome_extraction_error.txt")
  )
  stop("BBJ outcome 抓取失败: ", msg)
})

if (is.null(out_dat) || nrow(out_dat) == 0) {
  stop("BBJ outcome 抓取返回空结果。")
}

fwrite(
  out_dat,
  file.path(result_dir, "43_bbj_outcome_raw.tsv"),
  sep = "\t", na = "NA"
)

# ---------- outcome 抓取审计 ----------
outcome_audit <- data.frame(
  item = c(
    "bbj_outcome_id",
    "n_requested_snps",
    "n_returned_rows",
    "n_unique_returned_snps"
  ),
  value = c(
    bbj_outcome_id,
    length(all_snps),
    nrow(out_dat),
    length(unique(out_dat$SNP))
  ),
  stringsAsFactors = FALSE
)

fwrite(
  outcome_audit,
  file.path(registry_dir, "43_outcome_extraction_audit.tsv"),
  sep = "\t", na = "NA"
)

# ---------- 按基因 harmonise ----------
harmonised_list <- list()
harmonise_summary_list <- list()

for (g in genes) {
  exp_g <- exp_dat %>% filter(id.exposure == g)
  
  # 只取该基因 instrument 对应的 outcome 记录
  out_g <- out_dat %>% filter(SNP %in% exp_g$SNP)
  
  summary_row <- data.frame(
    gene = g,
    n_exp_snps = nrow(exp_g),
    n_outcome_rows = nrow(out_g),
    n_harmonised_rows = 0L,
    n_mr_keep = 0L,
    outcome_id = bbj_outcome_id,
    stringsAsFactors = FALSE
  )
  
  if (nrow(exp_g) == 0 || nrow(out_g) == 0) {
    harmonise_summary_list[[g]] <- summary_row
    next
  }
  
  harm_g <- tryCatch({
    TwoSampleMR::harmonise_data(
      exposure_dat = exp_g,
      outcome_dat = out_g,
      action = harmonise_action
    )
  }, error = function(e) {
    msg <- conditionMessage(e)
    writeLines(
      c(
        paste0("Harmonisation failed for gene: ", g),
        paste0("Created at: ", Sys.time()),
        paste0("Error: ", msg)
      ),
      file.path(log_dir, paste0("43_harmonise_error_", g, ".txt"))
    )
    NULL
  })
  
  if (is.null(harm_g) || nrow(harm_g) == 0) {
    harmonise_summary_list[[g]] <- summary_row
    next
  }
  
  summary_row$n_harmonised_rows <- nrow(harm_g)
  summary_row$n_mr_keep <- sum(harm_g$mr_keep, na.rm = TRUE)
  
  harmonised_list[[g]] <- harm_g
  harmonise_summary_list[[g]] <- summary_row
  
  fwrite(
    harm_g,
    file.path(proc_dir, paste0("43_", g, "_harmonised.tsv")),
    sep = "\t", na = "NA"
  )
  
  fwrite(
    harm_g %>% filter(mr_keep),
    file.path(proc_dir, paste0("43_", g, "_harmonised_mr_keep.tsv")),
    sep = "\t", na = "NA"
  )
}

harmonised_all <- bind_rows(harmonised_list)
harmonise_summary_tab <- bind_rows(harmonise_summary_list)

if (nrow(harmonised_all) > 0) {
  fwrite(
    harmonised_all,
    file.path(result_dir, "43_harmonised_all_genes.tsv"),
    sep = "\t", na = "NA"
  )
  
  fwrite(
    harmonised_all %>% filter(mr_keep),
    file.path(result_dir, "43_harmonised_mr_keep_all_genes.tsv"),
    sep = "\t", na = "NA"
  )
}

# ---------- 基因分层建议 ----------
gene_tier <- harmonise_summary_tab %>%
  mutate(
    analysis_tier = case_when(
      n_mr_keep >= 3 ~ "multi_snp_priority",
      n_mr_keep == 2 ~ "few_snp_exploratory",
      n_mr_keep == 1 ~ "single_snp_wald_ratio",
      TRUE ~ "insufficient_after_harmonisation"
    )
  )

fwrite(
  gene_tier,
  file.path(registry_dir, "43_gene_analysis_tier.tsv"),
  sep = "\t", na = "NA"
)

# ---------- 总状态 ----------
status_summary <- data.frame(
  item = c(
    "bbj_outcome_id",
    "n_genes_with_primary_instruments",
    "n_genes_with_harmonised_data",
    "n_genes_with_mr_keep_ge_3",
    "recommended_next_action"
  ),
  value = c(
    bbj_outcome_id,
    length(genes),
    sum(harmonise_summary_tab$n_harmonised_rows > 0, na.rm = TRUE),
    sum(harmonise_summary_tab$n_mr_keep >= 3, na.rm = TRUE),
    ifelse(
      sum(harmonise_summary_tab$n_mr_keep > 0, na.rm = TRUE) > 0,
      "Proceed to 44: run MR main analysis and sensitivity analyses",
      "No usable harmonised pairs: inspect 43_outcome_extraction_audit.tsv and per-gene harmonised files"
    )
  ),
  stringsAsFactors = FALSE
)

fwrite(
  status_summary,
  file.path(registry_dir, "43_status_summary.tsv"),
  sep = "\t", na = "NA"
)

# ---------- 日志 ----------
log_lines <- c(
  "=== Code 43 Summary ===",
  paste0("Created at: ", Sys.time()),
  "",
  "[Goal]",
  "Extract BBJ outcome data for clumped primary instruments and harmonise exposure-outcome pairs.",
  "",
  "[Outcome settings]",
  paste0("bbj_outcome_id = ", bbj_outcome_id),
  paste0("proxies = ", use_proxies),
  paste0("proxy_rsq = ", proxy_rsq),
  paste0("harmonise_action = ", harmonise_action),
  "",
  "[Outputs]",
  "- 43_api_auth_check.tsv",
  "- 43_primary_exposure_for_harmonisation.tsv",
  "- 43_bbj_outcome_raw.tsv",
  "- 43_outcome_extraction_audit.tsv",
  "- 43_<GENE>_harmonised.tsv",
  "- 43_<GENE>_harmonised_mr_keep.tsv",
  "- 43_harmonised_all_genes.tsv",
  "- 43_harmonised_mr_keep_all_genes.tsv",
  "- 43_gene_analysis_tier.tsv",
  "- 43_status_summary.tsv"
)

writeLines(
  log_lines,
  file.path(log_dir, "43_extract_outcome_and_harmonise_notes.txt")
)

cat("43_extract_bbj_outcome_and_harmonise.R finished successfully.\n")
cat("Generated files:\n")
cat("- 00_registry/43_api_auth_check.tsv\n")
cat("- 00_registry/43_outcome_extraction_audit.tsv\n")
cat("- 00_registry/43_gene_analysis_tier.tsv\n")
cat("- 00_registry/43_status_summary.tsv\n")
cat("- 02_processed_data/43_<GENE>_harmonised.tsv\n")
cat("- 02_processed_data/43_<GENE>_harmonised_mr_keep.tsv\n")
cat("- 03_results/43_primary_exposure_for_harmonisation.tsv\n")
cat("- 03_results/43_bbj_outcome_raw.tsv\n")
cat("- 03_results/43_harmonised_all_genes.tsv\n")
cat("- 03_results/43_harmonised_mr_keep_all_genes.tsv\n")