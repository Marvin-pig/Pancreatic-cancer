Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiI2MTk1MDE5OTBAcXEuY29tIiwiaWF0IjoxNzczOTc4NjM5LCJleHAiOjE3NzUxODgyMzl9.exWVTcTBoD5wCJz2mwm0ikv9rmHkzPOjhauRkA5aj_Jbl0pQpavaXK7oRWTOYhQ7f0lQs1v-c8LCWyFwNtlJbvFimMksWeuPjusBig410hR5IePfLXL44gbYpdENhQAwaNVuK6ELFGj-lSH1jPSWiNxpgP6ZwiIVf3Epwm6ubx5f8lggrVHpZw-N7DzXHyg1IfSE0HV3YKGo-CLp24wtF7XlzpQ-v6Wls83pQXemKcMoBcwGFqkxE5uG1Kmo5YCHA0YrLXG_H0tYh34Nn2Z0RUZlt0M8bUhKP5_jMGJk397T3mJV_vtEQqf-4nUnbdOh_QNvw4gcyddEmobKTvJbxQ")

# =========================
# 42b_ld_clump_with_jwt_and_prepare_bbj_snp_list.R
# 目的：
# 1. 检查 OPENGWAS_JWT 是否存在
# 2. 显式使用 JWT 调用 OpenGWAS clumping API
# 3. 重跑 strict / relax LD clumping
# 4. 生成后续 BBJ outcome 抓取所需 SNP 清单
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(TwoSampleMR)
  library(ieugwasr)
  library(httr)
  library(jsonlite)
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

strict_file <- file.path(result_dir, "41_instrument_candidates_strict_all_genes.tsv")
relax_file  <- file.path(result_dir, "41_instrument_candidates_relax_all_genes.tsv")

if (!file.exists(strict_file)) stop("缺少 strict candidates 文件: ", strict_file)
if (!file.exists(relax_file))  stop("缺少 relax candidates 文件: ", relax_file)

strict_dt <- fread(strict_file)
relax_dt  <- fread(relax_file)

# ---------- JWT ----------
jwt <- Sys.getenv("OPENGWAS_JWT", unset = "")
if (jwt == "") {
  stop(
    "未检测到 OPENGWAS_JWT。\n",
    "请先运行：\n",
    'Sys.setenv(OPENGWAS_JWT = "你的新JWT")\n',
    "然后重跑 42b。"
  )
}

# ---------- API 连通性测试 ----------
api_test <- tryCatch({
  r <- httr::GET(
    url = "https://api.opengwas.io/api/user",
    httr::add_headers(Authorization = paste("Bearer", jwt))
  )
  txt <- httr::content(r, as = "text", encoding = "UTF-8")
  list(
    status_code = httr::status_code(r),
    ok = httr::status_code(r) == 200,
    response_text = txt
  )
}, error = function(e) {
  list(
    status_code = NA_integer_,
    ok = FALSE,
    response_text = conditionMessage(e)
  )
})

fwrite(
  data.frame(
    item = c("api_user_status_code", "api_user_ok"),
    value = c(api_test$status_code, api_test$ok),
    stringsAsFactors = FALSE
  ),
  file.path(registry_dir, "42b_api_auth_check.tsv"),
  sep = "\t", na = "NA"
)

if (!isTRUE(api_test$ok)) {
  writeLines(
    c(
      "OpenGWAS JWT auth check failed.",
      paste0("Created at: ", Sys.time()),
      paste0("status_code: ", api_test$status_code),
      paste0("response: ", api_test$response_text)
    ),
    file.path(log_dir, "42b_api_auth_check_error.txt")
  )
  stop("JWT 认证失败，请确认 token 是否已重置并正确写入 OPENGWAS_JWT。")
}

# ---------- 参数 ----------
clump_pop_primary <- "EUR"
clump_kb  <- 10000
clump_r2  <- 0.001
clump_p   <- 1
min_primary_iv_n <- 3

required_cols <- c("gene", "SNP", "pval")
stopifnot(all(required_cols %in% names(strict_dt)))
stopifnot(all(required_cols %in% names(relax_dt)))

# ---------- 工具函数 ----------
prep_for_clump <- function(dt) {
  x <- copy(dt)
  x <- x %>%
    mutate(
      rsid = as.character(SNP),
      pval = as.numeric(pval),
      id = as.character(gene)
    ) %>%
    filter(
      !is.na(rsid), rsid != "",
      !is.na(pval), pval > 0, pval <= 1,
      !is.na(id), id != ""
    ) %>%
    distinct(id, rsid, .keep_all = TRUE)
  x
}

run_clump_safe <- function(dat, jwt_token, pop = "EUR", kb = 10000, r2 = 0.001, p = 1) {
  out <- NULL
  msg <- ""
  ok <- FALSE
  
  tryCatch({
    out <- ieugwasr::ld_clump(
      d = dat[, c("rsid", "pval", "id")],
      clump_kb = kb,
      clump_r2 = r2,
      clump_p = p,
      pop = pop,
      opengwas_jwt = jwt_token
    )
    ok <- !is.null(out)
  }, error = function(e) {
    msg <<- conditionMessage(e)
  })
  
  list(ok = ok, out = out, msg = msg)
}

clump_one_gene <- function(gene_name, dt_source, source_label, jwt_token, pop = "EUR") {
  sub <- dt_source %>% filter(gene == gene_name)
  
  pre_n <- nrow(sub)
  if (pre_n == 0) {
    return(list(
      gene = gene_name,
      source = source_label,
      pre_n = 0L,
      post_n = 0L,
      ok = FALSE,
      msg = "no_candidates",
      clumped = NULL
    ))
  }
  
  clump_input <- prep_for_clump(sub)
  if (nrow(clump_input) == 0) {
    return(list(
      gene = gene_name,
      source = source_label,
      pre_n = pre_n,
      post_n = 0L,
      ok = FALSE,
      msg = "no_valid_clump_input",
      clumped = NULL
    ))
  }
  
  res <- run_clump_safe(
    dat = clump_input,
    jwt_token = jwt_token,
    pop = pop,
    kb = clump_kb,
    r2 = clump_r2,
    p = clump_p
  )
  
  if (!res$ok || is.null(res$out) || nrow(res$out) == 0) {
    return(list(
      gene = gene_name,
      source = source_label,
      pre_n = pre_n,
      post_n = 0L,
      ok = FALSE,
      msg = ifelse(res$msg == "", "clumping_returned_empty", res$msg),
      clumped = NULL
    ))
  }
  
  keep_rsid <- unique(as.character(res$out$rsid))
  clumped_full <- sub %>%
    filter(SNP %in% keep_rsid) %>%
    distinct(SNP, .keep_all = TRUE)
  
  list(
    gene = gene_name,
    source = source_label,
    pre_n = pre_n,
    post_n = nrow(clumped_full),
    ok = TRUE,
    msg = "",
    clumped = clumped_full
  )
}

# ---------- 主流程 ----------
genes <- sort(unique(c(as.character(strict_dt$gene), as.character(relax_dt$gene))))

strict_summary_list <- list()
relax_summary_list  <- list()
primary_choice_list <- list()

strict_clumped_list <- list()
relax_clumped_list  <- list()
final_primary_list  <- list()

for (g in genes) {
  # 1) strict 主 clumping
  strict_res <- clump_one_gene(
    gene_name = g,
    dt_source = strict_dt,
    source_label = "strict",
    jwt_token = jwt,
    pop = clump_pop_primary
  )
  
  strict_summary_list[[g]] <- data.frame(
    gene = strict_res$gene,
    source = strict_res$source,
    pre_n = strict_res$pre_n,
    post_n = strict_res$post_n,
    ok = strict_res$ok,
    message = strict_res$msg,
    stringsAsFactors = FALSE
  )
  
  if (!is.null(strict_res$clumped) && nrow(strict_res$clumped) > 0) {
    strict_clumped_list[[g]] <- strict_res$clumped
  }
  
  # 2) relax 备用 clumping
  relax_res <- clump_one_gene(
    gene_name = g,
    dt_source = relax_dt,
    source_label = "relax",
    jwt_token = jwt,
    pop = clump_pop_primary
  )
  
  relax_summary_list[[g]] <- data.frame(
    gene = relax_res$gene,
    source = relax_res$source,
    pre_n = relax_res$pre_n,
    post_n = relax_res$post_n,
    ok = relax_res$ok,
    message = relax_res$msg,
    stringsAsFactors = FALSE
  )
  
  if (!is.null(relax_res$clumped) && nrow(relax_res$clumped) > 0) {
    relax_clumped_list[[g]] <- relax_res$clumped
  }
  
  # 3) 选择 primary
  chosen_source <- NA_character_
  chosen_dt <- NULL
  chosen_n <- 0L
  
  if (!is.null(strict_res$clumped) && nrow(strict_res$clumped) > 0) {
    chosen_source <- "strict_primary"
    chosen_dt <- strict_res$clumped
    chosen_n <- nrow(chosen_dt)
  } else if (!is.null(relax_res$clumped) && nrow(relax_res$clumped) > 0) {
    chosen_source <- "relax_backup"
    chosen_dt <- relax_res$clumped
    chosen_n <- nrow(chosen_dt)
  } else {
    chosen_source <- "none"
    chosen_n <- 0L
  }
  
  primary_choice_list[[g]] <- data.frame(
    gene = g,
    chosen_source = chosen_source,
    chosen_n = chosen_n,
    stringsAsFactors = FALSE
  )
  
  if (!is.null(chosen_dt) && nrow(chosen_dt) > 0) {
    final_primary_list[[g]] <- chosen_dt %>%
      mutate(primary_source = chosen_source)
  }
}

strict_summary_tab <- bind_rows(strict_summary_list)
relax_summary_tab  <- bind_rows(relax_summary_list)
primary_choice_tab <- bind_rows(primary_choice_list)

strict_clumped_tab <- bind_rows(strict_clumped_list)
relax_clumped_tab  <- bind_rows(relax_clumped_list)
final_primary_tab  <- bind_rows(final_primary_list)

# ---------- 输出每基因 clumped 文件 ----------
if (nrow(strict_clumped_tab) > 0) {
  for (g in unique(strict_clumped_tab$gene)) {
    fwrite(
      strict_clumped_tab %>% filter(gene == g),
      file.path(proc_dir, paste0("42b_", g, "_clumped_strict.tsv")),
      sep = "\t", na = "NA"
    )
  }
}

if (nrow(relax_clumped_tab) > 0) {
  for (g in unique(relax_clumped_tab$gene)) {
    fwrite(
      relax_clumped_tab %>% filter(gene == g),
      file.path(proc_dir, paste0("42b_", g, "_clumped_relax.tsv")),
      sep = "\t", na = "NA"
    )
  }
}

if (nrow(final_primary_tab) > 0) {
  for (g in unique(final_primary_tab$gene)) {
    fwrite(
      final_primary_tab %>% filter(gene == g),
      file.path(proc_dir, paste0("42b_", g, "_clumped_primary_for_outcome.tsv")),
      sep = "\t", na = "NA"
    )
  }
}

# ---------- 输出总表 ----------
fwrite(
  strict_summary_tab,
  file.path(registry_dir, "42b_clumping_summary_strict.tsv"),
  sep = "\t", na = "NA"
)

fwrite(
  relax_summary_tab,
  file.path(registry_dir, "42b_clumping_summary_relax.tsv"),
  sep = "\t", na = "NA"
)

fwrite(
  primary_choice_tab,
  file.path(registry_dir, "42b_primary_instrument_choice.tsv"),
  sep = "\t", na = "NA"
)

if (nrow(strict_clumped_tab) > 0) {
  fwrite(
    strict_clumped_tab,
    file.path(result_dir, "42b_clumped_instruments_strict_all_genes.tsv"),
    sep = "\t", na = "NA"
  )
  
  fwrite(
    strict_clumped_tab %>% distinct(SNP),
    file.path(result_dir, "42b_bbj_outcome_snp_list_strict_clumped.tsv"),
    sep = "\t", na = "NA"
  )
}

if (nrow(relax_clumped_tab) > 0) {
  fwrite(
    relax_clumped_tab,
    file.path(result_dir, "42b_clumped_instruments_relax_all_genes.tsv"),
    sep = "\t", na = "NA"
  )
  
  fwrite(
    relax_clumped_tab %>% distinct(SNP),
    file.path(result_dir, "42b_bbj_outcome_snp_list_relax_clumped.tsv"),
    sep = "\t", na = "NA"
  )
}

if (nrow(final_primary_tab) > 0) {
  fwrite(
    final_primary_tab,
    file.path(result_dir, "42b_clumped_instruments_primary_all_genes.tsv"),
    sep = "\t", na = "NA"
  )
  
  fwrite(
    final_primary_tab %>%
      distinct(SNP) %>%
      select(SNP),
    file.path(result_dir, "42b_bbj_outcome_snp_list_primary.tsv"),
    sep = "\t", na = "NA"
  )
}

# ---------- 总状态 ----------
n_primary_genes <- sum(primary_choice_tab$chosen_n > 0, na.rm = TRUE)
n_primary_ge3   <- sum(primary_choice_tab$chosen_n >= min_primary_iv_n, na.rm = TRUE)

status_summary <- data.frame(
  item = c(
    "api_user_ok",
    "clump_pop_primary",
    "n_genes_total",
    "n_genes_with_primary_instruments",
    "n_genes_with_primary_iv_ge_3",
    "recommended_next_action"
  ),
  value = c(
    api_test$ok,
    clump_pop_primary,
    length(genes),
    n_primary_genes,
    n_primary_ge3,
    ifelse(
      n_primary_genes > 0,
      "Proceed to 43: extract BBJ outcome data and harmonise exposure-outcome pairs",
      "JWT auth passed but clumping still empty: inspect 42b_clumping_summary_strict.tsv / relax.tsv"
    )
  ),
  stringsAsFactors = FALSE
)

fwrite(
  status_summary,
  file.path(registry_dir, "42b_status_summary.tsv"),
  sep = "\t", na = "NA"
)

# ---------- 日志 ----------
log_lines <- c(
  "=== Code 42b Summary ===",
  paste0("Created at: ", Sys.time()),
  "",
  "[Goal]",
  "Re-run LD clumping with explicit OpenGWAS JWT authentication.",
  "",
  "[Parameters]",
  paste0("pop = ", clump_pop_primary),
  paste0("clump_kb = ", clump_kb),
  paste0("clump_r2 = ", clump_r2),
  paste0("min_primary_iv_n = ", min_primary_iv_n),
  "",
  "[API auth check]",
  paste0("status_code = ", api_test$status_code),
  paste0("ok = ", api_test$ok),
  "",
  "[Outputs]",
  "- 42b_api_auth_check.tsv",
  "- 42b_clumping_summary_strict.tsv",
  "- 42b_clumping_summary_relax.tsv",
  "- 42b_primary_instrument_choice.tsv",
  "- 42b_status_summary.tsv",
  "- 42b_bbj_outcome_snp_list_primary.tsv"
)

writeLines(
  log_lines,
  file.path(log_dir, "42b_ld_clump_notes.txt")
)

cat("42b_ld_clump_with_jwt_and_prepare_bbj_snp_list.R finished successfully.\n")
cat("Generated files:\n")
cat("- 00_registry/42b_api_auth_check.tsv\n")
cat("- 00_registry/42b_clumping_summary_strict.tsv\n")
cat("- 00_registry/42b_clumping_summary_relax.tsv\n")
cat("- 00_registry/42b_primary_instrument_choice.tsv\n")
cat("- 00_registry/42b_status_summary.tsv\n")
cat("- 02_processed_data/42b_<GENE>_clumped_strict.tsv\n")
cat("- 02_processed_data/42b_<GENE>_clumped_relax.tsv\n")
cat("- 02_processed_data/42b_<GENE>_clumped_primary_for_outcome.tsv\n")
cat("- 03_results/42b_clumped_instruments_primary_all_genes.tsv\n")
cat("- 03_results/42b_bbj_outcome_snp_list_primary.tsv\n")


