# =========================
# 46b_find_public_primary_outcome_via_gwas_catalog_api.R
# 目的：
# 1. 通过 GWAS Catalog public summary statistics API 搜索 pancreatic cancer 相关 trait
# 2. 列出这些 trait 下带 summary statistics 的 studies
# 3. 用你当前 primary instrument SNP 列表测试各 study 的实际 overlap
# 4. 选择最像 PanScan/PanC4 primary outcome 的 study
# 5. 生成可直接进入 47 的 outcome 表（仅覆盖当前 primary SNP）
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(httr)
  library(jsonlite)
  library(purrr)
})

# ---------- 路径 ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "05_MR"

result_dir   <- file.path(root_dir, "03_results")
registry_dir <- file.path(root_dir, "00_registry")
proc_dir     <- file.path(root_dir, "02_processed_data")
raw_outcome_dir <- file.path(root_dir, "01_raw_data/outcome/PanScan_PanC4")
log_dir      <- file.path(root_dir, "06_logs")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_outcome_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# 你当前已经有的 primary SNP 列表（42b）
primary_snp_file <- file.path(result_dir, "42b_bbj_outcome_snp_list_primary.tsv")
if (!file.exists(primary_snp_file)) {
  stop("缺少 42b primary SNP 列表: ", primary_snp_file)
}

primary_snps <- fread(primary_snp_file)
if (!"SNP" %in% names(primary_snps)) {
  stop("42b_bbj_outcome_snp_list_primary.tsv 缺少 SNP 列。")
}
query_snps <- unique(as.character(primary_snps$SNP))
query_snps <- query_snps[!is.na(query_snps) & query_snps != ""]

if (length(query_snps) == 0) {
  stop("当前 primary SNP 列表为空。")
}

# ---------- 参数 ----------
# 目标关键词：尽量贴近 PanScan / PanC4 pancreatic cancer risk GWAS
target_keywords <- c("pancreatic", "pancreas", "pancreatic cancer", "panscan", "panc4")
target_pmids <- c("29422604")  # 2018 Nat Commun PanScan/PanC4 meta-analysis

max_trait_pages <- 100
page_size <- 200

# ---------- API endpoint ----------
base_api <- "https://www.ebi.ac.uk/gwas/summary-statistics/api"
base_labs <- "https://www.ebi.ac.uk/gwas/labs/rest/api"

# ---------- 工具函数 ----------
safe_get_json <- function(url, pause_sec = 0.1) {
  Sys.sleep(pause_sec)
  out <- tryCatch({
    r <- httr::GET(url, httr::timeout(60))
    txt <- httr::content(r, as = "text", encoding = "UTF-8")
    list(
      ok = httr::status_code(r) == 200,
      status_code = httr::status_code(r),
      text = txt,
      json = if (httr::status_code(r) == 200) jsonlite::fromJSON(txt, simplifyVector = FALSE) else NULL
    )
  }, error = function(e) {
    list(ok = FALSE, status_code = NA_integer_, text = conditionMessage(e), json = NULL)
  })
  out
}

extract_embedded_df <- function(x, embedded_name = NULL) {
  if (is.null(x) || is.null(x$`_embedded`)) return(data.frame())
  emb <- x$`_embedded`
  
  if (!is.null(embedded_name) && embedded_name %in% names(emb)) {
    y <- emb[[embedded_name]]
  } else {
    y <- emb[[1]]
  }
  
  if (is.null(y)) return(data.frame())
  
  # y 可能是 list of records
  if (is.data.frame(y)) return(as.data.frame(y, stringsAsFactors = FALSE))
  
  if (is.list(y)) {
    rows <- lapply(y, function(z) {
      if (is.null(z)) return(NULL)
      
      # 尽量扁平化常见字段
      flat <- list()
      nm <- names(z)
      for (k in nm) {
        v <- z[[k]]
        if (is.atomic(v) && length(v) <= 1) {
          flat[[k]] <- ifelse(length(v) == 0, NA, as.character(v))
        } else if (is.list(v) && "_links" != k) {
          # 少量一层扁平
          if (!is.null(names(v)) && all(lengths(v) == 1)) {
            for (kk in names(v)) {
              flat[[paste0(k, ".", kk)]] <- as.character(v[[kk]])
            }
          }
        }
      }
      as.data.frame(flat, stringsAsFactors = FALSE)
    })
    rows <- rows[!sapply(rows, is.null)]
    if (length(rows) == 0) return(data.frame())
    return(bind_rows(rows))
  }
  
  data.frame()
}

follow_next_link <- function(json_obj) {
  if (is.null(json_obj) || is.null(json_obj[["_links"]]) || is.null(json_obj[["_links"]][["next"]])) {
    return(NA_character_)
  }
  nxt <- json_obj[["_links"]][["next"]][["href"]]
  if (is.null(nxt) || length(nxt) == 0) return(NA_character_)
  as.character(nxt)
}

# 从 summary stats API 列 traits
get_all_traits <- function(max_pages = 100, size = 200) {
  url <- paste0(base_api, "/traits?size=", size)
  all_rows <- list()
  i <- 1
  
  while (!is.na(url) && nzchar(url) && i <= max_pages) {
    res <- safe_get_json(url)
    if (!res$ok || is.null(res$json)) break
    
    df <- extract_embedded_df(res$json, embedded_name = "traits")
    if (nrow(df) > 0) {
      all_rows[[length(all_rows) + 1]] <- df
    }
    
    url <- follow_next_link(res$json)
    i <- i + 1
  }
  
  bind_rows(all_rows)
}

# 给 trait 列 studies
get_trait_studies <- function(trait_id, size = 200) {
  url <- paste0(base_api, "/traits/", trait_id, "/studies?size=", size)
  res <- safe_get_json(url)
  if (!res$ok || is.null(res$json)) return(data.frame())
  df <- extract_embedded_df(res$json, embedded_name = "studies")
  df
}

# 补 study metadata（如果 labs api 可用）
get_study_metadata <- function(study_acc) {
  url <- paste0(base_labs, "/studies/", study_acc)
  res <- safe_get_json(url)
  if (!res$ok || is.null(res$json)) {
    return(data.frame(
      study_accession = study_acc,
      title = NA_character_,
      first_author = NA_character_,
      pubmed_id = NA_character_,
      reported_trait = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  
  j <- res$json
  data.frame(
    study_accession = study_acc,
    title = if (!is.null(j$title)) as.character(j$title) else NA_character_,
    first_author = if (!is.null(j$firstAuthor$fullname)) as.character(j$firstAuthor$fullname) else if (!is.null(j$author)) as.character(j$author) else NA_character_,
    pubmed_id = if (!is.null(j$publicationInfo$pubmedId)) as.character(j$publicationInfo$pubmedId) else if (!is.null(j$pubmedId)) as.character(j$pubmedId) else NA_character_,
    reported_trait = if (!is.null(j$reportedTrait)) as.character(j$reportedTrait) else NA_character_,
    stringsAsFactors = FALSE
  )
}

# 查询单个 SNP 在指定 study 的 association
get_variant_assoc_for_study <- function(rsid, study_acc) {
  url <- paste0(
    base_api, "/associations/", rsid,
    "?study_accession=", study_acc,
    "&reveal=all"
  )
  res <- safe_get_json(url, pause_sec = 0.05)
  if (!res$ok || is.null(res$json)) return(data.frame())
  
  df <- extract_embedded_df(res$json, embedded_name = "associations")
  if (nrow(df) == 0) return(data.frame())
  
  # 只保留目标 study
  if ("study_accession" %in% names(df)) {
    df <- df %>% filter(study_accession == study_acc)
  }
  df
}

# 从 association 资源标准化 outcome
standardise_outcome <- function(df, outcome_label) {
  if (nrow(df) == 0) return(data.frame())
  
  out <- df %>%
    transmute(
      SNP = if ("variant_id" %in% names(.)) as.character(variant_id) else NA_character_,
      effect_allele = if ("effect_allele" %in% names(.)) as.character(effect_allele) else NA_character_,
      other_allele = if ("other_allele" %in% names(.)) as.character(other_allele) else NA_character_,
      beta_raw = if ("beta" %in% names(.)) suppressWarnings(as.numeric(beta)) else NA_real_,
      se_raw = if ("se" %in% names(.)) suppressWarnings(as.numeric(se)) else NA_real_,
      or_raw = if ("odds_ratio" %in% names(.)) suppressWarnings(as.numeric(odds_ratio)) else NA_real_,
      ci_lower_raw = if ("ci_lower" %in% names(.)) suppressWarnings(as.numeric(ci_lower)) else NA_real_,
      ci_upper_raw = if ("ci_upper" %in% names(.)) suppressWarnings(as.numeric(ci_upper)) else NA_real_,
      pval = if ("p_value" %in% names(.)) suppressWarnings(as.numeric(p_value)) else NA_real_,
      eaf = if ("effect_allele_frequency" %in% names(.)) suppressWarnings(as.numeric(effect_allele_frequency)) else NA_real_,
      chr = if ("chromosome" %in% names(.)) suppressWarnings(as.integer(chromosome)) else NA_integer_,
      pos = if ("base_pair_location" %in% names(.)) suppressWarnings(as.integer(base_pair_location)) else NA_integer_,
      study_accession = if ("study_accession" %in% names(.)) as.character(study_accession) else NA_character_,
      trait = if ("trait" %in% names(.)) as.character(trait) else NA_character_,
      outcome = outcome_label
    ) %>%
    mutate(
      beta = case_when(
        !is.na(beta_raw) ~ beta_raw,
        is.na(beta_raw) & !is.na(or_raw) & or_raw > 0 ~ log(or_raw),
        TRUE ~ NA_real_
      ),
      se = case_when(
        !is.na(se_raw) ~ se_raw,
        is.na(se_raw) & !is.na(or_raw) & !is.na(ci_lower_raw) & !is.na(ci_upper_raw) &
          or_raw > 0 & ci_lower_raw > 0 & ci_upper_raw > 0 ~
          (log(ci_upper_raw) - log(ci_lower_raw)) / (2 * 1.96),
        TRUE ~ NA_real_
      )
    ) %>%
    select(SNP, effect_allele, other_allele, beta, se, pval, eaf, outcome, study_accession, trait, chr, pos)
  
  out
}

# ---------- Step 1: 列出所有 traits，并筛 pancreatic ----------
traits_all <- get_all_traits(max_pages = max_trait_pages, size = page_size)

# 尽量兼容不同字段名
trait_label_col <- names(traits_all)[str_detect(names(traits_all), regex("^trait$|trait_name|label|name", ignore_case = TRUE))][1]
trait_id_col <- names(traits_all)[str_detect(names(traits_all), regex("efo|trait", ignore_case = TRUE))][1]

if (is.na(trait_label_col)) trait_label_col <- names(traits_all)[1]
if (is.na(trait_id_col)) trait_id_col <- names(traits_all)[1]

traits_all <- traits_all %>%
  mutate(
    trait_label_guess = as.character(.data[[trait_label_col]]),
    trait_id_guess = as.character(.data[[trait_id_col]])
  )

traits_panc <- traits_all %>%
  filter(str_detect(trait_label_guess, regex("pancre", ignore_case = TRUE)) |
           str_detect(trait_id_guess, regex("pancre|EFO_", ignore_case = TRUE))) %>%
  distinct(trait_id_guess, trait_label_guess, .keep_all = TRUE)

fwrite(
  traits_all,
  file.path(registry_dir, "46b_gwas_catalog_traits_all.tsv"),
  sep = "\t", na = "NA"
)

fwrite(
  traits_panc,
  file.path(registry_dir, "46b_gwas_catalog_pancreatic_traits.tsv"),
  sep = "\t", na = "NA"
)

if (nrow(traits_panc) == 0) {
  stop("GWAS Catalog traits 中未检索到 pancreatic 相关 trait。请检查 46b_gwas_catalog_traits_all.tsv。")
}

# ---------- Step 2: 枚举 pancreatic traits 下的 studies ----------
study_rows <- list()

for (i in seq_len(nrow(traits_panc))) {
  trait_id <- traits_panc$trait_id_guess[i]
  trait_label <- traits_panc$trait_label_guess[i]
  
  st <- get_trait_studies(trait_id)
  if (nrow(st) == 0) next
  
  # 尝试识别 study accession
  study_acc_col <- names(st)[str_detect(names(st), regex("study_accession|accession", ignore_case = TRUE))][1]
  if (is.na(study_acc_col)) next
  
  st2 <- st %>%
    mutate(
      trait_id = trait_id,
      trait_label = trait_label,
      study_accession = as.character(.data[[study_acc_col]])
    )
  
  study_rows[[length(study_rows) + 1]] <- st2
}

studies_raw <- bind_rows(study_rows) %>% distinct()

if (nrow(studies_raw) == 0) {
  stop("pancreatic traits 下未检索到带 summary statistics 的 studies。")
}

# ---------- Step 3: 补 study metadata ----------
study_accessions <- unique(as.character(studies_raw$study_accession))
study_meta <- bind_rows(lapply(study_accessions, get_study_metadata))

studies_annot <- studies_raw %>%
  left_join(study_meta, by = "study_accession") %>%
  mutate(
    keyword_score =
      2 * str_detect(paste(title, reported_trait, trait_label, sep = " | "), regex("pancreatic cancer|pancreas|pancreatic", ignore_case = TRUE)) +
      2 * str_detect(paste(title, reported_trait, sep = " | "), regex("panscan|panc4", ignore_case = TRUE)) +
      3 * (pubmed_id %in% target_pmids)
  ) %>%
  arrange(desc(keyword_score), study_accession)

fwrite(
  studies_annot,
  file.path(registry_dir, "46b_pancreatic_trait_studies_annotated.tsv"),
  sep = "\t", na = "NA"
)

# ---------- Step 4: 用当前 primary SNP 测试每个 study 的 overlap ----------
# 为了避免太慢，优先测试 keyword_score 最高的前 20 个 studies
top_n_studies_to_test <- min(20, nrow(studies_annot))
studies_to_test <- studies_annot %>% slice(1:top_n_studies_to_test)

overlap_rows <- list()
assoc_cache_rows <- list()

for (i in seq_len(nrow(studies_to_test))) {
  study_acc <- as.character(studies_to_test$study_accession[i])
  
  hit_rows <- list()
  for (rs in query_snps) {
    aa <- get_variant_assoc_for_study(rs, study_acc)
    if (nrow(aa) > 0) {
      aa$queried_snp <- rs
      hit_rows[[length(hit_rows) + 1]] <- aa
    }
  }
  
  hit_df <- bind_rows(hit_rows)
  
  if (nrow(hit_df) > 0) {
    assoc_cache_rows[[length(assoc_cache_rows) + 1]] <- hit_df
  }
  
  usable_n <- 0L
  if (nrow(hit_df) > 0) {
    tmp_std <- standardise_outcome(hit_df, outcome_label = paste0("PanScan_PanC4_candidate_", study_acc))
    usable_n <- sum(
      !is.na(tmp_std$SNP) &
        !is.na(tmp_std$effect_allele) &
        !is.na(tmp_std$other_allele) &
        !is.na(tmp_std$beta) &
        !is.na(tmp_std$se) &
        !is.na(tmp_std$pval)
    )
  }
  
  overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
    study_accession = study_acc,
    title = studies_to_test$title[i],
    reported_trait = studies_to_test$reported_trait[i],
    trait_label = studies_to_test$trait_label[i],
    pubmed_id = studies_to_test$pubmed_id[i],
    keyword_score = studies_to_test$keyword_score[i],
    queried_snp_n = length(query_snps),
    returned_assoc_rows = nrow(hit_df),
    unique_returned_snps = ifelse(nrow(hit_df) > 0 && "variant_id" %in% names(hit_df), length(unique(hit_df$variant_id)), 0L),
    usable_outcome_rows = usable_n,
    stringsAsFactors = FALSE
  )
}

overlap_audit <- bind_rows(overlap_rows) %>%
  arrange(desc(usable_outcome_rows), desc(unique_returned_snps), desc(keyword_score))

fwrite(
  overlap_audit,
  file.path(registry_dir, "46b_study_overlap_audit.tsv"),
  sep = "\t", na = "NA"
)

assoc_cache <- bind_rows(assoc_cache_rows)
if (nrow(assoc_cache) > 0) {
  fwrite(
    assoc_cache,
    file.path(proc_dir, "46b_association_cache_raw.tsv"),
    sep = "\t", na = "NA"
  )
}

if (nrow(overlap_audit) == 0) {
  stop("没有任何 study 完成 overlap 测试。")
}

# ---------- Step 5: 选择最佳候选 study ----------
best_study <- overlap_audit %>% slice(1)

fwrite(
  best_study,
  file.path(registry_dir, "46b_selected_primary_outcome_candidate.tsv"),
  sep = "\t", na = "NA"
)

selected_acc <- as.character(best_study$study_accession[1])

# ---------- Step 6: 标准化选中 study 的 outcome（仅当前 primary SNP 覆盖） ----------
selected_raw <- assoc_cache %>%
  filter(study_accession == selected_acc)

selected_outcome <- standardise_outcome(
  selected_raw,
  outcome_label = paste0("PanScan_PanC4_", selected_acc)
) %>%
  distinct(SNP, .keep_all = TRUE)

fwrite(
  selected_outcome,
  file.path(result_dir, "46b_selected_primary_outcome_overlap.tsv"),
  sep = "\t", na = "NA"
)

# 同时保存为 47 可直接读取的 outcome 文件（当前仅覆盖 primary SNP）
fwrite(
  selected_outcome,
  file.path(raw_outcome_dir, "46b_pancreatic_cancer_sumstats_for_primary_snps.tsv"),
  sep = "\t", na = "NA"
)

# ---------- Step 7: outcome 结构检查 ----------
required_outcome_cols <- c("SNP", "effect_allele", "other_allele", "beta", "se", "pval", "eaf", "outcome")

attachment_check <- data.frame(
  column = required_outcome_cols,
  present = required_outcome_cols %in% names(selected_outcome),
  stringsAsFactors = FALSE
)

fwrite(
  attachment_check,
  file.path(registry_dir, "46b_primary_outcome_attachment_check.tsv"),
  sep = "\t", na = "NA"
)

# ---------- Step 8: 总状态 ----------
status_summary <- data.frame(
  item = c(
    "n_primary_snps_query",
    "n_pancreatic_traits_found",
    "n_studies_tested",
    "selected_study_accession",
    "selected_study_pubmed_id",
    "selected_study_keyword_score",
    "selected_outcome_rows",
    "recommended_next_action"
  ),
  value = c(
    length(query_snps),
    nrow(traits_panc),
    nrow(studies_to_test),
    selected_acc,
    best_study$pubmed_id[1],
    best_study$keyword_score[1],
    nrow(selected_outcome),
    ifelse(
      nrow(selected_outcome) > 0,
      "Proceed to 47 using 46b_pancreatic_cancer_sumstats_for_primary_snps.tsv",
      "No usable public primary-outcome overlap found; inspect 46b_study_overlap_audit.tsv"
    )
  ),
  stringsAsFactors = FALSE
)

fwrite(
  status_summary,
  file.path(registry_dir, "46b_status_summary.tsv"),
  sep = "\t", na = "NA"
)

# ---------- 日志 ----------
writeLines(
  c(
    "46b completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Find a usable public pancreatic cancer primary outcome study via GWAS Catalog summary statistics API.",
    "",
    "[Key outputs]",
    "- 46b_gwas_catalog_pancreatic_traits.tsv",
    "- 46b_pancreatic_trait_studies_annotated.tsv",
    "- 46b_study_overlap_audit.tsv",
    "- 46b_selected_primary_outcome_candidate.tsv",
    "- 46b_selected_primary_outcome_overlap.tsv",
    "- 46b_primary_outcome_attachment_check.tsv",
    "- 46b_status_summary.tsv"
  ),
  file.path(log_dir, "46b_find_public_primary_outcome_notes.txt")
)

cat("46b_find_public_primary_outcome_via_gwas_catalog_api.R finished successfully.\n")
cat("Generated files:\n")
cat("- 00_registry/46b_gwas_catalog_traits_all.tsv\n")
cat("- 00_registry/46b_gwas_catalog_pancreatic_traits.tsv\n")
cat("- 00_registry/46b_pancreatic_trait_studies_annotated.tsv\n")
cat("- 00_registry/46b_study_overlap_audit.tsv\n")
cat("- 00_registry/46b_selected_primary_outcome_candidate.tsv\n")
cat("- 00_registry/46b_primary_outcome_attachment_check.tsv\n")
cat("- 00_registry/46b_status_summary.tsv\n")
cat("- 02_processed_data/46b_association_cache_raw.tsv\n")
cat("- 03_results/46b_selected_primary_outcome_overlap.tsv\n")
cat("- 01_raw_data/outcome/PanScan_PanC4/46b_pancreatic_cancer_sumstats_for_primary_snps.tsv\n")