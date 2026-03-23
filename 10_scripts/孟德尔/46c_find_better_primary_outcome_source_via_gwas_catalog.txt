# =========================
# 46c_find_better_primary_outcome_source_via_gwas_catalog.R
# 目的：
# 1. 从 GWAS Catalog REST API 中严格检索 pancreatic cancer risk studies
# 2. 优先锁定 PanScan / PanC4 / PMID 29422604 / pancreatic cancer / European 风格 study
# 3. 用当前 primary instrument SNP 列表测试各候选 study 的真实 overlap
# 4. 只在 study metadata 与 overlap 都达到最低标准时，才写出可进入 47 的 outcome 文件
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

result_dir      <- file.path(root_dir, "03_results")
registry_dir    <- file.path(root_dir, "00_registry")
proc_dir        <- file.path(root_dir, "02_processed_data")
raw_outcome_dir <- file.path(root_dir, "01_raw_data/outcome/PanScan_PanC4")
log_dir         <- file.path(root_dir, "06_logs")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_outcome_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 输入：当前 primary SNP ----------
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

cat("已加载 ", length(query_snps), " 个 primary SNPs\n")

# ---------- API endpoints ----------
# 修正：统一使用 GWAS Catalog 主 REST API
base_rest <- "https://www.ebi.ac.uk/gwas/rest/api"

# ---------- 参数 ----------
target_keywords_strong <- c(
  "pancreatic cancer",
  "pancreas cancer",
  "pancreatic carcinoma",
  "pancreatic ductal adenocarcinoma",
  "pancreatic neoplasm"
)

target_keywords_support <- c(
  "panscan",
  "panc4",
  "pancreas",
  "pancreatic"
)

target_pmids <- c("29422604")
min_usable_overlap <- 3L
max_studies_to_test <- 30L
page_size <- 200L

# =============================================
# 工具函数
# =============================================

safe_get_json <- function(url, pause_sec = 0.15, max_retries = 3) {
  for (attempt in seq_len(max_retries)) {
    Sys.sleep(pause_sec)
    tryCatch({
      r <- httr::GET(
        url,
        httr::timeout(60),
        httr::add_headers(Accept = "application/json")
      )
      sc <- httr::status_code(r)
      txt <- httr::content(r, as = "text", encoding = "UTF-8")
      
      if (sc == 200) {
        return(list(
          ok = TRUE,
          status_code = sc,
          text = txt,
          json = jsonlite::fromJSON(txt, simplifyVector = FALSE)
        ))
      }
      
      # 429 = rate limit -> 等待后重试
      if (sc == 429 && attempt < max_retries) {
        cat("  [Rate limited] attempt ", attempt, " / ", max_retries, " — waiting 5s...\n")
        Sys.sleep(5)
        next
      }
      
      return(list(ok = FALSE, status_code = sc, text = txt, json = NULL))
    }, error = function(e) {
      if (attempt == max_retries) {
        return(list(ok = FALSE, status_code = NA_integer_,
                    text = conditionMessage(e), json = NULL))
      }
      Sys.sleep(2)
    })
  }
  list(ok = FALSE, status_code = NA_integer_, text = "max retries exceeded", json = NULL)
}

# 从 HAL+JSON 的 _embedded 中提取 data.frame
extract_embedded_df <- function(x, embedded_name = NULL) {
  if (is.null(x) || is.null(x[["_embedded"]])) return(data.frame())
  emb <- x[["_embedded"]]
  
  if (!is.null(embedded_name) && embedded_name %in% names(emb)) {
    y <- emb[[embedded_name]]
  } else {
    y <- emb[[1]]
  }
  
  if (is.null(y)) return(data.frame())
  if (is.data.frame(y)) return(as.data.frame(y, stringsAsFactors = FALSE))
  
  if (is.list(y)) {
    rows <- lapply(y, function(z) {
      if (is.null(z)) return(NULL)
      flat <- list()
      for (k in names(z)) {
        v <- z[[k]]
        if (is.atomic(v) && length(v) <= 1) {
          flat[[k]] <- if (length(v) == 0) NA_character_ else as.character(v)
        } else if (is.list(v) && !is.null(names(v)) && k != "_links") {
          # 展平一层嵌套
          for (kk in names(v)) {
            vv <- v[[kk]]
            if (is.atomic(vv) && length(vv) <= 1) {
              flat[[paste0(k, ".", kk)]] <- if (length(vv) == 0) NA_character_ else as.character(vv)
            }
          }
        }
      }
      if (length(flat) == 0) return(NULL)
      as.data.frame(flat, stringsAsFactors = FALSE)
    })
    rows <- rows[!sapply(rows, is.null)]
    if (length(rows) == 0) return(data.frame())
    return(bind_rows(rows))
  }
  
  data.frame()
}

# 获取分页的 next link
follow_next_link <- function(json_obj) {
  if (is.null(json_obj) ||
      is.null(json_obj[["_links"]]) ||
      is.null(json_obj[["_links"]][["next"]])) {
    return(NA_character_)
  }
  nxt <- json_obj[["_links"]][["next"]][["href"]]
  if (is.null(nxt) || length(nxt) == 0) return(NA_character_)
  as.character(nxt)
}

safe_text <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x
}

score_text <- function(text_vec, keywords, weight = 1L) {
  txt <- tolower(paste(text_vec, collapse = " | "))
  sum(sapply(keywords, function(k) str_detect(txt, fixed(tolower(k))))) * weight
}

# =============================================
# Step 1: 获取所有 EFO traits，筛选 pancreatic 相关
# =============================================
cat("Step 1: 获取 EFO traits...\n")

get_all_efo_traits <- function(size = 200, max_pages = 100) {
  # 修正：使用主 REST API 的 /efoTraits 端点，分页参数为 page + size
  url <- paste0(base_rest, "/efoTraits?page=0&size=", size)
  all_rows <- list()
  page_i <- 1
  
  while (!is.na(url) && nzchar(url) && page_i <= max_pages) {
    cat("  Trait page ", page_i, " ...\n")
    res <- safe_get_json(url)
    if (!res$ok || is.null(res$json)) {
      cat("  [WARN] Trait page fetch failed (status: ", res$status_code, ")\n")
      break
    }
    
    df <- extract_embedded_df(res$json, embedded_name = "efoTraits")
    if (nrow(df) > 0) {
      all_rows[[length(all_rows) + 1]] <- df
      cat("    -> got ", nrow(df), " traits\n")
    } else {
      break
    }
    
    url <- follow_next_link(res$json)
    page_i <- page_i + 1
  }
  
  bind_rows(all_rows)
}

traits_all <- get_all_efo_traits(size = page_size, max_pages = 100)
cat("  Total EFO traits fetched: ", nrow(traits_all), "\n")

# 确定 trait 标签列和 ID 列
# GWAS REST API /efoTraits 返回字段一般为: trait, shortForm, uri
trait_label_col <- if ("trait" %in% names(traits_all)) {
  "trait"
} else {
  cands <- names(traits_all)[str_detect(names(traits_all),
                                        regex("trait|label|name", ignore_case = TRUE))]
  if (length(cands) > 0) cands[1] else names(traits_all)[1]
}

# shortForm 通常为 EFO_XXXX 格式
trait_id_col <- if ("shortForm" %in% names(traits_all)) {
  "shortForm"
} else {
  cands <- names(traits_all)[str_detect(names(traits_all),
                                        regex("efo|short|id", ignore_case = TRUE))]
  if (length(cands) > 0) cands[1] else trait_label_col
}

traits_all <- traits_all %>%
  mutate(
    trait_label_guess = safe_text(.data[[trait_label_col]]),
    trait_id_guess    = safe_text(.data[[trait_id_col]])
  )

# 严格筛选 pancreatic 相关
traits_panc_strict <- traits_all %>%
  filter(str_detect(trait_label_guess, regex("pancre", ignore_case = TRUE))) %>%
  distinct(trait_id_guess, trait_label_guess, .keep_all = TRUE)

cat("  Pancreatic-related traits: ", nrow(traits_panc_strict), "\n")

fwrite(traits_all,
       file.path(registry_dir, "46c_gwas_catalog_traits_all.tsv"),
       sep = "\t", na = "NA")

fwrite(traits_panc_strict,
       file.path(registry_dir, "46c_gwas_catalog_pancreatic_traits_strict.tsv"),
       sep = "\t", na = "NA")

if (nrow(traits_panc_strict) == 0) {
  stop("严格 trait 过滤后未检索到 pancreatic 相关 trait。")
}

# =============================================
# Step 2: 获取每个 pancreatic trait 关联的 studies
# =============================================
cat("Step 2: 获取 pancreatic trait 对应的 studies...\n")

get_trait_studies <- function(trait_short_form, size = 200) {
  # 修正：使用 /efoTraits/{shortForm}/studies
  # 注意 shortForm 可能含有冒号等特殊字符，需要 URL encode
  encoded_id <- URLencode(trait_short_form, reserved = TRUE)
  url <- paste0(base_rest, "/efoTraits/", encoded_id, "/studies?size=", size)
  res <- safe_get_json(url)
  if (!res$ok || is.null(res$json)) return(data.frame())
  extract_embedded_df(res$json, embedded_name = "studies")
}

study_rows <- list()

for (i in seq_len(nrow(traits_panc_strict))) {
  trait_id    <- traits_panc_strict$trait_id_guess[i]
  trait_label <- traits_panc_strict$trait_label_guess[i]
  cat("  Trait [", i, "/", nrow(traits_panc_strict), "]: ", trait_label, " (", trait_id, ")\n")
  
  st <- get_trait_studies(trait_id)
  if (nrow(st) == 0) {
    cat("    -> 0 studies\n")
    next
  }
  
  # 找 study accession 列
  study_acc_col <- if ("accessionId" %in% names(st)) {
    "accessionId"
  } else {
    cands <- names(st)[str_detect(names(st),
                                  regex("accession", ignore_case = TRUE))]
    if (length(cands) > 0) cands[1] else NA_character_
  }
  
  if (is.na(study_acc_col)) {
    cat("    -> no accession column found, skipping\n")
    next
  }
  
  st2 <- st %>%
    mutate(
      trait_id = trait_id,
      trait_label = trait_label,
      study_accession = as.character(.data[[study_acc_col]])
    )
  
  cat("    -> ", nrow(st2), " studies\n")
  study_rows[[length(study_rows) + 1]] <- st2
}

studies_raw <- bind_rows(study_rows) %>%
  distinct(study_accession, .keep_all = TRUE)

cat("  Unique candidate studies: ", nrow(studies_raw), "\n")

if (nrow(studies_raw) == 0) {
  stop("严格 pancreatic traits 下未检索到 studies。")
}

# =============================================
# Step 3: 获取 study metadata 并打分排序
# =============================================
cat("Step 3: 获取 study metadata 并排序...\n")

get_study_metadata <- function(study_acc) {
  url <- paste0(base_rest, "/studies/", URLencode(study_acc, reserved = TRUE))
  res <- safe_get_json(url)
  
  default_row <- data.frame(
    study_accession    = study_acc,
    title              = NA_character_,
    first_author       = NA_character_,
    pubmed_id          = NA_character_,
    reported_trait      = NA_character_,
    initial_sample     = NA_character_,
    replication_sample = NA_character_,
    stringsAsFactors   = FALSE
  )
  
  if (!res$ok || is.null(res$json)) return(default_row)
  
  j <- res$json
  
  # 安全提取嵌套字段
  safe_extract <- function(obj, ...) {
    keys <- list(...)
    curr <- obj
    for (k in keys) {
      if (is.null(curr) || !k %in% names(curr)) return(NA_character_)
      curr <- curr[[k]]
    }
    if (is.null(curr) || length(curr) == 0) return(NA_character_)
    as.character(curr)
  }
  
  # diseaseTrait 有时是嵌套对象
  reported <- safe_extract(j, "diseaseTrait", "trait")
  if (is.na(reported)) reported <- safe_extract(j, "diseaseTrait")
  
  # initialSampleSize / replicateSampleSize 为字符串
  initial_sample <- safe_extract(j, "initialSampleSize")
  replication_sample <- safe_extract(j, "replicateSampleSize")
  
  # pubmedId 在 publicationInfo 或直接字段
  pmid <- safe_extract(j, "publicationInfo", "pubmedId")
  if (is.na(pmid)) pmid <- safe_extract(j, "pubmedId")
  
  # author
  author <- safe_extract(j, "publicationInfo", "author", "fullname")
  if (is.na(author)) author <- safe_extract(j, "author")
  
  # title
  title <- safe_extract(j, "publicationInfo", "title")
  if (is.na(title)) title <- safe_extract(j, "title")
  
  data.frame(
    study_accession    = study_acc,
    title              = title,
    first_author       = author,
    pubmed_id          = pmid,
    reported_trait      = reported,
    initial_sample     = initial_sample,
    replication_sample = replication_sample,
    stringsAsFactors   = FALSE
  )
}

study_accessions <- unique(as.character(studies_raw$study_accession))
cat("  Fetching metadata for ", length(study_accessions), " studies...\n")

study_meta <- bind_rows(lapply(seq_along(study_accessions), function(idx) {
  if (idx %% 10 == 0) cat("    metadata ", idx, "/", length(study_accessions), "\n")
  get_study_metadata(study_accessions[idx])
}))

# 合并并打分
studies_annot <- studies_raw %>%
  left_join(study_meta, by = "study_accession") %>%
  mutate(
    searchable_text = paste(
      safe_text(title),
      safe_text(reported_trait),
      safe_text(trait_label),
      safe_text(first_author),
      safe_text(initial_sample),
      safe_text(replication_sample),
      sep = " | "
    ),
    score_strong  = sapply(searchable_text, function(t) score_text(t, target_keywords_strong, weight = 4L)),
    score_support = sapply(searchable_text, function(t) score_text(t, target_keywords_support, weight = 2L)),
    score_pmid    = ifelse(!is.na(pubmed_id) & pubmed_id %in% target_pmids, 8L, 0L),
    score_europe  = ifelse(str_detect(tolower(searchable_text), "europe"), 1L, 0L),
    score_total   = score_strong + score_support + score_pmid + score_europe,
    metadata_complete = (!is.na(title) & title != "") |
      (!is.na(reported_trait) & reported_trait != "") |
      (!is.na(pubmed_id) & pubmed_id != "")
  ) %>%
  arrange(desc(score_total), desc(metadata_complete), study_accession)

fwrite(studies_annot,
       file.path(registry_dir, "46c_pancreatic_trait_studies_annotated.tsv"),
       sep = "\t", na = "NA")

cat("  Top 5 studies by score:\n")
top5 <- studies_annot %>%
  select(study_accession, score_total, pubmed_id, reported_trait, title) %>%
  slice(1:min(5, n()))
print(as.data.frame(top5))

# 筛选值得测试的 studies
studies_to_test <- studies_annot %>%
  filter(score_total > 0 | metadata_complete) %>%
  slice(1:min(max_studies_to_test, n()))

if (nrow(studies_to_test) == 0) {
  stop("没有任何值得测试的 pancreatic cancer study。")
}

cat("  Will test overlap for ", nrow(studies_to_test), " studies\n")

# =============================================
# Step 4: SNP overlap 测试
# 修正：通过 /singleNucleotidePolymorphisms/{rsid}/associations 获取关联，
#        然后按 study_accession 过滤
# =============================================
cat("Step 4: 测试 SNP overlap...\n")

# 先为每个 SNP 一次性获取所有 associations（避免 SNP × study 的组合查询）
cat("  Fetching associations for all ", length(query_snps), " query SNPs...\n")

get_snp_associations <- function(rsid) {
  url <- paste0(base_rest, "/singleNucleotidePolymorphisms/", rsid, "/associations")
  res <- safe_get_json(url, pause_sec = 0.2)
  if (!res$ok || is.null(res$json)) return(data.frame())
  
  df <- extract_embedded_df(res$json, embedded_name = "associations")
  if (nrow(df) == 0) return(data.frame())
  
  df$queried_snp <- rsid
  df
}

all_assoc_list <- list()
for (idx in seq_along(query_snps)) {
  rs <- query_snps[idx]
  if (idx %% 5 == 0 || idx == 1) {
    cat("  SNP ", idx, "/", length(query_snps), ": ", rs, "\n")
  }
  aa <- get_snp_associations(rs)
  if (nrow(aa) > 0) {
    all_assoc_list[[length(all_assoc_list) + 1]] <- aa
  }
}

all_assoc_raw <- bind_rows(all_assoc_list)
cat("  Total association rows fetched: ", nrow(all_assoc_raw), "\n")

if (nrow(all_assoc_raw) > 0) {
  fwrite(all_assoc_raw,
         file.path(proc_dir, "46c_association_cache_raw.tsv"),
         sep = "\t", na = "NA")
}

# 识别 study_accession 列 — associations 响应中可能叫不同名字
# 通常在 _links 或嵌套中; 如果 extract_embedded_df 展平了就会在列里
identify_study_col <- function(df) {
  cands <- c("study_accession", "studyId", "study.accessionId",
             "accessionId", "study")
  for (c in cands) {
    if (c %in% names(df)) return(c)
  }
  # 尝试模糊匹配
  matched <- names(df)[str_detect(names(df), regex("accession|study", ignore_case = TRUE))]
  if (length(matched) > 0) return(matched[1])
  NA_character_
}

# 从 association 响应中提取 study accession
# 如果直接列中没有，尝试从 _links 中的 study href 提取
extract_study_from_links <- function(df) {
  if ("study_accession" %in% names(df) &&
      !all(is.na(df$study_accession))) {
    return(df)
  }
  
  # 检查 _links.study.href 样式的列
  link_cols <- names(df)[str_detect(names(df), regex("link.*study.*href|study.*link.*href", ignore_case = TRUE))]
  if (length(link_cols) > 0) {
    # href 格式通常为 .../studies/GCSTXXXXXX
    df$study_accession <- str_extract(df[[link_cols[1]]], "GCST\\d+")
    return(df)
  }
  
  # 尝试其他列
  study_col <- identify_study_col(df)
  if (!is.na(study_col) && study_col != "study_accession") {
    df$study_accession <- as.character(df[[study_col]])
    # 尝试提取 GCST ID
    has_gcst <- str_detect(df$study_accession, "GCST\\d+")
    df$study_accession <- ifelse(has_gcst,
                                 str_extract(df$study_accession, "GCST\\d+"),
                                 df$study_accession)
    return(df)
  }
  
  df$study_accession <- NA_character_
  df
}

if (nrow(all_assoc_raw) > 0) {
  all_assoc_raw <- extract_study_from_links(all_assoc_raw)
}

# 标准化 outcome 的函数
standardise_outcome <- function(df, outcome_label) {
  empty_out <- data.frame(
    SNP = character(), effect_allele = character(), other_allele = character(),
    beta = numeric(), se = numeric(), pval = numeric(), eaf = numeric(),
    chr = integer(), pos = integer(), study_accession = character(),
    outcome = character(), stringsAsFactors = FALSE
  )
  if (nrow(df) == 0) return(empty_out)
  
  # 灵活匹配列名
  find_col <- function(patterns, data_names) {
    for (p in patterns) {
      m <- data_names[str_detect(data_names, regex(p, ignore_case = TRUE))]
      if (length(m) > 0) return(m[1])
    }
    NA_character_
  }
  
  col_snp <- find_col(c("^queried_snp$", "^variant_id$", "^rsId$", "^snp$", "rsid"), names(df))
  col_ea  <- find_col(c("effect_allele", "riskAllele", "risk.allele", "strongestRiskAlleles"), names(df))
  col_oa  <- find_col(c("other_allele"), names(df))
  col_beta <- find_col(c("^beta$", "betaNum"), names(df))
  col_se  <- find_col(c("^se$", "standardError"), names(df))
  col_or  <- find_col(c("^orPerCopyNum$", "^odds_ratio$", "^or$"), names(df))
  col_pval <- find_col(c("^pvalue$", "^p_value$", "^pval$", "pvalueMantissa"), names(df))
  col_eaf <- find_col(c("^eaf$", "riskFrequency", "effect_allele_frequency"), names(df))
  col_chr <- find_col(c("^chrom", "^chr$"), names(df))
  col_pos <- find_col(c("^base_pair", "^pos$", "^position$"), names(df))
  
  # 如果有 pvalueMantissa + pvalueExponent，需要组合
  col_pval_mantissa <- find_col(c("pvalueMantissa"), names(df))
  col_pval_exponent <- find_col(c("pvalueExponent"), names(df))
  
  out <- df %>%
    mutate(
      .snp_val = if (!is.na(col_snp)) as.character(.data[[col_snp]]) else NA_character_,
      .ea_val  = if (!is.na(col_ea)) as.character(.data[[col_ea]]) else NA_character_,
      .oa_val  = if (!is.na(col_oa)) as.character(.data[[col_oa]]) else NA_character_,
      .beta_val = if (!is.na(col_beta)) suppressWarnings(as.numeric(.data[[col_beta]])) else NA_real_,
      .se_val   = if (!is.na(col_se)) suppressWarnings(as.numeric(.data[[col_se]])) else NA_real_,
      .or_val   = if (!is.na(col_or)) suppressWarnings(as.numeric(.data[[col_or]])) else NA_real_,
      .eaf_val  = if (!is.na(col_eaf)) suppressWarnings(as.numeric(.data[[col_eaf]])) else NA_real_,
      .chr_val  = if (!is.na(col_chr)) suppressWarnings(as.integer(.data[[col_chr]])) else NA_integer_,
      .pos_val  = if (!is.na(col_pos)) suppressWarnings(as.integer(.data[[col_pos]])) else NA_integer_
    )
  
  # 处理 p-value
  if (!is.na(col_pval_mantissa) && !is.na(col_pval_exponent)) {
    out <- out %>%
      mutate(
        .pval_val = suppressWarnings(
          as.numeric(.data[[col_pval_mantissa]]) *
            10^as.numeric(.data[[col_pval_exponent]])
        )
      )
  } else if (!is.na(col_pval)) {
    out <- out %>%
      mutate(.pval_val = suppressWarnings(as.numeric(.data[[col_pval]])))
  } else {
    out <- out %>% mutate(.pval_val = NA_real_)
  }
  
  # effect_allele 有时是 "rs12345-A" 格式（strongestRiskAlleles），需要解析
  out <- out %>%
    mutate(
      .ea_val = ifelse(str_detect(.ea_val, "-"),
                       str_extract(.ea_val, "[A-Za-z]+$"),
                       .ea_val)
    )
  
  out <- out %>%
    transmute(
      SNP = .snp_val,
      effect_allele = toupper(.ea_val),
      other_allele  = toupper(.oa_val),
      beta = case_when(
        !is.na(.beta_val) ~ .beta_val,
        is.na(.beta_val) & !is.na(.or_val) & .or_val > 0 ~ log(.or_val),
        TRUE ~ NA_real_
      ),
      se = .se_val,
      pval = .pval_val,
      eaf  = .eaf_val,
      chr  = .chr_val,
      pos  = .pos_val,
      study_accession = if ("study_accession" %in% names(df)) as.character(df$study_accession) else NA_character_,
      outcome = outcome_label
    )
  
  out
}

# 对每个候选 study 评估 overlap
test_study_accessions <- unique(as.character(studies_to_test$study_accession))

overlap_rows <- list()

for (i in seq_len(nrow(studies_to_test))) {
  study_acc <- as.character(studies_to_test$study_accession[i])
  
  # 从全部 association 缓存中过滤当前 study
  if (nrow(all_assoc_raw) > 0 && "study_accession" %in% names(all_assoc_raw)) {
    hit_df <- all_assoc_raw %>%
      filter(study_accession == study_acc)
  } else {
    hit_df <- data.frame()
  }
  
  usable_n <- 0L
  if (nrow(hit_df) > 0) {
    tmp_std <- standardise_outcome(hit_df,
                                   outcome_label = paste0("Pancreatic_cancer_candidate_", study_acc))
    usable_n <- sum(
      !is.na(tmp_std$SNP) &
        !is.na(tmp_std$beta) &
        !is.na(tmp_std$pval),
      na.rm = TRUE
    )
  }
  
  overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
    study_accession     = study_acc,
    title               = studies_to_test$title[i],
    reported_trait       = studies_to_test$reported_trait[i],
    trait_label         = studies_to_test$trait_label[i],
    pubmed_id           = studies_to_test$pubmed_id[i],
    score_total         = studies_to_test$score_total[i],
    metadata_complete   = studies_to_test$metadata_complete[i],
    queried_snp_n       = length(query_snps),
    returned_assoc_rows = nrow(hit_df),
    unique_returned_snps = if (nrow(hit_df) > 0 && "queried_snp" %in% names(hit_df)) {
      length(unique(hit_df$queried_snp))
    } else { 0L },
    usable_outcome_rows = usable_n,
    stringsAsFactors    = FALSE
  )
}

overlap_audit <- bind_rows(overlap_rows) %>%
  mutate(
    overlap_ratio = ifelse(queried_snp_n > 0, usable_outcome_rows / queried_snp_n, 0),
    final_rank    = score_total * 1000 + usable_outcome_rows * 10 + unique_returned_snps
  ) %>%
  arrange(desc(final_rank), desc(usable_outcome_rows), desc(score_total), desc(metadata_complete))

fwrite(overlap_audit,
       file.path(registry_dir, "46c_study_overlap_audit.tsv"),
       sep = "\t", na = "NA")

cat("  Overlap audit (top 5):\n")
print(as.data.frame(overlap_audit %>%
                      select(study_accession, score_total, usable_outcome_rows, final_rank, reported_trait) %>%
                      slice(1:min(5, n()))))

if (nrow(overlap_audit) == 0) {
  stop("没有任何 study 完成 overlap 测试。")
}

# =============================================
# Step 5: 选择最佳候选 study
# =============================================
cat("Step 5: 选择最佳候选 study...\n")

best_study <- overlap_audit %>% slice(1)
selected_acc <- as.character(best_study$study_accession[1])

cat("  Selected: ", selected_acc,
    " (score=", best_study$score_total[1],
    ", usable=", best_study$usable_outcome_rows[1], ")\n")

fwrite(best_study,
       file.path(registry_dir, "46c_selected_primary_outcome_candidate.tsv"),
       sep = "\t", na = "NA")

# 提取并标准化该 study 的 outcome 数据
selected_raw <- if (nrow(all_assoc_raw) > 0 && "study_accession" %in% names(all_assoc_raw)) {
  all_assoc_raw %>% filter(study_accession == selected_acc)
} else {
  data.frame()
}

selected_outcome <- standardise_outcome(
  selected_raw,
  outcome_label = paste0("Pancreatic_cancer_", selected_acc)
) %>%
  distinct(SNP, .keep_all = TRUE)

fwrite(selected_outcome,
       file.path(result_dir, "46c_selected_primary_outcome_overlap.tsv"),
       sep = "\t", na = "NA")

# =============================================
# Step 6: 结构检查
# =============================================
cat("Step 6: 结构检查...\n")

required_outcome_cols <- c("SNP", "effect_allele", "other_allele", "beta", "se", "pval", "eaf", "outcome")
attachment_check <- data.frame(
  column  = required_outcome_cols,
  present = required_outcome_cols %in% names(selected_outcome),
  n_valid = sapply(required_outcome_cols, function(col) {
    if (col %in% names(selected_outcome)) sum(!is.na(selected_outcome[[col]])) else 0L
  }),
  stringsAsFactors = FALSE
)

fwrite(attachment_check,
       file.path(registry_dir, "46c_primary_outcome_attachment_check.tsv"),
       sep = "\t", na = "NA")

print(attachment_check)

# =============================================
# Step 7: 只有满足最低标准才写入 47 可用文件
# =============================================
cat("Step 7: 判断是否满足最低标准...\n")

candidate_is_acceptable <- (
  nrow(selected_outcome) >= min_usable_overlap &&
    isTRUE(best_study$score_total[1] > 0) &&
    isTRUE(best_study$metadata_complete[1])
)

target_file <- file.path(raw_outcome_dir, "46c_pancreatic_cancer_sumstats_for_primary_snps.tsv")

if (candidate_is_acceptable) {
  fwrite(selected_outcome, target_file, sep = "\t", na = "NA")
  cat("  [PASS] 候选 study 可用，已写入: ", target_file, "\n")
} else {
  cat("  [FAIL] 候选 study 未达到最低标准，不写入 47 可用文件\n")
  cat("    - nrow(selected_outcome) = ", nrow(selected_outcome), " (需要 >= ", min_usable_overlap, ")\n")
  cat("    - score_total = ", best_study$score_total[1], " (需要 > 0)\n")
  cat("    - metadata_complete = ", best_study$metadata_complete[1], " (需要 TRUE)\n")
}

# =============================================
# Step 8: 总状态
# =============================================
cat("Step 8: 输出总状态...\n")

status_summary <- data.frame(
  item = c(
    "n_primary_snps_query",
    "n_strict_pancreatic_traits_found",
    "n_candidate_studies_tested",
    "selected_study_accession",
    "selected_study_pubmed_id",
    "selected_study_score_total",
    "selected_study_metadata_complete",
    "selected_outcome_rows",
    "candidate_is_acceptable_for_47",
    "recommended_next_action"
  ),
  value = c(
    as.character(length(query_snps)),
    as.character(nrow(traits_panc_strict)),
    as.character(nrow(studies_to_test)),
    selected_acc,
    as.character(best_study$pubmed_id[1]),
    as.character(best_study$score_total[1]),
    as.character(best_study$metadata_complete[1]),
    as.character(nrow(selected_outcome)),
    as.character(candidate_is_acceptable),
    ifelse(
      candidate_is_acceptable,
      "Proceed to 47 using 46c_pancreatic_cancer_sumstats_for_primary_snps.tsv",
      "Do NOT proceed to 47 yet; current public candidate is not close enough to PanScan/PanC4 primary outcome"
    )
  ),
  stringsAsFactors = FALSE
)

fwrite(status_summary,
       file.path(registry_dir, "46c_status_summary.tsv"),
       sep = "\t", na = "NA")

# ---------- 日志 ----------
writeLines(
  c(
    "46c completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Search for a scientifically closer pancreatic cancer primary outcome source",
    "with stricter GWAS Catalog filtering via the official REST API.",
    "",
    "[API used]",
    paste0("Base: ", base_rest),
    "",
    "[Key outputs]",
    "- 46c_gwas_catalog_traits_all.tsv",
    "- 46c_gwas_catalog_pancreatic_traits_strict.tsv",
    "- 46c_pancreatic_trait_studies_annotated.tsv",
    "- 46c_study_overlap_audit.tsv",
    "- 46c_selected_primary_outcome_candidate.tsv",
    "- 46c_selected_primary_outcome_overlap.tsv",
    "- 46c_primary_outcome_attachment_check.tsv",
    "- 46c_status_summary.tsv",
    "- 46c_association_cache_raw.tsv",
    "",
    "[Result]",
    paste0("Candidate acceptable: ", candidate_is_acceptable),
    paste0("Selected study: ", selected_acc),
    paste0("Usable outcome rows: ", nrow(selected_outcome))
  ),
  file.path(log_dir, "46c_find_better_primary_outcome_notes.txt")
)

cat("\n========================================\n")
cat("46c_find_better_primary_outcome_source_via_gwas_catalog.R 运行完成\n")
cat("========================================\n")
cat("Generated files:\n")
cat("- 00_registry/46c_gwas_catalog_traits_all.tsv\n")
cat("- 00_registry/46c_gwas_catalog_pancreatic_traits_strict.tsv\n")
cat("- 00_registry/46c_pancreatic_trait_studies_annotated.tsv\n")
cat("- 00_registry/46c_study_overlap_audit.tsv\n")
cat("- 00_registry/46c_selected_primary_outcome_candidate.tsv\n")
cat("- 00_registry/46c_primary_outcome_attachment_check.tsv\n")
cat("- 00_registry/46c_status_summary.tsv\n")
cat("- 02_processed_data/46c_association_cache_raw.tsv\n")
cat("- 03_results/46c_selected_primary_outcome_overlap.tsv\n")
if (candidate_is_acceptable) {
  cat("- 01_raw_data/outcome/PanScan_PanC4/46c_pancreatic_cancer_sumstats_for_primary_snps.tsv [WRITTEN]\n")
} else {
  cat("- 01_raw_data/outcome/PanScan_PanC4/46c_pancreatic_cancer_sumstats_for_primary_snps.tsv [NOT WRITTEN - below threshold]\n")
}
cat("========================================\n")