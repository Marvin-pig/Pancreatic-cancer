# =========================================================
# 07F_mcd_subtype_clinical_comparison.R
# 作用：
# 1. 读取 cohort 与 K3 subtype annotation
# 2. 生成正式 phenotype 表
# 3. 比较 MCD subtype 与临床特征差异
# =========================================================

library(data.table)

project_root <- Sys.getenv("PROJECT_ROOT", unset = "/Users/wmz_mac/Desktop/胰腺癌")
proc_dir  <- file.path(project_root, "02_processed_data", "TCGA_PAAD")
clust_dir <- file.path(project_root, "04_bulk_analysis", "02_clustering")
cohort_file <- file.path(proc_dir, "tcga_paad_cohort_master_ordered.tsv")
anno_file   <- file.path(clust_dir, "tcga_paad_k3_cluster_annotation.tsv")

# ---- detect sample column ----
detect_sample_col <- function(df) {
  cand <- c("sample_id", "sample_barcode", "SampleID", "sample", "barcode", "submitter_id")
  hit  <- cand[cand %in% colnames(df)]
  if (length(hit) == 0) stop("未找到样本列，候选列名：", paste(cand, collapse = ", "))
  hit[1]
}

cohort <- fread(cohort_file, data.table = FALSE)
anno   <- fread(anno_file,   data.table = FALSE)

# FIX #1：anno 的样本列也动态检测，不再硬编码 "sample_id"
sample_col_cohort <- detect_sample_col(cohort)
sample_col_anno   <- detect_sample_col(anno)

# ---- merge phenotype ----
dat <- merge(
  cohort, anno,
  by.x = sample_col_cohort,
  by.y = sample_col_anno,
  all.x = TRUE,
  sort  = FALSE
)

# 保持原始样本顺序（match 可返回 NA，先检查）
row_idx <- match(cohort[[sample_col_cohort]], dat[[sample_col_cohort]])
if (any(is.na(row_idx))) {
  stop("merge 后 match 返回 NA，cohort 中存在不在 dat 中的样本，请检查样本列匹配逻辑。")
}
dat <- dat[row_idx, , drop = FALSE]
stopifnot(identical(dat[[sample_col_cohort]], cohort[[sample_col_cohort]]))

# FIX #2：subtype_label 缺失改为警告 + 过滤，而非直接 crash
n_na_subtype <- sum(is.na(dat$subtype_label))
if (n_na_subtype > 0) {
  warning(sprintf(
    "%d 个样本在 anno 中无 subtype 注释，已从后续分析中排除。",
    n_na_subtype
  ))
  dat <- dat[!is.na(dat$subtype_label), , drop = FALSE]
}
stopifnot(nrow(dat) > 0)

# ---- save phenotype ----
fwrite(
  dat,
  file.path(clust_dir, "tcga_paad_mcd_subtype_phenotype.tsv"),
  sep   = "\t",
  quote = FALSE
)

# ---- helper functions ----
safe_chr <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "Unknown", "unknown", "UNK")] <- NA
  x
}

# FIX #3 + #6：Fisher/Chi-square 选择逻辑
#   - 先计算期望频数（而非观测值）判断是否需要 Fisher
#   - 对 >2×2 列联表，Fisher 改用 Monte Carlo 模拟（simulate.p.value），避免 workspace 溢出
run_contingency_test <- function(tab) {
  # 期望频数 = 行边际 * 列边际 / 总数
  expected   <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  use_fisher <- any(expected < 5)
  is_2x2     <- all(dim(tab) == 2)
  
  if (use_fisher) {
    if (is_2x2) {
      res       <- fisher.test(tab)
      test_name <- "Fisher"
    } else {
      # 大于 2×2 时用蒙特卡洛模拟，B=2000 次
      res       <- fisher.test(tab, simulate.p.value = TRUE, B = 2000)
      test_name <- "Fisher(simulated)"
    }
  } else {
    res       <- chisq.test(tab)
    test_name <- "Chi-square"
  }
  list(p_value = res$p.value, test_name = test_name)
}

# ---- statistical comparisons ----
clinical_results <- list()

# 1) age: Kruskal-Wallis
if ("age_at_diagnosis_years" %in% colnames(dat)) {
  # FIX #7：确保为数值型，且至少 2 个亚型有数据
  age_col <- suppressWarnings(as.numeric(dat$age_at_diagnosis_years))
  tmp <- dat[!is.na(age_col) & !is.na(dat$subtype_label), , drop = FALSE]
  tmp$age_at_diagnosis_years <- age_col[!is.na(age_col) & !is.na(dat$subtype_label)]
  
  n_groups <- length(unique(tmp$subtype_label))
  if (n_groups >= 2 && nrow(tmp) > 0) {
    fit <- kruskal.test(age_at_diagnosis_years ~ subtype_label, data = tmp)
    clinical_results[[length(clinical_results) + 1]] <- data.frame(
      variable      = "age_at_diagnosis_years",
      variable_type = "continuous",
      test          = "Kruskal-Wallis",
      p_value       = fit$p.value,
      stringsAsFactors = FALSE
    )
  } else {
    warning("age_at_diagnosis_years: 有效亚型组数 < 2，跳过 Kruskal-Wallis 检验。")
  }
}

# 2) categorical variables
cat_vars <- c("gender_final", "ajcc_stage_final", "tumor_grade_final")

for (v in cat_vars) {
  if (!v %in% colnames(dat)) next
  
  tmp     <- dat[, c("subtype_label", v), drop = FALSE]
  tmp[[v]] <- safe_chr(tmp[[v]])
  tmp      <- tmp[!is.na(tmp$subtype_label) & !is.na(tmp[[v]]), , drop = FALSE]
  
  if (nrow(tmp) == 0) {
    warning(sprintf("%s: 过滤后无有效数据，跳过。", v))
    next
  }
  
  tab <- table(tmp$subtype_label, tmp[[v]])
  
  # 去除全为 0 的行/列（避免 degenerate table）
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  
  if (any(dim(tab) < 2)) {
    warning(sprintf("%s: 列联表维度不足（行或列 < 2），跳过。", v))
    next
  }
  
  res <- tryCatch(
    run_contingency_test(tab),
    error = function(e) {
      warning(sprintf("%s 检验失败：%s", v, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(res)) next
  
  clinical_results[[length(clinical_results) + 1]] <- data.frame(
    variable      = v,
    variable_type = "categorical",
    test          = res$test_name,
    p_value       = res$p_value,
    stringsAsFactors = FALSE
  )
  
  # 保存列联表（FIX #4：移除 fwrite 不支持的 row.names 参数）
  tab_df <- as.data.frame.matrix(tab)
  tab_df$subtype_label <- rownames(tab_df)
  
  fwrite(
    tab_df,
    file.path(clust_dir, paste0("tcga_paad_", v, "_by_mcd_subtype.tsv")),
    sep   = "\t",
    quote = FALSE
    # 注意：fwrite 不支持 row.names，已删除；subtype_label 列已显式加入 tab_df
  )
}

# FIX #5：clinical_results 为空时给出明确提示，不 crash
if (length(clinical_results) == 0) {
  warning("未找到任何可分析的临床变量，clinical_results 为空。")
  clinical_results <- data.frame(
    variable = character(0), variable_type = character(0),
    test = character(0), p_value = numeric(0), p_adj_BH = numeric(0)
  )
} else {
  clinical_results        <- do.call(rbind, clinical_results)
  clinical_results$p_adj_BH <- p.adjust(clinical_results$p_value, method = "BH")
}

fwrite(
  clinical_results,
  file.path(clust_dir, "tcga_paad_mcd_subtype_clinical_comparison.tsv"),
  sep   = "\t",
  quote = FALSE
)

# ---- subtype summary（FIX #8：先排除 NA，再统计）----
subtype_tbl <- table(dat$subtype_label)
summary_df  <- data.frame(
  subtype_label = names(subtype_tbl),
  n_samples     = as.integer(subtype_tbl),
  stringsAsFactors = FALSE
)
summary_df <- summary_df[order(summary_df$subtype_label), ]

fwrite(
  summary_df,
  file.path(clust_dir, "tcga_paad_mcd_subtype_clinical_summary.tsv"),
  sep   = "\t",
  quote = FALSE
)

print(clinical_results)
print(summary_df)
message("MCD subtype clinical comparison finished.")
