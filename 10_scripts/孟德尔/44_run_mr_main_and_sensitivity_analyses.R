# =========================
# 44_run_mr_main_and_sensitivity_analyses.R
# 目的：
# 1. 读取 43 的 harmonised mr_keep 数据
# 2. 按基因分层运行 MR
# 3. 输出主结果、异质性、多效性、单 SNP、leave-one-out 等结果
# 4. 生成"可写入主文 vs 仅补充材料"的审计结论
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(TwoSampleMR)
})

# ---------- 路径 ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "05_MR"

result_dir   <- file.path(root_dir, "03_results")
registry_dir <- file.path(root_dir, "00_registry")
proc_dir     <- file.path(root_dir, "02_processed_data")
fig_dir      <- file.path(root_dir, "04_figures")
log_dir      <- file.path(root_dir, "06_logs")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

harm_file <- file.path(result_dir, "43_harmonised_mr_keep_all_genes.tsv")
tier_file <- file.path(registry_dir, "43_gene_analysis_tier.tsv")

if (!file.exists(harm_file)) stop("缺少 harmonised mr_keep 文件: ", harm_file)
if (!file.exists(tier_file)) stop("缺少 gene analysis tier 文件: ", tier_file)

harm <- fread(harm_file)
tier <- fread(tier_file)

if (!"gene" %in% names(harm) && "id.exposure" %in% names(harm)) {
  harm$gene <- harm$id.exposure
}

if (!"gene" %in% names(harm)) {
  stop("43_harmonised_mr_keep_all_genes.tsv 中缺少 gene / id.exposure 信息。")
}

# ---------- 基础清理 ----------
harm <- harm %>%
  mutate(
    gene = as.character(gene),
    mr_keep = as.logical(mr_keep)
  ) %>%
  filter(mr_keep)

genes <- unique(as.character(harm$gene))
if (length(genes) == 0) stop("43 结果中没有 mr_keep=TRUE 的记录，不能进入 44。")

# ---------- 方法集合 ----------
# 多 SNP 主分析（>= 3 个 IV）
mr_methods_multi <- c(
  "mr_ivw",
  "mr_weighted_median",
  "mr_egger_regression"
)

# ---------- 工具函数 ----------
safe_run <- function(fun, ..., .label = "task") {
  out <- NULL
  err <- ""
  ok <- TRUE
  tryCatch({
    out <- fun(...)
  }, error = function(e) {
    ok  <<- FALSE
    err <<- conditionMessage(e)
  })
  list(ok = ok, out = out, err = err, label = .label)
}

# 给 OR 与 95%CI
add_or_ci <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (!all(c("b", "se") %in% names(df))) return(df)
  
  df %>%
    mutate(
      OR       = exp(b),
      OR_lci95 = exp(b - 1.96 * se),
      OR_uci95 = exp(b + 1.96 * se)
    )
}

# 结果解释分层
classify_gene_for_reporting <- function(n_mr_keep) {
  case_when(
    n_mr_keep >= 3 ~ "main_text_candidate",
    n_mr_keep == 2 ~ "supplementary_exploratory",
    n_mr_keep == 1 ~ "supplementary_wald_ratio_only",
    TRUE           ~ "insufficient"
  )
}

# ---------- 审计表 ----------
gene_input_summary <- harm %>%
  group_by(gene) %>%
  summarise(
    n_mr_keep  = n(),
    exposure   = first(exposure),
    outcome    = first(outcome),
    outcome_id = first(id.outcome),
    .groups    = "drop"
  ) %>%
  mutate(
    reporting_tier = classify_gene_for_reporting(n_mr_keep)
  )

fwrite(
  gene_input_summary,
  file.path(registry_dir, "44_gene_input_summary.tsv"),
  sep = "\t", na = "NA"
)

# ---------- 循环每个基因 ----------
mr_res_list             <- list()
heterogeneity_list      <- list()
pleiotropy_list         <- list()
singlesnp_list          <- list()
loo_list                <- list()
method_audit_list       <- list()
reporting_decision_list <- list()

for (g in genes) {
  
  # ====== 修复 1：每轮循环开头重置，避免上轮残留 ======
  mr_fit   <- NULL
  het_fit  <- NULL
  pleio_fit <- NULL
  ss_fit   <- NULL
  loo_fit  <- NULL
  
  dat_g <- harm %>% filter(gene == g)
  n_iv  <- nrow(dat_g)
  report_tier <- classify_gene_for_reporting(n_iv)
  
  # 保存每基因输入
  fwrite(
    dat_g,
    file.path(proc_dir, paste0("44_", g, "_harmonised_input.tsv")),
    sep = "\t", na = "NA"
  )
  
  # 方法选择
  if (n_iv >= 3) {
    methods_use <- mr_methods_multi
  } else if (n_iv == 2) {
    methods_use <- c("mr_ivw")
  } else if (n_iv == 1) {
    methods_use <- c("mr_wald_ratio")
  } else {
    methods_use <- character(0)
  }
  
  method_audit_list[[g]] <- data.frame(
    gene             = g,
    n_mr_keep        = n_iv,
    methods_selected = paste(methods_use, collapse = ";"),
    reporting_tier   = report_tier,
    stringsAsFactors = FALSE
  )
  
  if (length(methods_use) == 0) next
  
  # ---------- 主 MR ----------
  mr_fit <- safe_run(
    TwoSampleMR::mr,
    dat         = dat_g,
    method_list = methods_use,
    .label      = paste0("mr_", g)
  )
  
  if (mr_fit$ok && !is.null(mr_fit$out) && nrow(mr_fit$out) > 0) {
    x <- mr_fit$out %>%
      mutate(gene = g) %>%
      add_or_ci()
    mr_res_list[[g]] <- x
    
    fwrite(
      x,
      file.path(proc_dir, paste0("44_", g, "_mr_results.tsv")),
      sep = "\t", na = "NA"
    )
  } else {
    writeLines(
      c(
        paste0("MR failed for gene: ", g),
        paste0("Created at: ", Sys.time()),
        paste0("Error: ", mr_fit$err)
      ),
      file.path(log_dir, paste0("44_mr_error_", g, ".txt"))
    )
  }
  
  # ---------- 灵敏度分析（仅 >= 3 IV）----------
  if (n_iv >= 3) {
    
    # --- 异质性 ---
    het_fit <- safe_run(
      TwoSampleMR::mr_heterogeneity,
      dat         = dat_g,
      method_list = c("mr_ivw", "mr_egger_regression"),
      .label      = paste0("heterogeneity_", g)
    )
    
    if (het_fit$ok && !is.null(het_fit$out) && nrow(het_fit$out) > 0) {
      x <- het_fit$out %>% mutate(gene = g)
      heterogeneity_list[[g]] <- x
      fwrite(
        x,
        file.path(proc_dir, paste0("44_", g, "_heterogeneity.tsv")),
        sep = "\t", na = "NA"
      )
    }
    
    # --- 多效性 ---
    pleio_fit <- safe_run(
      TwoSampleMR::mr_pleiotropy_test,
      dat    = dat_g,
      .label = paste0("pleiotropy_", g)
    )
    
    if (pleio_fit$ok && !is.null(pleio_fit$out) && nrow(pleio_fit$out) > 0) {
      x <- pleio_fit$out %>% mutate(gene = g)
      pleiotropy_list[[g]] <- x
      fwrite(
        x,
        file.path(proc_dir, paste0("44_", g, "_pleiotropy.tsv")),
        sep = "\t", na = "NA"
      )
    }
    
    # --- single SNP ---
    ss_fit <- safe_run(
      TwoSampleMR::mr_singlesnp,
      dat    = dat_g,
      .label = paste0("singlesnp_", g)
    )
    
    if (ss_fit$ok && !is.null(ss_fit$out) && nrow(ss_fit$out) > 0) {
      x <- ss_fit$out %>%
        mutate(gene = g) %>%
        add_or_ci()
      singlesnp_list[[g]] <- x
      fwrite(
        x,
        file.path(proc_dir, paste0("44_", g, "_singlesnp.tsv")),
        sep = "\t", na = "NA"
      )
    }
    
    # --- leave-one-out ---
    loo_fit <- safe_run(
      TwoSampleMR::mr_leaveoneout,
      dat    = dat_g,
      .label = paste0("leaveoneout_", g)
    )
    
    if (loo_fit$ok && !is.null(loo_fit$out) && nrow(loo_fit$out) > 0) {
      x <- loo_fit$out %>%
        mutate(gene = g) %>%
        add_or_ci()
      loo_list[[g]] <- x
      fwrite(
        x,
        file.path(proc_dir, paste0("44_", g, "_leaveoneout.tsv")),
        sep = "\t", na = "NA"
      )
    }
    
    # ========== 图 ==========
    
    # ====== 修复 2：scatter plot 第二参数应为 harmonised data (dat_g) ======
    if (!is.null(mr_fit) && mr_fit$ok &&
        !is.null(mr_fit$out) && nrow(mr_fit$out) > 0) {
      sc <- tryCatch(
        TwoSampleMR::mr_scatter_plot(mr_fit$out, dat_g),
        error = function(e) NULL
      )
      if (!is.null(sc) && length(sc) >= 1) {
        pdf(file.path(fig_dir, paste0("44_", g, "_scatter.pdf")),
            width = 7, height = 6)
        print(sc[[1]])
        dev.off()
      }
    }
    
    # forest（依赖 single SNP 结果）
    if (!is.null(ss_fit) && ss_fit$ok &&
        !is.null(ss_fit$out) && nrow(ss_fit$out) > 0) {
      fp <- tryCatch(
        TwoSampleMR::mr_forest_plot(ss_fit$out),
        error = function(e) NULL
      )
      if (!is.null(fp) && length(fp) >= 1) {
        pdf(file.path(fig_dir, paste0("44_", g, "_forest.pdf")),
            width = 7, height = 8)
        print(fp[[1]])
        dev.off()
      }
    }
    
    # leave-one-out plot
    if (!is.null(loo_fit) && loo_fit$ok &&
        !is.null(loo_fit$out) && nrow(loo_fit$out) > 0) {
      lp <- tryCatch(
        TwoSampleMR::mr_leaveoneout_plot(loo_fit$out),
        error = function(e) NULL
      )
      if (!is.null(lp) && length(lp) >= 1) {
        pdf(file.path(fig_dir, paste0("44_", g, "_leaveoneout.pdf")),
            width = 7, height = 8)
        print(lp[[1]])
        dev.off()
      }
    }
    
    # funnel（依赖 single SNP 结果）
    if (!is.null(ss_fit) && ss_fit$ok &&
        !is.null(ss_fit$out) && nrow(ss_fit$out) > 0) {
      funp <- tryCatch(
        TwoSampleMR::mr_funnel_plot(ss_fit$out),
        error = function(e) NULL
      )
      if (!is.null(funp) && length(funp) >= 1) {
        pdf(file.path(fig_dir, paste0("44_", g, "_funnel.pdf")),
            width = 7, height = 8)
        print(funp[[1]])
        dev.off()
      }
    }
    
  } # end if (n_iv >= 3)
  
  # ---------- 可写入主文判断 ----------
  reporting_decision_list[[g]] <- data.frame(
    gene                = g,
    n_mr_keep           = n_iv,
    reporting_tier      = report_tier,
    recommended_section = case_when(
      n_iv >= 3 ~ "main_results_or_main_supplement",
      n_iv == 2 ~ "supplementary_results",
      n_iv == 1 ~ "supplementary_wald_ratio_only",
      TRUE      ~ "exclude"
    ),
    stringsAsFactors = FALSE
  )
  
} # end for (g in genes)

# ---------- 汇总输出 ----------
mr_res_tab        <- bind_rows(mr_res_list)
heterogeneity_tab <- bind_rows(heterogeneity_list)
pleiotropy_tab    <- bind_rows(pleiotropy_list)
singlesnp_tab     <- bind_rows(singlesnp_list)
loo_tab           <- bind_rows(loo_list)
method_audit_tab  <- bind_rows(method_audit_list)
reporting_tab     <- bind_rows(reporting_decision_list)

if (nrow(mr_res_tab) > 0) {
  fwrite(
    mr_res_tab,
    file.path(result_dir, "44_mr_main_results.tsv"),
    sep = "\t", na = "NA"
  )
}

if (nrow(heterogeneity_tab) > 0) {
  fwrite(
    heterogeneity_tab,
    file.path(result_dir, "44_mr_heterogeneity.tsv"),
    sep = "\t", na = "NA"
  )
}

if (nrow(pleiotropy_tab) > 0) {
  fwrite(
    pleiotropy_tab,
    file.path(result_dir, "44_mr_pleiotropy.tsv"),
    sep = "\t", na = "NA"
  )
}

if (nrow(singlesnp_tab) > 0) {
  fwrite(
    singlesnp_tab,
    file.path(result_dir, "44_mr_singlesnp.tsv"),
    sep = "\t", na = "NA"
  )
}

if (nrow(loo_tab) > 0) {
  fwrite(
    loo_tab,
    file.path(result_dir, "44_mr_leaveoneout.tsv"),
    sep = "\t", na = "NA"
  )
}

fwrite(
  method_audit_tab,
  file.path(registry_dir, "44_method_audit.tsv"),
  sep = "\t", na = "NA"
)

fwrite(
  reporting_tab,
  file.path(registry_dir, "44_reporting_decision.tsv"),
  sep = "\t", na = "NA"
)

# ---------- 主文优先表 ----------
main_text_candidates <- mr_res_tab %>%
  filter(gene %in% reporting_tab$gene[reporting_tab$reporting_tier == "main_text_candidate"])

if (nrow(main_text_candidates) > 0) {
  fwrite(
    main_text_candidates,
    file.path(result_dir, "44_main_text_candidate_results.tsv"),
    sep = "\t", na = "NA"
  )
}

# ---------- 总状态 ----------
status_summary <- data.frame(
  item = c(
    "n_genes_input",
    "n_genes_with_mr_results",
    "n_genes_multi_snp_priority",
    "n_genes_with_heterogeneity_output",
    "n_genes_with_pleiotropy_output",
    "recommended_next_action"
  ),
  value = c(
    length(genes),
    length(unique(mr_res_tab$gene)),
    sum(reporting_tab$reporting_tier == "main_text_candidate", na.rm = TRUE),
    length(unique(heterogeneity_tab$gene)),
    length(unique(pleiotropy_tab$gene)),
    ifelse(
      nrow(mr_res_tab) > 0,
      "Proceed to 45: interpret MR results, identify lead causal candidates, and prepare result text/figures",
      "No MR results generated: inspect per-gene logs in 06_logs"
    )
  ),
  stringsAsFactors = FALSE
)

fwrite(
  status_summary,
  file.path(registry_dir, "44_status_summary.tsv"),
  sep = "\t", na = "NA"
)

# ---------- 日志 ----------
log_lines <- c(
  "=== Code 44 Summary ===",
  paste0("Created at: ", Sys.time()),
  "",
  "[Goal]",
  "Run MR main analysis and sensitivity analyses based on harmonised mr_keep data.",
  "",
  "[Outputs]",
  "- 44_mr_main_results.tsv",
  "- 44_mr_heterogeneity.tsv",
  "- 44_mr_pleiotropy.tsv",
  "- 44_mr_singlesnp.tsv",
  "- 44_mr_leaveoneout.tsv",
  "- 44_method_audit.tsv",
  "- 44_reporting_decision.tsv",
  "- 44_main_text_candidate_results.tsv",
  "- 44_status_summary.tsv",
  "",
  "[Figures for multi-SNP genes]",
  "- 44_<GENE>_scatter.pdf",
  "- 44_<GENE>_forest.pdf",
  "- 44_<GENE>_leaveoneout.pdf",
  "- 44_<GENE>_funnel.pdf"
)

writeLines(
  log_lines,
  file.path(log_dir, "44_run_mr_main_and_sensitivity_notes.txt")
)

cat("44_run_mr_main_and_sensitivity_analyses.R finished successfully.\n")
cat("Generated files:\n")
cat("- 00_registry/44_method_audit.tsv\n")
cat("- 00_registry/44_reporting_decision.tsv\n")
cat("- 00_registry/44_status_summary.tsv\n")
cat("- 03_results/44_mr_main_results.tsv\n")
cat("- 03_results/44_mr_heterogeneity.tsv\n")
cat("- 03_results/44_mr_pleiotropy.tsv\n")
cat("- 03_results/44_mr_singlesnp.tsv\n")
cat("- 03_results/44_mr_leaveoneout.tsv\n")
cat("- 03_results/44_main_text_candidate_results.tsv\n")
cat("- 04_figures/44_<GENE>_scatter.pdf\n")
cat("- 04_figures/44_<GENE>_forest.pdf\n")
cat("- 04_figures/44_<GENE>_leaveoneout.pdf\n")
cat("- 04_figures/44_<GENE>_funnel.pdf\n")