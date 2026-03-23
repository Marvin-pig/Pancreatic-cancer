# =========================
# 12_prognostic_gene_screening_fix16.R
# 目的：
# 1. 统一 expr 与 survival 的 TCGA sample_id（前16位）
# 2. 整合 MCD 主线候选基因池
# 3. 进行单因素 Cox 批量筛选
# 4. 输出 shortlist，供 13_lasso_cox_model_train.R 使用
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(survival)
  library(broom)
  library(ggplot2)
  library(stringr)
})

# ---------- 路径 ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")
surv_file <- "04_bulk_analysis/03_survival_model/01_tcga_strict_pdac_survival_ready.tsv"
expr_file <- "02_processed_data/TCGA_PAAD/tcga_paad_expr_strict_pdac_symbol_counts_filtered.rds"

mech_file1 <- "04_bulk_analysis/09_maintext/tcga_paad_mcd_maintext_mechanism_gene_table.tsv"
mech_file2 <- "04_bulk_analysis/08_mechanism/tcga_paad_mcd_key_mechanism_gene_candidates.tsv"

out_dir <- "04_bulk_analysis/03_survival_model"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 工具函数 ----------
tcga16 <- function(x) substr(as.character(x), 1, 16)

extract_gene_symbols <- function(df) {
  cn <- colnames(df)
  gene_col_candidates <- c(
    "gene", "Gene", "gene_symbol", "GeneSymbol", "symbol", "SYMBOL",
    "external_gene_name", "hgnc_symbol", "Gene.symbol", "gene.symbol"
  )
  gene_col <- intersect(gene_col_candidates, cn)
  if (length(gene_col) == 0) return(character(0))
  
  genes <- as.character(df[[gene_col[1]]])
  genes <- str_trim(genes)
  genes <- genes[genes != "" & !is.na(genes)]
  unique(genes)
}

run_univcox_one_gene <- function(gene, expr_mat, clin_df) {
  vec <- as.numeric(expr_mat[gene, ])
  
  tmp <- clin_df %>%
    transmute(
      sample_id = sample_id,
      OS_time = as.numeric(OS_time),
      OS_event = as.numeric(OS_event),
      gene_expr = vec
    ) %>%
    filter(!is.na(OS_time), !is.na(OS_event), !is.na(gene_expr))
  
  if (nrow(tmp) < 30) return(NULL)
  if (sd(tmp$gene_expr) == 0) return(NULL)
  
  fit <- tryCatch(
    coxph(Surv(OS_time, OS_event) ~ gene_expr, data = tmp),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  
  res <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  if (nrow(res) == 0) return(NULL)
  
  data.frame(
    gene = gene,
    HR = res$estimate[1],
    CI_low = res$conf.low[1],
    CI_high = res$conf.high[1],
    p.value = res$p.value[1],
    statistic = res$statistic[1],
    stringsAsFactors = FALSE
  )
}

# ---------- 读取 ----------
surv_dat <- fread(surv_file)
expr <- readRDS(expr_file)

if (!is.matrix(expr)) {
  expr <- as.matrix(expr)
}

cat("surv_dat dim:", dim(surv_dat), "\n")
cat("expr dim:", dim(expr), "\n")

# ---------- 检查表达矩阵 ----------
if (is.null(rownames(expr))) {
  stop("expr 没有 rownames，无法识别基因名。")
}
if (is.null(colnames(expr))) {
  stop("expr 没有 colnames，无法识别样本名。")
}

# 行应该是基因，列应该是样本
if (sum(grepl("^TCGA-", rownames(expr))) > sum(grepl("^TCGA-", colnames(expr)))) {
  stop("expr 方向异常：看起来样本在行、基因在列。请先转置。")
}

# ---------- 统一 sample ID 到前16位 ----------
surv_dat$sample_id_16 <- tcga16(surv_dat$sample_id)
expr_sample_16 <- tcga16(colnames(expr))

common_samples <- intersect(expr_sample_16, surv_dat$sample_id_16)

cat("matched samples after 16-char harmonization:", length(common_samples), "\n")

if (length(common_samples) < 30) {
  stop(paste0("统一到前16位后仍匹配样本过少：n=", length(common_samples)))
}

# survival 去重并重排
surv_dat <- surv_dat %>%
  filter(sample_id_16 %in% common_samples) %>%
  distinct(sample_id_16, .keep_all = TRUE)

surv_dat <- surv_dat[match(common_samples, surv_dat$sample_id_16), ]

# expr 处理重复 sample
keep_expr <- !duplicated(expr_sample_16) & expr_sample_16 %in% common_samples
expr <- expr[, keep_expr, drop = FALSE]
expr_sample_16 <- tcga16(colnames(expr))

expr <- expr[, match(common_samples, expr_sample_16), drop = FALSE]
expr_sample_16 <- tcga16(colnames(expr))

stopifnot(all(expr_sample_16 == surv_dat$sample_id_16))

cat("final matched expr dim:", dim(expr), "\n")
cat("final survival dim:", dim(surv_dat), "\n")

# ---------- 读取候选基因 ----------
candidate_genes <- character(0)

if (file.exists(mech_file1)) {
  x1 <- fread(mech_file1)
  candidate_genes <- c(candidate_genes, extract_gene_symbols(x1))
}

if (file.exists(mech_file2)) {
  x2 <- fread(mech_file2)
  candidate_genes <- c(candidate_genes, extract_gene_symbols(x2))
}

# 你当前主线可手动补充的核心基因
manual_core_genes <- c(
  "GPX4", "SLC31A1", "ABI1", "FLNA",
  "CXCL9", "CXCL10", "CXCL11", "FGL2",
  "IL7R", "CMKLR1", "NOTCH2", "INHBA", "TGFBR3"
)

candidate_genes <- unique(c(candidate_genes, manual_core_genes))
candidate_genes <- intersect(candidate_genes, rownames(expr))

cat("candidate genes found in expr:", length(candidate_genes), "\n")

if (length(candidate_genes) < 5) {
  stop("候选基因太少。请检查机制表中的基因列名，或把机制表列名贴给我。")
}

candidate_gene_pool <- data.frame(
  gene = candidate_genes,
  stringsAsFactors = FALSE
)

fwrite(candidate_gene_pool,
       file.path(out_dir, "12_candidate_gene_pool.tsv"),
       sep = "\t")

# ---------- 单因素 Cox ----------
cox_list <- lapply(candidate_genes, run_univcox_one_gene, expr_mat = expr, clin_df = surv_dat)
cox_res <- bind_rows(cox_list)

if (nrow(cox_res) == 0) {
  stop("没有成功得到单因素 Cox 结果。")
}

cox_res <- cox_res %>%
  arrange(p.value) %>%
  mutate(FDR = p.adjust(p.value, method = "BH"))

fwrite(cox_res,
       file.path(out_dir, "12_univcox_candidate_genes.tsv"),
       sep = "\t")

# ---------- shortlist ----------
shortlist <- cox_res %>%
  filter(p.value < 0.1)

if (nrow(shortlist) < 8) {
  shortlist <- cox_res %>% slice_head(n = min(15, nrow(cox_res)))
}

if (nrow(shortlist) > 25) {
  shortlist <- shortlist %>% slice_head(n = 25)
}

fwrite(shortlist,
       file.path(out_dir, "12_prognostic_gene_shortlist.tsv"),
       sep = "\t")

print(shortlist)

# ---------- 绘图 ----------
plot_df <- cox_res %>%
  slice_head(n = min(15, nrow(cox_res))) %>%
  mutate(gene = factor(gene, levels = rev(gene)))

pdf(file.path(out_dir, "12_univcox_forest_top.pdf"), width = 7, height = 5.5)
ggplot(plot_df, aes(x = gene, y = HR, ymin = CI_low, ymax = CI_high)) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = 2) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    x = NULL,
    y = "Hazard ratio (univariable Cox)",
    title = "Top candidate prognostic genes in MCD axis"
  )
dev.off()

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "12_sessionInfo.txt"))

cat("\n12_prognostic_gene_screening_fix16.R finished successfully.\n")