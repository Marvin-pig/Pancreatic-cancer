# =========================
# 11_mcd_subtype_km_cox.R
# 目的：
# 1. 读取 survival-ready 数据
# 2. 统一 MCD subtype 命名
# 3. 绘制 OS Kaplan-Meier 曲线（含两两 pairwise log-rank）
# 4. 进行 log-rank 检验
# 5. 进行单因素与多因素 Cox 回归（含 forest plot）
# 6. 输出投稿级基础结果表
# =========================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(broom)
})

# ---------- 文件路径 ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")
infile  <- "04_bulk_analysis/03_survival_model/01_tcga_strict_pdac_survival_ready.tsv"
out_dir <- "04_bulk_analysis/03_survival_model"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 读取数据 ----------
dat <- fread(infile) %>% as.data.frame()   # 统一为 data.frame，避免后续 data.table/dplyr 混用问题
cat("Input dim:", dim(dat), "\n")

# ---------- 基本检查 ----------
required_cols <- c("sample_id", "OS_time", "OS_event", "mcd_subtype")
missing_cols <- setdiff(required_cols, colnames(dat))
if (length(missing_cols) > 0) {
  stop("缺少必需列: ", paste(missing_cols, collapse = ", "))
}

# ---------- 检查 age 列是否存在 ----------
age_col <- "age_at_diagnosis_years"
if (!age_col %in% colnames(dat)) {
  # 尝试常见备用列名
  alt_age <- intersect(c("age", "age_at_diagnosis", "age_years"), colnames(dat))
  if (length(alt_age) > 0) {
    age_col <- alt_age[1]
    warning("列 'age_at_diagnosis_years' 不存在，改用列: ", age_col)
  } else {
    stop("找不到年龄列，请检查列名。当前列名: ", paste(colnames(dat), collapse = ", "))
  }
}

# ---------- subtype 命名映射 ----------
# 根据聚类结果确认 C1/C2/C3 对应 MCD-low / intermediate / high 后修改此处
subtype_map <- c(
  "C1" = "C1",
  "C2" = "C2",
  "C3" = "C3"
)
dat <- dat %>%
  mutate(
    mcd_subtype_raw = mcd_subtype,
    mcd_subtype     = dplyr::recode(as.character(mcd_subtype), !!!subtype_map),
    mcd_subtype     = factor(mcd_subtype, levels = c("C1", "C2", "C3"))
  )

# ---------- 临床变量清洗 ----------
na_vals <- c("", "NA", "N/A", "Unknown", "not reported", "Not Reported", "[Not Available]")

clean_var <- function(x) {
  x <- str_trim(as.character(x))
  x[x %in% na_vals] <- NA
  x
}

dat <- dat %>%
  mutate(
    age      = suppressWarnings(as.numeric(.data[[age_col]])),
    gender   = clean_var(gender_final),
    stage    = clean_var(ajcc_stage_final),
    grade    = clean_var(tumor_grade_final),
    OS_time  = as.numeric(OS_time),
    OS_event = as.numeric(OS_event)
  )

# 性别标准化
dat$gender <- dplyr::case_when(
  tolower(dat$gender) %in% c("male", "m")     ~ "Male",
  tolower(dat$gender) %in% c("female", "f")   ~ "Female",
  TRUE ~ NA_character_
)
dat$gender <- factor(dat$gender, levels = c("Male", "Female"))

# stage 有序因子（罗马数字 + 字母亚期均覆盖）
stage_levels <- c(
  "Stage I", "Stage IA", "Stage IB",
  "Stage II", "Stage IIA", "Stage IIB",
  "Stage III", "Stage IIIA", "Stage IIIB",
  "Stage IV"
)
present_stages <- intersect(stage_levels, unique(na.omit(dat$stage)))
dat$stage <- factor(dat$stage, levels = present_stages, ordered = TRUE)

# grade 有序因子
grade_levels <- c("G1", "G2", "G3", "G4",
                  "Well differentiated", "Moderately differentiated",
                  "Poorly differentiated", "Undifferentiated")
present_grades <- intersect(grade_levels, unique(na.omit(dat$grade)))
dat$grade <- factor(dat$grade, levels = present_grades, ordered = TRUE)

# ---------- 生存摘要 ----------
surv_summary <- dat %>%
  group_by(mcd_subtype) %>%
  summarise(
    n              = n(),
    events         = sum(OS_event, na.rm = TRUE),
    event_rate     = round(mean(OS_event, na.rm = TRUE), 3),
    median_OS_days = median(OS_time, na.rm = TRUE),
    .groups        = "drop"
  )
fwrite(surv_summary,
       file.path(out_dir, "11_survival_summary_by_subtype.tsv"),
       sep = "\t")
print(surv_summary)

# ---------- KM 曲线 ----------
fit_km <- survfit(Surv(OS_time, OS_event) ~ mcd_subtype, data = dat)

# 色盲友好配色（okabe-ito 前3色）
subtype_colors <- c("#E69F00", "#56B4E9", "#009E73")

p_km <- ggsurvplot(
  fit_km,
  data             = dat,
  pval             = TRUE,
  pval.method      = TRUE,
  risk.table       = TRUE,
  conf.int         = FALSE,
  censor           = TRUE,
  surv.median.line = "hv",
  xlab             = "Time (days)",
  ylab             = "Overall survival probability",
  legend.title     = "MCD subtype",
  legend.labs      = levels(dat$mcd_subtype),
  palette          = subtype_colors,
  risk.table.height = 0.28,
  ggtheme          = theme_bw(base_size = 12)
)

# Bug fix: ggsurvplot 返回列表，用 arrange_ggsurvplots 或手动拼图保存
km_combined <- cowplot::plot_grid(
  p_km$plot,
  p_km$table,
  ncol   = 1,
  rel_heights = c(3, 1)
)
ggsave(
  file.path(out_dir, "11_km_os_by_mcd_subtype.pdf"),
  plot   = km_combined,
  width  = 7.5,
  height = 7,
  device = "pdf"
)

# ---------- log-rank 检验 ----------
logrank_fit <- survdiff(Surv(OS_time, OS_event) ~ mcd_subtype, data = dat)
logrank_p   <- 1 - pchisq(logrank_fit$chisq, df = length(logrank_fit$n) - 1)
logrank_res <- data.frame(
  chisq   = logrank_fit$chisq,
  df      = length(logrank_fit$n) - 1,
  p_value = logrank_p
)
fwrite(logrank_res,
       file.path(out_dir, "11_logrank_os_by_mcd_subtype.tsv"),
       sep = "\t")
print(logrank_res)

# ---------- 两两 pairwise log-rank（3组时必要） ----------
pairwise_res <- pairwise_survdiff(
  Surv(OS_time, OS_event) ~ mcd_subtype,
  data     = dat,
  p.adjust.method = "BH"   # Benjamini-Hochberg 校正
)
pairwise_pmat <- as.data.frame(pairwise_res$p.value)
pairwise_pmat$comparison <- rownames(pairwise_pmat)
fwrite(pairwise_pmat,
       file.path(out_dir, "11_pairwise_logrank_bh.tsv"),
       sep = "\t")
cat("\nPairwise log-rank (BH adjusted):\n")
print(pairwise_res$p.value)

# ---------- 单因素 Cox ----------
cox_uni     <- coxph(Surv(OS_time, OS_event) ~ mcd_subtype, data = dat)
cox_uni_res <- broom::tidy(cox_uni, exponentiate = TRUE, conf.int = TRUE) %>%
  dplyr::rename(HR = estimate, CI_low = conf.low, CI_high = conf.high)
fwrite(cox_uni_res,
       file.path(out_dir, "11_univcox_mcd_subtype.tsv"),
       sep = "\t")
print(cox_uni_res)

# ---------- 多因素 Cox ----------
vars_for_model <- c("OS_time", "OS_event", "mcd_subtype", "age", "gender", "stage", "grade")

# Bug fix: 使用 dplyr::select 代替 data.table 的 [, ..vars] 语法
dat_model <- dat %>% dplyr::select(dplyr::all_of(vars_for_model))

na_summary <- sapply(dat_model, function(x) sum(is.na(x)))
cat("\nNA summary per variable:\n")
print(na_summary)

dat_cc <- na.omit(dat_model)
cat("Complete-case n for multivariable Cox:", nrow(dat_cc), "\n")

if (nrow(dat_cc) >= 80) {
  cox_multi         <- coxph(Surv(OS_time, OS_event) ~ mcd_subtype + age + gender + stage + grade,
                              data = dat_cc)
  multi_formula_used <- "mcd_subtype + age + gender + stage + grade"
} else {
  dat_cc <- dat %>%
    dplyr::select(OS_time, OS_event, mcd_subtype, age, gender) %>%
    na.omit()
  cox_multi         <- coxph(Surv(OS_time, OS_event) ~ mcd_subtype + age + gender,
                              data = dat_cc)
  multi_formula_used <- "mcd_subtype + age + gender"
  warning("完整案例 <80，已降阶为简化模型: ", multi_formula_used)
}

cox_multi_res <- broom::tidy(cox_multi, exponentiate = TRUE, conf.int = TRUE) %>%
  dplyr::rename(HR = estimate, CI_low = conf.low, CI_high = conf.high) %>%
  dplyr::mutate(model_formula = multi_formula_used)
fwrite(cox_multi_res,
       file.path(out_dir, "11_multivcox_mcd_subtype.tsv"),
       sep = "\t")
print(cox_multi_res)

# ---------- Forest plot（多因素 Cox） ----------
forest_p <- ggforest(
  cox_multi,
  data      = dat_cc,
  main      = paste("Multivariable Cox:", multi_formula_used),
  cpositions = c(0.02, 0.22, 0.4),
  fontsize  = 0.8
)
ggsave(
  file.path(out_dir, "11_forestplot_multivcox.pdf"),
  plot   = forest_p,
  width  = 10,
  height = max(4, 0.4 * nrow(cox_multi_res) + 2),
  device = "pdf"
)

# ---------- 保存 session 信息 ----------
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "11_sessionInfo.txt"))

cat("\n11_mcd_subtype_km_cox.R finished successfully.\n")
