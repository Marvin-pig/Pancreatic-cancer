# =========================
# 14_internal_validation_and_clinical_independence.R
# 目的：
# 1. 基于 13_risk_score_training.tsv 做内部验证
# 2. 计算 time-dependent ROC 和 C-index
# 3. 做 risk score 的单因素/多因素 Cox
# 4. 规范 stage / grade 编码
# 5. 检查比例风险假设
# 6. 输出最终模型固定信息
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(broom)
  library(timeROC)
})

# ---------- 路径 ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")
risk_file <- "04_bulk_analysis/03_survival_model/13_risk_score_training.tsv"
surv_file <- "04_bulk_analysis/03_survival_model/01_tcga_strict_pdac_survival_ready.tsv"
nz_file   <- "04_bulk_analysis/03_survival_model/13_lasso_nonzero_genes.tsv"
out_dir   <- "04_bulk_analysis/03_survival_model"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 工具函数 ----------
clean_stage <- function(x) {
  x <- str_trim(as.character(x))
  x[x %in% c("", "NA", "N/A", "Unknown", "not reported")] <- NA
  x
}

clean_grade <- function(x) {
  x <- str_trim(as.character(x))
  x[x %in% c("", "NA", "N/A", "Unknown", "not reported")] <- NA
  x
}

# FIX #2：改用更精确的 regex，先匹配 III/IV，再匹配 I/II，
#         避免 "STAGE I" 作为前缀误匹配 "STAGE III"/"STAGE IV"
collapse_stage <- function(x) {
  x   <- toupper(clean_stage(x))
  out <- rep(NA_character_, length(x))
  # 先赋高期别，避免被后续 "STAGE I" 模式误覆盖
  out[grepl("STAGE\\s*III|STAGE\\s*IV|^III$|^IV$", x)] <- "III-IV"
  out[grepl("^STAGE\\s*I[AB]?$|^STAGE\\s*II[AB]?$|^I$|^II$", x)] <- "I-II"
  out
}

collapse_grade <- function(x) {
  x   <- toupper(clean_grade(x))
  out <- rep(NA_character_, length(x))
  out[grepl("G1|G2|LOW|MODERATE|WELL|MODERATELY", x)] <- "G1-2"
  out[grepl("G3|G4|HIGH|POOR|UNDIFFERENTIATED",   x)] <- "G3-4"
  out
}

# ---------- 读取 ----------
risk_dat <- fread(risk_file)
surv_dat <- fread(surv_file)
nz_dat   <- fread(nz_file)

stopifnot(all(c("sample_id", "OS_time", "OS_event", "risk_score", "risk_group") %in% colnames(risk_dat)))
stopifnot("sample_id" %in% colnames(surv_dat))

# ---------- 合并临床信息 ----------
dat <- risk_dat %>%
  left_join(
    surv_dat %>% transmute(
      sample_id   = sample_id,
      age         = suppressWarnings(as.numeric(age_at_diagnosis_years)),
      gender      = as.character(gender_final),
      ajcc_stage  = as.character(ajcc_stage_final),
      tumor_grade = as.character(tumor_grade_final)
    ),
    by = "sample_id"
  ) %>%
  as.data.frame()   # FIX #1 前置：明确转为 data.frame，消除后续 data.table 语法依赖

# ---------- 清洗临床变量 ----------
dat <- dat %>%
  mutate(
    gender = case_when(
      tolower(gender) %in% c("male",   "m") ~ "Male",
      tolower(gender) %in% c("female", "f") ~ "Female",
      TRUE ~ NA_character_
    ),
    stage_collapsed = collapse_stage(ajcc_stage),
    grade_collapsed = collapse_grade(tumor_grade)
  )

dat$gender          <- factor(dat$gender,          levels = c("Male",  "Female"))
dat$stage_collapsed <- factor(dat$stage_collapsed, levels = c("I-II",  "III-IV"))
dat$grade_collapsed <- factor(dat$grade_collapsed, levels = c("G1-2",  "G3-4"))
dat$risk_group      <- factor(dat$risk_group,      levels = c("Low",   "High"))

# ---------- 1. time-dependent ROC ----------
roc_times <- c(365, 730, 1095)

roc_fit <- timeROC(
  T         = dat$OS_time,
  delta     = dat$OS_event,
  marker    = dat$risk_score,
  cause     = 1,
  weighting = "marginal",
  times     = roc_times,
  iid       = TRUE
)

# timeROC 返回的 AUC 是具名数值向量（长度 = length(times)），直接使用即可
auc_vals <- roc_fit$AUC

roc_tab <- data.frame(
  time_days = roc_times,
  AUC       = auc_vals
)
fwrite(roc_tab, file.path(out_dir, "14_timeROC_training.tsv"), sep = "\t")

# FIX #3：去掉非法参数 title=FALSE；用独立的 title() 添加主标题
pdf(file.path(out_dir, "14_timeROC_training.pdf"), width = 6.5, height = 5.5)
plot(roc_fit, time = 365,  col = 1)
plot(roc_fit, time = 730,  add = TRUE, col = 2)
plot(roc_fit, time = 1095, add = TRUE, col = 3)
legend(
  "bottomright",
  legend = paste0(c("1-year", "2-year", "3-year"),
                  " AUC = ", sprintf("%.3f", auc_vals)),
  col = 1:3, lwd = 2, bty = "n"
)
title(main = "Time-dependent ROC for MCD-derived risk score")
dev.off()

# ---------- 2. C-index ----------
cox_cindex <- coxph(Surv(OS_time, OS_event) ~ risk_score, data = dat, x = TRUE, y = TRUE)
sum_cindex <- summary(cox_cindex)

# FIX #4：用具名索引代替位置索引，更健壮
cindex_tab <- data.frame(
  c_index = sum_cindex$concordance["C"],
  se      = sum_cindex$concordance["se(C)"]
)
fwrite(cindex_tab, file.path(out_dir, "14_cindex_training.tsv"), sep = "\t")

# ---------- 3. 单因素 Cox ----------
uni_results <- list()
vars_uni    <- c("risk_score", "risk_group", "age", "gender", "stage_collapsed", "grade_collapsed")

for (v in vars_uni) {
  fml <- as.formula(paste0("Surv(OS_time, OS_event) ~ ", v))
  
  # FIX #1：dat 已是 data.frame，直接用标准列选择，不能用 data.table 的 with=FALSE
  tmp <- dat[, c("OS_time", "OS_event", v)]
  tmp <- na.omit(tmp)
  if (nrow(tmp) < 30) next
  
  fit <- tryCatch(coxph(fml, data = tmp), error = function(e) NULL)
  if (is.null(fit)) next
  
  res <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  res$variable   <- v
  res$n_complete <- nrow(tmp)
  uni_results[[v]] <- res
}

uni_tab <- bind_rows(uni_results)
fwrite(uni_tab, file.path(out_dir, "14_univcox_clinical_and_risk.tsv"), sep = "\t")

# ---------- 4. 多因素 Cox ----------
# 主模型：risk_score + age + gender + stage + grade
vars_full    <- c("OS_time", "OS_event", "risk_score", "age", "gender", "stage_collapsed", "grade_collapsed")
vars_reduced <- c("OS_time", "OS_event", "risk_score", "age", "gender")

multi_dat <- na.omit(dat[, vars_full])
cat("Complete-case n for multivariable Cox (full model):", nrow(multi_dat), "\n")

if (nrow(multi_dat) >= 80) {
  fit_multi <- coxph(
    Surv(OS_time, OS_event) ~ risk_score + age + gender + stage_collapsed + grade_collapsed,
    data = multi_dat
  )
  multi_formula_used <- "Surv(OS_time, OS_event) ~ risk_score + age + gender + stage_collapsed + grade_collapsed"
} else {
  multi_dat <- na.omit(dat[, vars_reduced])
  # FIX #5：补充回退时的 n 提示，便于调试
  cat("Full model n < 80; falling back to reduced model, n =", nrow(multi_dat), "\n")
  fit_multi <- coxph(
    Surv(OS_time, OS_event) ~ risk_score + age + gender,
    data = multi_dat
  )
  multi_formula_used <- "Surv(OS_time, OS_event) ~ risk_score + age + gender"
}

multi_tab <- broom::tidy(fit_multi, exponentiate = TRUE, conf.int = TRUE)
multi_tab$model_formula <- multi_formula_used
multi_tab$n_complete    <- nrow(multi_dat)
fwrite(multi_tab, file.path(out_dir, "14_multivcox_clinical_and_risk.tsv"), sep = "\t")

# ---------- 5. PH 假设 ----------
ph_test <- cox.zph(fit_multi)
ph_tab  <- as.data.frame(ph_test$table)
ph_tab$term <- rownames(ph_tab)
rownames(ph_tab) <- NULL
fwrite(ph_tab, file.path(out_dir, "14_ph_assumption_test.tsv"), sep = "\t")

# ---------- 6. 固定最终模型公式 ----------
# 先默认固定为 13 中 lambda.min 版本，后续若比较 1se 更优再替换
model_formula_tab <- nz_dat %>%
  filter(model == "lambda.min") %>%
  arrange(desc(abs(coefficient))) %>%
  mutate(component = paste0(sprintf("%.6f", coefficient), " * ", gene))

formula_text <- paste(model_formula_tab$component, collapse = " + ")

final_formula_tab <- data.frame(
  model_version      = "lambda.min",
  n_genes            = nrow(model_formula_tab),
  risk_score_formula = paste0("risk_score = ", formula_text),
  cutoff_rule        = "training median",
  stringsAsFactors   = FALSE
)
fwrite(final_formula_tab, file.path(out_dir, "14_final_model_formula.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "14_sessionInfo.txt"))
cat("\n14_internal_validation_and_clinical_independence.R finished successfully.\n")
