# =========================
# 20_external_validation_refinement_and_publication_ready_interpretation.R
# 目的：
# 1. 对 GSE21501 外部验证做优化解释
# 2. 保留训练集 cutoff 结果
# 3. 增加外部中位数 cutoff 敏感性分析
# 4. 计算 12/24/36 月 timeROC
# 5. 比较训练集与外部 risk score 分布
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(survival)
  library(survminer)
  library(broom)
  library(timeROC)
  library(ggplot2)
})
setwd("/Users/wmz_mac/Desktop/胰腺癌")
risk_ext_file <- "04_bulk_analysis/04_external_validation/04_signature_validation/GSE21501/19_GSE21501_external_risk_score.tsv"
train_risk_file <- "04_bulk_analysis/03_survival_model/13_risk_score_training.tsv"

out_dir <- "04_bulk_analysis/04_external_validation/04_signature_validation/GSE21501"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

risk_ext <- fread(risk_ext_file)
train_risk <- fread(train_risk_file)

stopifnot(all(c("OS_time", "OS_event", "risk_score") %in% colnames(risk_ext)))
stopifnot("risk_score" %in% colnames(train_risk))

# ---------- 1. 训练集 vs 外部 risk score 分布 ----------
dist_tab <- data.frame(
  cohort = c("TCGA_train", "GSE21501_external"),
  n = c(nrow(train_risk), nrow(risk_ext)),
  mean = c(mean(train_risk$risk_score, na.rm = TRUE), mean(risk_ext$risk_score, na.rm = TRUE)),
  sd = c(sd(train_risk$risk_score, na.rm = TRUE), sd(risk_ext$risk_score, na.rm = TRUE)),
  median = c(median(train_risk$risk_score, na.rm = TRUE), median(risk_ext$risk_score, na.rm = TRUE)),
  q25 = c(quantile(train_risk$risk_score, 0.25, na.rm = TRUE), quantile(risk_ext$risk_score, 0.25, na.rm = TRUE)),
  q75 = c(quantile(train_risk$risk_score, 0.75, na.rm = TRUE), quantile(risk_ext$risk_score, 0.75, na.rm = TRUE))
)

fwrite(dist_tab, file.path(out_dir, "20_risk_score_distribution_train_vs_external.tsv"), sep = "\t")

plot_df <- bind_rows(
  data.frame(cohort = "TCGA_train", risk_score = train_risk$risk_score),
  data.frame(cohort = "GSE21501_external", risk_score = risk_ext$risk_score)
)

pdf(file.path(out_dir, "20_risk_score_distribution_train_vs_external.pdf"), width = 6.5, height = 5.2)
ggplot(plot_df, aes(x = cohort, y = risk_score)) +
  geom_boxplot() +
  theme_bw(base_size = 12) +
  labs(title = "Risk score distribution: training vs external")
dev.off()

# ---------- 2. 外部中位数 cutoff 敏感性分析 ----------
ext_cutoff <- median(risk_ext$risk_score, na.rm = TRUE)

risk_ext2 <- risk_ext %>%
  mutate(
    risk_group_ext_median = ifelse(risk_score >= ext_cutoff, "High", "Low")
  )

risk_ext2$risk_group_ext_median <- factor(risk_ext2$risk_group_ext_median, levels = c("Low", "High"))

fwrite(risk_ext2,
       file.path(out_dir, "20_GSE21501_external_risk_score_with_ext_median.tsv"),
       sep = "\t")

# KM: external median cutoff
fit_km_ext <- survfit(Surv(OS_time, OS_event) ~ risk_group_ext_median, data = risk_ext2)

pdf(file.path(out_dir, "20_GSE21501_km_external_median_cutoff.pdf"), width = 7.2, height = 6.2)
p <- ggsurvplot(
  fit_km_ext,
  data = risk_ext2,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  xlab = "Overall survival time (months)",
  ylab = "Overall survival probability",
  legend.title = "Risk group",
  legend.labs = c("Low", "High"),
  risk.table.height = 0.25,
  ggtheme = theme_bw(base_size = 12)
)
print(p)
dev.off()

cox_group_ext <- coxph(Surv(OS_time, OS_event) ~ risk_group_ext_median, data = risk_ext2)
cox_group_ext_res <- broom::tidy(cox_group_ext, exponentiate = TRUE, conf.int = TRUE)

fwrite(cox_group_ext_res,
       file.path(out_dir, "20_GSE21501_univcox_external_median_group.tsv"),
       sep = "\t")

# ---------- 3. 连续型 risk_score 再次确认 ----------
cox_score <- coxph(Surv(OS_time, OS_event) ~ risk_score, data = risk_ext2)
cox_score_res <- broom::tidy(cox_score, exponentiate = TRUE, conf.int = TRUE)

fwrite(cox_score_res,
       file.path(out_dir, "20_GSE21501_univcox_risk_score_confirm.tsv"),
       sep = "\t")

# ---------- 4. timeROC：按月 ----------
roc_times <- c(12, 24, 36)

roc_fit <- timeROC(
  T = risk_ext2$OS_time,
  delta = risk_ext2$OS_event,
  marker = risk_ext2$risk_score,
  cause = 1,
  weighting = "marginal",
  times = roc_times,
  iid = TRUE
)

roc_tab <- data.frame(
  time_months = roc_times,
  AUC = roc_fit$AUC
)

fwrite(roc_tab,
       file.path(out_dir, "20_GSE21501_timeROC_months.tsv"),
       sep = "\t")

pdf(file.path(out_dir, "20_GSE21501_timeROC_months.pdf"), width = 6.5, height = 5.5)
plot(roc_fit, time = 12, col = 1, title = FALSE)
plot(roc_fit, time = 24, add = TRUE, col = 2)
plot(roc_fit, time = 36, add = TRUE, col = 3)
legend(
  "bottomright",
  legend = paste0(c("12-month", "24-month", "36-month"), " AUC = ", sprintf("%.3f", roc_fit$AUC)),
  col = 1:3,
  lwd = 2,
  bty = "n"
)
title(main = "Time-dependent ROC in GSE21501")
dev.off()

# ---------- 5. 外部验证解释摘要 ----------
summary_tab <- data.frame(
  metric = c(
    "n_external_samples",
    "n_external_events",
    "training_cutoff",
    "external_median_cutoff",
    "continuous_HR",
    "continuous_p"
  ),
  value = c(
    nrow(risk_ext2),
    sum(risk_ext2$OS_event, na.rm = TRUE),
    median(train_risk$risk_score, na.rm = TRUE),
    ext_cutoff,
    cox_score_res$estimate[1],
    cox_score_res$p.value[1]
  )
)

fwrite(summary_tab,
       file.path(out_dir, "20_GSE21501_refined_validation_summary.tsv"),
       sep = "\t")

cat("20_external_validation_refinement_and_publication_ready_interpretation.R finished successfully.\n")