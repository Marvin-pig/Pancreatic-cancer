# =========================
# 31_begin_manual_mr_data_attachment_strategy.R
# 目的：
# 1. 读取 30 步生成的预期文件清单
# 2. 自动区分"最小可运行必需文件"和"支持性文件"
# 3. 生成第一批手动挂载优先顺序
# 4. 生成 MR 数据挂载状态表
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ---------- 路径设置 ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "05_MR"

expected_file   <- file.path(root_dir, "02_processed_data/30_batch1_expected_input_files.tsv")
manual_plan_file <- file.path(root_dir, "00_registry/30_manual_data_acquisition_plan.tsv")

# 输出目录（若不存在则创建）
out_registry <- file.path(root_dir, "00_registry")
out_logs     <- file.path(root_dir, "06_logs")
dir.create(out_registry, recursive = TRUE, showWarnings = FALSE)
dir.create(out_logs,     recursive = TRUE, showWarnings = FALSE)

# ---------- 安全读取输入 ----------
if (!file.exists(expected_file)) {
  stop("[ERROR] 30 步输出文件不存在: ", expected_file,
       "\n请先运行 30_*.R 脚本生成该文件。")
}

exp <- fread(expected_file)

# 校验必需列
required_cols <- c("gene", "exposure_file_expected",
                   "gtex_file_expected", "outcome_file_expected",
                   "run_priority")
missing_cols  <- setdiff(required_cols, names(exp))
if (length(missing_cols) > 0) {
  stop("[ERROR] 输入文件缺少必需列: ",
       paste(missing_cols, collapse = ", "))
}

# manual_plan_file 目前仅做存在性记录，暂不参与逻辑
plan_exists <- file.exists(manual_plan_file)
if (!plan_exists) {
  message("[WARN] 手动获取计划文件不存在，跳过: ", manual_plan_file)
}

# ---------- 分类：区分文件角色 ----------
attach_plan <- exp %>%
  mutate(
    file_role = case_when(
      str_detect(exposure_file_expected, "eQTLGen")       ~ "primary_exposure",
      str_detect(gtex_file_expected,     "GTEx_pancreas")  ~ "supporting_exposure",
      TRUE                                                  ~ "other"
    )
  )

# ---------- 整理为长表 ----------
attach_long <- bind_rows(
  # 每个基因的主暴露文件（eQTLGen）
  attach_plan %>%
    transmute(
      gene      = gene,
      file_type = "primary_exposure",
      expected_path = exposure_file_expected,
      priority  = ifelse(run_priority == "HIGH", "HIGH", "MEDIUM")
    ),
  # 每个基因的辅助暴露文件（GTEx）
  attach_plan %>%
    transmute(
      gene      = gene,
      file_type = "supporting_exposure",
      expected_path = gtex_file_expected,
      priority  = "LOW"
    )
)

# 加入唯一 outcome 文件（优先级最高）
outcome_paths <- unique(exp$outcome_file_expected)
outcome_rows  <- data.frame(
  gene          = "ALL_BATCH1",
  file_type     = "primary_outcome",
  expected_path = outcome_paths,
  priority      = "TOP",
  stringsAsFactors = FALSE
)
attach_long <- bind_rows(outcome_rows, attach_long)

# ---------- 检查文件存在性 ----------
attach_long <- attach_long %>%
  mutate(
    exists = file.exists(file.path(root_dir, expected_path))
  )

# ---------- 最小可运行条件 ----------
primary_exposure_ready <- any(
  attach_long$file_type == "primary_exposure" & attach_long$exists
)
primary_outcome_ready <- any(
  attach_long$file_type == "primary_outcome" & attach_long$exists
)
minimum_runnable <- primary_exposure_ready & primary_outcome_ready

status_tab <- data.frame(
  item = c(
    "primary_exposure_any_ready",
    "primary_outcome_ready",
    "minimum_MR_runnable"
  ),
  status = c(
    primary_exposure_ready,
    primary_outcome_ready,
    minimum_runnable
  ),
  stringsAsFactors = FALSE
)

# ---------- 推荐优先顺序 ----------
priority_order <- attach_long %>%
  mutate(
    priority_rank = case_when(
      priority == "TOP"    ~ 1L,
      priority == "HIGH"   ~ 2L,
      priority == "MEDIUM" ~ 3L,
      TRUE                 ~ 4L
    )
  ) %>%
  arrange(priority_rank, gene, file_type)

# ---------- 统计摘要（写入日志） ----------
n_total    <- nrow(attach_long)
n_exists   <- sum(attach_long$exists)
n_missing  <- n_total - n_exists
n_high     <- sum(attach_long$priority %in% c("TOP", "HIGH") & !attach_long$exists)

# ---------- 输出 ----------
fwrite(
  attach_long,
  file.path(out_registry, "31_mr_data_attachment_status.tsv"),
  sep = "\t"
)

fwrite(
  priority_order,
  file.path(out_registry, "31_mr_data_attachment_priority.tsv"),
  sep = "\t"
)

fwrite(
  status_tab,
  file.path(out_registry, "31_mr_minimum_runnable_status.tsv"),
  sep = "\t"
)

# ---------- 日志 ----------
notes <- c(
  "=== MR Manual Data Attachment Strategy ===",
  paste0("Created at: ", Sys.time()),
  "",
  "[Summary]",
  paste0("  Total file slots:   ", n_total),
  paste0("  Already present:    ", n_exists),
  paste0("  Still missing:      ", n_missing),
  paste0("  HIGH/TOP & missing: ", n_high),
  "",
  "[Runnable Status]",
  paste0("  Primary exposure ready : ", primary_exposure_ready),
  paste0("  Primary outcome ready  : ", primary_outcome_ready),
  paste0("  Minimum MR runnable    : ", minimum_runnable),
  "",
  "[Rule]",
  "  Do NOT start main MR until minimum_MR_runnable becomes TRUE.",
  "",
  "[Priority Order for Manual Data Placement]",
  "  1. Place PanScan/PanC4 outcome GWAS summary file",
  "  2. Place eQTLGen files for HIGH-priority batch1 genes",
  "  3. Add GTEx pancreas support files later",
  "",
  "[Input Files]",
  paste0("  expected_file    : ", expected_file),
  paste0("  manual_plan_file : ", manual_plan_file,
         ifelse(plan_exists, " (found)", " (NOT found)"))
)

writeLines(notes, file.path(out_logs, "31_mr_attachment_strategy_notes.txt"))

# ---------- 控制台反馈 ----------
cat("31_begin_manual_mr_data_attachment_strategy.R finished successfully.\n")
cat("  Files present / total:", n_exists, "/", n_total, "\n")
cat("  Minimum MR runnable  :", minimum_runnable, "\n")
if (!minimum_runnable) {
  cat("  [ACTION REQUIRED] 请按优先顺序手动放置数据文件，详见:\n")
  cat("    ", file.path(out_registry, "31_mr_data_attachment_priority.tsv"), "\n")
}