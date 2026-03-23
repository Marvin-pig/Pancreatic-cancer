# =========================
# 54g_finalize_stage54_outputs.R
# 目的：
# 1. 汇总 54b / 54c / 54f 的导入结果
# 2. 生成 54 阶段最终冻结版输出文件
# 3. 为 55 阶段 QC 分析提供统一输入
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "06_scRNA"

registry_dir <- file.path(root_dir, "00_registry")
proc_dir     <- file.path(root_dir, "02_processed_data")
table_dir    <- file.path(root_dir, "05_tables")
log_dir      <- file.path(root_dir, "06_logs")

dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------
# Step 1. 定义候选输入文件
# ---------------------------------------------------------
f_54b_plan   <- file.path(registry_dir, "54b_gse155698_import_plan.tsv")
f_54c_plan   <- file.path(registry_dir, "54c_gse155698_final_import_plan.tsv")
f_54f_h5plan <- file.path(registry_dir, "54f_gsm4710703_h5_plan.tsv")
f_54f_merge  <- file.path(registry_dir, "54f_gsm4710703_merge_status.tsv")

f_54b_qc <- file.path(table_dir, "54b_gse155698_pre_qc_summary.tsv")
f_54c_qc <- file.path(table_dir, "54c_gse155698_pre_qc_summary.tsv")
f_54f_qc <- file.path(table_dir, "54f_gsm4710703_pre_qc_summary.tsv")

rds_54b <- file.path(proc_dir, "54b_gse155698_tumor_tissue_raw_seurat.rds")
rds_54c <- file.path(proc_dir, "54c_gse155698_tumor_tissue_raw_seurat.rds")
rds_54f <- file.path(proc_dir, "54f_gse155698_tumor_tissue_raw_seurat_16samples.rds")

# ---------------------------------------------------------
# Step 2. 读取 merge 状态，确认 sample14 是否已救回
# ---------------------------------------------------------
sample14_rescued <- FALSE
merge_status_value <- NA_character_

if (file.exists(f_54f_merge)) {
  merge_tab <- fread(f_54f_merge)
  if (all(c("item", "value") %in% names(merge_tab))) {
    merge_status_value <- merge_tab$value[merge_tab$item == "merge_status"][1]
    sample14_rescued <- isTRUE(merge_status_value == "merged_success") ||
      isTRUE(merge_status_value == "sample_already_present")
  }
}

# ---------------------------------------------------------
# Step 3. 选择最终 Seurat 对象
# 优先级：54f merged > 54c > 54b
# ---------------------------------------------------------
final_rds_source <- NA_character_

if (file.exists(rds_54f)) {
  final_rds_source <- rds_54f
} else if (file.exists(rds_54c)) {
  final_rds_source <- rds_54c
} else if (file.exists(rds_54b)) {
  final_rds_source <- rds_54b
} else {
  stop("未找到任何可用的 54 阶段 Seurat 对象。")
}

# 复制为统一最终文件名
final_rds_target <- file.path(proc_dir, "54_gse155698_tumor_tissue_raw_seurat_final.rds")
file.copy(final_rds_source, final_rds_target, overwrite = TRUE)

# ---------------------------------------------------------
# Step 4. 生成最终 import plan
# 逻辑：
# - 优先用 54c final import plan
# - 若 sample14 已 rescue，则将其 method 标记为 Read10X_h5_rescued
# ---------------------------------------------------------
plan_tab <- NULL

if (file.exists(f_54c_plan)) {
  plan_tab <- fread(f_54c_plan)
} else if (file.exists(f_54b_plan)) {
  plan_tab <- fread(f_54b_plan)
} else {
  stop("缺少 54b/54c 导入计划文件。")
}

# 若存在 sample14 rescue 信息，则更新导入方式
if (sample14_rescued && file.exists(f_54f_h5plan)) {
  h5_tab <- fread(f_54f_h5plan)
  sample14_name <- h5_tab$sample_name[1]
  sample14_h5   <- h5_tab$h5_file[1]
  
  if ("sample_name" %in% names(plan_tab)) {
    if (!sample14_name %in% plan_tab$sample_name) {
      # 极少数情况下 54c 没写进去，则补一行
      plan_tab <- bind_rows(
        plan_tab,
        data.frame(
          sample_name = sample14_name,
          original_import_method = "Read10X_h5",
          final_import_method = "Read10X_h5_rescued",
          final_matrix_dir = NA_character_,
          final_matrix_file = sample14_h5,
          stringsAsFactors = FALSE
        )
      )
    } else {
      # 已存在则更新
      if ("original_import_method" %in% names(plan_tab)) {
        plan_tab$original_import_method[plan_tab$sample_name == sample14_name] <- "Read10X_h5"
      }
      if ("final_import_method" %in% names(plan_tab)) {
        plan_tab$final_import_method[plan_tab$sample_name == sample14_name] <- "Read10X_h5_rescued"
      }
      if ("final_matrix_file" %in% names(plan_tab)) {
        plan_tab$final_matrix_file[plan_tab$sample_name == sample14_name] <- sample14_h5
      }
    }
  }
}

# 补充阶段字段
if (!"stage54_final_include" %in% names(plan_tab) && "sample_name" %in% names(plan_tab)) {
  plan_tab$stage54_final_include <- "YES"
}

fwrite(
  plan_tab,
  file.path(registry_dir, "54_final_import_plan.tsv"),
  sep = "\t",
  na = "NA"
)

# ---------------------------------------------------------
# Step 5. 生成最终 pre-QC summary
# 逻辑：
# - 先取 54c（或 54b）主体
# - 再用 54f 的 sample14 QC 覆盖/补入
# ---------------------------------------------------------
qc_tab <- NULL

if (file.exists(f_54c_qc)) {
  qc_tab <- fread(f_54c_qc)
} else if (file.exists(f_54b_qc)) {
  qc_tab <- fread(f_54b_qc)
} else {
  stop("缺少 54b/54c QC 摘要文件。")
}

if (file.exists(f_54f_qc)) {
  qc14 <- fread(f_54f_qc)
  
  if ("sample_name" %in% names(qc_tab) && "sample_name" %in% names(qc14)) {
    qc_tab <- qc_tab %>% filter(sample_name != qc14$sample_name[1])
    qc_tab <- bind_rows(qc_tab, qc14)
  }
}

# 排序：按 PDAC_TISSUE 编号
if ("sample_name" %in% names(qc_tab)) {
  qc_tab <- qc_tab %>%
    mutate(sample_order = suppressWarnings(as.integer(str_extract(sample_name, "(?<=PDAC_TISSUE_)\\d+$")))) %>%
    arrange(sample_order, sample_name) %>%
    select(-sample_order)
}

fwrite(
  qc_tab,
  file.path(table_dir, "54_final_pre_qc_summary.tsv"),
  sep = "\t",
  na = "NA"
)

# ---------------------------------------------------------
# Step 6. 生成最终阶段状态摘要
# ---------------------------------------------------------
n_samples_final <- if ("sample_name" %in% names(qc_tab)) length(unique(qc_tab$sample_name)) else NA_integer_

stage_summary <- data.frame(
  item = c(
    "stage54_completed",
    "main_dataset",
    "analysis_scope",
    "sample14_rescued",
    "sample14_merge_status",
    "final_seurat_object_written",
    "final_import_plan_written",
    "final_pre_qc_summary_written",
    "n_tumor_tissue_samples_final",
    "next_stage"
  ),
  value = c(
    TRUE,
    "GSE155698",
    "PDAC_TISSUE_* only",
    sample14_rescued,
    ifelse(is.na(merge_status_value), "NA", merge_status_value),
    file.exists(final_rds_target),
    file.exists(file.path(registry_dir, "54_final_import_plan.tsv")),
    file.exists(file.path(table_dir, "54_final_pre_qc_summary.tsv")),
    n_samples_final,
    "55_qc_metrics_and_filtering"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  stage_summary,
  file.path(registry_dir, "54_final_stage_summary.tsv"),
  sep = "\t",
  na = "NA"
)

# ---------------------------------------------------------
# Step 7. 记录日志
# ---------------------------------------------------------
writeLines(
  c(
    "54g completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Finalize all stage-54 outputs after rescuing sample14 if successful.",
    "",
    "[Final outputs]",
    "- 02_processed_data/54_gse155698_tumor_tissue_raw_seurat_final.rds",
    "- 00_registry/54_final_import_plan.tsv",
    "- 05_tables/54_final_pre_qc_summary.tsv",
    "- 00_registry/54_final_stage_summary.tsv",
    "",
    "[Next step]",
    "Start stage 55: QC metrics, mitochondrial percentage calculation, threshold setting, and pre-filter QC plots."
  ),
  file.path(log_dir, "54g_finalize_stage54_outputs_notes.txt")
)

cat("54g_finalize_stage54_outputs.R finished successfully.\n")
cat("Generated files:\n")
cat("- 02_processed_data/54_gse155698_tumor_tissue_raw_seurat_final.rds\n")
cat("- 00_registry/54_final_import_plan.tsv\n")
cat("- 05_tables/54_final_pre_qc_summary.tsv\n")
cat("- 00_registry/54_final_stage_summary.tsv\n")
cat("- 06_logs/54g_finalize_stage54_outputs_notes.txt\n")