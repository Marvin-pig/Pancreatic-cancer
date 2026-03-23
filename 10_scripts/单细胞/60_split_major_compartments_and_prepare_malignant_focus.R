# =========================
# 60_split_major_compartments_and_prepare_malignant_focus.R
# 目的：
# 1. 基于 59 已注释对象拆分 major compartments
# 2. 单独输出 epithelial / immune / fibroblast / unresolved 子对象
# 3. 为后续 malignant-focused refinement 做准备
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(Seurat)
})

setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "06_scRNA"

registry_dir <- file.path(root_dir, "00_registry")
proc_dir     <- file.path(root_dir, "02_processed_data")
result_dir   <- file.path(root_dir, "03_results")
log_dir      <- file.path(root_dir, "06_logs")

dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

obj_file <- file.path(proc_dir, "59_gse155698_tumor_preliminary_annotated.rds")
if (!file.exists(obj_file)) {
  stop("未找到 59_gse155698_tumor_preliminary_annotated.rds")
}
obj <- readRDS(obj_file)

meta <- obj@meta.data
if (!"major_compartment" %in% colnames(meta)) {
  stop("meta.data 缺少 major_compartment")
}

# =========================================================
# Step 1. 统计 major compartments
# =========================================================
comp_summary <- meta %>%
  mutate(cell_id = rownames(meta)) %>%
  count(major_compartment, name = "n_cells") %>%
  mutate(pct = n_cells / sum(n_cells))

fwrite(
  comp_summary,
  file.path(result_dir, "60_major_compartment_split_summary.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 2. 拆分对象
# =========================================================
obj_epithelial <- subset(obj, subset = major_compartment %in% c("epithelial", "acinar_like"))
obj_immune     <- subset(obj, subset = major_compartment %in% c("myeloid", "T_NK", "B_plasma"))
obj_fibroblast <- subset(obj, subset = major_compartment %in% c("fibroblast"))
obj_unresolved <- subset(obj, subset = major_compartment %in% c("unresolved"))

saveRDS(obj_epithelial, file.path(proc_dir, "60_epithelial_acinar_focus.rds"))
saveRDS(obj_immune,     file.path(proc_dir, "60_immune_focus.rds"))
saveRDS(obj_fibroblast, file.path(proc_dir, "60_fibroblast_focus.rds"))
saveRDS(obj_unresolved, file.path(proc_dir, "60_unresolved_focus.rds"))

# =========================================================
# =========================================================
# Step 3. unresolved 细分登记
# =========================================================
unresolved_meta <- obj_unresolved@meta.data
unresolved_meta$cell_id <- rownames(unresolved_meta)

sample_col <- if ("orig.sample" %in% colnames(unresolved_meta)) {
  "orig.sample"
} else if ("orig.ident" %in% colnames(unresolved_meta)) {
  "orig.ident"
} else {
  NA_character_
}

required_cols <- c(
  "cell_id",
  "cluster_id",
  "refined_subtype",
  "major_compartment",
  "confidence",
  "suspicious_flag",
  "suspicious_reason"
)

existing_required_cols <- required_cols[required_cols %in% colnames(unresolved_meta)]

if (!is.na(sample_col)) {
  unresolved_tab <- unresolved_meta %>%
    dplyr::select(all_of(c("cell_id", sample_col, existing_required_cols[existing_required_cols != "cell_id"])))
  
  colnames(unresolved_tab)[colnames(unresolved_tab) == sample_col] <- "sample_id"
} else {
  unresolved_tab <- unresolved_meta %>%
    dplyr::select(all_of(existing_required_cols))
  
  unresolved_tab$sample_id <- NA_character_
  unresolved_tab <- unresolved_tab %>%
    dplyr::select(sample_id, dplyr::everything())
}

fwrite(
  unresolved_tab,
  file.path(result_dir, "60_unresolved_marker_recheck.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 4. malignant-focused 策略登记
# =========================================================
strategy_tab <- data.frame(
  item = c(
    "primary_downstream_focus",
    "epithelial_focus_definition",
    "immune_focus_definition",
    "fibroblast_focus_definition",
    "unresolved_handling",
    "next_stage"
  ),
  value = c(
    "malignant_refinement",
    "major_compartment in epithelial or acinar_like",
    "major_compartment in myeloid/T_NK/B_plasma",
    "major_compartment == fibroblast",
    "keep separate and re-check before absorption into lineage-specific objects",
    "61_refine_epithelial_malignant_compartment"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  strategy_tab,
  file.path(registry_dir, "60_malignant_focus_strategy.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 5. 阶段摘要
# =========================================================
stage_summary <- data.frame(
  item = c(
    "stage60_completed",
    "major_compartment_summary_written",
    "epithelial_object_written",
    "immune_object_written",
    "fibroblast_object_written",
    "unresolved_object_written",
    "next_stage"
  ),
  value = c(
    TRUE,
    file.exists(file.path(result_dir, "60_major_compartment_split_summary.tsv")),
    file.exists(file.path(proc_dir, "60_epithelial_acinar_focus.rds")),
    file.exists(file.path(proc_dir, "60_immune_focus.rds")),
    file.exists(file.path(proc_dir, "60_fibroblast_focus.rds")),
    file.exists(file.path(proc_dir, "60_unresolved_focus.rds")),
    "61_refine_epithelial_malignant_compartment"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  stage_summary,
  file.path(registry_dir, "60_stage_summary.tsv"),
  sep = "\t", na = "NA"
)

writeLines(
  c(
    "60 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Split major compartments and prepare malignant-focused downstream analysis.",
    "",
    "[Next step]",
    "Refine epithelial/malignant compartment and revisit epithelial-adjacent unresolved clusters."
  ),
  file.path(log_dir, "60_split_major_compartments_notes.txt")
)

cat("60_split_major_compartments_and_prepare_malignant_focus.R finished successfully.\n")