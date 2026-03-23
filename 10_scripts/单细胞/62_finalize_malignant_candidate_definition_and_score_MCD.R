# =========================
# 62_finalize_malignant_candidate_definition_and_score_MCD.R
# 目的：
# 1. 基于 61 的 refinement 决策，最终冻结 malignant candidate pool
# 2. 将 17/26 吸收进 epithelial，24 单独保留为 acinar-like，10/21 排除
# 3. 在 malignant candidate pool 上计算 MCD genes / 12-gene signature 单细胞分数
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
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

full_file <- file.path(proc_dir, "59_gse155698_tumor_preliminary_annotated.rds")
if (!file.exists(full_file)) {
  stop("缺少 59_gse155698_tumor_preliminary_annotated.rds")
}
obj <- readRDS(full_file)

meta <- obj@meta.data
meta$cell_id <- rownames(meta)

if (!"cluster_id" %in% colnames(meta)) {
  meta$cluster_id <- as.character(Idents(obj))
} else {
  meta$cluster_id <- as.character(meta$cluster_id)
}

if (!"major_compartment" %in% colnames(meta)) {
  stop("meta.data 中缺少 major_compartment")
}

# =========================================================
# Step 1. 固定 61 的 refinement 决策
# =========================================================
decision_tab <- data.frame(
  cluster_id = c("10", "17", "21", "24", "26"),
  proposed_direction = c(
    "immune_contamination",
    "epithelial_like_candidate",
    "immune_contamination",
    "acinar_like_candidate",
    "epithelial_like_candidate"
  ),
  confidence = c("high", "high", "high", "high", "high"),
  final_action = c(
    "exclude_from_malignant_pool",
    "absorb_into_epithelial",
    "exclude_from_malignant_pool",
    "absorb_into_acinar",
    "absorb_into_epithelial"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  decision_tab,
  file.path(registry_dir, "62_final_malignant_refinement_decision.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 2. 生成最终分组标签
# =========================================================
meta$malignant_pool_label <- "other"

# 原始 epithelial
meta$malignant_pool_label[
  meta$major_compartment %in% c("epithelial")
] <- "malignant_candidate"

# 吸收的 epithelial-like unresolved
meta$malignant_pool_label[
  meta$cluster_id %in% c("17", "26")
] <- "malignant_candidate"

# acinar-like 单独保留
meta$malignant_pool_label[
  meta$major_compartment %in% c("acinar_like") | meta$cluster_id %in% c("24")
] <- "acinar_like_pool"

# 明确排除
meta$malignant_pool_label[
  meta$cluster_id %in% c("10", "21")
] <- "excluded_contamination"

# 其余 unresolved 保留待后审
meta$malignant_pool_label[
  meta$major_compartment %in% c("unresolved") &
    !meta$cluster_id %in% c("10", "17", "21", "24", "26")
] <- "residual_review_pool"

rownames(meta) <- colnames(obj)
obj@meta.data <- meta

# =========================================================
# Step 3. 拆分对象
# =========================================================
obj_malignant <- subset(obj, subset = malignant_pool_label == "malignant_candidate")
obj_acinar    <- subset(obj, subset = malignant_pool_label == "acinar_like_pool")
obj_excluded  <- subset(obj, subset = malignant_pool_label == "excluded_contamination")
obj_review    <- subset(obj, subset = malignant_pool_label == "residual_review_pool")

saveRDS(obj_malignant, file.path(proc_dir, "62_malignant_candidate_pool.rds"))
saveRDS(obj_acinar,    file.path(proc_dir, "62_acinar_like_pool.rds"))
saveRDS(obj_excluded,  file.path(proc_dir, "62_excluded_contamination_pool.rds"))
saveRDS(obj_review,    file.path(proc_dir, "62_residual_review_pool.rds"))

# =========================================================
# Step 4. 冻结 12-gene signature
# 注意方向先统一按 signature genes 做 module score
# 后续如需方向化 risk score，可再单独加权
# =========================================================
signature_12 <- c(
  "LCP1", "GBP4", "MET", "SNAPC4", "TGFBI", "SLC25A11",
  "RFC4", "IDO1", "NMUR1", "NDUFA7", "INHBA", "IL15RA"
)

# 先给一个 MCD panel 占位：用当前主线里最稳定、在 scRNA 中最可能关注的 genes
mcd_focus_genes <- c(
  "GPX4", "SLC31A1", "TGFBI", "IDO1", "IL15RA", "INHBA", "MET"
)

sig_present <- signature_12[signature_12 %in% rownames(obj_malignant)]
mcd_present <- mcd_focus_genes[mcd_focus_genes %in% rownames(obj_malignant)]

score_registry <- data.frame(
  score_name = c("signature12_module", "mcd_focus_module"),
  genes_present = c(
    paste(sig_present, collapse = ","),
    paste(mcd_present, collapse = ",")
  ),
  n_genes_present = c(length(sig_present), length(mcd_present)),
  stringsAsFactors = FALSE
)

fwrite(
  score_registry,
  file.path(registry_dir, "62_score_gene_registry.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 5. 在 malignant candidate pool 上计算分数
# =========================================================
if (length(sig_present) >= 3) {
  obj_malignant <- AddModuleScore(
    obj_malignant,
    features = list(sig_present),
    name = "signature12_module"
  )
} else {
  obj_malignant$signature12_module1 <- NA_real_
}

if (length(mcd_present) >= 3) {
  obj_malignant <- AddModuleScore(
    obj_malignant,
    features = list(mcd_present),
    name = "mcd_focus_module"
  )
} else {
  obj_malignant$mcd_focus_module1 <- NA_real_
}

saveRDS(
  obj_malignant,
  file.path(proc_dir, "62_malignant_candidate_pool_scored.rds")
)

# =========================================================
# Step 6. 输出 malignant pool 摘要
# =========================================================
mal_meta <- obj_malignant@meta.data
mal_meta$cell_id <- rownames(mal_meta)

sample_col <- if ("orig.sample" %in% colnames(mal_meta)) {
  "orig.sample"
} else if ("orig.ident" %in% colnames(mal_meta)) {
  "orig.ident"
} else {
  NA_character_
}

if (!is.na(sample_col)) {
  mal_meta$sample_name <- mal_meta[[sample_col]]
} else {
  mal_meta$sample_name <- "unknown_sample"
}

if (!"cluster_id" %in% colnames(mal_meta)) {
  mal_meta$cluster_id <- as.character(Idents(obj_malignant))
}

malignant_summary <- mal_meta %>%
  group_by(cluster_id, sample_name) %>%
  summarise(
    n_cells = n(),
    median_signature12_module = median(signature12_module1, na.rm = TRUE),
    median_mcd_focus_module = median(mcd_focus_module1, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(
  malignant_summary,
  file.path(result_dir, "62_malignant_candidate_score_summary.tsv"),
  sep = "\t", na = "NA"
)

pool_summary <- meta %>%
  count(malignant_pool_label, name = "n_cells") %>%
  mutate(pct = n_cells / sum(n_cells))

fwrite(
  pool_summary,
  file.path(result_dir, "62_malignant_pool_partition_summary.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 7. 阶段摘要
# =========================================================
stage_summary <- data.frame(
  item = c(
    "stage62_completed",
    "malignant_pool_written",
    "acinar_pool_written",
    "excluded_pool_written",
    "review_pool_written",
    "score_registry_written",
    "malignant_score_summary_written",
    "next_stage"
  ),
  value = c(
    TRUE,
    file.exists(file.path(proc_dir, "62_malignant_candidate_pool_scored.rds")),
    file.exists(file.path(proc_dir, "62_acinar_like_pool.rds")),
    file.exists(file.path(proc_dir, "62_excluded_contamination_pool.rds")),
    file.exists(file.path(proc_dir, "62_residual_review_pool.rds")),
    file.exists(file.path(registry_dir, "62_score_gene_registry.tsv")),
    file.exists(file.path(result_dir, "62_malignant_candidate_score_summary.tsv")),
    "63_visualize_MCD_and_signature_in_malignant_pool"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  stage_summary,
  file.path(registry_dir, "62_stage_summary.tsv"),
  sep = "\t", na = "NA"
)

writeLines(
  c(
    "62 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Finalize malignant candidate pool and compute MCD / 12-gene signature module scores.",
    "",
    "[Next step]",
    "Visualize MCD focus genes and signature score distribution in malignant candidate pool."
  ),
  file.path(log_dir, "62_finalize_malignant_candidate_definition_notes.txt")
)

cat("62_finalize_malignant_candidate_definition_and_score_MCD.R finished successfully.\n")