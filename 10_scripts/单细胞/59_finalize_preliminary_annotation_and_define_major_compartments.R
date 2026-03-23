# =========================
# 59_finalize_preliminary_annotation_and_define_major_compartments.R
# 目的：
# 1. 将 58 的 preliminary annotation 模板正式填充
# 2. 定义 cluster -> refined subtype -> major compartment 映射
# 3. 输出 major compartment UMAP，为后续 malignant/immune/stromal 拆分做准备
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(Seurat)
  library(ggplot2)
})

setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "06_scRNA"

registry_dir <- file.path(root_dir, "00_registry")
proc_dir     <- file.path(root_dir, "02_processed_data")
result_dir   <- file.path(root_dir, "03_results")
figure_dir   <- file.path(root_dir, "04_figures")
log_dir      <- file.path(root_dir, "06_logs")

dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Step 1. 读取 57 clustered 对象
# =========================================================
obj_file <- file.path(proc_dir, "57_gse155698_tumor_preintegration_umap_clustered.rds")
if (!file.exists(obj_file)) {
  stop("未找到 57_gse155698_tumor_preintegration_umap_clustered.rds")
}
obj <- readRDS(obj_file)

cluster_ids <- sort(unique(as.character(Idents(obj))))

# =========================================================
# Step 2. 手动填充当前第一版 annotation
# 可根据你后续复核继续改
# =========================================================
anno_map <- data.frame(
  cluster_id = cluster_ids,
  refined_subtype = NA_character_,
  major_compartment = NA_character_,
  confidence = "low",
  suspicious_flag = "NO",
  suspicious_reason = NA_character_,
  stringsAsFactors = FALSE
)

fill_cluster <- function(df, cid, subtype, comp, conf = "high", suspicious = "NO", reason = NA_character_) {
  df[df$cluster_id == cid, "refined_subtype"] <- subtype
  df[df$cluster_id == cid, "major_compartment"] <- comp
  df[df$cluster_id == cid, "confidence"] <- conf
  df[df$cluster_id == cid, "suspicious_flag"] <- suspicious
  df[df$cluster_id == cid, "suspicious_reason"] <- reason
  df
}

anno_map <- fill_cluster(anno_map, "0",  "CD4_Treg_like_T",          "T_NK",       "high")
anno_map <- fill_cluster(anno_map, "1",  "CD8_effector_like_T",      "T_NK",       "high")
anno_map <- fill_cluster(anno_map, "2",  "FOLR2_C1QC_macrophage",    "myeloid",    "high")
anno_map <- fill_cluster(anno_map, "3",  "neutrophil_like",          "myeloid",    "medium")
anno_map <- fill_cluster(anno_map, "4",  "ductal_like_epithelial",   "epithelial", "high")
anno_map <- fill_cluster(anno_map, "6",  "inflammatory_neutrophil",  "myeloid",    "medium")
anno_map <- fill_cluster(anno_map, "7",  "SPP1_macrophage",          "myeloid",    "high")
anno_map <- fill_cluster(anno_map, "9",  "mast_cell",                "myeloid",    "high")
anno_map <- fill_cluster(anno_map, "11", "FCN1_monocyte",            "myeloid",    "high")
anno_map <- fill_cluster(anno_map, "12", "NK_cell",                  "T_NK",       "high")
anno_map <- fill_cluster(anno_map, "14", "plasma_cell",              "B_plasma",   "high")
anno_map <- fill_cluster(anno_map, "15", "B_cell",                   "B_plasma",   "high")
anno_map <- fill_cluster(anno_map, "18", "APOC1_TREM2_macrophage",   "myeloid",    "high")
anno_map <- fill_cluster(anno_map, "23", "CAF_matrix",               "fibroblast", "high")
anno_map <- fill_cluster(anno_map, "24", "acinar_like",              "acinar_like","high")
anno_map <- fill_cluster(anno_map, "25", "cytotoxic_T_NK",           "T_NK",       "high")
anno_map <- fill_cluster(anno_map, "28", "plasmablast_like",         "B_plasma",   "high")

# 未明确 cluster 统一标记
anno_map$refined_subtype[is.na(anno_map$refined_subtype)] <- "needs_marker_review"
anno_map$major_compartment[is.na(anno_map$major_compartment)] <- "unresolved"
anno_map$confidence[anno_map$refined_subtype == "needs_marker_review"] <- "low"
anno_map$suspicious_flag[anno_map$refined_subtype == "needs_marker_review"] <- "YES"
anno_map$suspicious_reason[anno_map$refined_subtype == "needs_marker_review"] <- "marker_not_clear_in_current_digest"

fwrite(
  anno_map,
  file.path(registry_dir, "59_final_preliminary_annotation.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 3. 写入对象 metadata
# =========================================================
meta <- obj@meta.data
meta$cluster_id <- as.character(Idents(obj))

meta <- meta %>%
  left_join(anno_map, by = "cluster_id")

rownames(meta) <- colnames(obj)
obj@meta.data <- meta

# =========================================================
# Step 4. 输出 cluster -> compartment 摘要
# =========================================================
cluster_compartment_summary <- meta %>%
  count(cluster_id, refined_subtype, major_compartment, confidence, suspicious_flag, name = "n_cells") %>%
  arrange(major_compartment, cluster_id)

fwrite(
  cluster_compartment_summary,
  file.path(result_dir, "59_cluster_compartment_summary.tsv"),
  sep = "\t", na = "NA"
)

major_compartment_summary <- meta %>%
  count(major_compartment, name = "n_cells") %>%
  mutate(pct = n_cells / sum(n_cells))

fwrite(
  major_compartment_summary,
  file.path(result_dir, "59_major_compartment_summary.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 5. UMAP by major compartment
# =========================================================
pdf(file.path(figure_dir, "59_umap_by_major_compartment.pdf"), width = 8, height = 6)
print(DimPlot(obj, reduction = "umap", group.by = "major_compartment", label = TRUE))
dev.off()

pdf(file.path(figure_dir, "59_umap_by_refined_subtype.pdf"), width = 10, height = 7)
print(DimPlot(obj, reduction = "umap", group.by = "refined_subtype", label = TRUE))
dev.off()

# =========================================================
# Step 6. 保存带注释对象
# =========================================================
saveRDS(
  obj,
  file.path(proc_dir, "59_gse155698_tumor_preliminary_annotated.rds")
)

# =========================================================
# Step 7. 阶段摘要
# =========================================================
stage_summary <- data.frame(
  item = c(
    "stage59_completed",
    "annotation_written",
    "cluster_compartment_summary_written",
    "major_compartment_summary_written",
    "annotated_object_written",
    "next_stage"
  ),
  value = c(
    TRUE,
    file.exists(file.path(registry_dir, "59_final_preliminary_annotation.tsv")),
    file.exists(file.path(result_dir, "59_cluster_compartment_summary.tsv")),
    file.exists(file.path(result_dir, "59_major_compartment_summary.tsv")),
    file.exists(file.path(proc_dir, "59_gse155698_tumor_preliminary_annotated.rds")),
    "60_split_major_compartments_and_start_malignant_focus"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  stage_summary,
  file.path(registry_dir, "59_stage_summary.tsv"),
  sep = "\t", na = "NA"
)

writeLines(
  c(
    "59 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Finalize preliminary annotation and define major compartments.",
    "",
    "[Next step]",
    "Split epithelial / myeloid / T_NK / B_plasma / fibroblast compartments and begin malignant-focused analysis."
  ),
  file.path(log_dir, "59_finalize_preliminary_annotation_notes.txt")
)

cat("59_finalize_preliminary_annotation_and_define_major_compartments.R finished successfully.\n")