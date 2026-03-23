# =========================
# 61_refine_epithelial_malignant_compartment.R
# 目的：
# 1. 基于 60 拆分结果，优先细化 epithelial / malignant compartment
# 2. 重新审查 unresolved 中与 epithelial 相邻、且最常见的 cluster
# 3. 生成 malignant refinement 候选池
# 4. 强化兼容性：自动识别 sample 列、cluster 列、FC 列，并记录 FindMarkers 错误
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
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

epi_file   <- file.path(proc_dir, "60_epithelial_acinar_focus.rds")
unres_file <- file.path(proc_dir, "60_unresolved_focus.rds")
full_file  <- file.path(proc_dir, "59_gse155698_tumor_preliminary_annotated.rds")

if (!file.exists(epi_file)) stop("缺少 60_epithelial_acinar_focus.rds")
if (!file.exists(unres_file)) stop("缺少 60_unresolved_focus.rds")
if (!file.exists(full_file)) stop("缺少 59_gse155698_tumor_preliminary_annotated.rds")

obj_epi   <- readRDS(epi_file)
obj_unres <- readRDS(unres_file)
obj_full  <- readRDS(full_file)

# =========================================================
# Step 0. 自动识别列名
# =========================================================
detect_sample_col <- function(meta_df) {
  candidate_cols <- c("orig.sample", "orig.ident", "sample_name", "sample", "orig_sample")
  hit <- candidate_cols[candidate_cols %in% colnames(meta_df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

detect_cluster_col <- function(meta_df) {
  candidate_cols <- c("cluster_id", "seurat_clusters", "cluster", "Cluster")
  hit <- candidate_cols[candidate_cols %in% colnames(meta_df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

sample_col_unres <- detect_sample_col(obj_unres@meta.data)
sample_col_full  <- detect_sample_col(obj_full@meta.data)

if (is.na(sample_col_unres)) {
  stop("obj_unres@meta.data 中未找到样本列（orig.sample / orig.ident / sample_name 等）")
}
if (is.na(sample_col_full)) {
  stop("obj_full@meta.data 中未找到样本列（orig.sample / orig.ident / sample_name 等）")
}

cluster_col_unres <- detect_cluster_col(obj_unres@meta.data)
cluster_col_full  <- detect_cluster_col(obj_full@meta.data)

# 补 cluster_id
if (is.na(cluster_col_unres)) {
  obj_unres@meta.data$cluster_id <- as.character(Idents(obj_unres))
  cluster_col_unres <- "cluster_id"
}
if (is.na(cluster_col_full)) {
  obj_full@meta.data$cluster_id <- as.character(Idents(obj_full))
  cluster_col_full <- "cluster_id"
}

obj_unres@meta.data$cluster_id <- as.character(obj_unres@meta.data[[cluster_col_unres]])
obj_full@meta.data$cluster_id  <- as.character(obj_full@meta.data[[cluster_col_full]])

obj_unres@meta.data$sample_name <- obj_unres@meta.data[[sample_col_unres]]
obj_full@meta.data$sample_name  <- obj_full@meta.data[[sample_col_full]]

# =========================================================
# Step 0b. 设置 assay / JoinLayers（兼容 Seurat v5）
# =========================================================
if ("RNA" %in% Assays(obj_full)) {
  DefaultAssay(obj_full) <- "RNA"
}
if ("RNA" %in% Assays(obj_unres)) {
  DefaultAssay(obj_unres) <- "RNA"
}
if ("RNA" %in% Assays(obj_epi)) {
  DefaultAssay(obj_epi) <- "RNA"
}

# 尝试 JoinLayers；失败也不中断
obj_full <- tryCatch({
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    JoinLayers(obj_full)
  } else {
    obj_full
  }
}, error = function(e) obj_full)

# =========================================================
# Step 1. unresolved 中优先关注的 cluster
# =========================================================
priority_clusters <- c("10", "17", "21", "26")

available_clusters <- sort(unique(obj_full@meta.data$cluster_id))
priority_present <- priority_clusters[priority_clusters %in% available_clusters]
priority_missing <- setdiff(priority_clusters, priority_present)

priority_presence_tab <- data.frame(
  cluster_id = priority_clusters,
  present_in_obj_full = priority_clusters %in% available_clusters,
  stringsAsFactors = FALSE
)

fwrite(
  priority_presence_tab,
  file.path(result_dir, "61_priority_cluster_presence_check.tsv"),
  sep = "\t", na = "NA"
)

if (length(priority_present) == 0) {
  stop(
    paste0(
      "priority clusters 均不存在于对象中。当前可用 cluster_id 包括：",
      paste(head(available_clusters, 50), collapse = ", ")
    )
  )
}

meta_unres <- obj_unres@meta.data
meta_unres$cell_id <- rownames(meta_unres)

priority_tab <- meta_unres %>%
  filter(cluster_id %in% priority_present) %>%
  count(cluster_id, sample_name, name = "n_cells") %>%
  group_by(cluster_id) %>%
  mutate(
    cluster_total = sum(n_cells),
    pct_in_cluster = n_cells / cluster_total
  ) %>%
  ungroup() %>%
  arrange(cluster_id, desc(n_cells))

fwrite(
  priority_tab,
  file.path(result_dir, "61_unresolved_priority_cluster_sample_summary.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 2. 抽取 epithelial + priority unresolved
# =========================================================
meta_full <- obj_full@meta.data
meta_full$cell_id <- rownames(meta_full)

if (!"major_compartment" %in% colnames(meta_full)) {
  stop("obj_full@meta.data 中缺少 major_compartment")
}

cells_epi <- rownames(meta_full)[meta_full$major_compartment %in% c("epithelial", "acinar_like")]
cells_unres_priority <- rownames(meta_full)[
  meta_full$major_compartment %in% c("unresolved") &
    meta_full$cluster_id %in% priority_present
]

candidate_cells <- unique(c(cells_epi, cells_unres_priority))

if (length(candidate_cells) == 0) {
  stop("未构建出 malignant candidate pool，请检查 major_compartment 和 priority clusters。")
}

obj_candidate <- subset(obj_full, cells = candidate_cells)

saveRDS(
  obj_candidate,
  file.path(proc_dir, "61_epithelial_malignant_candidate_pool.rds")
)

candidate_meta <- obj_candidate@meta.data
candidate_meta$cell_id <- rownames(candidate_meta)

if (!"cluster_id" %in% colnames(candidate_meta)) {
  candidate_meta$cluster_id <- as.character(Idents(obj_candidate))
} else {
  candidate_meta$cluster_id <- as.character(candidate_meta$cluster_id)
}

candidate_summary <- candidate_meta %>%
  count(major_compartment, cluster_id, name = "n_cells") %>%
  arrange(major_compartment, cluster_id)

fwrite(
  candidate_summary,
  file.path(result_dir, "61_candidate_pool_cluster_summary.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 3. marker panel + dotplot
# =========================================================
marker_panel <- list(
  epithelial_ductal = c("EPCAM", "KRT19", "KRT18", "KRT8", "KRT7", "MUC1", "MSLN", "CEACAM6", "KRT17"),
  acinar_like = c("PRSS1", "PRSS3", "CPA1", "CPB1", "REG1A", "REG1B", "CELA2A", "CELA2B", "CTRC", "SYCN"),
  cycling = c("MKI67", "TOP2A", "UBE2C", "TYMS", "CENPF"),
  myeloid_exclusion = c("LST1", "C1QA", "C1QB", "APOE", "FCER1A", "HLA-DRA"),
  tnk_exclusion = c("CD3D", "CD3E", "NKG7", "GNLY", "GZMB"),
  fibroblast_exclusion = c("COL1A1", "COL1A2", "DCN", "LUM", "SPARC", "TAGLN")
)

marker_table <- stack(marker_panel)
colnames(marker_table) <- c("gene", "panel")

fwrite(
  marker_table,
  file.path(registry_dir, "61_malignant_refinement_marker_panel.tsv"),
  sep = "\t", na = "NA"
)

all_markers <- unique(marker_table$gene)
present_markers <- all_markers[all_markers %in% rownames(obj_candidate)]

if (length(present_markers) > 0) {
  pdf(file.path(figure_dir, "61_candidate_pool_dotplot.pdf"), width = 14, height = 8)
  print(DotPlot(obj_candidate, features = rev(present_markers)) + RotatedAxis())
  dev.off()
}

# =========================================================
# Step 4. 逐 cluster 跑 FindMarkers，并记录错误
# =========================================================
Idents(obj_full) <- obj_full@meta.data$cluster_id

marker_error_log <- list()
priority_marker_list <- list()

for (clu in priority_present) {
  if (!clu %in% levels(Idents(obj_full))) {
    marker_error_log[[length(marker_error_log) + 1]] <- data.frame(
      cluster_id = clu,
      status = "failed",
      message = "cluster not found in Idents(obj_full)",
      stringsAsFactors = FALSE
    )
    next
  }
  
  mk <- tryCatch(
    {
      res <- FindMarkers(
        obj_full,
        ident.1 = clu,
        only.pos = TRUE,
        min.pct = 0.2,
        logfc.threshold = 0.25,
        verbose = FALSE
      )
      res$gene <- rownames(res)
      res$cluster_id <- clu
      res
    },
    error = function(e) {
      marker_error_log[[length(marker_error_log) + 1]] <<- data.frame(
        cluster_id = clu,
        status = "failed",
        message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      NULL
    }
  )
  
  if (!is.null(mk) && nrow(mk) > 0) {
    priority_marker_list[[length(priority_marker_list) + 1]] <- mk
    marker_error_log[[length(marker_error_log) + 1]] <- data.frame(
      cluster_id = clu,
      status = "success",
      message = NA_character_,
      stringsAsFactors = FALSE
    )
  } else if (is.null(mk)) {
    # error already recorded
  } else {
    marker_error_log[[length(marker_error_log) + 1]] <- data.frame(
      cluster_id = clu,
      status = "empty_result",
      message = "FindMarkers returned 0 rows",
      stringsAsFactors = FALSE
    )
  }
}

marker_error_log <- bind_rows(marker_error_log)
fwrite(
  marker_error_log,
  file.path(result_dir, "61_priority_unresolved_marker_error_log.tsv"),
  sep = "\t", na = "NA"
)

priority_markers <- bind_rows(priority_marker_list)

if (nrow(priority_markers) == 0) {
  priority_markers <- data.frame(
    cluster_id = character(0),
    gene = character(0),
    stringsAsFactors = FALSE
  )
}

fwrite(
  as.data.frame(priority_markers),
  file.path(result_dir, "61_priority_unresolved_markers.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 5. 兼容不同 Seurat 版本的 FC 列名
# =========================================================
fc_candidates <- c("avg_log2FC", "avg_logFC", "avg_diff")
fc_col <- fc_candidates[fc_candidates %in% colnames(priority_markers)]

if (length(fc_col) > 0 && nrow(priority_markers) > 0) {
  fc_col <- fc_col[1]
  
  priority_marker_digest <- priority_markers %>%
    group_by(cluster_id) %>%
    arrange(desc(.data[[fc_col]]), .by_group = TRUE) %>%
    slice_head(n = 15) %>%
    summarise(
      top_markers = paste(gene, collapse = ", "),
      fc_column_used = fc_col,
      .groups = "drop"
    )
} else {
  priority_marker_digest <- data.frame(
    cluster_id = priority_present,
    top_markers = NA_character_,
    fc_column_used = NA_character_,
    stringsAsFactors = FALSE
  )
}

fwrite(
  priority_marker_digest,
  file.path(result_dir, "61_priority_unresolved_marker_digest.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 6. 决策模板
# =========================================================
decision_template <- data.frame(
  cluster_id = priority_present,
  proposed_direction = NA_character_,   # epithelial_like_candidate / acinar_like_candidate / immune_contamination / stromal_contamination / low_quality_or_doublet_like
  confidence = NA_character_,
  supporting_markers = NA_character_,
  final_action = NA_character_,         # absorb_into_epithelial / absorb_into_acinar / exclude_from_malignant_pool / keep_for_review
  stringsAsFactors = FALSE
)

if (nrow(priority_marker_digest) > 0) {
  decision_template <- decision_template %>%
    left_join(priority_marker_digest %>% select(cluster_id, top_markers), by = "cluster_id") %>%
    mutate(
      supporting_markers = ifelse(is.na(top_markers), supporting_markers, top_markers)
    ) %>%
    select(-top_markers)
}

fwrite(
  decision_template,
  file.path(registry_dir, "61_malignant_refinement_decision_template.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 7. 策略登记
# =========================================================
strategy <- data.frame(
  item = c(
    "primary_focus",
    "candidate_pool_definition",
    "priority_unresolved_clusters_requested",
    "priority_unresolved_clusters_present",
    "sample_column_unresolved",
    "sample_column_full",
    "fc_column_detected",
    "refinement_goal",
    "next_stage"
  ),
  value = c(
    "epithelial_malignant_refinement",
    "epithelial + acinar_like + unresolved priority clusters",
    paste(priority_clusters, collapse = ","),
    paste(priority_present, collapse = ","),
    sample_col_unres,
    sample_col_full,
    ifelse(length(fc_col) == 0, "none_detected", fc_col[1]),
    "Decide which unresolved clusters should be absorbed into epithelial/acinar pool versus excluded",
    "62_finalize_malignant_candidate_definition"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  strategy,
  file.path(registry_dir, "61_malignant_refinement_strategy.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 8. 阶段摘要
# =========================================================
stage_summary <- data.frame(
  item = c(
    "stage61_completed",
    "candidate_pool_written",
    "priority_presence_check_written",
    "priority_cluster_sample_summary_written",
    "candidate_cluster_summary_written",
    "marker_panel_written",
    "priority_marker_table_written",
    "priority_marker_digest_written",
    "priority_marker_error_log_written",
    "decision_template_written",
    "candidate_dotplot_written",
    "next_stage"
  ),
  value = c(
    TRUE,
    file.exists(file.path(proc_dir, "61_epithelial_malignant_candidate_pool.rds")),
    file.exists(file.path(result_dir, "61_priority_cluster_presence_check.tsv")),
    file.exists(file.path(result_dir, "61_unresolved_priority_cluster_sample_summary.tsv")),
    file.exists(file.path(result_dir, "61_candidate_pool_cluster_summary.tsv")),
    file.exists(file.path(registry_dir, "61_malignant_refinement_marker_panel.tsv")),
    file.exists(file.path(result_dir, "61_priority_unresolved_markers.tsv")),
    file.exists(file.path(result_dir, "61_priority_unresolved_marker_digest.tsv")),
    file.exists(file.path(result_dir, "61_priority_unresolved_marker_error_log.tsv")),
    file.exists(file.path(registry_dir, "61_malignant_refinement_decision_template.tsv")),
    file.exists(file.path(figure_dir, "61_candidate_pool_dotplot.pdf")),
    "62_finalize_malignant_candidate_definition"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  stage_summary,
  file.path(registry_dir, "61_stage_summary.tsv"),
  sep = "\t", na = "NA"
)

writeLines(
  c(
    "61 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Refine epithelial/malignant compartment and re-check epithelial-adjacent unresolved clusters.",
    "",
    "[Important debug files]",
    "61_priority_cluster_presence_check.tsv",
    "61_priority_unresolved_marker_error_log.tsv",
    "",
    "[Next step]",
    "Finalize malignant candidate definition and exclude obvious non-epithelial contaminants."
  ),
  file.path(log_dir, "61_refine_epithelial_malignant_notes.txt")
)

cat("61_refine_epithelial_malignant_compartment.R finished successfully.\n")