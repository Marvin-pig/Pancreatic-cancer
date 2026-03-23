# =========================
# 58_preliminary_annotation_gse155698_tumor.R
# 目的：
# 1. 基于 57 的 pre-integration clustered Seurat 对象
# 2. 用 canonical markers 做 preliminary annotation
# 3. 输出 cluster-marker digest 和 cluster-celltype mapping
# 4. 修复 cluster_id / cluster 类型不一致导致的 left_join 报错
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

# =========================================================
# Step 1. 读取 57 对象
# =========================================================
obj_file <- file.path(proc_dir, "57_gse155698_tumor_preintegration_umap_clustered.rds")
if (!file.exists(obj_file)) {
  stop("未找到 57_gse155698_tumor_preintegration_umap_clustered.rds")
}
obj <- readRDS(obj_file)

# 确保 cluster identity 存在
if (length(Idents(obj)) == 0) {
  stop("Seurat 对象中未找到 cluster identities。")
}

# =========================================================
# Step 2. 定义 canonical markers
# =========================================================
marker_panel <- list(
  epithelial_malignant = c("EPCAM", "KRT19", "KRT18", "KRT8", "MUC1", "KRT17", "CEACAM6", "MSLN", "KRT7"),
  myeloid = c("LST1", "C1QA", "C1QB", "C1QC", "APOE", "FCER1A", "HLA-DRA", "TYMP", "CTSB"),
  t_nk = c("CD3D", "CD3E", "TRBC1", "NKG7", "GNLY", "GZMB", "IL7R"),
  b_plasma = c("MS4A1", "CD79A", "CD79B", "MZB1", "SDC1", "IGKC", "JCHAIN", "IGJ"),
  fibroblast_caf = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "TAGLN", "COL6A3", "SPARC", "RGS5"),
  endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "RAMP2", "PLVAP"),
  acinar_like = c("PRSS1", "PRSS3", "CPA1", "CPB1", "REG1A", "REG1B"),
  cycling = c("MKI67", "TOP2A", "UBE2C", "TYMS", "CENPF")
)

marker_table <- stack(marker_panel)
colnames(marker_table) <- c("gene", "panel")

fwrite(
  marker_table,
  file.path(registry_dir, "58_marker_panel_registry.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 3. DotPlot
# =========================================================
all_markers <- unique(marker_table$gene)
present_markers <- all_markers[all_markers %in% rownames(obj)]

if (length(present_markers) == 0) {
  warning("没有 canonical markers 出现在当前 Seurat 对象的基因列表中。")
} else {
  pdf(file.path(figure_dir, "58_dotplot_canonical_markers.pdf"), width = 14, height = 8)
  print(
    DotPlot(obj, features = rev(present_markers)) +
      RotatedAxis()
  )
  dev.off()
}

# =========================================================
# Step 4. 读取 57 的 top markers
# =========================================================
markers_file <- file.path(result_dir, "57_top10_markers_per_cluster.tsv")
if (!file.exists(markers_file)) {
  stop("未找到 57_top10_markers_per_cluster.tsv")
}
top10_markers <- fread(markers_file)

if (!"cluster" %in% colnames(top10_markers)) {
  stop("57_top10_markers_per_cluster.tsv 中未找到 'cluster' 列。")
}
if (!"gene" %in% colnames(top10_markers)) {
  stop("57_top10_markers_per_cluster.tsv 中未找到 'gene' 列。")
}

# 关键修复：统一 cluster 类型为 character
top10_markers$cluster <- as.character(top10_markers$cluster)

# =========================================================
# Step 5. 生成 cluster marker digest
# =========================================================
marker_digest <- top10_markers %>%
  group_by(cluster) %>%
  summarise(
    top_markers = paste(unique(gene)[1:min(10, length(unique(gene)))], collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(cluster = as.character(cluster))

fwrite(
  marker_digest,
  file.path(result_dir, "58_cluster_marker_digest.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 6. 生成 preliminary annotation 模板
# 关键修复：cluster_id 也统一转 character，再 left_join
# =========================================================
cluster_ids <- sort(unique(as.character(Idents(obj))))

anno_template <- data.frame(
  cluster_id = cluster_ids,
  preliminary_celltype = NA_character_,
  confidence = NA_character_,
  suspicious_flag = NA_character_,
  suspicious_reason = NA_character_,
  stringsAsFactors = FALSE
) %>%
  mutate(cluster_id = as.character(cluster_id)) %>%
  left_join(
    marker_digest %>% mutate(cluster = as.character(cluster)),
    by = c("cluster_id" = "cluster")
  )

fwrite(
  anno_template,
  file.path(registry_dir, "58_preliminary_cluster_annotation.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 7. FeaturePlot（按大类 marker 分组输出）
# =========================================================
plot_gene_group <- function(genes, out_pdf) {
  genes_use <- genes[genes %in% rownames(obj)]
  if (length(genes_use) == 0) return(NULL)
  
  pdf(out_pdf, width = 12, height = 8)
  print(FeaturePlot(obj, features = genes_use, ncol = 3))
  dev.off()
}

plot_gene_group(
  marker_panel$epithelial_malignant,
  file.path(figure_dir, "58_featureplot_epithelial_malignant.pdf")
)

plot_gene_group(
  marker_panel$myeloid,
  file.path(figure_dir, "58_featureplot_myeloid.pdf")
)

plot_gene_group(
  marker_panel$t_nk,
  file.path(figure_dir, "58_featureplot_t_nk.pdf")
)

plot_gene_group(
  marker_panel$b_plasma,
  file.path(figure_dir, "58_featureplot_b_plasma.pdf")
)

plot_gene_group(
  marker_panel$fibroblast_caf,
  file.path(figure_dir, "58_featureplot_fibroblast_caf.pdf")
)

plot_gene_group(
  marker_panel$endothelial,
  file.path(figure_dir, "58_featureplot_endothelial.pdf")
)

plot_gene_group(
  marker_panel$acinar_like,
  file.path(figure_dir, "58_featureplot_acinar_like.pdf")
)

plot_gene_group(
  marker_panel$cycling,
  file.path(figure_dir, "58_featureplot_cycling.pdf")
)

# =========================================================
# Step 8. 阶段摘要
# =========================================================
stage_summary <- data.frame(
  item = c(
    "stage58_completed",
    "marker_panel_written",
    "dotplot_written",
    "cluster_marker_digest_written",
    "preliminary_annotation_template_written",
    "next_stage"
  ),
  value = c(
    TRUE,
    file.exists(file.path(registry_dir, "58_marker_panel_registry.tsv")),
    file.exists(file.path(figure_dir, "58_dotplot_canonical_markers.pdf")),
    file.exists(file.path(result_dir, "58_cluster_marker_digest.tsv")),
    file.exists(file.path(registry_dir, "58_preliminary_cluster_annotation.tsv")),
    "59_finalize_preliminary_annotation_and_define_major_compartments"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  stage_summary,
  file.path(registry_dir, "58_stage_summary.tsv"),
  sep = "\t",
  na = "NA"
)

# =========================================================
# Step 9. 日志
# =========================================================
writeLines(
  c(
    "58 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Perform preliminary annotation using canonical markers.",
    "",
    "[Key fix]",
    "Converted cluster_id and marker cluster columns to character before left_join.",
    "",
    "[Next step]",
    "Fill 58_preliminary_cluster_annotation.tsv and define major cell compartments."
  ),
  file.path(log_dir, "58_preliminary_annotation_notes.txt")
)

cat("58_preliminary_annotation_gse155698_tumor.R finished successfully.\n")
cat("Generated files:\n")
cat("- 00_registry/58_marker_panel_registry.tsv\n")
cat("- 03_results/58_cluster_marker_digest.tsv\n")
cat("- 00_registry/58_preliminary_cluster_annotation.tsv\n")
cat("- 00_registry/58_stage_summary.tsv\n")
cat("- 04_figures/58_dotplot_canonical_markers.pdf\n")
cat("- 04_figures/58_featureplot_epithelial_malignant.pdf\n")
cat("- 04_figures/58_featureplot_myeloid.pdf\n")
cat("- 04_figures/58_featureplot_t_nk.pdf\n")
cat("- 04_figures/58_featureplot_b_plasma.pdf\n")
cat("- 04_figures/58_featureplot_fibroblast_caf.pdf\n")
cat("- 04_figures/58_featureplot_endothelial.pdf\n")
cat("- 04_figures/58_featureplot_acinar_like.pdf\n")
cat("- 04_figures/58_featureplot_cycling.pdf\n")