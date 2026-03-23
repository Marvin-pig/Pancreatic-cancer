# =========================
# 63_visualize_MCD_and_signature_in_malignant_pool.R
# 目的：
# 1. 读取 62 阶段已打分的 malignant candidate pool
# 2. 可视化 12-gene signature 与 MCD focus module
# 3. 输出 cluster/sample 层面的 score 摘要
# 4. 生成单基因表达图，为结果段写作做准备
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

obj_file <- file.path(proc_dir, "62_malignant_candidate_pool_scored.rds")
if (!file.exists(obj_file)) {
  stop("未找到 62_malignant_candidate_pool_scored.rds")
}
obj <- readRDS(obj_file)

meta <- obj@meta.data
meta$cell_id <- rownames(meta)

if (!"cluster_id" %in% colnames(meta)) {
  meta$cluster_id <- as.character(Idents(obj))
} else {
  meta$cluster_id <- as.character(meta$cluster_id)
}

sample_col <- if ("orig.sample" %in% colnames(meta)) {
  "orig.sample"
} else if ("orig.ident" %in% colnames(meta)) {
  "orig.ident"
} else {
  NA_character_
}

if (is.na(sample_col)) {
  meta$sample_name <- "unknown_sample"
} else {
  meta$sample_name <- meta[[sample_col]]
}

# =========================================================
# Step 1. 冻结当前用于可视化的 gene panel
# =========================================================
genes_to_plot <- c(
  "MET", "TGFBI", "IDO1", "INHBA", "IL15RA",
  "GBP4", "SLC25A11", "RFC4", "SNAPC4", "LCP1"
)

genes_present <- genes_to_plot[genes_to_plot %in% rownames(obj)]

gene_plot_registry <- data.frame(
  requested_gene = genes_to_plot,
  present_in_object = genes_to_plot %in% rownames(obj),
  stringsAsFactors = FALSE
)

fwrite(
  gene_plot_registry,
  file.path(registry_dir, "63_gene_plot_registry.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 2. score UMAP
# =========================================================
if ("signature12_module1" %in% colnames(obj@meta.data)) {
  pdf(file.path(figure_dir, "63_featureplot_signature12_module.pdf"), width = 7, height = 6)
  print(FeaturePlot(obj, features = "signature12_module1"))
  dev.off()
}

if ("mcd_focus_module1" %in% colnames(obj@meta.data)) {
  pdf(file.path(figure_dir, "63_featureplot_mcd_focus_module.pdf"), width = 7, height = 6)
  print(FeaturePlot(obj, features = "mcd_focus_module1"))
  dev.off()
}

# =========================================================
# Step 3. 单基因 FeaturePlot
# =========================================================
if (length(genes_present) > 0) {
  pdf(file.path(figure_dir, "63_featureplot_key_genes.pdf"), width = 12, height = 10)
  print(FeaturePlot(obj, features = genes_present, ncol = 3))
  dev.off()
}

# =========================================================
# Step 4. cluster-level score 摘要
# =========================================================
cluster_score_summary <- meta %>%
  group_by(cluster_id) %>%
  summarise(
    n_cells = n(),
    median_signature12_module = median(signature12_module1, na.rm = TRUE),
    mean_signature12_module = mean(signature12_module1, na.rm = TRUE),
    median_mcd_focus_module = median(mcd_focus_module1, na.rm = TRUE),
    mean_mcd_focus_module = mean(mcd_focus_module1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_signature12_module))

fwrite(
  cluster_score_summary,
  file.path(result_dir, "63_cluster_score_summary.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 5. sample-level score 摘要
# =========================================================
sample_score_summary <- meta %>%
  group_by(sample_name) %>%
  summarise(
    n_cells = n(),
    median_signature12_module = median(signature12_module1, na.rm = TRUE),
    mean_signature12_module = mean(signature12_module1, na.rm = TRUE),
    median_mcd_focus_module = median(mcd_focus_module1, na.rm = TRUE),
    mean_mcd_focus_module = mean(mcd_focus_module1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_signature12_module))

fwrite(
  sample_score_summary,
  file.path(result_dir, "63_sample_score_summary.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 6. 分数 violin 图
# =========================================================
pdf(file.path(figure_dir, "63_violin_signature12_by_cluster.pdf"), width = 8, height = 6)
print(VlnPlot(obj, features = "signature12_module1", group.by = "cluster_id", pt.size = 0))
dev.off()

pdf(file.path(figure_dir, "63_violin_mcd_focus_by_cluster.pdf"), width = 8, height = 6)
print(VlnPlot(obj, features = "mcd_focus_module1", group.by = "cluster_id", pt.size = 0))
dev.off()

# =========================================================
# Step 7. 结果写作 digest
# =========================================================
top_sig_cluster <- if (nrow(cluster_score_summary) > 0) cluster_score_summary$cluster_id[1] else NA_character_
top_mcd_cluster <- cluster_score_summary %>%
  arrange(desc(median_mcd_focus_module)) %>%
  slice(1)

writing_digest <- data.frame(
  item = c(
    "malignant_pool_cell_n",
    "top_signature_cluster",
    "top_signature_median",
    "top_mcd_cluster",
    "top_mcd_median",
    "n_signature_genes_present",
    "n_visualized_key_genes_present"
  ),
  value = c(
    ncol(obj),
    as.character(top_sig_cluster),
    ifelse(nrow(cluster_score_summary) > 0, as.character(cluster_score_summary$median_signature12_module[1]), NA),
    ifelse(nrow(top_mcd_cluster) > 0, as.character(top_mcd_cluster$cluster_id[1]), NA),
    ifelse(nrow(top_mcd_cluster) > 0, as.character(top_mcd_cluster$median_mcd_focus_module[1]), NA),
    sum(obj@meta.data$signature12_module1 %>% is.na() %>% `!`()),
    length(genes_present)
  ),
  stringsAsFactors = FALSE
)

fwrite(
  writing_digest,
  file.path(result_dir, "63_malignant_score_writing_digest.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 8. 阶段摘要
# =========================================================
stage_summary <- data.frame(
  item = c(
    "stage63_completed",
    "cluster_score_summary_written",
    "sample_score_summary_written",
    "gene_plot_registry_written",
    "signature_featureplot_written",
    "mcd_featureplot_written",
    "writing_digest_written",
    "next_stage"
  ),
  value = c(
    TRUE,
    file.exists(file.path(result_dir, "63_cluster_score_summary.tsv")),
    file.exists(file.path(result_dir, "63_sample_score_summary.tsv")),
    file.exists(file.path(registry_dir, "63_gene_plot_registry.tsv")),
    file.exists(file.path(figure_dir, "63_featureplot_signature12_module.pdf")),
    file.exists(file.path(figure_dir, "63_featureplot_mcd_focus_module.pdf")),
    file.exists(file.path(result_dir, "63_malignant_score_writing_digest.tsv")),
    "64_write_single_cell_result_paragraph_and_prepare_figure_layout"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  stage_summary,
  file.path(registry_dir, "63_stage_summary.tsv"),
  sep = "\t", na = "NA"
)

writeLines(
  c(
    "63 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Visualize MCD-related signals and 12-gene signature in malignant candidate pool.",
    "",
    "[Next step]",
    "Write the single-cell result paragraph and plan figure layout."
  ),
  file.path(log_dir, "63_visualize_MCD_and_signature_notes.txt")
)

cat("63_visualize_MCD_and_signature_in_malignant_pool.R finished successfully.\n")