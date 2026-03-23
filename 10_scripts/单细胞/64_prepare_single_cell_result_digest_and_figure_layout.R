# =========================
# 64_prepare_single_cell_result_digest_and_figure_layout.R
# 目的：
# 1. 汇总 63 阶段 malignant pool 的关键结果
# 2. 生成单细胞结果写作 digest
# 3. 输出主文 Figure 布局建议
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "06_scRNA"

registry_dir <- file.path(root_dir, "00_registry")
result_dir   <- file.path(root_dir, "03_results")
log_dir      <- file.path(root_dir, "06_logs")

dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

cluster_score_file <- file.path(result_dir, "63_cluster_score_summary.tsv")
sample_score_file  <- file.path(result_dir, "63_sample_score_summary.tsv")

if (!file.exists(cluster_score_file)) stop("缺少 63_cluster_score_summary.tsv")
if (!file.exists(sample_score_file)) stop("缺少 63_sample_score_summary.tsv")

cluster_score <- fread(cluster_score_file)
sample_score  <- fread(sample_score_file)

# =========================================================
# Step 1. 提取最关键摘要
# =========================================================
top_sig <- cluster_score %>%
  arrange(desc(median_signature12_module)) %>%
  slice(1)

top_mcd <- cluster_score %>%
  arrange(desc(median_mcd_focus_module)) %>%
  slice(1)

summary_digest <- data.frame(
  item = c(
    "n_clusters_in_malignant_pool",
    "top_signature_cluster",
    "top_signature_cluster_median",
    "top_mcd_cluster",
    "top_mcd_cluster_median",
    "n_samples_in_malignant_pool"
  ),
  value = c(
    nrow(cluster_score),
    top_sig$cluster_id[1],
    top_sig$median_signature12_module[1],
    top_mcd$cluster_id[1],
    top_mcd$median_mcd_focus_module[1],
    nrow(sample_score)
  ),
  stringsAsFactors = FALSE
)

fwrite(
  summary_digest,
  file.path(result_dir, "64_single_cell_result_digest.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 2. Figure layout
# =========================================================
figure_layout <- data.frame(
  panel = c("A", "B", "C", "D", "E", "F"),
  content = c(
    "UMAP of major compartments",
    "UMAP of refined subtypes / malignant candidate definition",
    "DotPlot of canonical markers",
    "FeaturePlot of signature12_module1",
    "FeaturePlot of mcd_focus_module1",
    "FeaturePlot of key genes (MET/TGFBI/IDO1/INHBA/IL15RA/GBP4/SLC25A11/SNAPC4)"
  ),
  interpretation = c(
    "Define global cellular ecosystem",
    "Show epithelial-focused malignant refinement",
    "Support lineage annotation",
    "Show bulk-derived prognostic signature localization",
    "Show MCD-related transcriptional program localization",
    "Show representative gene-level heterogeneity within malignant pool"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  figure_layout,
  file.path(registry_dir, "64_single_cell_figure_layout.tsv"),
  sep = "\t", na = "NA"
)

# =========================================================
# Step 3. 写作模板（中文）
# =========================================================
cn_template <- c(
  "单细胞分析显示，GSE155698 肿瘤组织样本可分为上皮/恶性样细胞、髓系细胞、T/NK 细胞、B/plasma 细胞、成纤维细胞以及少量腺泡样细胞等主要细胞群。",
  "在此基础上，我们进一步对上皮相关 cluster 及其相邻未解析 cluster 进行了针对性复核，并最终将 cluster 17 和 26 纳入 malignant candidate pool，而将 cluster 24 单独定义为 acinar-like，并排除 cluster 10 和 21 所代表的免疫污染群。",
  "值得注意的是，bulk 队列中构建的 12-gene prognostic signature 在 malignant candidate pool 中呈现连续而异质的空间富集，而非均匀分布。",
  "与此同时，MCD-focused module score 与 signature12 module score 在空间上表现出相近的富集趋势，提示 bulk 中识别出的风险相关转录程序主要定位于一部分 ductal-like malignant epithelial cells。",
  "在单基因层面，MET、TGFBI、GBP4、SLC25A11 和 SNAPC4 在 malignant pool 中均可检测到，而 IDO1、INHBA 与 IL15RA 则呈现更局灶的表达模式，提示不同 signature gene 在恶性细胞内部具有不同程度的程序异质性。"
)

writeLines(
  cn_template,
  file.path(result_dir, "64_single_cell_result_paragraph_cn.txt")
)

# =========================================================
# Step 4. 写作模板（英文）
# =========================================================
en_template <- c(
  "Single-cell transcriptomic analysis of the GSE155698 tumor tissue cohort resolved major cellular compartments, including epithelial/malignant-like cells, myeloid cells, T/NK cells, B/plasma cells, fibroblasts, and a minor acinar-like population.",
  "We then refined the epithelial-focused compartment by re-evaluating unresolved clusters adjacent to the epithelial region, ultimately incorporating clusters 17 and 26 into the malignant candidate pool, while retaining cluster 24 as acinar-like and excluding clusters 10 and 21 as immune contaminants.",
  "Notably, the bulk-derived 12-gene prognostic signature showed spatially coherent yet heterogeneous enrichment within the malignant candidate pool rather than uniform expression across all malignant-like cells.",
  "Moreover, the MCD-focused module score exhibited a distribution pattern broadly concordant with the 12-gene signature score, suggesting that the risk-associated transcriptional program identified in bulk tumors was primarily localized to a subset of ductal-like malignant epithelial cells.",
  "At the single-gene level, MET, TGFBI, GBP4, SLC25A11, and SNAPC4 were readily detectable in the malignant pool, whereas IDO1, INHBA, and IL15RA displayed more focal expression patterns, highlighting intratumoral heterogeneity of the signature-associated program."
)

writeLines(
  en_template,
  file.path(result_dir, "64_single_cell_result_paragraph_en.txt")
)

# =========================================================
# Step 5. 阶段摘要
# =========================================================
stage_summary <- data.frame(
  item = c(
    "stage64_completed",
    "result_digest_written",
    "figure_layout_written",
    "cn_paragraph_written",
    "en_paragraph_written",
    "next_stage"
  ),
  value = c(
    TRUE,
    file.exists(file.path(result_dir, "64_single_cell_result_digest.tsv")),
    file.exists(file.path(registry_dir, "64_single_cell_figure_layout.tsv")),
    file.exists(file.path(result_dir, "64_single_cell_result_paragraph_cn.txt")),
    file.exists(file.path(result_dir, "64_single_cell_result_paragraph_en.txt")),
    "65_prepare_single_cell_figure_legend_and_main_text_polish"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  stage_summary,
  file.path(registry_dir, "64_stage_summary.tsv"),
  sep = "\t", na = "NA"
)

writeLines(
  c(
    "64 completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Goal]",
    "Prepare single-cell result digest and figure layout for manuscript writing.",
    "",
    "[Next step]",
    "Draft figure legend and refine manuscript-ready wording."
  ),
  file.path(log_dir, "64_prepare_single_cell_result_digest_notes.txt")
)

cat("64_prepare_single_cell_result_digest_and_figure_layout.R finished successfully.\n")