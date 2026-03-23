# =========================================================
# 09D_mcd_maintext_mechanism_integration_fixed.R
# 作用：
# 1. 整合 09B Hallmark GSEA 与 09C mechanism candidates
# 2. 输出主文候选基因表
# 3. 绘制 Hallmark 机制条形图
# 4. 绘制 direct MCD-linked genes 箱线图
# 5. 绘制 immune/stromal bridge genes 箱线图
# 6. 生成 09D 主文整合图
# 说明：
# - 本版不再使用 cairo_pdf / X11
# - 全部改为 grDevices::pdf() 输出，规避 macOS cairo/X11 依赖问题
# =========================================================

project_root <- Sys.getenv("PROJECT_ROOT", unset = "/Users/wmz_mac/Desktop/胰腺癌")

proc_dir   <- file.path(project_root, "02_processed_data", "TCGA_PAAD")
clust_dir  <- file.path(project_root, "04_bulk_analysis", "02_clustering")
enrich_dir <- file.path(project_root, "04_bulk_analysis", "07_enrichment")
mech_dir   <- file.path(project_root, "04_bulk_analysis", "08_mechanism")
main_dir   <- file.path(project_root, "04_bulk_analysis", "09_maintext")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# 输入文件
# -------------------------------
expr_file  <- file.path(proc_dir,  "tcga_paad_expr_all_tp_symbol_counts_filtered.rds")
pheno_file <- file.path(clust_dir, "tcga_paad_mcd_subtype_phenotype.tsv")

gsea_hl_file <- file.path(enrich_dir, "tcga_paad_mcd_high_vs_low_hallmark_gsea.tsv")
gsea_hi_file <- file.path(enrich_dir, "tcga_paad_mcd_high_vs_intermediate_hallmark_gsea.tsv")

candidate_file <- file.path(mech_dir, "tcga_paad_mcd_key_mechanism_gene_candidates.tsv")
overlap_file   <- file.path(mech_dir, "tcga_paad_mcd_hallmark_leading_edge_mcd_overlap.tsv")

mcd_gene_set_file <- file.path(
  proc_dir,
  "tcga_paad_metabolic_cell_death_gene_set_analysis_ready.tsv"
)

# -------------------------------
# 参数
# -------------------------------
hallmark_focus <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_DNA_REPAIR"
)

direct_gene_priority <- c("GPX4", "ABI1", "FLNA", "SLC31A1")
bridge_gene_priority <- c("CXCL9", "CXCL10", "CXCL11", "FGL2", "IL7R", "CMKLR1")

# -------------------------------
# 载入包
# -------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)

# -------------------------------
# 工具函数
# -------------------------------
read_tsv_safe <- function(file) {
  if (!file.exists(file)) stop("找不到文件: ", file)
  fread(file, data.table = FALSE)
}

std_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  toupper(x)
}

infer_first_col <- function(df, candidates, required = TRUE, object_name = "data.frame") {
  hit <- intersect(candidates, colnames(df))[1]
  if (is.na(hit) && required) {
    stop(object_name, " 中未找到以下任一列：", paste(candidates, collapse = ", "))
  }
  hit
}

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid = element_blank(),
      legend.position = "right"
    )
}

# 关键修复：不再使用 ggsave(device = cairo_pdf)
save_pdf_plot <- function(plot_obj, out_file, width = 8, height = 6) {
  grDevices::pdf(file = out_file, width = width, height = height, onefile = FALSE)
  print(plot_obj)
  grDevices::dev.off()
}

make_boxplot_long <- function(expr_mat, pheno, gene_vec, sample_col, subtype_col) {
  gene_vec <- unique(gene_vec[gene_vec %in% rownames(expr_mat)])
  if (length(gene_vec) == 0) return(data.frame())
  
  logcpm <- edgeR::cpm(expr_mat, log = TRUE, prior.count = 1)
  
  df <- as.data.frame(
    t(logcpm[gene_vec, , drop = FALSE]),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  df$sample_id <- rownames(df)
  rownames(df) <- NULL
  
  df <- merge(
    df,
    pheno[, c(sample_col, subtype_col), drop = FALSE],
    by.x = "sample_id",
    by.y = sample_col,
    all.x = TRUE,
    sort = FALSE
  )
  colnames(df)[ncol(df)] <- "subtype_label"
  
  long_df <- tidyr::pivot_longer(
    df,
    cols = all_of(gene_vec),
    names_to = "gene",
    values_to = "logCPM"
  )
  
  long_df$subtype_label <- factor(
    long_df$subtype_label,
    levels = c("MCD-low", "MCD-intermediate", "MCD-high")
  )
  
  long_df
}

extract_hallmark_subset <- function(file, contrast_label, hallmark_focus) {
  df <- read_tsv_safe(file)
  
  id_col   <- infer_first_col(df, c("ID", "Description"), object_name = basename(file))
  nes_col  <- infer_first_col(df, c("NES"), object_name = basename(file))
  padj_col <- infer_first_col(
    df,
    c("p.adjust", "p_adj", "adj.P.Val"),
    required = FALSE,
    object_name = basename(file)
  )
  
  out <- df[df[[id_col]] %in% hallmark_focus, , drop = FALSE]
  if (nrow(out) == 0) return(data.frame())
  
  out$pathway_id <- out[[id_col]]
  out$NES2       <- as.numeric(out[[nes_col]])
  out$p_adjust2  <- if (!is.na(padj_col)) as.numeric(out[[padj_col]]) else NA_real_
  out$contrast   <- contrast_label
  
  out[, c("contrast", "pathway_id", "NES2", "p_adjust2"), drop = FALSE]
}

# -------------------------------
# 1. 读取输入数据
# -------------------------------
expr <- readRDS(expr_file)
if (!is.matrix(expr)) expr <- as.matrix(expr)

pheno <- read_tsv_safe(pheno_file)
sample_col  <- infer_first_col(pheno, c("sample_id", "sample_barcode"), object_name = basename(pheno_file))
subtype_col <- infer_first_col(pheno, c("subtype_label", "mcd_subtype", "subtype"), object_name = basename(pheno_file))

matched_idx <- match(colnames(expr), pheno[[sample_col]])
if (any(is.na(matched_idx))) {
  keep        <- !is.na(matched_idx)
  expr        <- expr[, keep, drop = FALSE]
  matched_idx <- matched_idx[keep]
}
pheno <- pheno[matched_idx, , drop = FALSE]
stopifnot(identical(colnames(expr), pheno[[sample_col]]))

pheno[[subtype_col]] <- factor(
  pheno[[subtype_col]],
  levels = c("MCD-low", "MCD-intermediate", "MCD-high")
)

candidate_df <- read_tsv_safe(candidate_file)
overlap_df   <- read_tsv_safe(overlap_file)

mcd_gs <- read_tsv_safe(mcd_gene_set_file)
mcd_symbol_col <- infer_first_col(
  mcd_gs,
  c("gene_symbol", "symbol", "SYMBOL", "GeneSymbol"),
  object_name = basename(mcd_gene_set_file)
)
mcd_gs$gene_symbol_std <- std_gene(mcd_gs[[mcd_symbol_col]])
mcd_gs <- mcd_gs[!is.na(mcd_gs$gene_symbol_std), , drop = FALSE]

# -------------------------------
# 2. 生成主文候选基因表
# -------------------------------
candidate_df$gene_symbol_std <- std_gene(candidate_df$gene_symbol_std)

if (!"gene_symbol" %in% colnames(candidate_df)) {
  candidate_df$gene_symbol <- candidate_df$gene_symbol_std
} else {
  candidate_df$gene_symbol <- ifelse(
    !is.na(candidate_df$gene_symbol) & candidate_df$gene_symbol != "",
    candidate_df$gene_symbol,
    candidate_df$gene_symbol_std
  )
}

direct_genes <- unique(std_gene(direct_gene_priority))
bridge_genes <- unique(std_gene(bridge_gene_priority))

sort_cols_present <- intersect(c("score", "n_contrast", "n_pathway"), colnames(candidate_df))

main_gene_table <- candidate_df %>%
  mutate(
    gene_symbol_std = std_gene(gene_symbol_std),
    mechanism_group = case_when(
      gene_symbol_std %in% direct_genes ~ "Direct MCD-linked",
      gene_symbol_std %in% bridge_genes ~ "Immune/stromal bridge",
      TRUE ~ "Other"
    )
  ) %>%
  filter(mechanism_group != "Other") %>%
  {
    tbl <- .
    tbl <- tbl[order(
      match(tbl$mechanism_group, c("Direct MCD-linked", "Immune/stromal bridge")),
      if ("score" %in% sort_cols_present) -tbl$score else rep(0, nrow(tbl)),
      if ("n_contrast" %in% sort_cols_present) -tbl$n_contrast else rep(0, nrow(tbl)),
      if ("n_pathway" %in% sort_cols_present) -tbl$n_pathway else rep(0, nrow(tbl))
    ), ]
    tbl
  }

if (length(sort_cols_present) < 3) {
  missing_cols <- setdiff(c("score", "n_contrast", "n_pathway"), sort_cols_present)
  message("注意：candidate_df 缺少排序列（", paste(missing_cols, collapse = ", "),
          "），已按现有列排序。")
}

fwrite(
  main_gene_table,
  file.path(main_dir, "tcga_paad_mcd_maintext_mechanism_gene_table.tsv"),
  sep = "\t",
  quote = FALSE
)

# -------------------------------
# 3. Hallmark 主文机制条形图
# -------------------------------
hl1 <- extract_hallmark_subset(
  gsea_hl_file,
  contrast_label = "MCD-high vs MCD-low",
  hallmark_focus = hallmark_focus
)

hl2 <- extract_hallmark_subset(
  gsea_hi_file,
  contrast_label = "MCD-high vs MCD-intermediate",
  hallmark_focus = hallmark_focus
)

hallmark_plot_df <- rbind(hl1, hl2)

if (nrow(hallmark_plot_df) == 0) {
  warning("hallmark_plot_df 为空，跳过 Hallmark 条形图。")
  p_hallmark <- ggplot() +
    labs(title = "No Hallmark data available") +
    theme_pub()
} else {
  hallmark_plot_df$pathway_label <- gsub("^HALLMARK_", "", hallmark_plot_df$pathway_id)
  hallmark_plot_df$pathway_label <- gsub("_", " ", hallmark_plot_df$pathway_label)
  
  hallmark_order <- c(
    "INFLAMMATORY RESPONSE",
    "INTERFERON GAMMA RESPONSE",
    "EPITHELIAL MESENCHYMAL TRANSITION",
    "OXIDATIVE PHOSPHORYLATION",
    "DNA REPAIR"
  )
  
  hallmark_plot_df$pathway_label <- factor(
    hallmark_plot_df$pathway_label,
    levels = rev(hallmark_order)
  )
  
  p_hallmark <- ggplot(hallmark_plot_df, aes(x = NES2, y = pathway_label, fill = contrast)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    geom_vline(xintercept = 0, linetype = 2) +
    labs(
      title = "Key Hallmark programs across major MCD contrasts",
      x = "Normalized enrichment score (NES)",
      y = NULL
    ) +
    theme_pub()
}

save_pdf_plot(
  p_hallmark,
  file.path(main_dir, "tcga_paad_mcd_maintext_hallmark_barplot.pdf"),
  width = 8.5,
  height = 4.8
)

# -------------------------------
# 4. Direct MCD-linked genes 箱线图
# -------------------------------
direct_present <- direct_gene_priority[direct_gene_priority %in% rownames(expr)]

if (length(direct_present) == 0) {
  message("注意：direct_gene_priority 中的基因均不在表达矩阵中，跳过 direct 箱线图。")
  p_direct <- ggplot() +
    labs(title = "No direct MCD-linked genes found in expression matrix") +
    theme_pub()
} else {
  direct_long <- make_boxplot_long(
    expr_mat = expr,
    pheno = pheno,
    gene_vec = direct_present,
    sample_col = sample_col,
    subtype_col = subtype_col
  )
  
  if (nrow(direct_long) == 0) {
    message("注意：direct_long 数据框为空，跳过 direct 箱线图。")
    p_direct <- ggplot() +
      labs(title = "No data for direct MCD-linked gene boxplot") +
      theme_pub()
  } else {
    p_direct <- ggplot(direct_long, aes(x = subtype_label, y = logCPM, fill = subtype_label)) +
      geom_boxplot(outlier.size = 0.4, width = 0.62, alpha = 0.85) +
      geom_jitter(width = 0.15, size = 0.35, alpha = 0.25, color = "grey30") +
      facet_wrap(~ gene, scales = "free_y", ncol = 2) +
      labs(
        title = "Direct MCD-linked genes across MCD subtypes",
        x = "MCD subtype",
        y = "logCPM"
      ) +
      theme_pub() +
      theme(legend.position = "none")
  }
}

save_pdf_plot(
  p_direct,
  file.path(main_dir, "tcga_paad_mcd_direct_linked_gene_boxplot.pdf"),
  width = 8.4,
  height = 6.2
)

# -------------------------------
# 5. Immune/stromal bridge genes 箱线图
# -------------------------------
bridge_present <- bridge_gene_priority[bridge_gene_priority %in% rownames(expr)]

if (length(bridge_present) == 0) {
  message("注意：bridge_gene_priority 中的基因均不在表达矩阵中，跳过 bridge 箱线图。")
  p_bridge <- ggplot() +
    labs(title = "No bridge genes found in expression matrix") +
    theme_pub()
} else {
  bridge_long <- make_boxplot_long(
    expr_mat = expr,
    pheno = pheno,
    gene_vec = bridge_present,
    sample_col = sample_col,
    subtype_col = subtype_col
  )
  
  if (nrow(bridge_long) == 0) {
    message("注意：bridge_long 数据框为空，跳过 bridge 箱线图。")
    p_bridge <- ggplot() +
      labs(title = "No data for bridge gene boxplot") +
      theme_pub()
  } else {
    p_bridge <- ggplot(bridge_long, aes(x = subtype_label, y = logCPM, fill = subtype_label)) +
      geom_boxplot(outlier.size = 0.4, width = 0.62, alpha = 0.85) +
      geom_jitter(width = 0.15, size = 0.35, alpha = 0.25, color = "grey30") +
      facet_wrap(~ gene, scales = "free_y", ncol = 3) +
      labs(
        title = "Immune/stromal bridge genes across MCD subtypes",
        x = "MCD subtype",
        y = "logCPM"
      ) +
      theme_pub() +
      theme(legend.position = "none")
  }
}

save_pdf_plot(
  p_bridge,
  file.path(main_dir, "tcga_paad_mcd_bridge_gene_boxplot.pdf"),
  width = 10.5,
  height = 6.2
)

# -------------------------------
# 6. 结果摘要表
# -------------------------------
result_summary <- data.frame(
  section = c(
    "Direct MCD-linked genes",
    "Immune/stromal bridge genes",
    "Key positive Hallmark programs"
  ),
  content = c(
    paste(unique(direct_present), collapse = "; "),
    paste(unique(bridge_present), collapse = "; "),
    paste(c("INFLAMMATORY RESPONSE", "INTERFERON GAMMA RESPONSE", "EPITHELIAL MESENCHYMAL TRANSITION"), collapse = "; ")
  ),
  stringsAsFactors = FALSE
)

fwrite(
  result_summary,
  file.path(main_dir, "tcga_paad_mcd_maintext_result_summary.tsv"),
  sep = "\t",
  quote = FALSE
)

# -------------------------------
# 7. 综合主图
# -------------------------------
if (has_patchwork) {
  library(patchwork)
  
  p_main <- (p_hallmark / p_direct) | p_bridge
  
  grDevices::pdf(
    file = file.path(main_dir, "tcga_paad_mcd_maintext_mechanism_integrated_figure.pdf"),
    width = 16,
    height = 10,
    onefile = FALSE
  )
  print(p_main)
  grDevices::dev.off()
} else {
  message("patchwork 未安装，已分别输出 Hallmark / direct genes / bridge genes 三张分图。")
}

# -------------------------------
# 8. 输出提示
# -------------------------------
message("09D main-text mechanism integration finished.")
message("输出文件：")
message("1) tcga_paad_mcd_maintext_mechanism_gene_table.tsv")
message("2) tcga_paad_mcd_maintext_hallmark_barplot.pdf")
message("3) tcga_paad_mcd_direct_linked_gene_boxplot.pdf")
message("4) tcga_paad_mcd_bridge_gene_boxplot.pdf")
message("5) tcga_paad_mcd_maintext_result_summary.tsv")
message("6) tcga_paad_mcd_maintext_mechanism_integrated_figure.pdf (if patchwork installed)")