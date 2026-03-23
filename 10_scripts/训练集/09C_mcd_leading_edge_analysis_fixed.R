# =========================================================
# 09C_mcd_leading_edge_analysis_fixed.R
# 作用：
# 1. 读取 09A DEG 结果与 09B Hallmark GSEA 结果
# 2. 解析 core_enrichment，提取 leading-edge genes
# 3. 与 MCD analysis-ready gene set 做 overlap
# 4. 生成关键机制基因候选表
# 5. 输出机制总结图
# 6. 增加显式诊断信息，避免 overlap 结果静默为空
# =========================================================

project_root <- Sys.getenv("PROJECT_ROOT", unset = "/Users/wmz_mac/Desktop/胰腺癌")

deg_dir      <- file.path(project_root, "04_bulk_analysis", "06_deg")
enrich_dir   <- file.path(project_root, "04_bulk_analysis", "07_enrichment")
out_dir      <- file.path(project_root, "04_bulk_analysis", "08_mechanism")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# 输入文件
# -------------------------------
deg_files <- list(
  high_vs_low          = file.path(deg_dir, "tcga_paad_mcd_high_vs_low_deg.tsv"),
  high_vs_intermediate = file.path(deg_dir, "tcga_paad_mcd_high_vs_intermediate_deg.tsv")
)

gsea_files <- list(
  high_vs_low          = file.path(enrich_dir, "tcga_paad_mcd_high_vs_low_hallmark_gsea.tsv"),
  high_vs_intermediate = file.path(enrich_dir, "tcga_paad_mcd_high_vs_intermediate_hallmark_gsea.tsv")
)

# 这里已经修正为你确认过的真实路径
mcd_gene_set_file <- file.path(
  project_root,
  "02_processed_data", "TCGA_PAAD",
  "tcga_paad_metabolic_cell_death_gene_set_analysis_ready.tsv"
)

# -------------------------------
# 参数
# -------------------------------
priority_hallmark <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_DNA_REPAIR"
)

fallback_top_positive <- 3
fallback_top_negative <- 2
candidate_top_n       <- 20

# -------------------------------
# 载入包
# -------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

has_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)

# -------------------------------
# 工具函数
# -------------------------------
read_tsv_safe <- function(file) {
  if (!file.exists(file)) stop("找不到文件: ", file)
  fread(file, data.table = FALSE)
}

infer_first_col <- function(df, candidates, required = TRUE, object_name = "data.frame") {
  hit <- intersect(candidates, colnames(df))[1]
  if (is.na(hit) && required) {
    stop(
      object_name, " 中未找到以下任一列名：",
      paste(candidates, collapse = ", ")
    )
  }
  hit
}

std_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  toupper(x)
}

split_core_enrichment <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character(0))
  parts <- unlist(strsplit(x, "[/;,|]"))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  unique(parts)
}

collapse_unique <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = "; ")
}

sign_preserving_max_abs <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  x[which.max(abs(x))[1]]
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

pick_hallmark_pathways <- function(gsea_df,
                                   priority_terms = priority_hallmark,
                                   top_pos = fallback_top_positive,
                                   top_neg = fallback_top_negative) {
  if (nrow(gsea_df) == 0) return(character(0))
  
  id_col   <- infer_first_col(gsea_df, c("ID", "Description"), object_name = "gsea_df")
  nes_col  <- infer_first_col(gsea_df, c("NES"), object_name = "gsea_df")
  padj_col <- infer_first_col(
    gsea_df,
    c("p.adjust", "p_adj", "adj.P.Val"),
    required = FALSE,
    object_name = "gsea_df"
  )
  
  gsea_df <- gsea_df[!is.na(gsea_df[[id_col]]) & !is.na(gsea_df[[nes_col]]), , drop = FALSE]
  
  selected <- intersect(priority_terms, gsea_df[[id_col]])
  
  tmp <- gsea_df
  if (!is.na(padj_col)) {
    sig_df <- tmp[!is.na(tmp[[padj_col]]) & tmp[[padj_col]] < 0.25, , drop = FALSE]
    if (nrow(sig_df) > 0) tmp <- sig_df
  }
  
  pos_terms <- character(0)
  neg_terms <- character(0)
  
  if (!is.na(padj_col)) {
    tmp_pos <- tmp[tmp[[nes_col]] > 0, , drop = FALSE]
    tmp_neg <- tmp[tmp[[nes_col]] < 0, , drop = FALSE]
    
    if (nrow(tmp_pos) > 0) {
      tmp_pos   <- tmp_pos[order(tmp_pos[[padj_col]], -abs(tmp_pos[[nes_col]])), , drop = FALSE]
      pos_terms <- tmp_pos[[id_col]]
    }
    if (nrow(tmp_neg) > 0) {
      tmp_neg   <- tmp_neg[order(tmp_neg[[padj_col]], -abs(tmp_neg[[nes_col]])), , drop = FALSE]
      neg_terms <- tmp_neg[[id_col]]
    }
  } else {
    tmp_pos <- tmp[tmp[[nes_col]] > 0, , drop = FALSE]
    tmp_neg <- tmp[tmp[[nes_col]] < 0, , drop = FALSE]
    
    if (nrow(tmp_pos) > 0) {
      tmp_pos   <- tmp_pos[order(-abs(tmp_pos[[nes_col]])), , drop = FALSE]
      pos_terms <- tmp_pos[[id_col]]
    }
    if (nrow(tmp_neg) > 0) {
      tmp_neg   <- tmp_neg[order(-abs(tmp_neg[[nes_col]])), , drop = FALSE]
      neg_terms <- tmp_neg[[id_col]]
    }
  }
  
  auto_terms <- c(head(pos_terms, top_pos), head(neg_terms, top_neg))
  auto_terms <- unique(auto_terms[!is.na(auto_terms)])
  
  selected <- unique(c(selected, auto_terms))
  selected
}

read_deg_with_std <- function(file) {
  df <- read_tsv_safe(file)
  
  gene_col <- infer_first_col(
    df,
    c("gene_symbol", "symbol", "SYMBOL", "GeneSymbol"),
    object_name = basename(file)
  )
  
  required_cols <- c("logFC", "adj.P.Val")
  miss <- setdiff(required_cols, colnames(df))
  if (length(miss) > 0) {
    stop(basename(file), " 缺少必要列: ", paste(miss, collapse = ", "))
  }
  
  if (!"t" %in% colnames(df)) {
    if ("P.Value" %in% colnames(df)) {
      df$t <- sign(df$logFC) * -log10(pmax(df$P.Value, 1e-300))
    } else {
      df$t <- df$logFC
    }
  }
  
  df$gene_symbol_raw <- as.character(df[[gene_col]])
  df$gene_symbol_std <- std_gene(df[[gene_col]])
  
  df <- df[!is.na(df$gene_symbol_std), , drop = FALSE]
  df <- df[!duplicated(df$gene_symbol_std), , drop = FALSE]
  rownames(df) <- NULL
  df
}

read_mcd_gene_set <- function(file) {
  if (!file.exists(file)) {
    warning("未找到 MCD analysis-ready gene set 文件：", file)
    return(NULL)
  }
  
  gs <- read_tsv_safe(file)
  
  symbol_col <- infer_first_col(
    gs,
    c("gene_symbol", "symbol", "SYMBOL", "GeneSymbol"),
    object_name = basename(file)
  )
  
  type_col <- infer_first_col(
    gs,
    c("death_type", "cell_death_type", "category", "process", "mcd_class"),
    required = FALSE,
    object_name = basename(file)
  )
  
  tier_col <- infer_first_col(
    gs,
    c("tier", "set", "main_or_extended", "gene_set_tier", "evidence_level"),
    required = FALSE,
    object_name = basename(file)
  )
  
  out <- data.frame(
    gene_symbol_raw = as.character(gs[[symbol_col]]),
    gene_symbol_std = std_gene(gs[[symbol_col]]),
    stringsAsFactors = FALSE
  )
  
  if (!is.na(type_col)) out$death_type <- as.character(gs[[type_col]])
  if (!is.na(tier_col)) out$tier <- as.character(gs[[tier_col]])
  
  out <- out[!is.na(out$gene_symbol_std), , drop = FALSE]
  out <- out[!duplicated(out$gene_symbol_std), , drop = FALSE]
  rownames(out) <- NULL
  out
}

extract_leading_edge_long <- function(gsea_file, deg_df, contrast_name, selected_terms) {
  gsea_df  <- read_tsv_safe(gsea_file)
  id_col   <- infer_first_col(gsea_df, c("ID", "Description"), object_name = basename(gsea_file))
  desc_col <- infer_first_col(gsea_df, c("Description", "ID"), object_name = basename(gsea_file))
  nes_col  <- infer_first_col(gsea_df, c("NES"), object_name = basename(gsea_file))
  padj_col <- infer_first_col(
    gsea_df,
    c("p.adjust", "p_adj", "adj.P.Val"),
    required = FALSE,
    object_name = basename(gsea_file)
  )
  core_col <- infer_first_col(gsea_df, c("core_enrichment"), object_name = basename(gsea_file))
  
  keep <- gsea_df[[id_col]] %in% selected_terms
  sub  <- gsea_df[keep, , drop = FALSE]
  
  if (nrow(sub) == 0) {
    warning("contrast ", contrast_name, " 未选中任何 Hallmark 通路。")
    return(data.frame())
  }
  
  out_list <- lapply(seq_len(nrow(sub)), function(i) {
    genes <- split_core_enrichment(sub[[core_col]][i])
    if (length(genes) == 0) return(NULL)
    
    data.frame(
      contrast        = contrast_name,
      pathway_id      = as.character(sub[[id_col]][i]),
      pathway_label   = as.character(sub[[desc_col]][i]),
      NES             = as.numeric(sub[[nes_col]][i]),
      p_adjust        = if (!is.na(padj_col)) as.numeric(sub[[padj_col]][i]) else NA_real_,
      gene_symbol_raw = genes,
      gene_symbol_std = std_gene(genes),
      stringsAsFactors = FALSE
    )
  })
  
  out <- data.table::rbindlist(out_list, fill = TRUE)
  out <- as.data.frame(out)
  out <- out[!is.na(out$gene_symbol_std), , drop = FALSE]
  
  out <- merge(
    out,
    deg_df[, c("gene_symbol_raw", "gene_symbol_std", "logFC", "adj.P.Val", "t"), drop = FALSE],
    by = "gene_symbol_std",
    all.x = TRUE,
    sort = FALSE,
    suffixes = c("", "_deg")
  )
  
  out$gene_symbol <- ifelse(
    !is.na(out$gene_symbol_raw_deg) & out$gene_symbol_raw_deg != "",
    out$gene_symbol_raw_deg,
    out$gene_symbol_raw
  )
  
  out$direction <- ifelse(out$NES >= 0, "positive_pathway", "negative_pathway")
  out$abs_logFC <- abs(out$logFC)
  
  out <- out[order(out$pathway_id, out$p_adjust, -out$abs_logFC, out$gene_symbol), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# -------------------------------
# 1. 读取输入
# -------------------------------
deg_list <- lapply(deg_files, read_deg_with_std)
mcd_gs   <- read_mcd_gene_set(mcd_gene_set_file)

# -------------------------------
# 2. 选定每个 contrast 的重点 Hallmark 通路
# -------------------------------
selected_pathways_list <- list()

for (nm in names(gsea_files)) {
  gsea_df <- read_tsv_safe(gsea_files[[nm]])
  selected_pathways_list[[nm]] <- pick_hallmark_pathways(gsea_df)
}

selected_pathways_df <- data.frame(
  contrast   = rep(names(selected_pathways_list), lengths(selected_pathways_list)),
  pathway_id = unlist(selected_pathways_list),
  stringsAsFactors = FALSE
)

fwrite(
  selected_pathways_df,
  file.path(out_dir, "tcga_paad_mcd_hallmark_selected_pathways.tsv"),
  sep = "\t",
  quote = FALSE
)

# -------------------------------
# 3. 提取 leading-edge
# -------------------------------
leading_edge_list <- list()

for (nm in names(gsea_files)) {
  leading_edge_list[[nm]] <- extract_leading_edge_long(
    gsea_file      = gsea_files[[nm]],
    deg_df         = deg_list[[nm]],
    contrast_name  = nm,
    selected_terms = selected_pathways_list[[nm]]
  )
}

leading_high_low <- leading_edge_list[["high_vs_low"]]
leading_high_mid <- leading_edge_list[["high_vs_intermediate"]]

if (nrow(leading_high_low) > 0) {
  fwrite(
    leading_high_low,
    file.path(out_dir, "tcga_paad_mcd_high_vs_low_hallmark_leading_edge.tsv"),
    sep = "\t",
    quote = FALSE
  )
}

if (nrow(leading_high_mid) > 0) {
  fwrite(
    leading_high_mid,
    file.path(out_dir, "tcga_paad_mcd_high_vs_intermediate_hallmark_leading_edge.tsv"),
    sep = "\t",
    quote = FALSE
  )
}

leading_all <- data.table::rbindlist(leading_edge_list, fill = TRUE)
leading_all <- as.data.frame(leading_all)

if (nrow(leading_all) == 0) {
  stop("未能从 Hallmark GSEA 中提取到任何 leading-edge genes。请检查 09B 输出。")
}

# -------------------------------
# 4. overlap with MCD gene set
# -------------------------------
if (!is.null(mcd_gs)) {
  leading_all <- merge(
    leading_all,
    mcd_gs,
    by = "gene_symbol_std",
    all.x = TRUE,
    sort = FALSE,
    suffixes = c("", "_mcd")
  )
  
  leading_all$is_mcd_gene <- !is.na(leading_all$gene_symbol_raw_mcd)
} else {
  leading_all$is_mcd_gene <- FALSE
  leading_all$death_type  <- NA_character_
  leading_all$tier        <- NA_character_
}

if (!"death_type" %in% colnames(leading_all)) leading_all$death_type <- NA_character_
if (!"tier" %in% colnames(leading_all)) leading_all$tier <- NA_character_

cat("MCD gene set loaded:", !is.null(mcd_gs), "\n")
if (!is.null(mcd_gs)) {
  cat("MCD unique genes:", nrow(mcd_gs), "\n")
}
cat("Leading-edge total rows:", nrow(leading_all), "\n")
cat("Leading-edge unique genes:", length(unique(leading_all$gene_symbol_std)), "\n")

direct_overlap <- if (!is.null(mcd_gs)) intersect(mcd_gs$gene_symbol_std, unique(leading_all$gene_symbol_std)) else character(0)
cat("Direct overlap count:", length(direct_overlap), "\n")
if (length(direct_overlap) > 0) {
  cat("Direct overlap genes:", paste(sort(direct_overlap), collapse = ", "), "\n")
}

overlap_df <- leading_all[leading_all$is_mcd_gene %in% TRUE, , drop = FALSE]

if (nrow(overlap_df) > 0) {
  overlap_summary <- do.call(rbind, lapply(
    split(overlap_df, overlap_df$gene_symbol_std),
    function(df) {
      data.frame(
        gene_symbol_std = df$gene_symbol_std[1],
        n_pathway_hit   = length(unique(df$pathway_id)),
        n_contrast_hit  = length(unique(df$contrast)),
        pathways        = collapse_unique(df$pathway_id),
        contrasts       = collapse_unique(df$contrast),
        death_types     = collapse_unique(df$death_type),
        tiers           = collapse_unique(df$tier),
        gene_symbols    = collapse_unique(df$gene_symbol),
        stringsAsFactors = FALSE
      )
    }
  ))
  
  rownames(overlap_summary) <- NULL
  
  fwrite(
    overlap_summary,
    file.path(out_dir, "tcga_paad_mcd_hallmark_leading_edge_mcd_overlap.tsv"),
    sep = "\t",
    quote = FALSE
  )
  
  cat("Written non-empty overlap file with rows:", nrow(overlap_summary), "\n")
} else {
  fwrite(
    data.frame(),
    file.path(out_dir, "tcga_paad_mcd_hallmark_leading_edge_mcd_overlap.tsv"),
    sep = "\t",
    quote = FALSE
  )
  
  cat("WARNING: overlap_df is empty, wrote empty overlap file.\n")
}

# -------------------------------
# 5. candidate detail / summary
# -------------------------------
candidate_detail <- leading_all[, c(
  "gene_symbol_std", "gene_symbol", "contrast", "pathway_id", "pathway_label",
  "NES", "p_adjust", "logFC", "adj.P.Val", "t",
  "is_mcd_gene", "death_type", "tier"
), drop = FALSE]

fwrite(
  candidate_detail,
  file.path(out_dir, "tcga_paad_mcd_key_mechanism_gene_candidates_detail.tsv"),
  sep = "\t",
  quote = FALSE
)

split_list <- split(candidate_detail, candidate_detail$gene_symbol_std)

candidate_summary <- do.call(rbind, lapply(split_list, function(df) {
  data.frame(
    gene_symbol_std = df$gene_symbol_std[1],
    gene_symbol     = collapse_unique(df$gene_symbol),
    n_record        = nrow(df),
    n_contrast      = length(unique(df$contrast)),
    n_pathway       = length(unique(df$pathway_id)),
    pathway_list    = collapse_unique(df$pathway_id),
    contrast_list   = collapse_unique(df$contrast),
    max_abs_logFC   = safe_max(abs(df$logFC)),
    signed_logFC    = sign_preserving_max_abs(df$logFC),
    min_adjP        = safe_min(df$adj.P.Val),
    max_abs_t       = safe_max(abs(df$t)),
    is_mcd_gene     = any(df$is_mcd_gene %in% TRUE),
    death_type      = collapse_unique(df$death_type),
    tier            = collapse_unique(df$tier),
    stringsAsFactors = FALSE
  )
}))

rownames(candidate_summary) <- NULL

candidate_summary$score <- 0
candidate_summary$score <- candidate_summary$score +
  ifelse(candidate_summary$is_mcd_gene, 5, 0) +
  ifelse(!is.na(candidate_summary$tier) &
           grepl("main", candidate_summary$tier, ignore.case = TRUE), 2, 0) +
  2 * pmin(candidate_summary$n_contrast, 2) +
  pmin(candidate_summary$n_pathway, 5) +
  ifelse(!is.na(candidate_summary$max_abs_logFC) & candidate_summary$max_abs_logFC >= 2, 3,
         ifelse(!is.na(candidate_summary$max_abs_logFC) & candidate_summary$max_abs_logFC >= 1.5, 2,
                ifelse(!is.na(candidate_summary$max_abs_logFC) & candidate_summary$max_abs_logFC >= 1, 1, 0))) +
  ifelse(!is.na(candidate_summary$min_adjP) & candidate_summary$min_adjP < 1e-10, 3,
         ifelse(!is.na(candidate_summary$min_adjP) & candidate_summary$min_adjP < 1e-5, 2,
                ifelse(!is.na(candidate_summary$min_adjP) & candidate_summary$min_adjP < 0.05, 1, 0)))

candidate_summary <- candidate_summary[
  order(
    -candidate_summary$score,
    -candidate_summary$n_contrast,
    -candidate_summary$n_pathway,
    -candidate_summary$max_abs_logFC,
    candidate_summary$min_adjP,
    candidate_summary$gene_symbol
  ),
  ,
  drop = FALSE
]

fwrite(
  candidate_summary,
  file.path(out_dir, "tcga_paad_mcd_key_mechanism_gene_candidates.tsv"),
  sep = "\t",
  quote = FALSE
)

# -------------------------------
# 6. summary figure
# -------------------------------
top_genes <- head(candidate_summary$gene_symbol_std, candidate_top_n)
plot_df   <- candidate_detail[candidate_detail$gene_symbol_std %in% top_genes, , drop = FALSE]

plot_split <- split(plot_df, interaction(plot_df$gene_symbol_std, plot_df$contrast, drop = TRUE))

plot_summary <- do.call(rbind, lapply(plot_split, function(df) {
  idx <- which.max(abs(df$logFC))[1]
  data.frame(
    gene_symbol_std = df$gene_symbol_std[idx],
    gene_symbol     = collapse_unique(df$gene_symbol),
    contrast        = df$contrast[idx],
    logFC           = df$logFC[idx],
    pathway_id      = df$pathway_id[idx],
    is_mcd_gene     = any(df$is_mcd_gene %in% TRUE),
    death_type      = collapse_unique(df$death_type),
    tier            = collapse_unique(df$tier),
    stringsAsFactors = FALSE
  )
}))

rownames(plot_summary) <- NULL

all_combo <- expand.grid(
  gene_symbol_std = top_genes,
  contrast        = names(deg_files),
  stringsAsFactors = FALSE
)

gene_label_map <- candidate_summary[
  candidate_summary$gene_symbol_std %in% top_genes,
  c("gene_symbol_std", "gene_symbol", "is_mcd_gene"),
  drop = FALSE
]

gene_label_map$display_gene <- ifelse(
  gene_label_map$is_mcd_gene,
  paste0(gene_label_map$gene_symbol, " *"),
  gene_label_map$gene_symbol
)

plot_summary <- merge(
  all_combo,
  plot_summary,
  by = c("gene_symbol_std", "contrast"),
  all.x = TRUE,
  sort = FALSE
)

plot_summary <- merge(
  plot_summary,
  gene_label_map[, c("gene_symbol_std", "display_gene")],
  by = "gene_symbol_std",
  all.x = TRUE
)

plot_summary$logFC[is.na(plot_summary$logFC)] <- 0
plot_summary$pathway_id[is.na(plot_summary$pathway_id)] <- ""

gene_order <- gene_label_map$display_gene[match(top_genes, gene_label_map$gene_symbol_std)]

plot_summary$display_gene <- factor(plot_summary$display_gene, levels = rev(gene_order))
plot_summary$contrast <- factor(
  plot_summary$contrast,
  levels = c("high_vs_low", "high_vs_intermediate"),
  labels = c("MCD-high vs MCD-low", "MCD-high vs MCD-intermediate")
)

p <- ggplot(plot_summary, aes(x = contrast, y = display_gene, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(pathway_id == "", "", "•")), size = 4) +
  scale_fill_gradient2(
    low = "#2c7bb6",
    mid = "white",
    high = "#d7191c",
    midpoint = 0,
    name = "log2FC"
  ) +
  labs(
    title = "Key mechanism candidate genes across major MCD contrasts",
    subtitle = "* indicates gene present in analysis-ready MCD gene set;\n• indicates the gene is a selected leading-edge member in that contrast",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 9),
    axis.text.x = element_text(angle = 20, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(
  file.path(out_dir, "tcga_paad_mcd_mechanism_summary_figure.pdf"),
  p,
  width = 8.6,
  height = max(5.5, 0.32 * length(top_genes) + 2.8)
)

if (has_pheatmap) {
  heat_wide <- reshape(
    plot_summary[, c("display_gene", "contrast", "logFC"), drop = FALSE],
    idvar = "display_gene",
    timevar = "contrast",
    direction = "wide"
  )
  
  rownames(heat_wide) <- heat_wide$display_gene
  heat_mat <- as.matrix(heat_wide[, -1, drop = FALSE])
  colnames(heat_mat) <- sub("^logFC\\.", "", colnames(heat_mat))
  
  pdf(
    file.path(out_dir, "tcga_paad_mcd_mechanism_candidate_heatmap.pdf"),
    width = 6.5,
    height = max(5.5, 0.3 * nrow(heat_mat) + 2.5)
  )
  pheatmap::pheatmap(
    heat_mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    main = "Key mechanism candidate genes across contrasts"
  )
  dev.off()
}

message("09C fixed analysis finished.")
message("输出文件：")
message("1) tcga_paad_mcd_high_vs_low_hallmark_leading_edge.tsv")
message("2) tcga_paad_mcd_high_vs_intermediate_hallmark_leading_edge.tsv")
message("3) tcga_paad_mcd_hallmark_leading_edge_mcd_overlap.tsv")
message("4) tcga_paad_mcd_key_mechanism_gene_candidates_detail.tsv")
message("5) tcga_paad_mcd_key_mechanism_gene_candidates.tsv")
message("6) tcga_paad_mcd_mechanism_summary_figure.pdf")
message("7) tcga_paad_mcd_mechanism_candidate_heatmap.pdf")
message("8) tcga_paad_mcd_hallmark_selected_pathways.tsv")