# =========================================================
# 09B_mcd_subtype_enrichment_analysis.R
# 修正说明：
#   [Fix1] msigdbr API 兼容新旧版本（collection / category）
#   [Fix2] msigdbr gene 列名兼容（gene_symbol / human_gene_symbol）
#   [Fix3] GSEA dotplot 加 split=".sign" + facet_grid（enrichplot >= 1.18）
#   [Fix4] PDF 设备用 finally 保护，防止设备泄漏
#   [Fix5] save_gsea_curve_pdf 实际过滤显著通路（fdr_thr）
#   [Fix6] GSEA 加 set.seed + seed=TRUE 保证可重复性
#   [Fix7] bitr 同时 suppressWarnings + suppressMessages
#   [Fix8] enrichKEGG 显式声明 use_internal_data + 有意义警告
# =========================================================

project_root <- Sys.getenv("PROJECT_ROOT", unset = "/Users/wmz_mac/Desktop/胰腺癌")
deg_dir      <- file.path(project_root, "04_bulk_analysis", "06_deg")
enrich_dir   <- file.path(project_root, "04_bulk_analysis", "07_enrichment")
dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# 输入文件
# -------------------------------
deg_files <- list(
  high_vs_low          = file.path(deg_dir, "tcga_paad_mcd_high_vs_low_deg.tsv"),
  high_vs_intermediate = file.path(deg_dir, "tcga_paad_mcd_high_vs_intermediate_deg.tsv"),
  intermediate_vs_low  = file.path(deg_dir, "tcga_paad_mcd_intermediate_vs_low_deg.tsv")
)

# -------------------------------
# 参数
# -------------------------------
fdr_cutoff   <- 0.05
logfc_cutoff <- 1

# -------------------------------
# 载入包
# -------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(enrichplot)
  library(ggplot2)
})

# -------------------------------
# 工具函数
# -------------------------------
read_deg_table <- function(file) {
  if (!file.exists(file)) stop("找不到文件: ", file)
  df <- fread(file, data.table = FALSE)
  required_cols <- c("gene_symbol", "logFC", "adj.P.Val")
  miss <- setdiff(required_cols, colnames(df))
  if (length(miss) > 0) {
    stop("DEG 文件缺少必要列: ", paste(miss, collapse = ", "), "\n文件: ", file)
  }
  if (!"t" %in% colnames(df)) {
    if ("P.Value" %in% colnames(df)) {
      df$t <- sign(df$logFC) * -log10(pmax(df$P.Value, 1e-300))
    } else {
      df$t <- df$logFC
    }
  }
  df <- df[!is.na(df$gene_symbol) & df$gene_symbol != "", , drop = FALSE]
  df <- df[!duplicated(df$gene_symbol), , drop = FALSE]
  return(df)
}

# [Fix7] bitr 同时压制 message 和 warning（unmapped gene 产生 warning）
symbol_to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])
  if (length(symbols) == 0) return(character(0))
  map <- suppressWarnings(suppressMessages(
    bitr(
      symbols,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    )
  ))
  if (is.null(map) || nrow(map) == 0) return(character(0))
  unique(map$ENTREZID)
}

safe_write_table <- function(df, out_file) {
  if (is.null(df)) return(invisible(NULL))
  fwrite(df, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

# [Fix3] 增加 is_gsea 参数；GSEA dotplot 需要 split=".sign" + facet_grid
# [Fix4] 用 finally 确保 dev.off() 始终执行
save_dotplot_pdf <- function(enrich_obj, out_file, title_text,
                             show_n = 15, is_gsea = FALSE) {
  if (is.null(enrich_obj)) return(invisible(NULL))
  df <- as.data.frame(enrich_obj)
  if (nrow(df) == 0) return(invisible(NULL))
  
  pdf(out_file, width = 8, height = 6)
  tryCatch({
    if (is_gsea) {
      # enrichplot >= 1.18：GSEA dotplot 必须传 split=".sign" 并追加 facet_grid
      p <- dotplot(
        enrich_obj,
        showCategory = min(show_n, nrow(df)),
        split        = ".sign"
      ) +
        facet_grid(. ~ .sign) +
        ggtitle(title_text) +
        theme(plot.title = element_text(hjust = 0.5))
    } else {
      p <- dotplot(enrich_obj, showCategory = min(show_n, nrow(df))) +
        ggtitle(title_text) +
        theme(plot.title = element_text(hjust = 0.5))
    }
    print(p)
  },
  error   = function(e) warning("dotplot 绘图失败 [", basename(out_file), "]: ",
                                conditionMessage(e)),
  finally = dev.off()   # [Fix4] 无论成败都关闭设备
  )
}

# [Fix5] 实际按 fdr_thr 过滤显著通路，无显著时 fallback 到全集
# [Fix4] 用 finally 确保 dev.off() 始终执行
save_gsea_curve_pdf <- function(gsea_obj, out_file, title_prefix,
                                top_n_each_side = 3, fdr_thr = 0.25) {
  if (is.null(gsea_obj)) return(invisible(NULL))
  gsea_df <- as.data.frame(gsea_obj)
  if (nrow(gsea_df) == 0) return(invisible(NULL))
  
  if ("p.adjust" %in% colnames(gsea_df)) {
    # 先筛显著通路；若无则 fallback 到全集
    sig_df <- gsea_df[gsea_df$p.adjust < fdr_thr, , drop = FALSE]
    if (nrow(sig_df) == 0) {
      message("  [GSEA curve] 无 p.adjust < ", fdr_thr, " 的通路，使用全部结果作为 fallback")
      sig_df <- gsea_df
    }
    sig_df <- sig_df[order(sig_df$p.adjust, -abs(sig_df$NES)), , drop = FALSE]
  } else {
    sig_df <- gsea_df[order(-abs(gsea_df$NES)), , drop = FALSE]
  }
  
  pos_ids <- head(sig_df$ID[sig_df$NES > 0], top_n_each_side)
  neg_ids <- head(sig_df$ID[sig_df$NES < 0], top_n_each_side)
  pathway_ids <- unique(c(pos_ids, neg_ids))
  pathway_ids <- pathway_ids[!is.na(pathway_ids)]
  if (length(pathway_ids) == 0) return(invisible(NULL))
  
  pdf(out_file, width = 8, height = max(4, 2.6 * length(pathway_ids)))
  tryCatch({
    for (pid in pathway_ids) {
      print(
        gseaplot2(
          gsea_obj,
          geneSetID = pid,
          title     = paste0(title_prefix, ": ", pid)
        )
      )
    }
  },
  error   = function(e) warning("gseaplot2 失败 [", basename(out_file), "]: ",
                                conditionMessage(e)),
  finally = dev.off()   # [Fix4]
  )
}

run_go_bp <- function(sig_symbols, universe_symbols) {
  if (length(sig_symbols) < 10) return(NULL)
  ego <- tryCatch(
    enrichGO(
      gene          = sig_symbols,
      universe      = universe_symbols,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.20,
      readable      = TRUE
    ),
    error = function(e) { warning("enrichGO 失败: ", conditionMessage(e)); NULL }
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  return(ego)
}

# [Fix8] 显式声明 use_internal_data，并给出有意义的网络错误提示
run_kegg <- function(sig_symbols, universe_symbols) {
  if (length(sig_symbols) < 10) return(NULL)
  sig_entrez      <- symbol_to_entrez(sig_symbols)
  universe_entrez <- symbol_to_entrez(universe_symbols)
  if (length(sig_entrez) < 10 || length(universe_entrez) < 10) return(NULL)
  
  ekegg <- tryCatch(
    enrichKEGG(
      gene              = sig_entrez,
      universe          = universe_entrez,
      organism          = "hsa",
      pAdjustMethod     = "BH",
      pvalueCutoff      = 0.05,
      qvalueCutoff      = 0.20,
      use_internal_data = FALSE   # [Fix8] 显式声明，需要网络访问 KEGG API
    ),
    error = function(e) {
      warning("enrichKEGG 失败（请检查网络连接 / KEGG API 可用性）: ",
              conditionMessage(e))
      NULL
    }
  )
  if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) return(NULL)
  
  ekegg <- tryCatch(
    setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
    error = function(e) ekegg
  )
  return(ekegg)
}

# [Fix6] 加 set.seed + seed=TRUE 保证 permutation 可重复
run_hallmark_gsea <- function(df_deg, hallmark_t2g) {
  rank_df <- df_deg[, c("gene_symbol", "t"), drop = FALSE]
  rank_df <- rank_df[!is.na(rank_df$gene_symbol) & !is.na(rank_df$t), , drop = FALSE]
  
  ord     <- order(abs(rank_df$t), decreasing = TRUE)
  rank_df <- rank_df[ord, , drop = FALSE]
  rank_df <- rank_df[!duplicated(rank_df$gene_symbol), , drop = FALSE]
  
  geneList        <- rank_df$t
  names(geneList) <- rank_df$gene_symbol
  geneList        <- sort(geneList, decreasing = TRUE)
  
  if (length(geneList) < 50) return(NULL)
  
  set.seed(42)   # [Fix6] 外层 seed 保证全局可重复
  gsea_res <- tryCatch(
    GSEA(
      geneList      = geneList,
      TERM2GENE     = hallmark_t2g,
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      minGSSize     = 10,
      maxGSSize     = 500,
      eps           = 0,
      seed          = TRUE,   # [Fix6] clusterProfiler 内部固定 permutation seed
      verbose       = FALSE
    ),
    error = function(e) { warning("GSEA 失败: ", conditionMessage(e)); NULL }
  )
  if (is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) return(NULL)
  return(gsea_res)
}

save_enrichment_result <- function(obj, out_tsv) {
  if (is.null(obj)) return(invisible(NULL))
  df <- as.data.frame(obj)
  if (nrow(df) == 0) return(invisible(NULL))
  safe_write_table(df, out_tsv)
}

# -------------------------------
# [Fix1+Fix2] Hallmark 基因集：兼容新旧 msigdbr API 及列名差异
# -------------------------------
hallmark_df <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "H"),  # msigdbr >= 2023.1
  error = function(e) {
    message("msigdbr: collection 参数不可用，回退到 category 参数（旧版 API）")
    msigdbr(species = "Homo sapiens", category = "H")   # msigdbr <= 7.x
  }
)

# [Fix2] 兼容列名：新版可能叫 human_gene_symbol
gene_col <- intersect(c("gene_symbol", "human_gene_symbol"), colnames(hallmark_df))[1]
if (is.na(gene_col)) {
  stop("无法识别 msigdbr 返回的基因列名，当前列名为: ",
       paste(colnames(hallmark_df), collapse = ", "))
}
hallmark_t2g <- unique(hallmark_df[, c("gs_name", gene_col)])
colnames(hallmark_t2g) <- c("term", "gene")

# -------------------------------
# 主循环
# -------------------------------
summary_list <- list()

for (contrast_name in names(deg_files)) {
  message("========== Running enrichment for: ", contrast_name, " ==========")
  
  deg_df           <- read_deg_table(deg_files[[contrast_name]])
  universe_symbols <- unique(deg_df$gene_symbol)
  
  up_symbols <- unique(
    deg_df$gene_symbol[
      !is.na(deg_df$adj.P.Val) & !is.na(deg_df$logFC) &
        deg_df$adj.P.Val < fdr_cutoff & deg_df$logFC >= logfc_cutoff
    ]
  )
  down_symbols <- unique(
    deg_df$gene_symbol[
      !is.na(deg_df$adj.P.Val) & !is.na(deg_df$logFC) &
        deg_df$adj.P.Val < fdr_cutoff & deg_df$logFC <= -logfc_cutoff
    ]
  )
  
  # ---- ORA: GO BP / KEGG ----
  go_up    <- run_go_bp(up_symbols,    universe_symbols)
  go_down  <- run_go_bp(down_symbols,  universe_symbols)
  kegg_up  <- run_kegg(up_symbols,    universe_symbols)
  kegg_down <- run_kegg(down_symbols, universe_symbols)
  
  save_enrichment_result(go_up,    file.path(enrich_dir, paste0("tcga_paad_mcd_", contrast_name, "_up_go_bp.tsv")))
  save_enrichment_result(go_down,  file.path(enrich_dir, paste0("tcga_paad_mcd_", contrast_name, "_down_go_bp.tsv")))
  save_enrichment_result(kegg_up,  file.path(enrich_dir, paste0("tcga_paad_mcd_", contrast_name, "_up_kegg.tsv")))
  save_enrichment_result(kegg_down, file.path(enrich_dir, paste0("tcga_paad_mcd_", contrast_name, "_down_kegg.tsv")))
  
  # ORA dotplot（is_gsea = FALSE，默认值）
  save_dotplot_pdf(go_up,    file.path(enrich_dir, paste0("tcga_paad_mcd_", contrast_name, "_up_go_bp_dotplot.pdf")),
                   paste0("GO BP (Up): ",   contrast_name))
  save_dotplot_pdf(go_down,  file.path(enrich_dir, paste0("tcga_paad_mcd_", contrast_name, "_down_go_bp_dotplot.pdf")),
                   paste0("GO BP (Down): ", contrast_name))
  save_dotplot_pdf(kegg_up,  file.path(enrich_dir, paste0("tcga_paad_mcd_", contrast_name, "_up_kegg_dotplot.pdf")),
                   paste0("KEGG (Up): ",   contrast_name))
  save_dotplot_pdf(kegg_down, file.path(enrich_dir, paste0("tcga_paad_mcd_", contrast_name, "_down_kegg_dotplot.pdf")),
                   paste0("KEGG (Down): ", contrast_name))
  
  # ---- GSEA: Hallmark ----
  hallmark_gsea <- run_hallmark_gsea(deg_df, hallmark_t2g)
  
  save_enrichment_result(
    hallmark_gsea,
    file.path(enrich_dir, paste0("tcga_paad_mcd_", contrast_name, "_hallmark_gsea.tsv"))
  )
  
  if (!is.null(hallmark_gsea) && nrow(as.data.frame(hallmark_gsea)) > 0) {
    # [Fix3] GSEA dotplot 传 is_gsea = TRUE
    save_dotplot_pdf(
      hallmark_gsea,
      file.path(enrich_dir, paste0("tcga_paad_mcd_", contrast_name, "_hallmark_gsea_dotplot.pdf")),
      paste0("Hallmark GSEA: ", contrast_name),
      show_n  = 20,
      is_gsea = TRUE   # [Fix3]
    )
    
    save_gsea_curve_pdf(
      hallmark_gsea,
      file.path(enrich_dir, paste0("tcga_paad_mcd_", contrast_name, "_hallmark_gsea_plot.pdf")),
      title_prefix = paste0("Hallmark GSEA (", contrast_name, ")")
    )
  }
  
  # ---- summary ----
  summary_list[[contrast_name]] <- data.frame(
    contrast         = contrast_name,
    n_universe_genes = length(universe_symbols),
    n_up_genes       = length(up_symbols),
    n_down_genes     = length(down_symbols),
    n_go_up          = ifelse(is.null(go_up),        0, nrow(as.data.frame(go_up))),
    n_go_down        = ifelse(is.null(go_down),      0, nrow(as.data.frame(go_down))),
    n_kegg_up        = ifelse(is.null(kegg_up),      0, nrow(as.data.frame(kegg_up))),
    n_kegg_down      = ifelse(is.null(kegg_down),    0, nrow(as.data.frame(kegg_down))),
    n_hallmark_gsea  = ifelse(is.null(hallmark_gsea), 0, nrow(as.data.frame(hallmark_gsea))),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_list)
safe_write_table(
  summary_df,
  file.path(enrich_dir, "tcga_paad_mcd_enrichment_summary.tsv")
)

message("09B enrichment analysis finished.")
print(summary_df)
