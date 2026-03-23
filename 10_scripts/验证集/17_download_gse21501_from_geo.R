# =========================
# 17_download_gse21501_from_geo.R
# 目的：
# 1. 下载 GEO GSE21501
# 2. 保存 Series Matrix / ExpressionSet
# 3. 导出原始表达表型对象，供后续整理表达矩阵和临床表
# =========================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(data.table)
})
setwd("/Users/wmz_mac/Desktop/胰腺癌")
out_dir <- "04_bulk_analysis/04_external_validation/01_raw_data/GEO/GSE21501"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Downloading GSE21501 from GEO...\n")

gset_list <- getGEO(
  "GSE21501",
  GSEMatrix = TRUE,
  getGPL = FALSE,
  destdir = out_dir
)

cat("Number of ExpressionSet objects:", length(gset_list), "\n")

# 如果有多个平台，先取第一个；后面再核查
gset <- gset_list[[1]]

saveRDS(gset, file.path(out_dir, "17_GSE21501_gset.rds"))

# 导出表达矩阵
expr_mat <- exprs(gset)
expr_df <- data.frame(
  probe_id = rownames(expr_mat),
  expr_mat,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

fwrite(
  expr_df,
  file.path(out_dir, "17_GSE21501_expression_probe_level.tsv"),
  sep = "\t"
)

# 导出样本注释
pheno_df <- pData(gset)
pheno_df <- data.frame(
  geo_accession = rownames(pheno_df),
  pheno_df,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

fwrite(
  pheno_df,
  file.path(out_dir, "17_GSE21501_pheno_raw.tsv"),
  sep = "\t"
)

# 导出 feature 注释（如果有）
feat_df <- fData(gset)
feat_df <- data.frame(
  probe_id = rownames(feat_df),
  feat_df,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

fwrite(
  feat_df,
  file.path(out_dir, "17_GSE21501_feature_raw.tsv"),
  sep = "\t"
)

cat("Download completed.\n")
cat("Files generated:\n")
cat("- 17_GSE21501_gset.rds\n")
cat("- 17_GSE21501_expression_probe_level.tsv\n")
cat("- 17_GSE21501_pheno_raw.tsv\n")
cat("- 17_GSE21501_feature_raw.tsv\n")
