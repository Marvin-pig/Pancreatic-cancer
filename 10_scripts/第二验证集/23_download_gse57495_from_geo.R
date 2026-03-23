# =========================
# 23_download_gse57495_from_geo.R
# 目的：
# 1. 下载 GSE57495
# 2. 保存 ExpressionSet / expression / pheno / feature
# 3. 为后续构建 gene-level expression 和 survival-ready clinical table 做准备
# =========================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(data.table)
})
setwd("/Users/wmz_mac/Desktop/胰腺癌")
out_dir <- "04_bulk_analysis/04_external_validation/01_raw_data/GEO/GSE57495"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Downloading GSE57495 from GEO...\n")

gset_list <- getGEO(
  "GSE57495",
  GSEMatrix = TRUE,
  getGPL = FALSE,
  destdir = out_dir
)

cat("Number of ExpressionSet objects:", length(gset_list), "\n")

gset <- gset_list[[1]]
saveRDS(gset, file.path(out_dir, "23_GSE57495_gset.rds"))

# expression
expr_mat <- exprs(gset)
expr_df <- data.frame(
  probe_id = rownames(expr_mat),
  expr_mat,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

fwrite(
  expr_df,
  file.path(out_dir, "23_GSE57495_expression_probe_level.tsv"),
  sep = "\t"
)

# pheno
pheno_df <- pData(gset)
pheno_df <- data.frame(
  geo_accession = rownames(pheno_df),
  pheno_df,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

fwrite(
  pheno_df,
  file.path(out_dir, "23_GSE57495_pheno_raw.tsv"),
  sep = "\t"
)

# feature
feat_df <- fData(gset)
feat_df <- data.frame(
  probe_id = rownames(feat_df),
  feat_df,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

fwrite(
  feat_df,
  file.path(out_dir, "23_GSE57495_feature_raw.tsv"),
  sep = "\t"
)

cat("Download completed.\n")
cat("Generated files:\n")
cat("- 23_GSE57495_gset.rds\n")
cat("- 23_GSE57495_expression_probe_level.tsv\n")
cat("- 23_GSE57495_pheno_raw.tsv\n")
cat("- 23_GSE57495_feature_raw.tsv\n")