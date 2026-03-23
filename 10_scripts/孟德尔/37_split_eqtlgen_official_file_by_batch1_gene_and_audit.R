# =========================
# 37_split_eqtlgen_official_file_by_batch1_gene_and_audit.R
# 目的：
# 1. 读取 eQTLGen 官方 significant cis-eQTL 主文件
# 2. 用 GeneSymbol（优先）检索 batch1 基因
# 3. 拆出每个基因的原始候选表
# 4. 生成全表命中审计、行数统计、结构审计
# 5. 暂不覆盖 01_raw_data/exposure/eQTLGen/<GENE>_cis_eqtl.tsv
#    因为当前还没完成 beta / se / eaf 的标准化
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ---------- 路径 ----------
setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "05_MR"

batch1_file <- file.path(root_dir, "02_processed_data/batch1_queries/29_batch1_priority_genes.tsv")
sig_file    <- file.path(root_dir, "01_raw_data/exposure/eQTLGen/36_eqtlgen_significant_cis_eqtls.txt.gz")

out_registry <- file.path(root_dir, "00_registry")
out_proc_dir <- file.path(root_dir, "02_processed_data/37_candidate_gene_tables")
out_logs     <- file.path(root_dir, "06_logs")

dir.create(out_registry, recursive = TRUE, showWarnings = FALSE)
dir.create(out_proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_logs, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(batch1_file)) {
  stop("未找到 batch1 基因文件: ", batch1_file)
}
if (!file.exists(sig_file)) {
  stop("未找到 eQTLGen 主文件: ", sig_file,
       "\n请先确认 36_eqtlgen_significant_cis_eqtls.txt.gz 已下载。")
}

# ---------- 读取 batch1 ----------
batch1 <- fread(batch1_file)
if (!"gene" %in% colnames(batch1)) {
  stop("29_batch1_priority_genes.tsv 缺少 gene 列。")
}
batch1_genes <- unique(as.character(batch1$gene))

# ---------- 读取表头并识别关键列 ----------
header <- names(fread(sig_file, nrows = 0))

detect_one <- function(cands, header_vec) {
  hit <- cands[cands %in% header_vec]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

col_gene_symbol <- detect_one(c("GeneSymbol", "gene_symbol", "SYMBOL", "symbol"), header)
col_gene_id     <- detect_one(c("Gene", "gene_id", "ENSEMBL", "EnsemblID"), header)
col_snp         <- detect_one(c("SNP", "rsid", "variant_id"), header)
col_chr         <- detect_one(c("SNPChr", "chr", "CHR"), header)
col_pos         <- detect_one(c("SNPPos", "pos", "BP"), header)
col_ea          <- detect_one(c("AssessedAllele", "effect_allele", "EA", "A1"), header)
col_oa          <- detect_one(c("OtherAllele", "other_allele", "OA", "A2"), header)
col_p           <- detect_one(c("Pvalue", "pval", "P", "p_value"), header)
col_z           <- detect_one(c("Zscore", "z", "Z"), header)
col_n           <- detect_one(c("NrSamples", "N", "samplesize"), header)

if (is.na(col_gene_symbol)) {
  stop("未识别到 GeneSymbol 列。当前表头为：\n", paste(header, collapse = ", "))
}

# ---------- 只读取必要列，降低内存 ----------
read_cols <- unique(na.omit(c(
  col_gene_symbol, col_gene_id, col_snp, col_chr, col_pos,
  col_ea, col_oa, col_p, col_z, col_n
)))

message("读取 eQTLGen 主文件必要列中...")
eqtl <- fread(sig_file, select = read_cols)

# 统一基因符号列
eqtl[[col_gene_symbol]] <- as.character(eqtl[[col_gene_symbol]])
eqtl$.__gene_symbol__.  <- eqtl[[col_gene_symbol]]

# ---------- 全表命中审计 ----------
gene_hit_audit <- data.frame(
  gene = batch1_genes,
  hit_in_full_table = batch1_genes %in% unique(eqtl$.__gene_symbol__.),
  stringsAsFactors = FALSE
)

# ---------- 按基因拆分原始候选表 ----------
gene_row_counts <- lapply(batch1_genes, function(g) {
  sub <- eqtl[.__gene_symbol__. == g]
  
  if (nrow(sub) > 0) {
    out_file <- file.path(out_proc_dir, paste0(g, "_eqtlgen_raw.tsv"))
    fwrite(sub, out_file, sep = "\t", na = "NA")
  }
  
  data.frame(
    gene = g,
    n_rows = nrow(sub),
    out_file = ifelse(nrow(sub) > 0,
                      file.path(out_proc_dir, paste0(g, "_eqtlgen_raw.tsv")),
                      NA_character_),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

# ---------- 结构审计 ----------
# 这里审计的是“能否向 30 号模板靠拢”，不是说现在已经合格
# 30 号模板需要：
# SNP, chr, pos, effect_allele, other_allele, beta, se, pval, eaf, gene, exposure
# 当前官方 significant 文件通常只有 SNP/chr/pos/alleles/pval/Zscore/Gene/GeneSymbol/NrSamples 等
# 因此这里明确标记缺口，不强行伪造 beta / se / eaf

structural_audit <- gene_row_counts %>%
  mutate(
    has_raw_table = n_rows > 0,
    has_SNP = !is.na(col_snp),
    has_chr = !is.na(col_chr),
    has_pos = !is.na(col_pos),
    has_effect_allele = !is.na(col_ea),
    has_other_allele  = !is.na(col_oa),
    has_pval = !is.na(col_p),
    has_gene_symbol = !is.na(col_gene_symbol),
    has_gene_id = !is.na(col_gene_id),
    has_zscore = !is.na(col_z),
    has_samplesize = !is.na(col_n),
    has_beta = FALSE,
    has_se   = FALSE,
    has_eaf  = FALSE,
    convertible_to_30_template_now = FALSE,
    reason = case_when(
      !has_raw_table ~ "gene_not_found_in_full_table",
      has_raw_table & (!has_SNP | !has_effect_allele | !has_other_allele | !has_pval) ~
        "raw_table_found_but_basic_variant_columns_incomplete",
      has_raw_table & !has_beta & !has_se & !has_eaf ~
        "raw_table_found_but_beta_se_eaf_missing_need_additional_reconstruction_or_support_files",
      TRUE ~ "needs_manual_review"
    )
  )

# ---------- 保存表头映射 ----------
schema_map <- data.frame(
  target_field = c(
    "gene_symbol", "gene_id", "SNP", "chr", "pos",
    "effect_allele", "other_allele", "pval", "zscore", "samplesize"
  ),
  detected_column = c(
    col_gene_symbol, col_gene_id, col_snp, col_chr, col_pos,
    col_ea, col_oa, col_p, col_z, col_n
  ),
  stringsAsFactors = FALSE
)

# ---------- 生成候选标准化预览 ----------
# 不写入正式 exposure 文件，只生成每个基因的“候选标准化预览”
# 目的：先让你看到哪些字段已有，哪些仍为空
candidate_preview_rows <- list()

for (g in batch1_genes) {
  f <- file.path(out_proc_dir, paste0(g, "_eqtlgen_raw.tsv"))
  if (!file.exists(f)) next
  
  x <- fread(f, nrows = 5)
  
  y <- data.table(
    SNP           = if (!is.na(col_snp) && col_snp %in% names(x)) as.character(x[[col_snp]]) else NA_character_,
    chr           = if (!is.na(col_chr) && col_chr %in% names(x)) x[[col_chr]] else NA,
    pos           = if (!is.na(col_pos) && col_pos %in% names(x)) x[[col_pos]] else NA,
    effect_allele = if (!is.na(col_ea)  && col_ea  %in% names(x)) as.character(x[[col_ea]]) else NA_character_,
    other_allele  = if (!is.na(col_oa)  && col_oa  %in% names(x)) as.character(x[[col_oa]]) else NA_character_,
    beta          = NA_real_,
    se            = NA_real_,
    pval          = if (!is.na(col_p)   && col_p   %in% names(x)) x[[col_p]] else NA,
    eaf           = NA_real_,
    gene          = g,
    exposure      = paste0("eQTLGen_", g)
  )
  
  fwrite(
    y,
    file.path(out_proc_dir, paste0(g, "_candidate_standardized_preview.tsv")),
    sep = "\t",
    na = "NA"
  )
  
  candidate_preview_rows[[g]] <- data.frame(
    gene = g,
    preview_file = file.path(out_proc_dir, paste0(g, "_candidate_standardized_preview.tsv")),
    stringsAsFactors = FALSE
  )
}

candidate_preview_index <- bind_rows(candidate_preview_rows)

# ---------- 输出 ----------
fwrite(
  gene_hit_audit,
  file.path(out_registry, "37_batch1_gene_full_hit_audit.tsv"),
  sep = "\t"
)

fwrite(
  gene_row_counts,
  file.path(out_registry, "37_eqtlgen_gene_row_counts.tsv"),
  sep = "\t"
)

fwrite(
  structural_audit,
  file.path(out_registry, "37_exposure_structural_audit.tsv"),
  sep = "\t"
)

fwrite(
  schema_map,
  file.path(out_registry, "37_eqtlgen_schema_map.tsv"),
  sep = "\t"
)

if (nrow(candidate_preview_index) > 0) {
  fwrite(
    candidate_preview_index,
    file.path(out_registry, "37_candidate_preview_index.tsv"),
    sep = "\t"
  )
}

# ---------- 日志 ----------
summary_lines <- c(
  "=== Code 37 Summary ===",
  paste0("Created at: ", Sys.time()),
  "",
  "[Input files]",
  paste0("batch1_file: ", batch1_file),
  paste0("sig_file   : ", sig_file),
  "",
  "[Detected columns]",
  paste0("GeneSymbol     : ", ifelse(is.na(col_gene_symbol), "NA", col_gene_symbol)),
  paste0("Gene ID        : ", ifelse(is.na(col_gene_id), "NA", col_gene_id)),
  paste0("SNP            : ", ifelse(is.na(col_snp), "NA", col_snp)),
  paste0("chr            : ", ifelse(is.na(col_chr), "NA", col_chr)),
  paste0("pos            : ", ifelse(is.na(col_pos), "NA", col_pos)),
  paste0("effect_allele  : ", ifelse(is.na(col_ea), "NA", col_ea)),
  paste0("other_allele   : ", ifelse(is.na(col_oa), "NA", col_oa)),
  paste0("pval           : ", ifelse(is.na(col_p), "NA", col_p)),
  paste0("zscore         : ", ifelse(is.na(col_z), "NA", col_z)),
  paste0("samplesize     : ", ifelse(is.na(col_n), "NA", col_n)),
  "",
  "[Hit summary]",
  paste0("Batch1 genes total     : ", length(batch1_genes)),
  paste0("Genes found in full tbl: ", sum(gene_hit_audit$hit_in_full_table)),
  paste0("Genes not found        : ", sum(!gene_hit_audit$hit_in_full_table)),
  "",
  "[Rule]",
  "Do NOT overwrite 01_raw_data/exposure/eQTLGen/<GENE>_cis_eqtl.tsv yet.",
  "Current outputs are raw candidate tables and previews only.",
  "",
  "[Next step]",
  "Use code 38 to decide whether beta/se can be reconstructed and whether EAF can be joined."
)

writeLines(
  summary_lines,
  file.path(out_logs, "37_split_eqtlgen_and_audit_notes.txt")
)

cat("37_split_eqtlgen_official_file_by_batch1_gene_and_audit.R finished successfully.\n")
cat("Generated files:\n")
cat("- 00_registry/37_batch1_gene_full_hit_audit.tsv\n")
cat("- 00_registry/37_eqtlgen_gene_row_counts.tsv\n")
cat("- 00_registry/37_exposure_structural_audit.tsv\n")
cat("- 00_registry/37_eqtlgen_schema_map.tsv\n")
cat("- 00_registry/37_candidate_preview_index.tsv (if any genes found)\n")
cat("- 02_processed_data/37_candidate_gene_tables/<GENE>_eqtlgen_raw.tsv\n")
cat("- 02_processed_data/37_candidate_gene_tables/<GENE>_candidate_standardized_preview.tsv\n")
cat("- 06_logs/37_split_eqtlgen_and_audit_notes.txt\n")