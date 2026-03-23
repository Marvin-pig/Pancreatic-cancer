# =========================================================
# 06A_create_metabolic_cell_death_starter_gene_lists.R
# 作用：
# 1. 创建 ferroptosis / cuproptosis / disulfidptosis starter gene list
# 2. 生成统一格式的 gene list TSV
# 3. 生成 master gene set table 及审计文件
# =========================================================

# --- 路径配置（优先读环境变量，回退到脚本所在目录）---
project_root <- Sys.getenv(
  "PANCREATIC_PROJECT_ROOT",
  unset = normalizePath(
    ifelse(nchar(getSrcDirectory(function() {})) > 0,
           getSrcDirectory(function() {}),
           ".")
  )
)

gene_dir <- file.path(project_root, "03_gene_sets")
dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

# --- 包含/排除标志常量 ---
INCLUDE <- 1L
EXCLUDE <- 0L

# =========================================================
# 辅助函数
# =========================================================

#' 构建标准化 gene 注释 data.frame
#'
#' @param genes       character vector，基因符号（自动转大写）
#' @param death_type  细胞死亡类型
#' @param subclass    基因功能子分类
#' @param evidence_tier  证据等级："A" 机制明确 / "B" 关联或扩展
#' @param source_class   证据来源类别
#' @param source_name    具体来源标签
#' @param include_main   是否纳入主分析（INCLUDE / EXCLUDE）
#' @param include_extended 是否纳入扩展分析
#' @param notes       备注
make_df <- function(genes, death_type, subclass, evidence_tier,
                    source_class, source_name,
                    include_main = INCLUDE,
                    include_extended = INCLUDE,
                    notes = "") {
  stopifnot(
    length(genes) > 0,
    evidence_tier %in% c("A", "B"),
    include_main %in% c(INCLUDE, EXCLUDE),
    include_extended %in% c(INCLUDE, EXCLUDE)
  )
  data.frame(
    gene_symbol      = toupper(trimws(genes)),
    death_type       = death_type,
    subclass         = subclass,
    evidence_tier    = evidence_tier,
    source_class     = source_class,
    source_name      = source_name,
    include_main     = include_main,
    include_extended = include_extended,
    notes            = notes,
    stringsAsFactors = FALSE
  )
}

#' 写出 TSV，统一 quote = FALSE
write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# =========================================================
# 1. Ferroptosis
#    Tier A: 机制明确的 canonical driver / suppressor
#    Tier B: 标志物或关联基因（扩展分析用）
#
# 调整说明：
#   - DPP4 : 仅在无 KRAS 突变结直肠癌中有直接证据，降为 Tier B
#   - SLC1A5: 谷氨酰胺转运，间接参与，降为 Tier B
#   - NOX1/NOX4: 背景 ROS 来源，降为 Tier B（context-specific）
# =========================================================

ferro_driver_A <- c(
  "ACSL4", "ALOX15", "LPCAT3", "NCOA4",
  "SAT1", "TFRC"
)

ferro_suppressor_A <- c(
  "GPX4", "SLC7A11", "AIFM2", "NFE2L2", "KEAP1",
  "FTH1", "FTL", "GCLC", "GCLM", "GSS",
  "CISD1", "NFS1", "FANCD2"
)

# Tier B：背景 ROS / 间接调控 / 仅部分背景有证据
ferro_extended_B <- c(
  "NOX1", "NOX4", "DPP4", "SLC1A5", "EMC2",
  "HMOX1", "ALOX5", "HSPB1", "GLS2", "MT1G", "FDFT1", "RPL8"
)

ferro_df <- rbind(
  make_df(ferro_driver_A, "ferroptosis", "driver", "A",
          "curated_database_or_seminal", "starter_curated_core",
          INCLUDE, INCLUDE,
          "canonical driver; recommended for main analysis"),
  make_df(ferro_suppressor_A, "ferroptosis", "suppressor", "A",
          "curated_database_or_seminal", "starter_curated_core",
          INCLUDE, INCLUDE,
          "canonical suppressor; recommended for main analysis"),
  make_df(ferro_extended_B, "ferroptosis", "associated_or_context_specific", "B",
          "review_or_extended_literature", "starter_extended",
          EXCLUDE, INCLUDE,
          "extended only; context-specific or indirect evidence")
)
ferro_df <- unique(ferro_df)

# =========================================================
# 2. Cuproptosis
#    Tier A: 铜离子-蛋白质毒性直接机制（FDX1 → 蛋白脂酰化 → 聚集）
#    Tier B: 铜转运 / 代谢扩展
#
# 调整说明：
#   - CDKN2A: 下游调控效应，非直接机制；降为 Tier B
#   - MTF1   : 铜响应转录因子，调控层，降为 Tier B
# =========================================================

cupro_core_A <- c(
  "FDX1", "LIAS", "LIPT1", "DLD", "DLAT",
  "PDHA1", "PDHB", "GLS", "DBT"
)

cupro_extended_B <- c(
  "SLC31A1", "ATP7A", "ATP7B",
  "GCSH", "DLST",
  "MTF1", "CDKN2A"
)

cupro_df <- rbind(
  make_df(cupro_core_A, "cuproptosis", "canonical_regulator", "A",
          "seminal_mechanism", "starter_curated_core",
          INCLUDE, INCLUDE,
          "direct cuproptosis mechanism; recommended for main analysis"),
  make_df(cupro_extended_B, "cuproptosis", "transport_or_regulatory", "B",
          "review_or_extended_literature", "starter_extended",
          EXCLUDE, INCLUDE,
          "extended only; copper transport or downstream regulatory layer")
)
cupro_df <- unique(cupro_df)

# =========================================================
# 3. Disulfidptosis
#    Tier A: SLC7A11 高表达 + 葡萄糖饥饿 → actin 塌陷直接参与
#    Tier B: Rac-WAVE 相关 / actin 骨架扩展
#
# 注：SLC7A11 同时出现于 ferroptosis（suppressor）和 disulfidptosis（core）
#     两者 death_type 不同，在 master 中保留两行，属于已知生物学重叠。
# =========================================================

disulf_core_A <- c(
  "SLC7A11", "RAC1", "NCKAP1", "ABI1", "BRK1", "CYFIP1", "WASF2"
)

disulf_extended_B <- c(
  "NCKAP1L", "WASF1", "ACTB", "FLNA", "TLN1"
)

disulf_df <- rbind(
  make_df(disulf_core_A, "disulfidptosis", "canonical_or_direct_mechanism", "A",
          "seminal_mechanism", "starter_curated_core",
          INCLUDE, INCLUDE,
          "direct disulfidptosis mechanism; recommended for main analysis"),
  make_df(disulf_extended_B, "disulfidptosis", "cytoskeleton_extended", "B",
          "review_or_extended_literature", "starter_extended",
          EXCLUDE, INCLUDE,
          "extended only; actin cytoskeleton related")
)
disulf_df <- unique(disulf_df)

# =========================================================
# 4. 导出各类 starter TSV（统一格式）
# =========================================================

write_tsv(ferro_df,
          file.path(gene_dir, "ferroptosis_genes_starter.tsv"))
write_tsv(cupro_df,
          file.path(gene_dir, "cuproptosis_genes_starter.tsv"))
write_tsv(disulf_df,
          file.path(gene_dir, "disulfidptosis_genes_starter.tsv"))

# =========================================================
# 5. Master gene set（保留跨类型重叠行，以 death_type 区分）
# =========================================================

master_df <- rbind(ferro_df, cupro_df, disulf_df)
master_df <- unique(master_df)
master_df <- master_df[
  order(master_df$death_type,
        master_df$include_main == EXCLUDE,   # include_main=1 排前
        master_df$gene_symbol),
]

write_tsv(
  master_df,
  file.path(gene_dir, "metabolic_cell_death_gene_sets_master_starter.tsv")
)

# =========================================================
# 6. 审计摘要
# =========================================================

make_audit <- function(df, type) {
  data.frame(
    death_type      = type,
    n_total         = nrow(df),
    n_main          = sum(df$include_main == INCLUDE),
    n_extended_only = sum(df$include_main == EXCLUDE &
                          df$include_extended == INCLUDE),
    n_tier_A        = sum(df$evidence_tier == "A"),
    n_tier_B        = sum(df$evidence_tier == "B")
  )
}

audit_df <- rbind(
  make_audit(ferro_df,  "ferroptosis"),
  make_audit(cupro_df,  "cuproptosis"),
  make_audit(disulf_df, "disulfidptosis")
)

write_tsv(
  audit_df,
  file.path(gene_dir, "metabolic_cell_death_gene_sets_starter_audit.tsv")
)

# =========================================================
# 7. 可重复性：记录 session 信息
# =========================================================

writeLines(
  capture.output(sessionInfo()),
  file.path(gene_dir, "session_info.txt")
)

# =========================================================
# 8. 控制台摘要
# =========================================================

cat("\n=== Starter gene lists created successfully ===\n")
cat("Output directory:", gene_dir, "\n\n")
print(audit_df, row.names = FALSE)

# 提示跨类型重叠基因
cross_genes <- intersect(
  ferro_df$gene_symbol,
  c(cupro_df$gene_symbol, disulf_df$gene_symbol)
)
if (length(cross_genes) > 0) {
  cat("\n[INFO] Cross-death-type overlapping genes (expected, kept separately):\n")
  cat(" ", paste(cross_genes, collapse = ", "), "\n")
}
