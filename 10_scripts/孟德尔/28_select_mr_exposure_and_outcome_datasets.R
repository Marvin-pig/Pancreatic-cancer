# =========================
# 28_select_mr_exposure_and_outcome_datasets.R
# 目的：
# 1. 冻结 MR 的 exposure / outcome 数据源
# 2. 确定第一批优先 MR 基因
# 3. 生成后续下载与执行清单
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

root_dir <- "05_MR"
dir.create(file.path(root_dir, "00_registry"), recursive = TRUE, showWarnings = FALSE)

target_file <- file.path(root_dir, "00_registry/27_mr_target_gene_registry.tsv")
targets <- fread(target_file)

# ---------- 第一批优先基因 ----------
priority_genes <- c(
  "MET", "TGFBI", "IDO1", "INHBA", "IL15RA",
  "GBP4", "GPX4", "SLC31A1", "TGFBR3", "CXCL10"
)

mr_gene_plan <- targets %>%
  mutate(
    batch = ifelse(gene %in% priority_genes, "batch1_priority", "batch2_extended"),
    run_priority = case_when(
      gene %in% priority_genes ~ "HIGH",
      target_class == "signature_gene" ~ "MEDIUM",
      TRUE ~ "LOW"
    )
  ) %>%
  arrange(match(run_priority, c("HIGH", "MEDIUM", "LOW")), gene)

fwrite(
  mr_gene_plan,
  file.path(root_dir, "00_registry/28_mr_gene_execution_plan.tsv"),
  sep = "\t"
)

# ---------- exposure 数据集计划 ----------
exposure_plan <- data.frame(
  exposure_dataset = c(
    "eQTLGen_whole_blood_cis_eQTL",
    "GTEx_v8_pancreas_cis_eQTL"
  ),
  role = c("primary_exposure_source", "tissue_specific_supporting_exposure"),
  ancestry = c("predominantly_European", "mixed_but_mainly_European"),
  priority = c("HIGH", "MEDIUM"),
  access_method = c(
    "download summary stats / OpenGWAS-compatible processing",
    "GTEx portal / processed cis-eQTL resources"
  ),
  note = c(
    "best first-pass broad coverage for candidate genes",
    "used for pancreas-specific biological support, likely fewer genes/instruments"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  exposure_plan,
  file.path(root_dir, "00_registry/28_mr_exposure_dataset_plan.tsv"),
  sep = "\t"
)

# ---------- outcome 数据集计划 ----------
outcome_plan <- data.frame(
  outcome_dataset = c(
    "PanScan_PanC4_European_pancreatic_cancer_GWAS",
    "BBJ_pancreatic_cancer_bbj-a-140"
  ),
  role = c("primary_outcome", "ancestry_sensitivity_outcome"),
  ancestry = c("European", "East_Asian"),
  priority = c("HIGH", "MEDIUM"),
  access_method = c(
    "public summary statistics / curated GWAS source",
    "OpenGWAS"
  ),
  note = c(
    "main pancreatic cancer risk outcome for primary MR",
    "do not use as primary result; use as sensitivity analysis"
  ),
  stringsAsFactors = FALSE
)

fwrite(
  outcome_plan,
  file.path(root_dir, "00_registry/28_mr_outcome_dataset_plan.tsv"),
  sep = "\t"
)

# ---------- 下载与执行清单 ----------
todo_plan <- data.frame(
  step_id = c("28A", "28B", "29A", "29B", "30"),
  task = c(
    "prepare OpenGWAS account/JWT if using API-based outcome access",
    "confirm source and file path for PanScan/PanC4 summary statistics",
    "obtain eQTLGen cis-eQTL summary statistics or gene-level query interface",
    "obtain GTEx pancreas cis-eQTL resource for backup/supporting analysis",
    "run primary MR for batch1_priority genes"
  ),
  status = c("TODO", "TODO", "TODO", "TODO", "TODO"),
  stringsAsFactors = FALSE
)

fwrite(
  todo_plan,
  file.path(root_dir, "00_registry/28_mr_download_and_execution_todo.tsv"),
  sep = "\t"
)

# ---------- 操作说明 ----------
notes <- c(
  "MR dataset selection notes",
  paste0("Created at: ", Sys.time()),
  "",
  "[Primary exposure]",
  "eQTLGen whole blood cis-eQTL",
  "",
  "[Supporting exposure]",
  "GTEx v8 pancreas cis-eQTL",
  "",
  "[Primary outcome]",
  "PanScan/PanC4 European pancreatic cancer GWAS",
  "",
  "[Sensitivity outcome]",
  "BBJ pancreatic cancer (bbj-a-140)",
  "",
  "[Execution principle]",
  "Run batch1_priority genes first, then expand if signal quality is acceptable."
)

writeLines(
  notes,
  file.path(root_dir, "06_logs/28_mr_dataset_selection_notes.txt")
)

cat("28_select_mr_exposure_and_outcome_datasets.R finished successfully.\n")