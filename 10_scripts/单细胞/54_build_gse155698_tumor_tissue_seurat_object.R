# =========================
# 54c_build_gse155698_tumor_tissue_seurat_object_force_triplet_fallback.R
# 目的：
# 1. 优先使用 10X 三件套
# 2. 若 import_plan 指向 h5 但无 hdf5r，则自动回退到同目录的三件套
# 3. 最大限度完成 GSE155698 tumor tissue 样本导入
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(Seurat)
})

setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "06_scRNA"

registry_dir <- file.path(root_dir, "00_registry")
proc_dir     <- file.path(root_dir, "02_processed_data")
table_dir    <- file.path(root_dir, "05_tables")
log_dir      <- file.path(root_dir, "06_logs")

dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

import_plan_file <- file.path(registry_dir, "54b_gse155698_import_plan.tsv")
if (!file.exists(import_plan_file)) stop("缺少 54b_gse155698_import_plan.tsv")

import_plan <- fread(import_plan_file)

# 只保留 tumor tissue
import_plan <- import_plan %>%
  filter(str_detect(sample_name, "PDAC_TISSUE"))

# ---------------------------------------------------------
# 工具函数：若 h5 不可读，则尝试同目录 fallback 到三件套
# ---------------------------------------------------------
resolve_triplet_from_h5 <- function(h5_path) {
  base_dir <- dirname(h5_path)
  alt_dir1 <- base_dir
  alt_dir2 <- file.path(dirname(base_dir), "filtered_feature_bc_matrix")
  alt_dir3 <- file.path(dirname(base_dir), "filtered_gene_bc_matrices")
  
  cand_dirs <- unique(c(alt_dir1, alt_dir2, alt_dir3))
  
  for (d in cand_dirs) {
    m1 <- file.path(d, "matrix.mtx.gz")
    b1 <- file.path(d, "barcodes.tsv.gz")
    f1 <- file.path(d, "features.tsv.gz")
    g1 <- file.path(d, "genes.tsv.gz")
    
    m2 <- file.path(d, "matrix.mtx")
    b2 <- file.path(d, "barcodes.tsv")
    f2 <- file.path(d, "features.tsv")
    g2 <- file.path(d, "genes.tsv")
    
    # gz 版本
    if (file.exists(m1) && file.exists(b1) && (file.exists(f1) || file.exists(g1))) {
      return(list(
        ok = TRUE,
        matrix_dir = d
      ))
    }
    
    # 非 gz 版本
    if (file.exists(m2) && file.exists(b2) && (file.exists(f2) || file.exists(g2))) {
      return(list(
        ok = TRUE,
        matrix_dir = d
      ))
    }
  }
  
  list(ok = FALSE, matrix_dir = NA_character_)
}

obj_list <- list()
qc_list <- list()
fail_list <- list()
final_plan <- list()

for (i in seq_len(nrow(import_plan))) {
  smp <- import_plan$sample_name[i]
  method <- import_plan$import_method[i]
  
  cat("Importing:", smp, "using", method, "...\n")
  
  mat <- NULL
  used_method <- method
  used_dir <- import_plan$matrix_dir[i]
  used_file <- import_plan$matrix_file[i]
  ok <- FALSE
  msg <- ""
  
  tryCatch({
    if (method == "Read10X") {
      mat <- Read10X(data.dir = used_dir)
      
    } else if (method == "Read10X_h5") {
      
      # 先尝试 h5
      if (requireNamespace("hdf5r", quietly = TRUE)) {
        mat <- Read10X_h5(used_file)
      } else {
        # 自动 fallback
        fb <- resolve_triplet_from_h5(used_file)
        if (!isTRUE(fb$ok)) {
          stop("hdf5r 未安装，且未找到可替代的三件套目录")
        }
        used_method <- "Read10X_fallback_from_h5"
        used_dir <- fb$matrix_dir
        mat <- Read10X(data.dir = used_dir)
      }
      
    } else {
      stop("Unknown import method")
    }
    
    seu <- CreateSeuratObject(
      counts = mat,
      project = "GSE155698_PDACCohort",
      min.cells = 3,
      min.features = 200
    )
    
    seu$orig.ident <- smp
    seu$sample_name <- smp
    seu$sample_scope <- "tumor_tissue"
    
    obj_list[[smp]] <- seu
    
    qc_list[[smp]] <- data.frame(
      sample_name = smp,
      import_method = used_method,
      import_status = "SUCCESS",
      n_cells = ncol(seu),
      n_genes = nrow(seu),
      median_nFeature_RNA = median(seu$nFeature_RNA),
      median_nCount_RNA = median(seu$nCount_RNA),
      stringsAsFactors = FALSE
    )
    
    final_plan[[smp]] <- data.frame(
      sample_name = smp,
      original_import_method = method,
      final_import_method = used_method,
      final_matrix_dir = used_dir,
      final_matrix_file = used_file,
      stringsAsFactors = FALSE
    )
    
    ok <- TRUE
    msg <- "success"
  }, error = function(e) {
    msg <<- conditionMessage(e)
  })
  
  if (!ok) {
    fail_list[[smp]] <- data.frame(
      sample_name = smp,
      import_method = method,
      import_status = "FAILED",
      error_message = msg,
      stringsAsFactors = FALSE
    )
  }
}

qc_tab <- bind_rows(qc_list)
fail_tab <- bind_rows(fail_list)
plan_tab <- bind_rows(final_plan)

fwrite(
  plan_tab,
  file.path(registry_dir, "54c_gse155698_final_import_plan.tsv"),
  sep = "\t", na = "NA"
)

fwrite(
  qc_tab,
  file.path(table_dir, "54c_gse155698_pre_qc_summary.tsv"),
  sep = "\t", na = "NA"
)

if (nrow(fail_tab) > 0) {
  fwrite(
    fail_tab,
    file.path(registry_dir, "54c_gse155698_failed_samples.tsv"),
    sep = "\t", na = "NA"
  )
}

if (length(obj_list) == 0) {
  stop("没有任何样本成功导入。")
}

merged_seu <- obj_list[[1]]
if (length(obj_list) > 1) {
  for (j in 2:length(obj_list)) {
    merged_seu <- merge(merged_seu, y = obj_list[[j]])
  }
}

saveRDS(
  merged_seu,
  file.path(proc_dir, "54c_gse155698_tumor_tissue_raw_seurat.rds")
)

writeLines(
  c(
    "54c completed",
    paste0("Created at: ", Sys.time()),
    "",
    "[Key fix]",
    "If h5 import is unavailable, automatically fall back to matrix/barcodes/features.",
    "",
    "[Next step]",
    "Proceed to QC metrics, mitochondrial percentage calculation, and filtering."
  ),
  file.path(log_dir, "54c_build_gse155698_tumor_tissue_seurat_notes.txt")
)

cat("54c_build_gse155698_tumor_tissue_seurat_object_force_triplet_fallback.R finished successfully.\n")