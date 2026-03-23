# =========================================================
# 07C_compare_k2_k3_k4_survival.R
# =========================================================
project_root <- Sys.getenv("PROJECT_ROOT", unset = "/Users/wmz_mac/Desktop/胰腺癌")
proc_dir  <- file.path(project_root, "02_processed_data", "TCGA_PAAD")
clust_dir <- file.path(project_root, "04_bulk_analysis", "02_clustering")

cohort_file <- file.path(proc_dir, "tcga_paad_cohort_master_ordered.tsv")
k2_file <- file.path(clust_dir, "tcga_paad_cluster_assignment_k2.tsv")
k3_file <- file.path(clust_dir, "tcga_paad_cluster_assignment_k3.tsv")
k4_file <- file.path(clust_dir, "tcga_paad_cluster_assignment_k4.tsv")

library(data.table)
library(survival)

cohort <- fread(cohort_file, data.table = FALSE)
k2     <- fread(k2_file,     data.table = FALSE)
k3     <- fread(k3_file,     data.table = FALSE)
k4     <- fread(k4_file,     data.table = FALSE)

# ---- detect sample id column ----
detect_sample_col <- function(df) {
  cand <- c("sample_id", "sample_barcode", "SampleID", "sample", "barcode", "submitter_id")
  hit  <- cand[cand %in% colnames(df)]
  if (length(hit) == 0) stop(
    sprintf("未找到 sample_id 列，现有列：%s", paste(colnames(df), collapse = ", "))
  )
  hit[1]
}

sample_col    <- detect_sample_col(cohort)
k2_sample_col <- detect_sample_col(k2)
k3_sample_col <- detect_sample_col(k3)
k4_sample_col <- detect_sample_col(k4)

# ---- detect endpoint columns ----
time_col  <- intersect(c("OS.time", "OS_time", "os_time"), colnames(cohort))
event_col <- intersect(c("OS.status", "OS", "OS_event", "os_event"), colnames(cohort))
if (length(time_col)  == 0) stop("未找到 OS time 列。")
if (length(event_col) == 0) stop("未找到 OS event 列。")
time_col  <- time_col[1]
event_col <- event_col[1]

# ---- 重命名 k-file 的 cluster 列，避免 merge 后靠位置猜列名 ----
rename_cluster_col <- function(kdf, sample_col_name, new_name) {
  other_cols <- setdiff(colnames(kdf), sample_col_name)
  if (length(other_cols) == 0) stop(sprintf("%s 中未找到 cluster 列", new_name))
  colnames(kdf)[colnames(kdf) == other_cols[1]] <- new_name
  kdf
}

k2 <- rename_cluster_col(k2, k2_sample_col, "cluster_k2")
k3 <- rename_cluster_col(k3, k3_sample_col, "cluster_k3")
k4 <- rename_cluster_col(k4, k4_sample_col, "cluster_k4")

# ---- merge ----
dat <- merge(cohort, k2[, c(k2_sample_col, "cluster_k2")],
             by.x = sample_col, by.y = k2_sample_col,
             all.x = TRUE, sort = FALSE)
dat <- merge(dat,   k3[, c(k3_sample_col, "cluster_k3")],
             by.x = sample_col, by.y = k3_sample_col,
             all.x = TRUE, sort = FALSE)
dat <- merge(dat,   k4[, c(k4_sample_col, "cluster_k4")],
             by.x = sample_col, by.y = k4_sample_col,
             all.x = TRUE, sort = FALSE)

# ---- restore original row order ----
dat <- dat[match(cohort[[sample_col]], dat[[sample_col]]), , drop = FALSE]

# ---- check ----
stopifnot(identical(dat[[sample_col]], cohort[[sample_col]]))
na_check <- c(
  cluster_k2 = sum(is.na(dat$cluster_k2)),
  cluster_k3 = sum(is.na(dat$cluster_k3)),
  cluster_k4 = sum(is.na(dat$cluster_k4))
)
if (any(na_check > 0)) {
  stop(sprintf("Merge 后存在 NA：%s",
               paste(names(na_check[na_check > 0]),
                     na_check[na_check > 0], sep = "=", collapse = "; ")))
}

# ---- write merged phenotype ----
fwrite(dat,
       file.path(clust_dir, "tcga_paad_cluster_assignment_merged.tsv"),
       sep = "\t", quote = FALSE)

# ---- summary helper（修复：用 tapply 替代 aggregate 嵌套结构）----
make_summary <- function(df, cluster_col, event_col) {
  grp      <- df[[cluster_col]]
  evt      <- as.numeric(df[[event_col]])   # 确保数值型
  clusters <- sort(unique(grp))
  data.frame(
    cluster   = clusters,
    n_samples = as.integer(tapply(grp, grp, length)[as.character(clusters)]),
    n_events  = as.integer(tapply(evt, grp, sum, na.rm = TRUE)[as.character(clusters)]),
    stringsAsFactors = FALSE
  )
}

sum_k2 <- make_summary(dat, "cluster_k2", event_col); sum_k2$k <- "K2"
sum_k3 <- make_summary(dat, "cluster_k3", event_col); sum_k3$k <- "K3"
sum_k4 <- make_summary(dat, "cluster_k4", event_col); sum_k4$k <- "K4"

summary_all <- rbind(sum_k2, sum_k3, sum_k4)
summary_all <- summary_all[, c("k", "cluster", "n_samples", "n_events")]

fwrite(summary_all,
       file.path(clust_dir, "tcga_paad_cluster_k_comparison_summary.tsv"),
       sep = "\t", quote = FALSE)

# ---- survival screen helper ----
run_logrank <- function(df, cluster_col, time_col, event_col) {
  t_vec <- as.numeric(df[[time_col]])
  e_vec <- as.numeric(df[[event_col]])   # 修复：确保 event 为数值型
  grp   <- df[[cluster_col]]
  fit   <- survdiff(Surv(t_vec, e_vec) ~ grp)
  p     <- 1 - pchisq(fit$chisq, length(fit$n) - 1)
  data.frame(
    k       = cluster_col,
    chisq   = fit$chisq,
    df      = length(fit$n) - 1,
    p_value = p,
    stringsAsFactors = FALSE
  )
}

res_k2 <- run_logrank(dat, "cluster_k2", time_col, event_col)
res_k3 <- run_logrank(dat, "cluster_k3", time_col, event_col)
res_k4 <- run_logrank(dat, "cluster_k4", time_col, event_col)
km_stats <- rbind(res_k2, res_k3, res_k4)

fwrite(km_stats,
       file.path(clust_dir, "tcga_paad_k2_k3_k4_logrank_summary.tsv"),
       sep = "\t", quote = FALSE)

print(summary_all)
print(km_stats)
message("K comparison survival screen finished.")
