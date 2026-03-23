# =========================
# 34_start_real_mr_data_download_outcome_first.R
# 目的：
# 1. 优先接收并标准化挂载手动下载的 PanScan/PanC4 outcome 文件
# 2. 检查 OpenGWAS 的 BBJ outcome (bbj-a-140) 访问是否可用
# 3. 生成 outcome 数据接入状态表
# 4. 为后续 MR 正式运行建立真实 outcome 入口
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(httr)
  library(jsonlite)
  library(stringr)
})
setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "05_MR"
outcome_dir_primary <- file.path(root_dir, "01_raw_data/outcome/PanScan_PanC4")
outcome_dir_bbj <- file.path(root_dir, "01_raw_data/outcome/OpenGWAS_BBJ")
dir.create(outcome_dir_primary, recursive = TRUE, showWarnings = FALSE)
dir.create(outcome_dir_bbj, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# A. 你需要按实际情况修改这里
# =========================================================

# 方案 1：如果你已经手动下载了 PanScan/PanC4 文件，
# 把它放到某个本地路径，然后填到这里
manual_panscan_file <- ""  
# 例如：
# manual_panscan_file <- "/Users/你的用户名/Downloads/pancreatic_cancer_gwas.txt.gz"

# 方案 2：如果你要测试 OpenGWAS 的 BBJ 访问，
# 请把 JWT 放在环境变量里，或直接写到这里（更建议用环境变量）
jwt_token <- Sys.getenv("OPENGWAS_JWT", unset = "")

# 目标标准文件名（主 outcome）
primary_outcome_target <- file.path(outcome_dir_primary, "pancreatic_cancer_sumstats.tsv")

# BBJ 元信息输出文件
bbj_meta_file <- file.path(outcome_dir_bbj, "34_bbj_a_140_gwasinfo.json")
bbj_check_file <- file.path(outcome_dir_bbj, "34_bbj_a_140_access_check.tsv")

# =========================================================
# B. 工具函数
# =========================================================

safe_read_head <- function(path, n = 5) {
  x <- tryCatch(fread(path, nrows = n), error = function(e) NULL)
  x
}

normalize_to_tsv <- function(src, dst) {
  # 支持 .tsv / .txt / .csv / .gz
  message("Reading manual outcome file: ", src)
  
  # 尝试自动读取
  x <- tryCatch(
    fread(src),
    error = function(e) {
      message("fread failed on auto mode, retry with sep='\\t' ...")
      tryCatch(fread(src, sep = "\t"), error = function(e2) NULL)
    }
  )
  
  if (is.null(x)) {
    stop("无法读取手动 outcome 文件，请检查格式：", src)
  }
  
  fwrite(x, dst, sep = "\t")
  message("Standardized outcome file written to: ", dst)
  invisible(x)
}

check_required_outcome_cols <- function(df) {
  required <- c("SNP", "effect_allele", "other_allele", "beta", "se", "pval", "eaf", "outcome")
  missing <- setdiff(required, colnames(df))
  list(
    ok = length(missing) == 0,
    missing = missing
  )
}

# =========================================================
# C. Part 1: 优先挂载手动下载的 PanScan/PanC4 outcome 文件
# =========================================================

manual_attach_status <- data.frame(
  step = "manual_panscan_attach",
  file_provided = nzchar(manual_panscan_file),
  source_exists = FALSE,
  copied = FALSE,
  structurally_valid = FALSE,
  stringsAsFactors = FALSE
)

if (nzchar(manual_panscan_file)) {
  manual_attach_status$source_exists <- file.exists(manual_panscan_file)
  
  if (!file.exists(manual_panscan_file)) {
    warning("手动指定的 PanScan/PanC4 文件不存在：", manual_panscan_file)
  } else {
    # 标准化为项目内固定 outcome 文件
    dat <- normalize_to_tsv(manual_panscan_file, primary_outcome_target)
    manual_attach_status$copied <- file.exists(primary_outcome_target)
    
    # 检查结构
    structure_check <- check_required_outcome_cols(dat)
    manual_attach_status$structurally_valid <- structure_check$ok
    
    fwrite(
      data.frame(
        check_item = c("primary_outcome_file_exists", "required_columns_ok", "n_rows", "n_cols", "missing_columns"),
        value = c(
          file.exists(primary_outcome_target),
          structure_check$ok,
          nrow(dat),
          ncol(dat),
          paste(structure_check$missing, collapse = ",")
        ),
        stringsAsFactors = FALSE
      ),
      file.path(outcome_dir_primary, "34_primary_outcome_attachment_check.tsv"),
      sep = "\t"
    )
    
    # 保存前几行预览
    head_df <- dat[1:min(10, nrow(dat)), ]
    fwrite(
      head_df,
      file.path(outcome_dir_primary, "34_primary_outcome_head_preview.tsv"),
      sep = "\t"
    )
  }
}

# =========================================================
# D. Part 2: 检查 OpenGWAS 的 BBJ outcome (bbj-a-140) 是否可访问
# =========================================================

bbj_status <- data.frame(
  step = "bbj_opengwas_access_check",
  jwt_provided = nzchar(jwt_token),
  http_status = NA_integer_,
  access_ok = FALSE,
  saved_metadata = FALSE,
  stringsAsFactors = FALSE
)

# 这里只做“可访问性检查”和“元数据保存”
# 不在这里尝试全量下载整个 outcome summary stats
# 因为 OpenGWAS API 大多数端点需要 JWT，且全量下载策略常依赖具体端点/额度
if (nzchar(jwt_token)) {
  # gwasinfo 接口检查
  url <- "https://api.opengwas.io/api/gwasinfo?id=bbj-a-140"
  
  resp <- tryCatch(
    GET(
      url,
      add_headers(Authorization = paste("Bearer", jwt_token))
    ),
    error = function(e) NULL
  )
  
  if (!is.null(resp)) {
    bbj_status$http_status <- status_code(resp)
    bbj_status$access_ok <- status_code(resp) == 200
    
    txt <- content(resp, as = "text", encoding = "UTF-8")
    writeLines(txt, bbj_meta_file)
    bbj_status$saved_metadata <- file.exists(bbj_meta_file)
  }
}

fwrite(
  bbj_status,
  bbj_check_file,
  sep = "\t"
)

# =========================================================
# E. 总状态汇总
# =========================================================

primary_outcome_exists <- file.exists(primary_outcome_target)

# 如果 primary 已就绪，优先使用 primary
# 否则，如果 BBJ 可访问，说明至少 outcome 方向可启动
overall_status <- data.frame(
  item = c(
    "primary_outcome_local_ready",
    "bbj_access_ready",
    "any_outcome_route_ready",
    "recommended_next_action"
  ),
  value = c(
    primary_outcome_exists,
    bbj_status$access_ok[1],
    primary_outcome_exists || bbj_status$access_ok[1],
    ifelse(
      primary_outcome_exists,
      "Proceed to exposure data attachment (preferred primary outcome already ready)",
      ifelse(
        bbj_status$access_ok[1],
        "Proceed to exposure data attachment and keep BBJ as runnable fallback outcome",
        "Still need at least one real outcome route; attach PanScan/PanC4 file or fix OpenGWAS JWT"
      )
    )
  ),
  stringsAsFactors = FALSE
)

fwrite(
  manual_attach_status,
  file.path(root_dir, "00_registry/34_manual_panscan_attachment_status.tsv"),
  sep = "\t"
)

fwrite(
  overall_status,
  file.path(root_dir, "00_registry/34_outcome_route_status.tsv"),
  sep = "\t"
)

notes <- c(
  "34 outcome-first real data attachment notes",
  paste0("Created at: ", Sys.time()),
  "",
  "[What this script does]",
  "1. Attach manual PanScan/PanC4 file if provided",
  "2. Check OpenGWAS BBJ access using JWT",
  "3. Decide whether any outcome route is ready",
  "",
  "[Preferred route]",
  "Primary outcome = local PanScan/PanC4 file",
  "",
  "[Fallback route]",
  "BBJ (bbj-a-140) if OpenGWAS access works",
  "",
  "[Next step]",
  "Once any outcome route is ready, move to exposure data attachment."
)

writeLines(
  notes,
  file.path(root_dir, "06_logs/34_outcome_real_data_notes.txt")
)

cat("34_start_real_mr_data_download_outcome_first.R finished successfully.\n")
cat("Primary outcome local ready:", primary_outcome_exists, "\n")
cat("BBJ access ready:", bbj_status$access_ok[1], "\n")
cat("Any outcome route ready:", primary_outcome_exists || bbj_status$access_ok[1], "\n")