Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiI2MTk1MDE5OTBAcXEuY29tIiwiaWF0IjoxNzczOTEyNjM5LCJleHAiOjE3NzUxMjIyMzl9.DAEZ_VO7HtPEm27N6ErJrAeVD1D7yF_V6P8gXLRylRF_sm4KAeLT5fdDPylZCDaf17pjAlsH58NPVlu6BOEkgIGulfbLOBNt7OaVJ0-DcVMOH-hWZ8U80mH1DuU9RsRzGuuaFpTuLZlozAfynIgzLWz94K4Plv-TvRmX4NR_lKQmUMFVTx0lyAVD2JuHCd-OKXrs3wRCGOcOlPTbjol1cpgj7x-a5ntaqjOqVhce0xfCdyk49nYxaWndbwgJ0fFh-s8TbzKPq9uK1MVp4G0bWk8Zuto6NgRNjKA6XA2zP6FlddAZ7T8YFlqFXk0BTOih5AnGaDxV5WyDgZoC4_cUgQ")
# =========================
# 34b_check_opengwas_bbj_access_with_jwt.R
# 目的：
# 1. 用 JWT 检查 OpenGWAS API 是否可访问
# 2. 检查 bbj-a-140 是否可访问
# 3. 保存访问状态和元数据
# =========================

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(data.table)
})
setwd("/Users/wmz_mac/Desktop/胰腺癌")
root_dir <- "05_MR"
out_dir <- file.path(root_dir, "01_raw_data/outcome/OpenGWAS_BBJ")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

jwt_token <- Sys.getenv("OPENGWAS_JWT", unset = "")

if (!nzchar(jwt_token)) {
  stop("未检测到 OPENGWAS_JWT。请先在 R 中运行：Sys.setenv(OPENGWAS_JWT = '你的JWT')")
}

# ---------- 1. 检查 /user ----------
user_resp <- GET(
  "https://api.opengwas.io/api/user",
  add_headers(Authorization = paste("Bearer", jwt_token))
)

user_status <- status_code(user_resp)
user_text <- content(user_resp, as = "text", encoding = "UTF-8")

writeLines(user_text, file.path(out_dir, "34b_opengwas_user_response.json"))

# ---------- 2. 检查 bbj-a-140 ----------
bbj_resp <- GET(
  "https://api.opengwas.io/api/gwasinfo?id=bbj-a-140",
  add_headers(Authorization = paste("Bearer", jwt_token))
)

bbj_status <- status_code(bbj_resp)
bbj_text <- content(bbj_resp, as = "text", encoding = "UTF-8")

writeLines(bbj_text, file.path(out_dir, "34b_bbj_a_140_gwasinfo.json"))

# ---------- 3. 输出状态表 ----------
status_tab <- data.frame(
  check_item = c(
    "jwt_present",
    "user_endpoint_status_200",
    "bbj_a_140_access_status_200",
    "any_outcome_route_ready"
  ),
  status = c(
    TRUE,
    user_status == 200,
    bbj_status == 200,
    (user_status == 200) && (bbj_status == 200)
  ),
  stringsAsFactors = FALSE
)

fwrite(
  status_tab,
  file.path(root_dir, "00_registry/34b_bbj_outcome_access_status.tsv"),
  sep = "\t"
)

cat("OpenGWAS /user status:", user_status, "\n")
cat("BBJ bbj-a-140 status:", bbj_status, "\n")
cat("Saved files:\n")
cat("- 34b_opengwas_user_response.json\n")
cat("- 34b_bbj_a_140_gwasinfo.json\n")
cat("- 34b_bbj_outcome_access_status.tsv\n")