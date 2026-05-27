############################################################
## Pipeline test: Part B only (HPA search) for AD, DLB, FTD
## Reuses existing gwas_list_*.RDS; saves HPA_data_{disease}_v24_{date}.RDS
############################################################

library(dplyr)
library(httr)
library(jsonlite)
library(tidyr)

HPA_VERSION <- "v25"
date_str    <- format(Sys.Date(), "%Y%m%d")
path        <- paste0("https://", HPA_VERSION, ".proteinatlas.org/api/search_download.php")

hpa_columns <- paste(c(
  "g", "eg", "gd", "pc", "upbp", "up_mf", "di",
  "rnats", "rnatd", "rnabrs", "rnabrd",
  "ecbrain", "ectissue", "rnascs", "rnascd", "ecsinglecell",
  "rtcte", "rnasnbs", "rnasnbd",
  "brain_RNA__tau", "Brain_sn_RNA__tau", "t_RNA__tau",
  "sc_RNA__tau",
  "rnatss", "rnascss", "rnasnbss", "rnacass", "rnabrss"
), collapse = ",")

diseases <- list(
  AD  = "data/search_result/gwas_list_AD.RDS",
  DLB = "data/search_result/gwas_list_DLB.RDS",
  FTD = "data/search_result/gwas_list_FTD.RDS"
)

for (disease_name in names(diseases)) {
  cat("\n============================\n")
  cat("Disease:", disease_name, "\n")
  cat("============================\n")

  gene_list      <- readRDS(diseases[[disease_name]])
  efo_meta       <- attr(gene_list, "search_efo_meta")
  search_queries <- unique(gene_list$ensembl_id)
  cat("Unique genes to query:", length(search_queries), "\n")

  all_data <- list()
  failed   <- character(0)

  for (i in seq_along(search_queries)) {
    sq <- search_queries[i]
    if (i %% 50 == 0) cat("  Progress:", i, "/", length(search_queries), "\n")

    req <- GET(
      url   = path,
      query = list(search = sq, format = "json",
                   columns = hpa_columns, compress = "no")
    )

    if (req$status_code == 200) {
      resp <- content(req, as = "text", encoding = "UTF-8")
      df   <- tryCatch(fromJSON(resp, flatten = TRUE) %>% data.frame(),
                       error = function(e) NULL)
      if (!is.null(df) && nrow(df) > 0) all_data[[sq]] <- df
    } else {
      failed <- c(failed, sq)
    }
  }

  if (length(failed) > 0)
    cat("  Failed queries (", length(failed), "):", paste(head(failed, 5), collapse=", "), "\n")

  search_result   <- bind_rows(all_data)
  filtered_result <- search_result %>%
    filter(Ensembl %in% search_queries) %>%
    distinct()

  attr(filtered_result, "search_efo_meta") <- efo_meta
  attr(filtered_result, "search_date")     <- Sys.time()
  attr(filtered_result, "HPA_version")     <- path

  out_file <- paste0("data/search_result/HPA_data_", disease_name, "_", HPA_VERSION, "_", date_str, ".RDS")
  saveRDS(filtered_result, out_file)
  cat("  Saved:", out_file, "\n")
  cat("  Rows in result:", nrow(filtered_result), "\n")
}

cat("\n============================\n")
cat("All done!\n")
cat("============================\n")
