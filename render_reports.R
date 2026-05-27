############################################################
## Render Generalized report.Rmd for AD, DLB, FTD (v24)
## Saves PDF to /Users/labrat/Downloads/new reports/
############################################################

library(rmarkdown)

date_str   <- format(Sys.Date(), "%Y%m%d")
output_dir <- "/Users/labrat/Downloads/new reports"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

rmd_original <- "Generalized report.Rmd"
rmd_content  <- readLines(rmd_original, warn = FALSE)

diseases <- list(
  AD  = list(gwas = "gwas_list_AD.RDS",
             hpa  = paste0("HPA_data_AD_v24_",  date_str, ".RDS")),
  DLB = list(gwas = "gwas_list_DLB.RDS",
             hpa  = paste0("HPA_data_DLB_v24_", date_str, ".RDS")),
  FTD = list(gwas = "gwas_list_FTD.RDS",
             hpa  = paste0("HPA_data_FTD_v24_", date_str, ".RDS"))
)

for (disease_name in names(diseases)) {
  cat("\n============================\n")
  cat("Rendering:", disease_name, "\n")
  cat("============================\n")

  # Substitute file names only — no other changes to the Rmd
  tmp_content <- rmd_content
  tmp_content <- gsub("gwas_list_FTD\\.RDS",
                      diseases[[disease_name]]$gwas, tmp_content)
  tmp_content <- gsub("HPA_data_FTD\\.RDS",
                      diseases[[disease_name]]$hpa,  tmp_content)

  tmp_rmd <- paste0("_render_", disease_name, ".Rmd")
  writeLines(tmp_content, tmp_rmd)

  out_pdf <- file.path(output_dir,
                       paste0("report_", disease_name, "_v24_", date_str, ".pdf"))

  tryCatch({
    render(tmp_rmd,
           output_format = "pdf_document",
           output_file   = out_pdf,
           quiet         = FALSE)
    cat("Saved:", out_pdf, "\n")
  }, error = function(e) {
    cat("ERROR rendering", disease_name, ":", conditionMessage(e), "\n")
  })

  file.remove(tmp_rmd)
}

cat("\n============================\n")
cat("All reports done!\n")
cat("============================\n")
