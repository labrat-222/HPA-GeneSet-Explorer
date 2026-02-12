############################################################
## GWAS → risk gene list → HPA query + data download
############################################################

## -------------------- Setup -------------------- ##

# install.packages("gwasrapidd")   # if needed
# BiocManager::install("hpar")     # if needed
library(dplyr)
library(tidyr)
library(gwasrapidd)
library(httr)
library(jsonlite)
library(ggplot2)

#set work directory to where the file is located
getwd()
source("functions.R")

############################################################
## PART A – GWAS search: get gene list and basic filtering
############################################################

# 1. Search traits (to inspect relevant EFO IDs if needed)
#all_traits <- get_traits() #see all traits and EFO ID available

#type in key words to search relevant traits (optional)
dplyr::filter(all_traits@traits, grepl("Alzheimer", trait, ignore.case = TRUE))

# 2. Search GWAS for selected traits / EFO IDs
efo_ids <-c("EFO_0006514","EFO_0006801","MONDO_0004975","EFO_1001870","EFO_0009268","OBA_2001000","EFO_0022957")
result_gwas <- search_gwasrapidd(efo_id = efo_ids)
nrow(result_gwas)

# Get annotation (EFO ID, trait, URI) from the traits table, for saving as metadata
efo_meta <- all_traits@traits %>%
  dplyr::filter(efo_id %in% efo_ids) %>%
  dplyr::distinct(efo_id, trait, uri)

# 3. Filter by p-value and remove NAs
result_gwas <- result_gwas %>%
  filter(pvalue < 1e-5) %>%
  tidyr::drop_na(pvalue)

nrow(result_gwas)
length(unique(result_gwas$ensembl_id))  # check for duplicated genes

# 4. Keep one row per gene (lowest p-value)
result_gwas_unique <- result_gwas %>%
  group_by(ensembl_id) %>%
  slice_min(order_by = pvalue, n = 1, with_ties = FALSE) %>%
  ungroup()
nrow(result_gwas_unique)

# 5. Save gene list with search term as metadata (optional, for record)
attr(result_gwas_unique, "search_efo_meta") <- efo_meta
attr(result_gwas_unique, "search_date")   <- Sys.time()
saveRDS(result_gwas_unique, "data/gwas_list.RDS")

############################################################
## PART B – Use risk gene list to query HPA and fetch data
############################################################

# Use the GWAS gene list directly (or import your own gene list)
gene_list <- result_gwas_unique
head(gene_list)
dim(gene_list)

########## Database search in HPA ##########
#choose which version of HPA to search

#path <- "https://www.proteinatlas.org/api/search_download.php"
path<-'https://v24.proteinatlas.org/api/search_download.php'
# Unique Ensembl IDs
search_queries <- unique(gene_list$ensembl_id)

all_data <- list()

for (search_query in search_queries) {
  req <- GET(
    url   = path,
    query = list(
      search   = search_query,
      format   = "json",
      columns  = paste(c(
        "g", "eg", "gd", "pc", "upbp", "up_mf", "di",
        "rnats", "rnatd", "rnabrs", "rnabrd",
        "ecbrain", "ectissue", "rnascs", "rnascd", "ecsinglecell",
        "rtcte", "rnasnbs", "rnasnbd",
        "brain_RNA__tau", "Brain_sn_RNA__tau", "t_RNA__tau",
        "sc_RNA__tau",
        "rnatss", "rnascss", "rnasnbss", "rnacass", "rnabrss"
      ), collapse = ","),
      compress = "no"
    )
  )
  
  if (req$status_code == 200) {
    resp <- content(req, as = "text", encoding = "UTF-8")
    df   <- fromJSON(resp, flatten = TRUE) %>% data.frame()
    all_data[[search_query]] <- df
  } else {
    message("Failed to retrieve data for: ", search_query)
  }
}

# Combine all data frames
search_result <- bind_rows(all_data)
head(search_result)
dim(search_result)

# Filter again to be safe
filtered_result <- search_result %>%
  filter(Ensembl %in% search_queries) %>%
  distinct()

head(filtered_result)
dim(filtered_result)

# Save final HPA data
# Attach metadata to the HPA result
attr(filtered_result, "search_efo_meta") <- efo_meta
attr(filtered_result, "search_date")     <- Sys.time()
attr(filtered_result, "HPA_version") <- path

# Save final HPA data with metadata included
saveRDS(filtered_result, "data/HPA_data.RDS")

# write.table(filtered_result, "data/HPA_data.txt",
#             sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

################ END ################
