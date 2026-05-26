############################################################
## Generate UMAP_and_cluster.{version}.txt from individual UMAP CSV files
## Run this once when adding a new HPA version to the project
############################################################

library(dplyr)

# Set the version to generate
HPA_VERSION <- "v25"
hpa_dir <- paste0("data/", HPA_VERSION, "/")

# Load the three UMAP files
brain   <- read.csv(paste0(hpa_dir, "umap_brain_data.",      HPA_VERSION, ".csv"))
sc      <- read.csv(paste0(hpa_dir, "umap_singlecell_data.", HPA_VERSION, ".csv"))
tissue  <- read.csv(paste0(hpa_dir, "umap_tissue_data.",     HPA_VERSION, ".csv"))

# Use brain as the base gene set (matches v24 behaviour)
# Left-join SC and tissue cluster numbers onto brain genes
umap_combined <- brain %>%
  select(
    gene             = ensembl_gene_id,
    UMAP_1_scaled    = x,
    UMAP_2_scaled    = y,
    `UMAP Cluster brain` = cluster_no
  ) %>%
  mutate(
    UMAP_1 = UMAP_1_scaled,   # raw coords not used by pipeline; mirror scaled
    UMAP_2 = UMAP_2_scaled
  ) %>%
  left_join(
    sc %>% select(ensembl_gene_id, `UMAP Cluster singlecell` = cluster_no),
    by = c("gene" = "ensembl_gene_id")
  ) %>%
  left_join(
    tissue %>% select(ensembl_gene_id, `UMAP Cluster tissue` = cluster_no),
    by = c("gene" = "ensembl_gene_id")
  ) %>%
  mutate(
    `UMAP Cluster blood`   = NA_integer_,  # not used by pipeline
    `UMAP Cluster celline` = NA_integer_   # not used by pipeline
  ) %>%
  select(
    gene, UMAP_1, UMAP_2, UMAP_1_scaled, UMAP_2_scaled,
    `UMAP Cluster blood`, `UMAP Cluster brain`,
    `UMAP Cluster singlecell`, `UMAP Cluster tissue`,
    `UMAP Cluster celline`
  )

out_path <- paste0(hpa_dir, "UMAP_and_cluster.", HPA_VERSION, ".txt")
write.table(umap_combined, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

############################################################
## Verification — compare structure against v24.1
############################################################
ref <- read.delim("data/v24.1/UMAP_and_cluster.v24.1.txt")

cat("\n===== VERIFICATION =====\n")
cat("Column names match v24.1:  ", identical(colnames(umap_combined), colnames(ref)), "\n")
cat("v24.1 gene count:          ", nrow(ref), "\n")
cat("v25  gene count:           ", nrow(umap_combined), "\n")
cat("Duplicate gene IDs in v25: ", sum(duplicated(umap_combined$gene)), "\n")
cat("NA in UMAP Cluster brain:  ", sum(is.na(umap_combined$`UMAP Cluster brain`)), "\n")
cat("NA in UMAP Cluster SC:     ", sum(is.na(umap_combined$`UMAP Cluster singlecell`)), "\n")
cat("NA in UMAP Cluster tissue: ", sum(is.na(umap_combined$`UMAP Cluster tissue`)), "\n")
cat("Brain cluster range:       ", min(umap_combined$`UMAP Cluster brain`),
    "-", max(umap_combined$`UMAP Cluster brain`), "\n")
cat("SC cluster range:          ", min(umap_combined$`UMAP Cluster singlecell`, na.rm=TRUE),
    "-", max(umap_combined$`UMAP Cluster singlecell`, na.rm=TRUE), "\n")
cat("Tissue cluster range:      ", min(umap_combined$`UMAP Cluster tissue`, na.rm=TRUE),
    "-", max(umap_combined$`UMAP Cluster tissue`, na.rm=TRUE), "\n")
cat("\nFile saved to:", out_path, "\n")
