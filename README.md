# HPA GeneSet Explorer

A reproducible R pipeline that maps GWAS-derived risk genes onto the [Human Protein Atlas (HPA)](https://www.proteinatlas.org/) to identify which tissues, cell types, and brain regions show significant enrichment for disease-associated gene expression.

Originally developed for Alzheimer's disease (AD), Lewy body dementia (DLB), and frontotemporal dementia (FTD), the framework is **fully generalisable** to any trait with a GWAS gene list.

---

## What this tool does

Starting from trait EFO identifiers, the pipeline:

1. Retrieves risk gene lists from the **GWAS Catalog** via `gwasrapidd`
2. Queries the **HPA API** (v24 or v25) for expression profiles across tissues, single-cell types, and brain regions
3. Maps risk genes onto HPA's UMAP-based expression clusters
4. Performs **over-representation analysis** using Fisher's exact test with Benjamini-Hochberg correction
5. Validates enrichment independently via **Monte Carlo simulation** (multinomial resampling, 1,000,000 iterations)
6. Renders a fully automated, versioned **PDF report** for each disease

The dual enrichment strategy — parametric (Fisher) plus simulation-based (Monte Carlo) — provides complementary statistical evidence and strengthens biological interpretation.

---

## Key features

- **Multi-level transcriptomic mapping**: tissue, single-cell, and brain region contexts in one run
- **Two complementary enrichment methods**: Fisher's exact test and Monte Carlo simulation
- **Automatic HPA version detection**: switching between v24 and v25 requires no code changes — version is read from file metadata
- **Disease-agnostic design**: add any new disease in two lines; the report template adapts automatically
- **Fully reproducible**: parameterised Rmd, versioned data directories, search dates recorded in file metadata
- **Graceful degradation**: sections that depend on optional classification files are automatically skipped with an informative note when files are absent

---

## Repository structure

```
HPA-GeneSet-Explorer/
├── Generalized report.Rmd    # Parameterised report template
├── Functions.R                # All analysis and plotting functions
├── run_hpa_search_v24.R       # HPA API batch search — v24
├── run_hpa_search_v25.R       # HPA API batch search — v25
├── render_reports.R           # Batch PDF renderer — AD, DLB, FTD (v24)
├── render_reports_v25.R       # Batch PDF renderer — AD, DLB, FTD (v25)
├── search_2_in_1.r            # Combined GWAS + HPA search for a single disease
├── generate_umap_cluster.R    # Utility: build cluster dictionaries from HPA data
└── data/
    ├── v24/                   # HPA v24 UMAP coordinates and cluster dictionaries
    ├── v25/                   # HPA v25 UMAP coordinates and cluster dictionaries
    └── search_result/         # GWAS gene lists and HPA search results (.RDS)
```

---

## Quick start

### Prerequisites

R >= 4.0 with the following packages:

```r
install.packages(c(
  "dplyr", "ggplot2", "rmarkdown", "httr", "jsonlite",
  "tidyr", "ggvenn", "VennDiagram", "ggVennDiagram",
  "viridis", "data.table", "gwasrapidd"
))
```

### Step 1 — Retrieve GWAS risk genes and HPA expression data

To run all three dementia diseases at once (recommended):

```bash
Rscript run_hpa_search_v24.R   # uses HPA v24
# or
Rscript run_hpa_search_v25.R   # uses HPA v25
```

This saves to `data/search_result/`:
- `gwas_list_<DISEASE>.RDS` — risk gene list with GWAS metadata
- `HPA_data_<DISEASE>_<VERSION>_<DATE>.RDS` — HPA expression data

For a single disease, open `search_2_in_1.r`, set your EFO IDs, and run.

### Step 2 — Generate the report

**Option A — Single disease in RStudio:**
Open `Generalized report.Rmd`, update the two file paths in the data-loading chunk, then click **Knit**.

**Option B — Batch render all three diseases:**

```r
source("render_reports.R")     # outputs PDFs for AD, DLB, FTD (v24)
source("render_reports_v25.R") # same for v25
```

Output PDFs are saved to the directory specified in the render script.

---

## Report contents

Each PDF report contains:

| Section | Content |
|---|---|
| **0. Search parameters** | Trait table with EFO IDs; GWAS and HPA search dates |
| **1. Summary statistics** | Gene counts: GWAS list → found in HPA → mappable to UMAP |
| **2a. Fisher enrichment** | Lollipop plots (raw p and BH-adjusted) for tissue, single-cell, and brain clusters |
| **2b. UMAP visualisation** | Enriched clusters highlighted on the full expression UMAP |
| **2c. Monte Carlo enrichment** | Z-score lollipop plots with empirical p-values; significant clusters highlighted in disease colour |
| **2d. Volcano plots** | Fold enrichment vs. −log(BH-adjusted p) with cluster name annotations |
| **3a. Specificity profiles** | Pie charts of tissue-specificity categories (HPA search data always shown; classification-file data shown when available) |
| **3b. Tissue distribution** | Bar plots of tissue-enriched genes per expression context |
| **3c. Venn diagrams** | Overlap of enriched gene sets across tissue, single-cell, and brain contexts |

---

## Methods

**Fisher's exact test**  
For each UMAP expression cluster, a one-sided 2×2 contingency table tests whether risk genes are over-represented relative to all genes in the HPA. Fold enrichment = observed overlap / expected overlap, where expected = (risk genes × cluster size) / total genes. P-values are corrected with the Benjamini-Hochberg procedure (FDR threshold: 0.1).

**Monte Carlo simulation**  
Risk genes are randomly redistributed across clusters 1,000,000 times, weighted proportionally to cluster size (multinomial sampling). A Z-score is computed as (observed − mean) / SD of the null distribution. Empirical p-values are derived from the proportion of simulations that meet or exceed the observed count. Clusters with empirical p < 0.05 and a positive Z-score are considered enriched.

---

## HPA version compatibility

| Version | Tissue clusters | Single-cell clusters | Brain clusters | Classification files |
|---|---|---|---|---|
| v24 | 83 | 80 | 56 | Included (`data/v24/`) |
| v25 | 83 | 111 | 57 | Not yet included |

The pipeline auto-detects the version from HPA search result metadata. Report sections that depend on classification files (3a bottom row, 3b) are skipped automatically when the files are absent, with an inline note.

---

## Adding a new disease

1. In `run_hpa_search_v24.R` (and/or v25), add an entry to the `diseases` list:

```r
diseases <- list(
  AD  = "data/search_result/gwas_list_AD.RDS",
  DLB = "data/search_result/gwas_list_DLB.RDS",
  FTD = "data/search_result/gwas_list_FTD.RDS",
  PD  = "data/search_result/gwas_list_PD.RDS"   # new disease
)
```

2. Add the same entry to `render_reports.R`. The report template requires no changes.

---

## Background

This pipeline was developed as part of a systems-level investigation of dementia genetics, integrating GWAS risk loci with transcriptomic data from the Human Protein Atlas to characterise the cell-type and tissue specificity of disease-associated genes. The framework prioritises reproducibility: all data are versioned, search provenance is stored in file metadata, and reports are generated from a single parameterised template without manual intervention.
