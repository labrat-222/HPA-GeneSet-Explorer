# HPA-GeneSet-Explorer

A scalable and user-friendly R framework integrating gene sets with the Human Protein Atlas (HPA) for multi-level transcriptomic mapping across tissues, brain regions, and single-cell types, enabling rapid identification of disease-relevant expression patterns.

## 🚀 Quick start

1.  Clone the repository.

``` bash
git clone https://github.com/labrat-222/HPA-GeneSet-Explorer.git
```

2.  Open `search_2_in_1.r`. Specify your gene list (or EFO IDs if using GWAS).Run the script to retrieve and process data from the Human Protein Atlas (HPA). This step will generate:

-   `data/search_result/gwas_list.RDS`
-   `dat/``search_result/``HPA_data.RDS`

⚠️ For convenience, it is recommended not to change the default file names or directory structure, as the report script expects these files.

3\. Open `Generalized report.Rmd` and click **"Knit"** on top of Rstudio.

4\. 🎉 That's it! A complete enrichment and visualization report will be automatically generated.

See *Generalized-report.pdf* for sample report.
