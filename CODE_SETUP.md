# Running the original R analysis

## Prepare the project

Copy the code folder into the root of your repository. Create a data folder there and place GSE43292_series_matrix.txt.gz inside it. Use your supplied GEO matrix. GPL6244.soft.gz is not required by these scripts.

Open an R session with the repository root as its working directory. For example, create an RStudio project in that folder. Do not run from inside code/.

## Install dependencies once

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "statmod",
                      "hugene10sttranscriptcluster.db", "org.Hs.eg.db",
                      "clusterProfiler", "enrichplot", "EnhancedVolcano"))
install.packages(c("pheatmap", "ggplot2", "png", "gridExtra"))
```

## Run in order

```r
source("code/DEG_Analysis.R")
source("code/DEG_Analysis_2.R")
source("code/GSEA_Full_Analysis.R")
source("code/ORA_DEG_Analysis.R")
```

Each script creates results/ if needed. Intermediate tables and explicitly saved figures go there. Interactive volcano, heatmap, and ORA plots appear in the R graphics device; export them manually if needed. Running the scripts again overwrites their output files.

The heatmap script repeats the model fit before plotting. The original heatmap selection takes 12 probe rows and then removes missing symbols; its title says top 10. This preparation preserves that selection logic rather than changing the analysis.

## Edits to the supplied scripts

- Removed the personal working-directory setting and duplicate GEO download alternatives.
- Used one local input path: data/GSE43292_series_matrix.txt.gz.
- Standardised the full ranked table filename across the producer and GSEA script.
- Directed intermediate tables and saved plots to results/.
- Removed the package installation command from the analysis script.
- Removed row.names = 1 when ORA imports Cleaned_DEGs.txt, because the producing script writes that file without row names.

The preparation retains the paired design, thresholds, annotation rules, enrichment calls, heatmap selection, and original plotting labels. The scripts still assign pairs from sample order, so use the supplied matrix without reordering columns. Plot labels that say Normal refer to macroscopically intact tissue, not healthy donors.

The GSEA script requires live KEGG access and does not set a random seed. Annotation updates, package versions, and stochastic estimation may change results. This package does not freeze the original software environment.

## Validation

File-path and producer/consumer checks completed. The R workflow has not run in this environment because R is unavailable. These prepared scripts do not generate the later portfolio summary figures or perform the Cytoscape network analysis. The original Cytoscape workflow remains separate.
