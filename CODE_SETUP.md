# Running the Analysis

## Requirements

- R and RStudio
- Internet access for package installation and queries
- The GSE43292 series matrix file

## Project setup

Download or clone this repository.

Create a `data` folder in the repository root and place `GSE43292_series_matrix.txt.gz` inside it. Keep the file compressed.

The scripts expect this structure:

```text
atherosclerosis-transcriptomics/
├── code/
│   ├── setup.R
│   ├── DEG_Analysis.R
│   ├── GSEA_Full_Analysis.R
│   └── ORA_DEG_Analysis.R
├── data/
│   └── GSE43292_series_matrix.txt.gz
├── figures/
├── tables/
├── CODE_SETUP.md
├── METHODS.md
└── README.md
```

Open RStudio with the repository root as the working directory.

## Dependencies

Install the required packages before running the analysis:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "GEOquery",
  "limma",
  "statmod",
  "hugene10sttranscriptcluster.db",
  "org.Hs.eg.db",
  "clusterProfiler",
  "enrichplot",
  "EnhancedVolcano"
))

install.packages(c(
  "pheatmap",
  "ggplot2",
  "png",
  "gridExtra"
))
```

## Run order

Run the scripts from the repository root in this order:

```r
source("code/DEG_Analysis.R")
source("code/GSEA_Full_Analysis.R")
source("code/ORA_DEG_Analysis.R")
```

| Script | Purpose |
|---|---|
| `DEG_Analysis.R` | Fits the paired limma model, annotates results, exports gene lists, and creates the volcano plot and heatmap |
| `GSEA_Full_Analysis.R` | Runs GO and KEGG GSEA using the full ranked gene list |
| `ORA_DEG_Analysis.R` | Runs GO and KEGG ORA using the cleaned significant gene list |
| `setup.R` | Checks the project folders and creates the output directory |

Each analysis script loads `setup.R` automatically.

## Outputs

The scripts write intermediate gene lists, enrichment tables, and exported plots to `results/`.

The volcano plot, heatmap, and ORA plots appear in the R graphics device. Export these separately to save image files.

Rerunning the scripts overwrites output files with matching names. The `tables/` and `figures/` folders contain the results shared in this repository.

## Data and execution details

Use the expression matrix in its original sample order. The paired model assigns patient pairs from consecutive sample columns.

The heatmap selects the 12 significant probe rows with the largest absolute log2 fold changes, then removes entries without gene symbols.

The enrichment scripts use installed annotation packages and online KEGG resources. Results may vary with package versions and database updates. The GSEA script does not set a random seed.

The R scripts do not perform the STRING and Cytoscape network analyses. See [METHODS.md](METHODS.md) for that workflow.
