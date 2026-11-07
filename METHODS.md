## Dataset

The analysis uses the GSE43292 microarray dataset from NCBI GEO, generated on the GPL6244 platform.

The dataset contains paired atheroma plaque and macroscopically intact carotid tissue samples from 32 hypertensive patients, giving 64 samples in total.

GEOquery imports the supplied RMA-normalised expression matrix into R. The workflow applies no additional normalisation.

## Differential expression analysis

The analysis uses limma with patient pairing included in the design:

`~ 0 + Group + Pairs`

The contrast compares atheroma plaque against intact tissue:

`Atheroma - Intact`

The workflow fits the model with `lmFit`, applies the contrast with `contrasts.fit`, and moderates the statistics with `eBayes(robust = TRUE)`.

Differential expression selection requires:

- Benjamini–Hochberg adjusted p-value < 0.01
- Absolute log2 fold change > 1

Positive log2 fold changes indicate higher expression in plaque. Negative values indicate lower expression in plaque.

## Gene annotation and filtering

The workflow maps probe identifiers to gene symbols using `hugene10sttranscriptcluster.db` and retains the first annotation when a probe has multiple symbol mappings.

For downstream analysis, it removes entries with missing or empty symbols and retains the first occurrence of each duplicate symbol.

For the full GSEA input, it resolves duplicate symbols separately by retaining the probe with the highest absolute log2 fold change.

## Gene set enrichment analysis

GSEA uses the full annotated gene list ranked by log2 fold change, without a differential expression significance filter.

The analysis uses clusterProfiler:

- `gseGO` tests Biological Process, Molecular Function, and Cellular Component terms with `ont = "ALL"` and gene symbols as identifiers.
- `gseKEGG` tests human KEGG pathways after symbol-to-Entrez mapping with `org.Hs.eg.db`.

Both analyses use:

| Parameter | Value |
|---|---|
| Minimum gene set size | 10 |
| Maximum gene set size | 500 |
| eps | 0 |
| pvalueCutoff | 0.05 |
| Multiple testing correction | Benjamini–Hochberg |

The KEGG analysis retrieves pathway annotations with `use_internal_data = FALSE`.

The scripts do not set a random seed.

## Overrepresentation analysis

ORA uses the cleaned list of significant differentially expressed genes, combining genes with higher and lower expression in plaque.

The workflow converts symbols to Entrez IDs with `org.Hs.eg.db` and removes missing mappings.

The analysis uses:

- `enrichGO` across all three GO ontologies, with Benjamini–Hochberg correction, `pvalueCutoff = 0.05`, and `qvalueCutoff = 0.05`.
- `enrichKEGG` for human pathways, with Benjamini–Hochberg correction and `pvalueCutoff = 0.05`.

The scripts do not supply a custom background universe. Each function therefore uses its default annotation background.

## STRING network construction

The network analysis submits genes with higher expression in plaque to STRING through the stringApp plugin in Cytoscape 3.10.4.

The STRING protein query uses a confidence threshold of 0.4.

## Network clustering and ranking

Network clustering uses the MCODE algorithm through clusterMaker2.

The report records the following settings:

| Setting | Value |
|---|---|
| Cluster score criterion | > 5 |
| Degree cutoff | 2 |
| Node score cutoff | 0.2 |
| k-core | 2 |
| Maximum depth | 100 |

CytoHubba ranks network nodes using Degree and Maximal Clique Centrality (MCC). The analysis exports the ten highest-ranking nodes from each method.

## Network functional enrichment

STRING functional enrichment within Cytoscape examines:

- The network of genes with higher expression in plaque
- The highest-scoring MCODE cluster
- The top-ten MCC module
- The top-ten Degree module

The analysis queries GO, KEGG, WikiPathways, and Reactome to support interpretation of the network modules.

## Visualisation

EnhancedVolcano displays log2 fold changes and adjusted p-values.

pheatmap displays selected genes with large absolute log2 fold changes using gene-wise z-scores. Samples follow their tissue groups.

enrichplot generates GSEA enrichment curves and ORA summaries. Cytoscape displays the interaction network, clusters, and node rankings.

## Interpretation

Macroscopically intact tissue does not represent healthy-donor tissue. Differences in cell composition may influence expression patterns.

Enrichment and network rankings identify associations and candidates for further study. They do not establish causality or validate clinical biomarkers.
