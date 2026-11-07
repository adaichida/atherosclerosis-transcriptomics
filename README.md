# Functional Enrichment And Network-Based Re-Analysis of Differentially Expressed Genes in Atherosclerosis: A Bioinformatics Extension of Yan et al. (2022)

## Overview

This project reproduces and extends the atherosclerosis analysis reported by Yan et al. (2022), using the public GSE43292 microarray dataset.

I completed this work during my MSc course in Network and Pathway Analysis. The project combines paired differential expression analysis with complementary enrichment approaches and network analysis to investigate molecular processes associated with atherosclerotic plaque.

## Aim

The project aimed to compare differential expression and network findings with the original study, then extend the analysis through:

- Gene Set Enrichment Analysis (GSEA) of the full ranked gene list.
- Overrepresentation analysis (ORA) of significant differentially expressed genes.
- Hub gene ranking using both Degree and Maximal Clique Centrality (MCC).
- Functional enrichment of network-derived clusters and hub gene modules.

These analyses examine disease-associated processes at the gene, pathway, and network levels.

## Dataset

| Feature | Description |
|---|---|
| GEO accession | GSE43292 |
| Platform | GPL6244 |
| Patients | 32 hypertensive patients |
| Samples | 32 atheroma plaque and 32 paired macroscopically intact tissue samples |
| Tissue source | Carotid endarterectomy |
| Expression data | Supplied RMA-normalised microarray matrix |

The workflow applies no additional normalisation.

## Analysis workflow

### 1. Differential expression

GEOquery imports the expression data into R. limma fits a paired model accounting for patient identity and compares plaque against intact tissue.

Selection requires a Benjamini–Hochberg adjusted p-value below 0.01 and an absolute log2 fold change above 1.

Probe annotation and removal of missing or duplicate gene symbols produce the gene list for downstream analysis.

### 2. Functional enrichment

clusterProfiler performs two complementary analyses:

- GSEA evaluates GO terms and KEGG pathways across the full gene list ranked by log2 fold change.
- ORA evaluates GO terms and KEGG pathways within the cleaned significant gene list.

GO analysis includes Biological Process, Molecular Function, and Cellular Component.

### 3. Network construction and clustering

STRING, accessed through stringApp in Cytoscape, constructs a functional association network from genes with higher expression in plaque, using a confidence threshold of 0.4.

The MCODE algorithm in clusterMaker2 identifies densely connected network clusters.

### 4. Hub gene ranking and module enrichment

CytoHubba ranks the top ten nodes using Degree and MCC.

STRING functional enrichment examines the full network, the leading MCODE cluster, and the Degree and MCC modules using GO, KEGG, WikiPathways, and Reactome.

## Main findings

- Differential expression analysis identified 147 significant probe sets, compared with 134 DEGs reported in the original study. Annotation and cleaning yielded 114 unique genes.
- GSEA highlighted immune-associated processes, including adaptive immune response and leukocyte-mediated immunity.
- ORA highlighted complement activation, humoral immune response, and pathways including PPAR signalling and tryptophan metabolism.
- The leading MCODE cluster contained CHI3L1, CD163, MMP7, CXCL10, CD36, DPP4, ANPEP, and MME.
- Degree and MCC rankings shared nine genes: ITGAX, CCR1, IL1RN, CXCL10, CD163, MMP9, MMP12, SELE, and CD36.
- Enrichment of the MCC module highlighted the fluid shear stress and atherosclerosis pathway, including SELE, HMOX1, and MMP9.

The findings connect immune activity with lipid metabolism, extracellular matrix remodelling, and vascular inflammatory processes.

## Tools

R, GEOquery, limma, clusterProfiler, hugene10sttranscriptcluster.db, org.Hs.eg.db, EnhancedVolcano, pheatmap, enrichplot, STRING, and Cytoscape.

Network analysis uses stringApp, clusterMaker2, and CytoHubba.

## Repository contents

| Location | Contents |
|---|---|
| [METHODS.md](METHODS.md) | Detailed analysis methods |
| [tables/](tables/) | Enrichment results and network rankings |
| [figures/](figures/) | Analysis visualisations |
| [code/](code/) | R scripts for differential expression, GSEA, ORA, and visualisation |
| [CODE_SETUP.md](CODE_SETUP.md) | R Packages, input files, and instructions for running the scripts |

## Interpretation

This project examines associations within one observational microarray dataset. Macroscopically intact tissue does not represent healthy-donor tissue.

The network analysis focuses on genes with higher expression in plaque. The findings identify candidates for further investigation rather than experimentally validated biomarkers or therapeutic targets.

## Data source

[NCBI GEO: GSE43292](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43292)

## Author

Jane Chidimma Ezeanowai

MSc Systems Biology with a specialisation in Bioinformatics  
University of Skövde
