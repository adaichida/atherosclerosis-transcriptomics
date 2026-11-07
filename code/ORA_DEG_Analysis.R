# Original coursework analysis with file-handling edits for this repository.
# Run from the project root; outputs go to results/.
source("code/setup.R")

# DEG OVER-REPRESENTATION ANALYSIS (ORA)

# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)


# Load cleaned DEG list
deg <- read.table("results/Cleaned_DEGs.txt",
                  header = TRUE, 
                  as.is = TRUE,
                  sep = "\t")


# Extract gene symbols
genes <- deg$symbol


# Convert gene symbols to ENTREZ IDs
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = genes,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")


# Remove missing ENTREZ IDs
entrez_ids <- na.omit(entrez_ids)


# GO OVER-REPRESENTATION ANALYSIS

ego <- enrichGO(gene = entrez_ids,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)


# View top GO results
head(ego, 50)


# Save GO results
write.csv(as.data.frame(ego),
          "results/GO_ORA_Results.csv",
          row.names = FALSE)


# GO dotplot
dotplot(ego,
        showCategory = 10,
        title = "GO ORA Enrichment")


# GO barplot
barplot(ego,
        showCategory = 10,
        title = "GO ORA Pathways")



# KEGG OVER-REPRESENTATION ANALYSIS

ekegg <- enrichKEGG(gene = entrez_ids,
                    organism = "hsa",
                    keyType = "kegg",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH")


# View top KEGG results
head(ekegg,50)


# Save KEGG results
write.csv(as.data.frame(ekegg),
          "results/KEGG_ORA_Results.csv",
          row.names = FALSE)


# KEGG dotplot
dotplot(ekegg,
        showCategory = 10,
        title = "KEGG ORA Enrichment")


# KEGG barplot
barplot(ekegg,
        showCategory = 10,
        title = "KEGG ORA Pathways")
