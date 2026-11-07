# Original coursework analysis with file-handling edits for this repository.
# Run from the project root; outputs go to results/.
source("code/setup.R")

# FULL DE GENE LIST GSEA ANALYSIS
# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)


# Load full gene list
full_gene_list <- read.table("results/Ranked_logFC_Full_Gene_List.txt", 
                             header = TRUE, 
                             as.is = TRUE,
                             row.names = 1,
                             sep = "\t")


# Check structure
head(full_gene_list)
dim(full_gene_list)


# Create ranked gene vector using logFC
gsea_input <- full_gene_list$logFC
names(gsea_input) <- full_gene_list$symbol

# Sort in descending order
gsea_input <- sort(gsea_input, decreasing = TRUE)

# Check
head(gsea_input)
tail(gsea_input)
length(gsea_input)



# Convert gene symbols to Entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = names(gsea_input),
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")


# Replace names with Entrez IDs
gsea_input_entrez <- gsea_input
names(gsea_input_entrez) <- entrez_ids


# Remove any NAs
gsea_input_entrez <- gsea_input_entrez[!is.na(names(gsea_input_entrez))]


# Sort again
gsea_input_entrez <- sort(gsea_input_entrez, decreasing = TRUE)


# Check
head(gsea_input_entrez)
length(gsea_input_entrez)



# using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).
# Run GSEA on KEGG pathways
# Run GSEA on KEGG pathways
gsea_kegg <- gseKEGG(geneList = gsea_input_entrez,
                     organism = "hsa",
                     keyType = "ncbi-geneid",
                     minGSSize = 10,
                     maxGSSize = 500,
                     eps = 0,
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     use_internal_data = FALSE,
                     verbose = TRUE)
# Check results
head(gsea_kegg)
dim(gsea_kegg)



# Run GSEA on ALL GO ontologies together
gsea_go <- gseGO(geneList = gsea_input,
                 OrgDb = org.Hs.eg.db,
                 ont = "ALL",
                 keyType = "SYMBOL",
                 minGSSize = 10,
                 maxGSSize = 500,
                 eps = 0,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 verbose = TRUE)

# Check results
dim(gsea_go)
head(gsea_go)


# Save KEGG results
write.csv(as.data.frame(gsea_kegg), 
          file = "results/GSEA_KEGG_results.csv", 
          row.names = FALSE)

# Save GO results
write.csv(as.data.frame(gsea_go), 
          file = "results/GSEA_GO_results.csv", 
          row.names = FALSE)


# Top 5 Enrichment plots for KEGG
for(i in 1:5){
  
  pdf(paste0("KEGG_Pathway_", i, ".pdf"),
      width = 8,
      height = 6)
  
  print(
    gseaplot2(gsea_kegg,
              geneSetID = i,
              title = gsea_kegg@result$Description[i])
  )
  
  dev.off()
  
}



# Top 5 Enrichment plots for GO 
for(i in 1:5){
  
  pdf(paste0("GO_Pathway_", i, ".pdf"),
      width = 8,
      height = 6)
  
  print(
    gseaplot2(gsea_go,
              geneSetID = i,
              title = gsea_go@result$Description[i])
  )
  
  dev.off()
  
}



# View Top five KEGG and GO terms
gsea_kegg@result[1:5, c("ID", "Description", "NES", "p.adjust")]
gsea_go@result[1:5, c("ID", "Description", "NES", "p.adjust")]






# OR





# COMBINE GO ANALYSIS PLOTS
# FOR-LOOP TO EXTRACT THE TOP 6 GO RESULTS
for(i in 1:6){
  
  png(paste0("GO_", i, ".png"),
      width = 2000,
      height = 1600,
      res = 300)
  
  print(
    gseaplot2(gsea_go,
              geneSetID = i,
              title = gsea_go@result$Description[i])
  )
  
  dev.off()
  
}



# FOR-LOOP TO EXTRACT THE TOP 6 KEGG RESULTS
for(i in 1:6){
  
  png(paste0("KEGG_", i, ".png"),
      width = 2000,
      height = 1600,
      res = 300)
  
  print(
    gseaplot2(gsea_kegg,
              geneSetID = i,
              title = gsea_kegg@result$Description[i])
  )
  
  dev.off()
  
}






# COMBINE PLOTS FOR GO
library(png)
library(grid)
library(gridExtra)

go_images <- lapply(1:6, function(i){
  
  img <- readPNG(paste0("GO_", i, ".png"))
  rasterGrob(img, interpolate = TRUE)
  
})

pdf("results/Combined_GO_GSEA_Plots.pdf",
    width = 14,
    height = 16)

grid.arrange(grobs = go_images,
             ncol = 2)

dev.off()




# COMBINE PLOTS FOR KEGG
kegg_images <- lapply(1:6, function(i){
  
  img <- readPNG(paste0("KEGG_", i, ".png"))
  rasterGrob(img, interpolate = TRUE)
  
})

pdf("results/Combined_KEGG_GSEA_Plots.pdf",
    width = 14,
    height = 16)

grid.arrange(grobs = kegg_images,
             ncol = 2)

dev.off()











