# Original coursework analysis with file-handling edits for this repository.
# Run from the project root; outputs go to results/.
source("code/setup.R")

# Run this script from the project root. See CODE_SETUP.md.
getwd()

list.files()

# IMPORT THE SUPPLIED RMA-NORMALISED GEO MATRIX
library(GEOquery)
eset <- getGEO(filename = "data/GSE43292_series_matrix.txt.gz")

# Check the dimension of the normalized gene expression matrix
dim(exprs(eset))


# Save the expression metadata in R as a variable (same as sdrf)
phenoData <- pData(eset)


# See the dimensions of pData from eset
dim(pData(eset))





# DEFINE EXPERIMENTAL GROUPS
# Define sample groups
group <- factor(phenoData$`tissue:ch1`)
pairs <- factor(rep(1:32, each = 2))

# Rename factor
levels(group) <- c("Atheroma", "Intact")





# QC VISUALIZATION
# BOXPLOT
# Increase bottom margin to make room for labels
par(mar = c(10, 4, 4, 2))

boxplot(exprs(eset), col = c(rep("red", 32), rep("blue", 32)),
        las = 2, main ="Boxplot of exprs(eset)", ylab ="Expression values")



# HISTOGRAM
# Divide plotting window
par(mfrow=c(1,2))

hist(exprs(eset), main = "Histogram of exprs(eset)", 
     xlab = "exprs(eset)", ylab = "Frequency")




# Q-Q PLOT
qqnorm(exprs(eset), main = "Normal Q-Q Plot", 
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles") 






# LINEAR MODELLING
# Load package
library(limma)


# Define design matrix
studyDesign <- data.frame(Group = group, Pairs = pairs)
designMatrix <- model.matrix(~ 0 + Group + Pairs, data = studyDesign)


# Change column 1 and 2 to preferred labels
colnames(designMatrix)[1:2] <- c("Atheroma", "Intact")

# Fit linear model
fit <- lmFit(exprs(eset), designMatrix)


# Define contrast matrix
contrastMatrix <- makeContrasts(AtheromaVsIntact = Atheroma - Intact,
                                levels = designMatrix)

# Specify contrasts
fit <- contrasts.fit(fit, contrastMatrix)


# Perform statistical testing
fit <- eBayes(fit, robust = TRUE)


# Get DEGs
results <- topTable(fit, coef = "AtheromaVsIntact", number = nrow(eset),
                    p.value = 0.01, lfc = 1)

# Check results
print(dim(results))
print(head(results, 30))






# DEGs RESULTS ANNOTATION
# Load annotation package for this platform
library(hugene10sttranscriptcluster.db)

# Add gene symbols to results
results$symbol <- mapIds(hugene10sttranscriptcluster.db,
                         keys = rownames(results),
                         column = "SYMBOL",
                         keytype = "PROBEID",
                         multiVals = "first")

# Check results with gene symbols
print(head(results, 60))


# Retain the significant annotated results for the heatmap.
deg_results <- results

# Save DEGs
write.table(results,
            file = "results/147_DEGs.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)









# DEGs PREPARATION FOR DOWNSTREAM ANALYSIS
# Import DEGs
deg <- read.table("results/147_DEGs.txt", header = TRUE, as.is = TRUE, sep = "\t")


# Clean the table "rows" for downstream analysis
deg <- deg[!is.na(deg$symbol), ]
deg <- deg[deg$symbol != "", ]
deg <- deg[!duplicated(deg$symbol), ]


# Check table
nrow(deg)


# Save cleaned DEGs
write.table(deg,
            "results/Cleaned_DEGs.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


# Extract full gene names for downstream analysis
g_symbols <- deg$symbol

head(g_symbols)
length(g_symbols)


# Save the gene names
write.table(g_symbols,
            file = "results/GeneSymbols_DEGs.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)





# Extract upregulated DEGs for downstream analysis
deg <- read.table("results/Cleaned_DEGs.txt",
                  header = TRUE,
                  sep = "\t",
                  as.is = TRUE)

up_deg <- deg[deg$logFC > 0, ]
nrow(up_deg)




# Extract gene names for downstream analysis
write.table(up_deg$symbol,
            "results/Gene_Symbols_Upregulated_DEGs.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)








# SAVE FULL GENE LIST FOR GSEA 
# Get ALL genes, not just DEGs
results <- topTable(fit, coef = "AtheromaVsIntact", 
                        number = Inf, # Pulls every row in the dataset
                        p.value = 1,   # Removes p-value filter
                        lfc = 0)      # Removes fold change filter

# Check the dimension - it should be thousands of rows
print(dim(results))


# Load annotation package for this platform
library(hugene10sttranscriptcluster.db)


# Add gene symbols to the full list
results$symbol <- mapIds(hugene10sttranscriptcluster.db,
                             keys = rownames(results),
                             column = "SYMBOL",
                             keytype = "PROBEID",
                             multiVals = "first")


# Remove rows without symbols
clean_results <- results[!is.na(results$symbol), ]


# If multiple probes map to the same symbol, keep the one with the highest absolute logFC
clean_results <- clean_results[order(abs(clean_results$logFC), decreasing = TRUE), ]
clean_results <- clean_results[!duplicated(clean_results$symbol), ]


# Create the final ranked vector
gsea_input <- clean_results$logFC
names(gsea_input) <- clean_results$symbol


# Sort in descending order for clusterProfiler
gsea_input <- sort(gsea_input, decreasing = TRUE)


# Save the full list with Probe IDs as row names
write.table(clean_results, 
            file = "results/Ranked_logFC_Full_Gene_List.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = TRUE)




# Reading the file back while recognizing the Probe IDs
full_data <- read.table("results/Ranked_logFC_Full_Gene_List.txt", 
                        header = TRUE, 
                        sep = "\t", 
                        row.names = 1)









# DEGs VISUALIZATION
# VOLCANO PLOT
# Install packages once using the instructions in CODE_SETUP.md.

library(EnhancedVolcano)

# Get full results without filtering
full_results <- topTable(fit, 
                         coef = "AtheromaVsIntact", 
                         number = Inf)



# Make sure symbol column exists
head(full_results)

# Check if symbol column exists
colnames(full_results)


library(hugene10sttranscriptcluster.db)

# Add gene symbols to full results
full_results$symbol <- mapIds(hugene10sttranscriptcluster.db,
                              keys = rownames(full_results),
                              column = "SYMBOL",
                              keytype = "PROBEID",
                              multiVals = "first")

# Replace NA symbols with probe ID
full_results$symbol[is.na(full_results$symbol)] <- rownames(full_results)[is.na(full_results$symbol)]

# Now run volcano plot
EnhancedVolcano(full_results,
                lab = full_results$symbol,
                x = "logFC",
                y = "adj.P.Val",
                pCutoff = 0.01,
                FCcutoff = 1,
                title = "Atherosclerosis vs Normal Tissue",
                subtitle = "GSE43292",
                pointSize = 1.5,
                labSize = 3)




# HEATMAP
# Preserve the original selection: take 12 probe rows, then remove missing symbols.
library(pheatmap)

# Rename groups for cleaner visualization
levels(group) <- c("AS", "Normal")

# Select probe rows by absolute logFC using the original heatmap settings.
top10 <- deg_results[order(abs(deg_results$logFC), decreasing = TRUE), ][1:12, ]

# Remove rows with missing gene symbols
top10 <- top10[!is.na(top10$symbol), ]

# Extract expression values using probe IDs
heatmap_data <- exprs(eset)[rownames(top10), ]

# Replace probe IDs with gene symbols
rownames(heatmap_data) <- top10$symbol

# Scale expression values by gene (Z-score)
heatmap_scaled <- t(scale(t(heatmap_data)))

# Create sample annotation
annotation_col <- data.frame(
  Group = group
)

rownames(annotation_col) <- colnames(heatmap_scaled)

# Order samples by group
sample_order <- order(group)

heatmap_scaled <- heatmap_scaled[, sample_order]

annotation_col <- annotation_col[sample_order, , drop = FALSE]

# Generate heatmap
pheatmap(heatmap_scaled,
         scale = "row",
         annotation_col = annotation_col,
         cluster_cols = FALSE,
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         border_color = NA,
         main = "Top 10 DEGs in Atherosclerosis and Normal Tissue Samples")
