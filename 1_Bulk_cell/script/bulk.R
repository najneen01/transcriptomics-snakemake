###BULK CELL TRANSCRIPTOMIC###
#Najneen Rejwana
#count data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247791
#metadata: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1040329&o=acc_s%3Aa

setwd("G:/Rstudio/TRANSCRIPTOMICS/bulk_trans/GSE247791")
getwd()


library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(ggrepel)

count_data <- read.csv("count_data_GSE247791_EDIT.csv", row.names = 1)
meta_data <- read.csv("meta_data_GSE247791_EDIT.csv", row.names = 1)
View(count_data)
View(meta_data)

str(count_data)
str(metadata)

# Check if sample names in count matrix match metadata row names
all(colnames(count_data) == rownames(meta_data))
#[1] FALSE


setdiff(colnames(count_data), rownames(meta_data))  # Samples in count_matrix but not in metadata
#[1] "ARK1_Abema_A"   "ARK1_Abema_B"   "ARK1_DMSO_A"    "ARK1_DMSO_B"    "HEC1B_Abema_A" 
#[6] "HEC1B_Abema_B"  "HEC1B_DMSO_A"   "HEC1B_DMSO_B"   "EMCA.A_DMSO_B"  "EMCA.A_DMSO_A" 
#[11] "EMCA.A_Abema_B" "EMCA.A_Abema_A" "CS.A_DMSO_B"    "CS.A_DMSO_A"    "CS.A_Abema_B"  
#[16] "CS.A_Abema_A"   "UPSC.A_DMSO_B"  "UPSC.A_DMSO_A"  "UPSC.A_Abema_B" "UPSC.A_Abema_A"
#[21] "UPSC.A_DMSO_1"  "UPSC.A_DMSO_2"  "CS.A_DMSO_1"    "CS.A_DMSO_2"    "EMCA.A_DMSO_1" 
#[26] "EMCA.A_DMSO_2"  "CS.B_DMSO_1"    "CS.B_DMSO_2"    "AN3_DMSO1"      "AN3_DMSO2"     
#[31] "ARK1_DMSO1"     "ARK1_DMSO2"     "HEC_DMSO1"      "HEC_DMSO2"      "KLE_DMSO1"     
#[36] "KLE_DMSO2"     
setdiff(rownames(meta_data), colnames(count_data))  # Samples in metadata but not in count_matrix
#[1] "SRR26830900" "SRR26830901" "SRR26830904" "SRR26830905" "SRR26830908" "SRR26830909"
#[7] "SRR26830912" "SRR26830913" "SRR26830916" "SRR26830917" "SRR26830920" "SRR26830921"
#[13] "SRR26830928" "SRR26830929" "SRR26830932" "SRR26830933" "SRR26830936" "SRR26830937"
#[19] "SRR26830938" "SRR26830939" "SRR26830940" "SRR26830941" "SRR26830942" "SRR26830943"
#[25] "SRR26830948" "SRR26830949" "SRR26830950" "SRR26830951" "SRR26830952" "SRR26830953"
#[31] "SRR26830954" "SRR26830955" "SRR26830956" "SRR26830957" "SRR26830958" "SRR26830959"

count_data_edit <- read.csv("count_data_GSE247791_EDIT.csv", row.names = 1)
meta_data_edit <- read.csv("meta_data_GSE247791_EDIT.csv", row.names = 1)
View(count_data_edit)
View(meta_data_edit)

str(count_data_edit)
str(meta_data_edit)


ncol(count_data_edit)  # Number of samples in count matrix
#[1] 36
nrow(meta_data_edit)   # Number of samples in metadata
#[1] 36
all(colnames(count_data_edit) == rownames(meta_data_edit))
#[1] FALSE


print(colnames(count_data_edit))
# [1] "SRR26830959" "SRR26830958" "SRR26830957" "SRR26830956" "SRR26830955" "SRR26830954"
#[7] "SRR26830953" "SRR26830952" "SRR26830951" "SRR26830950" "SRR26830949" "SRR26830948"
#[13] "SRR26830943" "SRR26830942" "SRR26830941" "SRR26830940" "SRR26830939" "SRR26830938"
#[19] "SRR26830937" "SRR26830936" "SRR26830933" "SRR26830932" "SRR26830929" "SRR26830928"
#[25] "SRR26830921" "SRR26830920" "SRR26830917" "SRR26830916" "SRR26830913" "SRR26830912"
#[31] "SRR26830909" "SRR26830908" "SRR26830905" "SRR26830904" "SRR26830901" "SRR26830900"
print(rownames(meta_data_edit))
#[1] "SRR26830900" "SRR26830901" "SRR26830904" "SRR26830905" "SRR26830908" "SRR26830909"
#[7] "SRR26830912" "SRR26830913" "SRR26830916" "SRR26830917" "SRR26830920" "SRR26830921"
#[13] "SRR26830928" "SRR26830929" "SRR26830932" "SRR26830933" "SRR26830936" "SRR26830937"
#[19] "SRR26830938" "SRR26830939" "SRR26830940" "SRR26830941" "SRR26830942" "SRR26830943"
#[25] "SRR26830948" "SRR26830949" "SRR26830950" "SRR26830951" "SRR26830952" "SRR26830953"
#[31] "SRR26830954" "SRR26830955" "SRR26830956" "SRR26830957" "SRR26830958" "SRR26830959"


###It seems that the order of the sample names between count_data_edit and meta_data_edit is reversed. The sample names are present in both, but they are in a different order.


#Align the Order:
count_data_edit <- count_data_edit[, order(colnames(count_data_edit))]
meta_data_edit <- meta_data_edit[order(rownames(meta_data_edit)), , drop=FALSE]

all(colnames(count_data_edit) == rownames(meta_data_edit))
#[1] TRUE



####DSEQ2 analysis---------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count_data_edit,
                              colData = meta_data_edit,
                              design = ~ treatment)
#Warning message:
#In DESeqDataSet(se, design = design, ignoreRank) :
# some variables in design formula are characters, converting to factor


#Check conversion
str(meta_data_edit$treatment)  #before conversion
#chr [1:36] "DMSO" "DMSO" "DMSO" "DMSO" "DMSO" "DMSO" "DMSO" "DMSO" "DMSO" "DMSO" "DMSO" ...

##After conversion
# Convert to factor if not already done
meta_data_edit$treatment <- factor(meta_data_edit$treatment)

# Show structure again
str(meta_data_edit$treatment)
#Factor w/ 2 levels "Abemaciclib",..: 2 2 2 2 2 2 2 2 2 2 ...
#It looks like the treatment column has been successfully converted into a factor with two levels: "Abemaciclib" and "DMSO". The numbers (e.g., 2 2 2 2 2 ...) represent the factor levels (where "Abemaciclib" is 1 and "DMSO" is 2 in the factor encoding).


table(meta_data_edit$treatment)   #Check how many samples fall under each treatment group,
#Abemaciclib        DMSO 
#10                 26 




############# Run DESeq2 analysis#####################################
dds <- DESeqDataSetFromMatrix(countData = count_data_edit,
                              colData = meta_data_edit,
                              design = ~ treatment)

# Pre-filter low count genes (optional but recommended)
dds <- dds[rowSums(counts(dds)) > 1, ]

# Run the DESeq function
dds <- DESeq(dds)

# Get results
res <- results(dds)

# View results
head(res)
View(res)



# Plot MA plot
plotMA(res, main="DESeq2 MA Plot")



# Load DESeq2 library
library(DESeq2)


# Create the MA plot
plotMA(res, main = "DESeq2 MA Plot")

# Define upregulated and downregulated genes
upregulated_genes <- rownames(res)[res$log2FoldChange > 1 & res$padj < 0.05]
downregulated_genes <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.05]

# Add labels for upregulated genes
text(x = log10(res[upregulated_genes, "baseMean"]), 
     y = res[upregulated_genes, "log2FoldChange"], 
     labels = upregulated_genes, 
     pos = 2,  # Position labels to the right of the points
     col = "green", 
     cex = 0.6)  # Adjust label size

# Add labels for downregulated genes
text(x = log10(res[downregulated_genes, "baseMean"]), 
     y = res[downregulated_genes, "log2FoldChange"], 
     labels = downregulated_genes, 
     pos = 2,  # Position labels to the left of the points
     col = "red", 
     cex = 0.6)  # Adjust label size

# Add a legend
legend("topright",                        # Position of the legend
       legend = c("Upregulated", "Downregulated"),  # Legend labels
       col = c("green", "red"),            # Colors for the legend
       pch = 16,                          # Use solid circles in the legend
       cex = 0.8,                         # Adjust legend text size
       bty = "n")                         # Remove legend box


cat("Upregulated Genes:\n", paste(upregulated_genes, collapse = ", "), "\n")
#NA, NA, AKR1B10, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, PDE4D, NA, NA, NA, NA, NA, NA, NA, NA, NA, WDR76, NA, NA, WDR62, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, KRT81, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, ABCB1, NA, NA, NA, NA, NA, NA, NA, NA, NA, FN1, MAGEA6, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, CKAP2L, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, DMBX1, NA, NA, NA, SPC25, NA, MAGEC2, NA, NA, NA, NA, MKI67, NA, NA, NA, HJURP, NA, MND1, NA, NA, NA, NA, NA, NA, NA, NA, SPINK6, NA, NA, NA, NA, NA, NA, NA, E2F8, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, RPS4Y1, CDKN3, PLAT, NA, NA, HSD17B6, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, RP11-169F17.1, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, CDC45, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, FLG, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, ARHGAP33, NA, NA, NA, ESPL1, NA, NA, NA, KIF2C, NA, KIF15, NA, NA, NA, NA, NA, NA, LMNB1, NA, NA, RUNX1T1, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, CCDC109B, CCNA2, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, GNG2, NA, NA, NA, NA, NA, NA, NA, HMGB2, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, SP8, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, E2F2, NA, CPA4, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, DDX3Y, NA, ALDH1A2, NA, NGFR, NA, NA, NA, GTSE1, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, PBK, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, MAGEC1, NA, NA, NA, NA, NA, NA, NA 

cat("Downregulated Genes:\n", paste(downregulated_genes, collapse = ", "), "\n")
# NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, GOLGA8Q, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, SLC6A13, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, IL1A, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA 



# Save the results to a file
write.csv(as.data.frame(res), file="deseq2_results.csv")


# Save full results for Upregulated and Downregulated DEGs
write.csv(as.data.frame(res_clean[res_clean$log2FoldChange > 1 & res_clean$padj < 0.05, ]), "result/upregulated_DEGs_full.csv")
write.csv(as.data.frame(res_clean[res_clean$log2FoldChange < -1 & res_clean$padj < 0.05, ]), "result/downregulated_DEGs_full.csv")







# Plot PCA
vsd <- vst(dds, blind = FALSE)  # Variance stabilizing transformation
pcaData <- plotPCA(vsd, intgroup = "treatment", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Samples")





# Remove rows with NA values in log2FoldChange or padj
res_clean <- res[!is.na(res$log2FoldChange) & !is.na(res$padj),]

# Filter for upregulated genes (log2 fold change > 1 and padj < 0.05)
upregulated_genes <- rownames(res_clean[res_clean$log2FoldChange > 1 & res_clean$padj < 0.05,])

# Filter for downregulated genes (log2 fold change < -1 and padj < 0.05)
downregulated_genes <- rownames(res_clean[res_clean$log2FoldChange < -1 & res_clean$padj < 0.05,])

# Print the names of the upregulated and downregulated genes
print("Upregulated genes:")
print(upregulated_genes)
# print(upregulated_genes)
#[1] "AKR1B10"       "PDE4D"         "WDR76"         "WDR62"         "KRT81"        
#[6] "ABCB1"         "FN1"           "MAGEA6"        "CKAP2L"        "DMBX1"        
#[11] "SPC25"         "MAGEC2"        "MKI67"         "HJURP"         "MND1"         
#[16] "SPINK6"        "E2F8"          "RPS4Y1"        "CDKN3"         "PLAT"         
#[21] "HSD17B6"       "RP11-169F17.1" "CDC45"         "FLG"           "ARHGAP33"     
#[26] "ESPL1"         "KIF2C"         "KIF15"         "LMNB1"         "RUNX1T1"      
#[31] "CCDC109B"      "CCNA2"         "GNG2"          "HMGB2"         "SP8"          
#[36] "E2F2"          "CPA4"          "DDX3Y"         "ALDH1A2"       "NGFR"         
#[41] "GTSE1"         "PBK"           "MAGEC1"       

print("Downregulated genes:")
print(downregulated_genes)
#[1] "GOLGA8Q" "SLC6A13" "IL1A" 

###Enhanced volcano plot----------------------------------------------------------
# Install EnhancedVolcano if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")

# Load EnhancedVolcano library
library(EnhancedVolcano)

# Clean data
res_clean <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]

dev.off()
graphics.off()

library(EnhancedVolcano)

res_clean <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]

keyvals <- ifelse(res_clean$log2FoldChange > 1 & res_clean$padj < 0.05, "Up",
                  ifelse(res_clean$log2FoldChange < -1 & res_clean$padj < 0.05, "Down", "NS"))

keyvals_col <- ifelse(keyvals == "Up", "green",
                      ifelse(keyvals == "Down", "red", "grey"))

names(keyvals_col) <- keyvals

labels <- ifelse(keyvals != "NS", rownames(res_clean), "")


# Label only top 15 genes (clean plot)
#top_genes <- head(rownames(res_clean[order(res_clean$padj), ]), 15)
#labels <- ifelse(rownames(res_clean) %in% top_genes,
 #                rownames(res_clean), "")

EnhancedVolcano(res_clean,
                lab = labels,
                x = "log2FoldChange",
                y = "padj",
                colCustom = keyvals_col,
                pCutoff = 0.05,
                FCcutoff = 1,
                drawConnectors = FALSE,
                max.overlaps = 20,
                title = "Volcano Plot")

#----------------------------------------------------------------------------------------

## Select top 20 differentially expressed genes by adjusted p-value
topGenes <- head(order(res$padj), 20)

# Plot heatmap
library(pheatmap)
pheatmap(assay(vsd)[topGenes,], 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         annotation_col = meta_data_edit[, "treatment", drop = FALSE],
         scale = "row", 
         main = "Heatmap of Top 20 DE Genes")













------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------
###pathway enrichment analysis
  
  # Install and load required packages
  if (!requireNamespace("BiocManager", quietly = TRUE, force = TRUE))
    install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

library(clusterProfiler)
library(org.Hs.eg.db)

# Define gene lists
upregulated_genes <- c("AKR1B10", "PDE4D", "WDR76", "WDR62", "KRT81", "ABCB1", "FN1", "MAGEA6", "CKAP2L", "DMBX1", 
                       "SPC25", "MAGEC2", "MKI67", "HJURP", "MND1", "SPINK6", "E2F8", "RPS4Y1", "CDKN3", "PLAT", 
                       "HSD17B6", "RP11-169F17.1", "CDC45", "FLG", "ARHGAP33", "ESPL1", "KIF2C", "KIF15", "LMNB1", 
                       "RUNX1T1", "CCDC109B", "CCNA2", "GNG2", "HMGB2", "SP8", "E2F2", "CPA4", "DDX3Y", "ALDH1A2", 
                       "NGFR", "GTSE1", "PBK", "MAGEC1")

downregulated_genes <- c("GOLGA8Q", "SLC6A13", "IL1A")

# Convert gene symbols to Entrez IDs
upregulated_entrez <- bitr(upregulated_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
downregulated_entrez <- bitr(downregulated_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Extract Entrez IDs
upregulated_entrez_ids <- upregulated_entrez$ENTREZID
downregulated_entrez_ids <- downregulated_entrez$ENTREZID

# View the results
print("Upregulated Entrez IDs:")
print(upregulated_entrez_ids)
#[1] "57016"  "5144"   "79968"  "284403" "3887"   "5243"   "2335"   "4105"   "150468"
#[10] "127343" "57405"  "51438"  "4288"   "55355"  "84057"  "404203" "79733"  "6192"  
#[19] "1033"   "5327"   "8630"   "8318"   "2312"   "115703" "9700"   "11004"  "56992" 
#[28] "4001"   "862"    "890"    "54331"  "3148"   "221833" "1870"   "51200"  "8653"  
#[37] "8854"   "4804"   "51512"  "55872"  "9947"  

print("Downregulated Entrez IDs:")
print(downregulated_entrez_ids)
#[1] "727909" "6540"   "3552"







## Perform GO enrichment for upregulated genes
go_up <- enrichGO(gene = upregulated_entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",  # Biological Process
                  pAdjustMethod = "BH",  # Benjamini-Hochberg correction
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)

# Perform GO enrichment for downregulated genes
go_down <- enrichGO(gene = downregulated_entrez_ids,
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP",  # Biological Process
                    pAdjustMethod = "BH",  # Benjamini-Hochberg correction
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

# Extract key columns from the GO enrichment results
go_results <- as.data.frame(go_up)[, c("ID", "Description", "p.adjust", "geneID")]

# Rename columns for clarity
colnames(go_results) <- c("GO Term ID", "GO Term Description", "Adjusted p-value", "Genes")

# Sort by adjusted p-value (most significant terms first)
go_results <- go_results[order(go_results$`Adjusted p-value`), ]

# Display the top 10 results
print(head(go_results, 10))

# Convert Entrez IDs to gene symbols
gene_symbols <- bitr(unlist(strsplit(go_results$Genes, "/")), fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

# Replace Entrez IDs with gene symbols in the results
go_results$Genes <- sapply(strsplit(go_results$Genes, "/"), function(ids) {
  symbols <- gene_symbols$SYMBOL[match(ids, gene_symbols$ENTREZID)]
  paste(symbols, collapse = ", ")
})

# Display the cleaned results
print(head(go_results, 10))


# Save the cleaned results to a CSV file
write.csv(go_results, file = "result/UPREGULATED_go_enrichment_cleaned_results.csv", row.names = FALSE)

----------------------------------------------------------------
  
  
  
  # Extract key columns from the GO enrichment results
  go_down_results <- as.data.frame(go_down)[, c("ID", "Description", "p.adjust", "geneID")]

# Rename columns for clarity
colnames(go_down_results) <- c("GO Term ID", "GO Term Description", "Adjusted p-value", "Genes")

# Sort by adjusted p-value (most significant terms first)
go_down_results <- go_down_results[order(go_down_results$`Adjusted p-value`), ]

# Convert Entrez IDs to gene symbols
gene_symbols_down <- bitr(unlist(strsplit(go_down_results$Genes, "/")), fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

# Replace Entrez IDs with gene symbols in the results
go_down_results$Genes <- sapply(strsplit(go_down_results$Genes, "/"), function(ids) {
  symbols <- gene_symbols_down$SYMBOL[match(ids, gene_symbols_down$ENTREZID)]
  paste(symbols, collapse = ", ")
})

# Display the cleaned results
print(head(go_down_results, 10))

# Save the cleaned results to a CSV file
write.csv(go_down_results, file = "result/DOWN_REGULATED_enrichment_cleaned_results.csv", row.names = FALSE)

