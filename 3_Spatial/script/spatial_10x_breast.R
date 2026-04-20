
##spatial transcriptomic- 10x visium##
#Najneen Rejwana
#Source : https://www.10xgenomics.com/datasets/fresh-frozen-visium-on-cytassist-human-breast-cancer-probe-based-whole-transcriptome-profiling-2-standard

#-------------folder arrangement should be like this-------------------------
# G:/Rstudio/TRANSCRIPTOMICS/spatial/10x_visium/breast_cancer/
#  │── filtered_feature_bc_matrix.h5
#  │── raw_feature_bc_matrix.h5
#│── spatial/
#    ├── tissue_lowres_image.png  (Ensure this is present)
#│   ├── scalefactors_json.json   (Required for correct scaling)
#│   ├── tissue_positions_list.csv (Required for spot positioning)





getwd()
setwd("G:/Rstudio/TRANSCRIPTOMICS/spatial/10x_visium/breast_cancer")

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(SeuratData)
library(patchwork)
library(tiff)
library(png)




  
#install.packages("Seurat")
#ls("package:Seurat")  

#install.packages("remotes")  
#remotes::install_github("satijalab/seurat-data", force = TRUE)

#install.packages("tiff")



setwd("G:/Rstudio/TRANSCRIPTOMICS/spatial/10x_visium/breast_cancer")
list.files()
#[1] "CytAssist_Fresh_Frozen_Human_Breast_Cancer_analysis.tar.gz"              
#[2] "CytAssist_Fresh_Frozen_Human_Breast_Cancer_cloupe.cloupe"                
#[3] "CytAssist_Fresh_Frozen_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
#[4] "CytAssist_Fresh_Frozen_Human_Breast_Cancer_image.png"                    
#[5] "CytAssist_Fresh_Frozen_Human_Breast_Cancer_image.tif"                    
#[6] "CytAssist_Fresh_Frozen_Human_Breast_Cancer_metrics_summary.csv"          
#[7] "CytAssist_Fresh_Frozen_Human_Breast_Cancer_raw_feature_bc_matrix.h5"     
#[8] "CytAssist_Fresh_Frozen_Human_Breast_Cancer_spatial.tar.gz"               
#[9] "spatial"                                                                 
#[10] "spatial_10x_breast.R"                                                    
#[11] "tissue_lowres_image.png" 

list.files("spatial")
#[1] "aligned_fiducials.jpg"     "aligned_tissue_image.jpg"  "cytassist_image.tiff"     
#[4] "detected_tissue_image.jpg" "scalefactors_json.json"    "spatial_enrichment.csv"   
#[7] "tissue_hires_image.png"    "tissue_lowres_image.png"   "tissue_positions.csv"     



##Loading Spatial data into Seurat-------------------------------------------------------------------
data_dir <- getwd()
h5_file <- "CytAssist_Fresh_Frozen_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"

# Load the spatial dataset
seurat_obj <- Load10X_Spatial(
  data.dir = data_dir,
  filename = h5_file,
  assay = "Spatial",
  slice = "breast_tissue",
  filter.matrix = TRUE
)

# Print summary of Seurat object
seurat_obj

#An object of class Seurat 
#18085 features across 1657 samples within 1 assay 
#Active assay: Spatial (18085 features, 0 variable features)
#1 layer present: counts
#1 spatial field of view present: breast_tissue






##Preprocessing data-------------------------------------------------------------------

# Normalize using SCTransform
seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = FALSE)

# Dimensionality reduction
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# Plot UMAP
DimPlot(seurat_obj, reduction = "umap", label = TRUE)



##Visualize Spatial data---------------------------------------------------------------

# View tissue image with spots
SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3)

# Plot marker gene expression in spatial context
SpatialFeaturePlot(seurat_obj, features = c("MS4A1", "EPCAM", "CD3D"))






##Differential gene expression------------------------------------------------------------

# Find differentially expressed genes
markers <- FindAllMarkers(seurat_obj, min.pct = 0.25, logfc.threshold = 0.25)

# View top marker genes
head(markers)
#     p_val           avg_log2F pct.1 pct.2 p_val_adj        cluster gene
#VIM    1.499643e-207   2.054761     1 0.917 2.251414e-203       0    VIM
#COL3A1 1.038675e-206   1.976601     1 0.979 1.559363e-202       0 COL3A1
#CD74   9.252261e-200   2.658123     1 0.923 1.389042e-195       0   CD74
#COL1A1 3.248479e-198   2.028522     1 1.000 4.876942e-194       0 COL1A1
#COL1A2 1.278452e-192   1.918285     1 0.949 1.919341e-188       0 COL1A2
#TIMP2  1.082596e-189   2.029384     1 0.795 1.625302e-185       0  TIMP2




# Define some marker genes for cell-type identification
marker_genes <- c("MS4A1", "EPCAM", "CD3D", "PTPRC", "CD19")

# Plot expression of these marker genes
SpatialFeaturePlot(seurat_obj, features = marker_genes)


# Visualizing expression of a gene in spatial context
SpatialFeaturePlot(seurat_obj, features = "MS4A1")








# Create a volcano plot
library(ggplot2)

# Add a column for significance based on p-value threshold
de_genes$significance <- ifelse(de_genes$p_val_adj < 0.05, "Significant", "Not Significant")

# Plot
ggplot(de_genes, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Not Significant" = "red", "Significant" = "blue")) +
  theme_minimal() +
  ggtitle("Volcano Plot: Cluster 0 vs Cluster 1") +
  xlab("Log2 Fold Change") +
  ylab("-log10 Adjusted P-value") +
  theme(legend.position = "top")




#Dotplot---------------------------------------------------
# Get the top 10 genes by p-value
top_genes <- head(de_genes[order(de_genes$p_val_adj), ], 10)

# Create a dot plot for the top 10 genes
top_gene_names <- rownames(top_genes)
DotPlot(seurat_obj, features = top_gene_names) + 
  ggtitle("Top 10 Differentially Expressed Genes: Cluster 0 vs Cluster 1")



#Heatmap----------------------------------------------------
# Select the top 10 genes based on p-value
top_genes_for_heatmap <- rownames(top_genes)

# Plot the heatmap
DoHeatmap(seurat_obj, features = top_genes_for_heatmap) + 
  ggtitle("Heatmap: Top 10 Differentially Expressed Genes")



#MA Plot-------------------------------------------------------
ggplot(de_genes, aes(x = avg_log2FC, y = avg_log2FC)) +
  geom_point(alpha = 0.6) +
  ggtitle("MA Plot: Cluster 0 vs Cluster 1") +
  xlab("Log2 Fold Change") +
  ylab("Log2 Fold Change") +
  theme_minimal()





##Number of genes----------------------------------------------------------------
#Total Number of Genes in de_genes
nrow(de_genes)
#[1] 11300


#Number of significant DEGs (Adjusted p-value < 0.05)
sum(de_genes$p_val_adj < 0.05, na.rm = TRUE)
#[1] 4701



#Number of upregulated and downregulated genes
# Count upregulated genes (log2FC > 0 & adj p-value < 0.05)
sum(de_genes$avg_log2FC > 0 & de_genes$p_val_adj < 0.05, na.rm = TRUE)
#[1] 1105

# Count downregulated genes (log2FC < 0 & adj p-value < 0.05)
sum(de_genes$avg_log2FC < 0 & de_genes$p_val_adj < 0.05, na.rm = TRUE)
#[1] 3596



#Filter Upregulated and Downregulated Genes-------------------
# Define thresholds
logFC_threshold <- 0  
pval_threshold <- 0.05  # Adjusted p-value threshold

# Filter upregulated genes (log2FC > 0 & adj p-value < 0.05)
upregulated_genes <- de_genes[de_genes$avg_log2FC > logFC_threshold & de_genes$p_val_adj < pval_threshold, ]

# Filter downregulated genes (log2FC < 0 & adj p-value < 0.05)
downregulated_genes <- de_genes[de_genes$avg_log2FC < logFC_threshold & de_genes$p_val_adj < pval_threshold, ]


# Save upregulated genes to CSV
write.csv(upregulated_genes, "upregulated_genes_breast_spatial.csv", row.names = TRUE)

# Save downregulated genes to CSV
write.csv(downregulated_genes, "downregulated_genes_breast_spatial.csv", row.names = TRUE)



# Save top upregulated genes
write.csv(top_10_up, "top_10_upregulated_genes.csv", row.names = TRUE)

# Save top downregulated genes
write.csv(top_10_down, "top_10_downregulated_genes.csv", row.names = TRUE)







#Visualize using Bar plot---------------------------------------------------------
library(ggplot2)

# Combine top 10 upregulated and downregulated genes into one dataframe
top_genes_combined <- rbind(
  data.frame(Gene = rownames(top_10_up), Log2FC = top_10_up$avg_log2FC, Direction = "Upregulated"),
  data.frame(Gene = rownames(top_10_down), Log2FC = top_10_down$avg_log2FC, Direction = "Downregulated")
)

# Plot
ggplot(top_genes_combined, aes(x = reorder(Gene, Log2FC), y = Log2FC, fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 Upregulated and Downregulated Genes",
       x = "Gene",
       y = "Log2 Fold Change") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue"))




#Dot plot---------------------------
DotPlot(seurat_obj, features = c(rownames(top_10_up), rownames(top_10_down))) +
  ggtitle("Expression of Top 10 Upregulated and Downregulated Genes")






##volcano plot for upregulated & downregulated genes---------

library(ggplot2)

# Define significance thresholds
logFC_threshold <- 0.5  # Change if needed
pval_threshold <- 0.05

# Add a new column to classify genes as Upregulated, Downregulated, or Not Significant
de_genes$Significance <- ifelse(
  de_genes$p_val_adj < pval_threshold & de_genes$avg_log2FC > logFC_threshold, "Upregulated",
  ifelse(de_genes$p_val_adj < pval_threshold & de_genes$avg_log2FC < -logFC_threshold, "Downregulated", "Not Significant")
)

# Create the volcano plot
ggplot(de_genes, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "purple", "Downregulated" = "orange", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Differential Expression",
       x = "Log2 Fold Change",
       y = "-log10 Adjusted P-value") +
  theme(legend.position = "top")







#Highlight most significant genes------------
library(ggplot2)
install.packages("ggrepel")
library(ggrepel)



library(ggplot2)
library(ggrepel)  

# Define thresholds
logFC_threshold <- 0.5  
pval_threshold <- 0.05  

# Classify genes based on significance
de_genes$Significance <- ifelse(
  de_genes$p_val_adj < pval_threshold & de_genes$avg_log2FC > logFC_threshold, "Upregulated",
  ifelse(de_genes$p_val_adj < pval_threshold & de_genes$avg_log2FC < -logFC_threshold, "Downregulated", "Not Significant")
)

# Select top 15 most significant genes from both upregulated and downregulated categories
top_up <- rownames(head(de_genes[order(-de_genes$avg_log2FC, de_genes$p_val_adj), ], 15))
top_down <- rownames(head(de_genes[order(de_genes$avg_log2FC, de_genes$p_val_adj), ], 15))
top_genes <- c(top_up, top_down)

# Create the volcano plot with black text labels
ggplot(de_genes, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "green", "Downregulated" = "yellow", "Not Significant" = "gray")) +
  geom_text_repel(data = subset(de_genes, rownames(de_genes) %in% top_genes), 
                  aes(label = rownames(subset(de_genes, rownames(de_genes) %in% top_genes))), 
                  color = "black",  # Set text color to black
                  size = 3, max.overlaps = 15) +
  theme_minimal() +
  labs(title = "Volcano Plot: Differential Expression",
       x = "Log2 Fold Change",
       y = "-log10 Adjusted P-value") +
  theme(legend.position = "top")







#Feature plot (for specific genes)------------------------------

library(Seurat)

# Visualize expression of specific genes in Seurat object------------------------
FeaturePlot(seurat_obj, features = c("ERBB2", "COL6A3", "VIM", "GRB7", "NDUFB9", "FOXD1"))






##Spatial Gene Expression Enrichment------------------------------------------------------------------------
#Enrichment of Specific Pathways: Use gene set enrichment analysis (GSEA) to identify if certain biological pathways are enriched in different spatial regions of the tissue.

install.packages("clusterProfiler")
library(clusterProfiler)

#Suppose for pyrimidine pathway-------------------
#To perform a Gene Set Enrichment Analysis (GSEA) for the pyrimidine pathway, you need to first identify the genes associated with this pathway. One common resource for this is the KEGG pathway database, which provides a collection of gene sets for various biological pathways, including the pyrimidine pathway.

# Check the available identities in this Seurat object
levels(seurat_obj)
#[1] "0" "1" "2" "3" "4" "5" "6" "7"



# Perform differential expression between Cluster 0 and Cluster 1-------------------------------
de_genes <- FindMarkers(seurat_obj, ident.1 = "0", ident.2 = "1")
print(de_genes)



#Prepare a ranked gene list-----------------------------------
# Create a ranked gene list based on avg_log2FC
ranked_genes <- de_genes[, "avg_log2FC"]
names(ranked_genes) <- rownames(de_genes)

# Sort genes from highest to lowest log fold change
ranked_genes <- sort(ranked_genes, decreasing = TRUE)


#Obtain pyrimidine pathway from KEGG---------------
if (!requireNamespace("clusterProfiler", quietly = TRUE)) install.packages("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) install.packages("org.Hs.eg.db")

library(clusterProfiler)
library(org.Hs.eg.db)  # Human gene annotation





#Retrieve KEGG Pathway Gene Sets----------------------------------------------------------------------------------------
# Download KEGG pathways for humans
kegg_pathways <- enrichKEGG(gene = NULL, organism = "hsa", keyType = "kegg")

# Get ranked gene list (logFC and p-value)
ranked_genes <- de_genes[, "avg_log2FC"]
names(ranked_genes) <- rownames(de_genes)
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Extract Entrez IDs from the ranked genes (make sure your data contains these)
gene_list <- bitr(names(ranked_genes), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Filter Entrez IDs and use for KEGG enrichment
entrez_gene_list <- gene_list$ENTREZID



# Run KEGG enrichment
kegg_pathways <- enrichKEGG(gene = entrez_gene_list, organism = "hsa")

# Check the results of KEGG enrichment
head(kegg_pathways@result)



#Extract Pyrimidine Pathway Genes----------------------------------------------------------------------------------------
# Extract the gene set for Pyrimidine metabolism (KEGG ID: hsa00240)
pyrimidine_genes <- kegg_pathways@result[kegg_pathways@result$ID == "hsa00240", "geneID"]

# Split and convert the gene IDs to gene symbols
pyrimidine_gene_list <- bitr(strsplit(pyrimidine_genes, "/")[[1]], fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

# Check extracted genes
head(pyrimidine_gene_list)
# ENTREZID SYMBOL
#1    51733   UPB1
#2   129607  CMPK2
#3     1806   DPYD
#4     1890   TYMP
#5     4907   NT5E
#6      953 ENTPD1




#Run GESA------------------------------------------------------------------------------------------------
gsea_results <- GSEA(
  geneList = ranked_genes,  # Ranked gene list from FindMarkers
  TERM2GENE = pyrimidine_gene_list,  
  pvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500,
  verbose = TRUE
)

# Check GSEA results
head(gsea_results@result)
#[1] ID              Description     setSize         enrichmentScore NES             pvalue         
#[7] p.adjust        qvalue         
#<0 rows> (or 0-length row.names)


##That means no significant pathways were identified by GSEA based on the pyrimidine pathway gene set.




#Save processed data------------------------------------------------------------
saveRDS(seurat_obj, file = "spatial_breast_cancer_seurat.rds")




























