
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
setwd("F:/GITHUB_1/Transcriptomics/3_Spatial/files_spatial")

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(SeuratData)
library(patchwork)
library(tiff)
library(png)





#install.packages("Seurat")
#ls("package:Seurat")  # Lists all functions in Seurat

#install.packages("remotes")  # Install remotes if not already installed
#remotes::install_github("satijalab/seurat-data", force = TRUE)

#install.packages("tiff")



setwd("F:/GITHUB_1/Transcriptomics/3_Spatial/files_spatial")
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



# Spatial Visualization-----------------------------------------------------
SpatialDimPlot( seurat_obj,label = TRUE,label.size = 3)



# Marker genes---------------------------------------------------------------
marker_genes <- c("MS4A1","EPCAM","CD3D","PTPRC","CD19")
SpatialFeaturePlot(seurat_obj,features = marker_genes)


# Single gene------------------------------------------------------------------
SpatialFeaturePlot(seurat_obj,features = "MS4A1")


# Differential Expression Analysis---------------------------------------------
# Cluster 0 vs Cluster 1

levels(seurat_obj)

de_genes <- FindMarkers(seurat_obj,ident.1 = "0",ident.2 = "1",logfc.threshold = 0.25,min.pct = 0.25)

# Add gene names
de_genes$gene <- rownames(de_genes)

# View results
head(de_genes)
#              p_val avg_log2FC pct.1 pct.2    p_val_adj   gene
#CANT1  4.731465e-101 -1.2944430     1     1 7.103348e-97  CANT1
#DDX5    8.718826e-99 -1.2988469     1     1 1.308957e-94   DDX5
#VIM     1.041647e-96  1.7565626     1     1 1.563824e-92    VIM
#PSMD3   3.756907e-96 -0.9890809     1     1 5.640244e-92  PSMD3
#ORMDL3  4.049394e-96 -1.1202142     1     1 6.079356e-92 ORMDL3
#ERBB2   1.171241e-95 -1.2484934     1     1 1.758384e-91  ERBB2
#> 
  



# Volcano Plot-------------------------------------------------------------------

de_genes$Significance <- ifelse(de_genes$p_val_adj < 0.05 & de_genes$avg_log2FC > 0,"Upregulated",
  
  ifelse(de_genes$p_val_adj < 0.05 & de_genes$avg_log2FC < 0,"Downregulated","Not Significant")
)

ggplot(de_genes,aes(x = avg_log2FC,y = -log10(p_val_adj),color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +scale_color_manual(values = c("Upregulated" = "red","Downregulated" = "blue","Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot",x = "Log2 Fold Change",y = "-log10 Adjusted P-value")



# Highlight Top Genes in Volcano Plot--------------------------------------------------------
top_up <- head(de_genes[order(-de_genes$avg_log2FC,de_genes$p_val_adj),],15)

top_down <- head(de_genes[order(de_genes$avg_log2FC,de_genes$p_val_adj),],15)

top_genes_labels <- c(rownames(top_up),rownames(top_down))

ggplot(de_genes,aes(x = avg_log2FC,y = -log10(p_val_adj),color = Significance)
) +
  geom_point(alpha = 0.6) +
  
  geom_text_repel(
    data = subset(
      de_genes,
      rownames(de_genes) %in% top_genes_labels
    ),
    aes(label = gene),
    size = 3,
    color = "black"
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Volcano Plot with Top Genes",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )



# Top 10 Genes--------------------------------------------------------------------

top_genes <- head(
  de_genes[
    order(de_genes$p_val_adj),
  ],
  10
)

top_gene_names <- top_genes$gene



# Dot Plot----------------------------------------------------------------------

DotPlot(
  seurat_obj,
  features = top_gene_names
) +
  RotatedAxis() +
  ggtitle("Top 10 Differentially Expressed Genes")


# Heatmap------------------------------------------------------------------------

DoHeatmap(
  seurat_obj,
  features = top_gene_names
) +
  ggtitle("Heatmap: Top 10 DE Genes")


# MA Plot-------------------------------------------------------------------------

ggplot(
  de_genes,
  aes(
    x = (pct.1 + pct.2) / 2,
    y = avg_log2FC
  )
) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "MA Plot",
    x = "Average Expression",
    y = "Log2 Fold Change"
  )



# DEG Statistics----------------------------------------------------------------

# Total genes
nrow(de_genes)
#[1] 6178

# Significant genes
sum(
  de_genes$p_val_adj < 0.05,
  na.rm = TRUE
)
#[1] 3914


# Upregulated genes
sum(
  de_genes$avg_log2FC > 0 &
    de_genes$p_val_adj < 0.05,
  na.rm = TRUE
)
#[1] 1043

# Downregulated genes
sum(
  de_genes$avg_log2FC < 0 &
    de_genes$p_val_adj < 0.05,
  na.rm = TRUE
)
# [1] 2871


# Filter Upregulated and Downregulated Genes

upregulated_genes <- subset(
  de_genes,
  avg_log2FC > 0 &
    p_val_adj < 0.05
)

downregulated_genes <- subset(
  de_genes,
  avg_log2FC < 0 &
    p_val_adj < 0.05
)


# Top 10 Upregulated Genes

top_10_up <- head(
  upregulated_genes[
    order(
      -upregulated_genes$avg_log2FC
    ),
  ],
  10
)


# Top 10 Downregulated Genes------------------------------------------------------

top_10_down <- head(
  downregulated_genes[
    order(
      downregulated_genes$avg_log2FC
    ),
  ],
  10
)

# Save CSV Files----------------------------------------------------------------

write.csv(
  upregulated_genes,
  "upregulated_genes.csv",
  row.names = TRUE
)

write.csv(
  downregulated_genes,
  "downregulated_genes.csv",
  row.names = TRUE
)

write.csv(
  top_10_up,
  "top_10_upregulated_genes.csv",
  row.names = TRUE
)

write.csv(
  top_10_down,
  "top_10_downregulated_genes.csv",
  row.names = TRUE
)


# Bar Plot-----------------------------------------------------------------------

top_genes_combined <- rbind(
  
  data.frame(
    Gene = rownames(top_10_up),
    Log2FC = top_10_up$avg_log2FC,
    Direction = "Upregulated"
  ),
  
  data.frame(
    Gene = rownames(top_10_down),
    Log2FC = top_10_down$avg_log2FC,
    Direction = "Downregulated"
  )
)

ggplot(
  top_genes_combined,
  aes(
    x = reorder(Gene, Log2FC),
    y = Log2FC,
    fill = Direction
  )
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top Upregulated and Downregulated Genes",
    x = "Gene",
    y = "Log2 Fold Change"
  )


# Dot Plot for Top Genes--------------------------------------------------------

DotPlot(
  seurat_obj,
  features = c(
    rownames(top_10_up),
    rownames(top_10_down)
  )
) +
  RotatedAxis() +
  ggtitle("Top Upregulated and Downregulated Genes")



# Feature Plot-------------------------------------------------------------------

FeaturePlot(
  seurat_obj,
  features = c(
    "ERBB2",
    "COL6A3",
    "VIM",
    "GRB7",
    "NDUFB9",
    "FOXD1"
  )
)


# GSEA Analysis-------------------------------------------------------------------

# Create ranked gene list
ranked_genes <- de_genes$avg_log2FC

names(ranked_genes) <- rownames(de_genes)

ranked_genes <- sort(
  ranked_genes,
  decreasing = TRUE
)

# Convert gene symbols to ENTREZ IDs

gene_conversion <- bitr(
  names(ranked_genes),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Merge with ranked genes

ranked_df <- merge(
  gene_conversion,
  data.frame(
    SYMBOL = names(ranked_genes),
    logFC = ranked_genes
  ),
  by = "SYMBOL"
)

geneList <- ranked_df$logFC

names(geneList) <- ranked_df$ENTREZID

geneList <- sort(
  geneList,
  decreasing = TRUE
)

# Run KEGG GSEA

gsea_results <- gseKEGG(
  geneList = geneList,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# View results

head(gsea_results@result)


# save GSEA results
gsea_df <- as.data.frame(gsea_results@result)

write.csv(gsea_df, "gsea_kegg_results.csv", row.names = FALSE)


# Save Seurat Object---------------------------------------------------------------

saveRDS(
  seurat_obj,
  file = "spatial_breast_cancer_seurat.rds"
)

# reload it:
seurat_obj <- readRDS("spatial_breast_cancer_seurat.rds")