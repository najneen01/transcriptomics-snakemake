Project Title

Spatial Transcriptomic Analysis of Human Breast Cancer (10x Visium)

📖 Overview

This project performs spatial transcriptomic analysis of human breast cancer tissue using 10x Genomics Visium data and the Seurat (R) pipeline.

The workflow includes:

Data loading and preprocessing
Dimensionality reduction and clustering
Spatial visualization
Differential gene expression analysis
Functional enrichment (GSEA)
📂 Dataset
Source: 10x Genomics
Dataset: Fresh Frozen Human Breast Cancer
Platform: Visium Spatial Gene Expression


### 1. Data Ingestion & Object Initialization
Loading the filtered H5 matrix and spatial imaging metadata into the Seurat environment.

seurat_obj <- Load10X_Spatial(
  data.dir = "F:/GITHUB_1/Transcriptomics/3_Spatial/files_spatial",
  filename = "CytAssist_Fresh_Frozen_Human_Breast_Cancer_filtered_feature_bc_matrix.h5",
  slice = "breast_tissue"
)

Description: Initializes the Seurat object by integrating the gene expression matrix with histological image coordinates.

### 2. Dimensionality Reduction (UMAP)

![umap](https://raw.githubusercontent.com/najneen01/transcriptomics-snakemake/main/3_Spatial/result/umap.jpeg)

Visualization of transcriptomic clusters in a 2D UMAP space.

DimPlot(seurat_obj, reduction = "umap", label = TRUE)


Description: Displays groups of spots with similar transcriptional profiles, identifying distinct cell populations within the breast cancer microenvironment

### 3. Spatial Tissue Mapping

![Spatial_Tissue_Mapping](https://raw.githubusercontent.com/najneen01/transcriptomics-snakemake/main/3_Spatial/result/Spatial_Tissue_Mapping.jpeg)

Projecting the identified clusters back onto the H&E stained tissue image.

SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3)

Description: Maps the cluster identities directly onto the tissue morphology, allowing for the observation of the spatial distribution of tumor and stroma.



### 5. Spatial Marker Gene Distribution

![spatial_marker_genes](https://raw.githubusercontent.com/najneen01/transcriptomics-snakemake/main/3_Spatial/result/spatial_marker_genes.jpeg)

Visualizing specific biological markers across the tissue section.

marker_genes <- c("MS4A1","EPCAM","CD3D","PTPRC","CD19")
SpatialFeaturePlot(seurat_obj,features = marker_genes)

Description: Shows the local expression intensity of key markers. For example, EPCAM helps identify epithelial tumor regions, while CD3D highlights immune infiltration.



### 4. Differential Expression & Volcano Plot

![volcano_plot](https://raw.githubusercontent.com/najneen01/transcriptomics-snakemake/main/3_Spatial/result/volcano_plot.jpeg)

Identifying significant biomarkers using FindMarkers and visualizing with a Volcano Plot.


de_genes$Significance <- ifelse(de_genes$p_val_adj < 0.05 & de_genes$avg_log2FC > 0,"Upregulated",
  
  ifelse(de_genes$p_val_adj < 0.05 & de_genes$avg_log2FC < 0,"Downregulated","Not Significant")
)


ggplot(de_genes,aes(x = avg_log2FC,y = -log10(p_val_adj),color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +scale_color_manual(values = c("Upregulated" = "red","Downregulated" = "blue","Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot",x = "Log2 Fold Change",y = "-log10 Adjusted P-value")

  

Description: Statistical summary of differentially expressed genes (DEGs). Highlights significant upregulated genes like VIM and ERBB2 across clusters.




### 7. Expression Profiling (Dot Plot- Top DEGs)

![dot_plot](https://raw.githubusercontent.com/najneen01/transcriptomics-snakemake/main/3_Spatial/result/dot_plot.jpeg)

A comparative view of the top 10 upregulated and downregulated genes.


DotPlot(
  seurat_obj,
  features = top_gene_names
) +
  RotatedAxis() +
  ggtitle("Top 10 Differentially Expressed Genes")


Description: Displays the average expression intensity and the percentage of spots expressing specific genes across identified clusters.










### 6. Expression Heatmap (Top DEGs)

![heatmap](https://raw.githubusercontent.com/najneen01/transcriptomics-snakemake/main/3_Spatial/result/heatmap.jpeg)


A high-resolution visualization of the expression patterns for the top 10 differentially expressed genes across every spot in the tissue.


DoHeatmap(
  seurat_obj,
  features = top_gene_names
) +
  ggtitle("Heatmap: Top 10 DE Genes")


Description: The heatmap provides a granular view of gene expression, where columns represent individual spatial spots and rows represent genes. This allows for a clear visual validation of how specifically the top markers (like CANT1, VIM, and ERBB2) are expressed within their respective clusters.







### 8. MA Plot (Log Fold Change vs. Expression)

![MA_plot](https://raw.githubusercontent.com/najneen01/transcriptomics-snakemake/main/3_Spatial/result/MA_plot.jpeg)

Visualizing the relationship between the intensity of gene expression and the magnitude of change between clusters.

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


Description: The MA plot helps detect biases in differential expression. It confirms that the identified DEGs are significant across various expression levels, ensuring that findings are not restricted only to highly expressed genes.


### 9. Top DEGs Bar Plot

A side-by-side comparison of the top 10 upregulated and top 10 downregulated genes.

![bar_plot](https://raw.githubusercontent.com/najneen01/transcriptomics-snakemake/main/3_Spatial/result/bar_plot.jpeg)


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

Description: This bar chart ranks the genes with the most dramatic expression changes. By showing both directions (Up/Down), it provides a clear summary of the transcriptional shifts occurring between the compared tissue regions.


### 10. Integrated Top Genes Dot Plot
A final comprehensive look at the expression frequency and intensity of all top markers.

![integrated_dot_plot](https://raw.githubusercontent.com/najneen01/transcriptomics-snakemake/main/3_Spatial/result/integrated_dot_plot.jpeg)


DotPlot(
  seurat_obj,
  features = c(
    rownames(top_10_up),
    rownames(top_10_down)
  )
) +
  RotatedAxis() +
  ggtitle("Top Upregulated and Downregulated Genes")



Description: This combined dot plot serves as a signature map for the clusters, showing how the top 20 genes (10 up, 10 down) define the transcriptional identity of each spatial zone.




### 11. Feature Plot (Gene Expression Density)
Visualizing the expression levels of key breast cancer drivers and structural markers across the UMAP projection.

![feature_plot](https://raw.githubusercontent.com/najneen01/transcriptomics-snakemake/main/3_Spatial/result/featuret_plot.jpeg)


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

Description: This plot identifies the density of specific gene transcripts within the UMAP clusters. It highlights markers like ERBB2 (HER2) and VIM (Vimentin), which are critical for understanding tumor subtype and epithelial-mesenchymal transition (EMT).



### 12. Pathway Enrichment (KEGG GSEA)
Functional interpretation of the gene expression changes using Gene Set Enrichment Analysis.

gsea_results <- gseKEGG(
  geneList = geneList,
  organism = "hsa",
  pvalueCutoff = 0.05
)


Description: Identifies enriched biological pathways (e.g., cell cycle or metabolic pathways) associated with the differentially expressed genes in the tumor tissue.


[View CSV](https://raw.githubusercontent.com/najneen01/transcriptomics-snakemake/main/3_Spatial/result/gsea_kegg_results.csv)


### Data Artifacts
Matrix: filtered_feature_bc_matrix.h5

Output: spatial_breast_cancer_seurat.rds (Processed Seurat Object)

Tables: upregulated_genes.csv, downregulated_genes.csv
