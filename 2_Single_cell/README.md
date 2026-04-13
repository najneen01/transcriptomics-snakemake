 # Single-Cell Transcriptomic Analysis using Seurat (R)

## Project Overview

This project performs single-cell RNA sequencing (scRNA-seq) analysis using the Seurat package in R. The workflow includes data preprocessing, quality control, dimensionality reduction, clustering, and marker visualization.

The analysis is based on publicly available data from GEO.

## Data Source
Dataset: GSE275330
Source: NCBI Gene Expression Omnibus (GEO)
Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275330

## Tools & Packages

The analysis was conducted in RStudio using the following packages:

Seurat
dplyr
Matrix
ggplot2
patchwork

## Workflow Summary
### 1. Data Loading
Imported matrix, barcodes, and features files
Constructed Seurat object
### 2. Quality Control (QC)
Filtered cells based on:
Number of features (genes)
Total counts
Mitochondrial gene percentage
### 3. Normalization & Feature Selection
Log normalization
Identification of highly variable genes
### 4. Dimensionality Reduction
Principal Component Analysis (PCA)
t-SNE and UMAP visualization
### 5. Clustering
Cell clustering based on selected principal components
### 6. Marker Analysis
Visualization of:
Neural markers (e.g., MKI67, NES, DCX)
Lymphoma markers (DLBCL, MCL)


## Results & Visualizations
Figure 1: QC Violin Plot
p1 <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("figures/Fig1_QC_violin.png", p1, width = 10, height = 5)

Quality control of single-cell RNA-seq data.
Violin plots showing the distribution of the number of detected genes (nFeature_RNA), total RNA counts (nCount_RNA), and mitochondrial gene percentage (percent.mt) across all cells.



Figure 2. Quality control scatter plots.
p2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt") +
      FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave("figures/Fig2_QC_scatter.png", p2, width = 10, height = 5)

Scatter plots illustrating relationships between total RNA counts and mitochondrial gene percentage, and between total RNA counts and number of detected genes, used to identify outliers and low-quality cells.



Figure 3. Identification of highly variable genes.
p3 <- VariableFeaturePlot(seurat)
ggsave("figures/Fig3_variable_features.png", p3, width = 6, height = 5)

Variable feature plot displaying genes with high cell-to-cell variation.





Figure 4. Top variable genes.
p4 <- LabelPoints(plot = p3, points = head(VariableFeatures(seurat), 20), repel = TRUE)
ggsave("figures/Fig4_top_variable_features.png", p4, width = 6, height = 5)

Top 20 highly variable genes highlighted on the variable feature plot.




Figure 5. Principal component analysis (PCA).
p5 <- DimPlot(seurat, reduction = "pca", dims = c(1, 2))
ggsave("figures/Fig5_PCA.png", p5, width = 6, height = 5)

Projection of cells onto the first two principal components, summarizing major sources of variation.



Figure 6. Determination of significant principal components.
p6 <- ElbowPlot(seurat)
ggsave("figures/Fig6_elbow.png", p6, width = 6, height = 5)

Elbow plot showing variance explained by each principal component to guide dimensionality selection.



Figure 7. Principal component heatmap.
png("figures/Fig7_PC_heatmap.png", width = 800, height = 600)
PCHeatmap(seurat, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

Heatmap displaying gene contributions across selected principal components.



Figure 8. Visualization of principal component structure.
p8 <- PCHeatmap(seurat, dims = 1:6, cells = 500, balanced = TRUE)
ggsave("figures/Fig8_PC_heatmap_layout.png", p8, width = 10, height = 6)

Arranged heatmaps showing gene loading patterns across top principal components.



Figure 9. Cell clustering using t-SNE and UMAP.
p9 <- TSNEPlot(seurat) + UMAPPlot(seurat)
ggsave("figures/Fig9_TSNE_UMAP.png", p9, width = 10, height = 5)

Low-dimensional visualization of cells using t-SNE and UMAP, illustrating clustering structure.




Figure 10. Neural marker gene expression.
p10 <- FeaturePlot(seurat, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"), ncol=3)
ggsave("figures/Fig10_neural_markers.png", p10, width = 12, height = 8)

Feature plots showing expression of neural-associated genes across cells.




Figure 11. Diffuse large B-cell lymphoma (DLBCL) markers.
p11 <- FeaturePlot(seurat, c("CD19","CD20","CD22","CD79a","BCL6","MYC","MUM1"), ncol=3)
ggsave("figures/Fig11_DLBCL_markers.png", p11, width = 12, height = 8)

Expression patterns of key B-cell lymphoma markers across cell populations.




Figure 12. Mantle cell lymphoma (MCL) markers.
p12 <- FeaturePlot(seurat, c("CD5","Cyclin D1","SOX11"), ncol=3)
ggsave("figures/Fig12_MCL_markers.png", p12, width = 10, height = 6)

Feature plots illustrating expression of MCL-associated genes.



Figure 13. Cluster annotation.
p13 <- DimPlot(seurat, reduction = "tsne", label = TRUE) +
       DimPlot(seurat, reduction = "umap", label = TRUE)

ggsave("figures/Fig13_clusters.png", p13, width = 10, height = 5)

t-SNE and UMAP plots with cluster labels indicating distinct cell populations.



Figure 14. Heatmap of marker genes.
png("figures/Fig14_heatmap.png", width = 800, height = 600)
DoHeatmap(seurat, features = cl_markers$gene) + NoLegend()
dev.off()

Heatmap showing expression of selected marker genes across clusters.


