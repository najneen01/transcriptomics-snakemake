 # Single-Cell Transcriptomic Analysis using Seurat (R)

## Project Overview

This project performs single-cell RNA sequencing (scRNA-seq) analysis using the Seurat package in R. The analysis is based on publicly available data from GEO.

## Data Source
Dataset: GSE275330
Source: NCBI Gene Expression Omnibus (GEO)
Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275330

## Project Structure
data/: Contains the raw 10X Genomics input files (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz).

script/: Contains single_snake_R.R, the core R script for Seurat analysis.

result/: Output directory (automatically managed by Snakemake). Contains the final RDS object, marker lists, and visualizations.

Snakefile: The workflow management file that defines the rules and dependencies.


## Workflow Overview:
The seurat_analysis rule executes the following steps:
### Data Ingestion: 
Loads sparse matrices and metadata.

### Quality Control: 
Filters cells based on feature counts ($500 < n < 5000$) and mitochondrial content ($< 5\%$).

### Normalization & Scaling: 
Standardizes the data and identifies the top 3000 variable features.

### Dimensional Reduction: 
Performs PCA, followed by UMAP and t-SNE for visualization.

### Clustering: 
Identifies cell clusters using SNN graphs.

### Marker Discovery: 
Calculates differentially expressed genes for every cluster.



## Getting Started
Prerequisites
Ensure you have R and Snakemake installed. The following R packages are required:

Seurat

dplyr

Matrix

ggplot2

## Running the Pipeline
To execute the pipeline using 4 CPU cores, run the following command in your terminal:

                 snakemake --cores 4




## Outputs
After a successful run, the result/ folder will contain:

single_snake_R.rds: The fully processed Seurat object (usable for downstream analysis).

markers.csv: A list of all cluster markers with p-values and fold-change logic.

umap_tsne.png: Side-by-side comparison of UMAP and t-SNE projections.

heatmap.png: Expression heatmap of the top 10 markers per cluster.

featureplot.png: Spatial expression plots for key genes (PRKCA, ITCH, ABL2, ALCAM).


## Troubleshooting
Deleted Results: If the script fails, Snakemake automatically deletes the result/ files to prevent downstream errors from corrupted data. Ensure all paths in the R script match the Snakefile.

Memory: For larger datasets, consider increasing the --cores or ensuring your system has sufficient RAM for ScaleData operations.






