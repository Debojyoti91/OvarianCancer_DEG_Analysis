# OvarianCancer_DEG_Analysis
Analysis of differential gene expression in ovarian cancer cells co-cultured under various conditions using DESeq2. Focuses on interactions between platinum-resistant and platinum-sensitive cells, highlighting the role of E2F1.


## Setup Instructions

### Prerequisites
- R (version 4.0 or later)
- RStudio (recommended)
- The following R packages:
  - DESeq2
  - AnnotationDbi
  - org.Hs.eg.db
  - ggplot2
  - dplyr
  - pheatmap
  - ComplexHeatmap
  - ggrepel
 

## Usage Instructions

### Prepare the Data
Place your `combined_counts.csv` and `metadata.csv` files in the `data/` directory.

### Run the Analysis
1. Open the `differential_expression_analysis.R` script in RStudio.
2. Set your working directory to the project folder, e.g., `setwd("/path/to/your_project/")`.
3. Source the script or run it line by line to perform the analysis.

### View Results
- Results will be saved in the `results/` directory.
- Generated plots will be saved in the `figures/` directory.

### Script Details

The `differential_expression_analysis.R` script performs the following tasks:
- Loads necessary libraries.
- Reads count data (`combined_counts.csv`) and metadata (`metadata.csv`).
- Prepares the DESeq2 dataset.
- Filters low-count genes.
- Runs DESeq2 analysis to identify differentially expressed genes (DEGs).
- Maps Ensembl IDs to gene names.
- Saves DEG results to CSV files.
- Generates PCA, MA, and Volcano plots.
- Creates a heatmap for the top 30 DEGs




### Installing Required Packages
```r
install.packages(c("ggplot2", "dplyr", "pheatmap", "ggrepel"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "AnnotationDbi", "org.Hs.eg.db", "ComplexHeatmap"))



