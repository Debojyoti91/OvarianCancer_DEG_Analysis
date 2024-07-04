# Load necessary libraries
setwd("/Users/debojyotidas/Dropbox/My Mac (Debojyotis-MacBook-Pro.local)/Desktop/NGS/bioinf_DEG/r_deg")
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)
library(ggrepel)

# Load the count data
cnt <- read.csv("data/combined_counts.csv")

# Load the metadata
met <- read.csv("data/metadata.csv")

# Set GeneID as row names for the count data
rownames(cnt) <- cnt$GeneID
cnt <- cnt[, -1]  # Remove the GeneID column

# Set Sample.ID as row names for the metadata
rownames(met) <- met$Sample.ID
met <- met[, -1]  # Remove the Sample.ID column

# Ensure that the sample names match
stopifnot(all(colnames(cnt) %in% rownames(met)))  # This should return TRUE

# Create DESeq2 dataset with the design formula ~ Condition
dds <- DESeqDataSetFromMatrix(
    countData = cnt,
    colData = met,
    design = ~ Condition
)

# Filter out low-count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run the DESeq2 analysis
dds <- DESeq(dds)

# Get the results for the default comparison with adjusted p-value threshold of 0.05
res <- results(dds, alpha = 0.05)

# Map Ensembl IDs to gene names for default comparison results
ensembl_ids <- rownames(res)
gene_names <- mapIds(org.Hs.eg.db, keys=ensembl_ids, column="SYMBOL", keytype="ENSEMBL", multiVals="first")

# Add gene names to the DESeq2 results
res$gene_name <- gene_names

# Save the default comparison results to a CSV file
write.csv(as.data.frame(res), "results/DESeq2_results_with_gene_names.csv")

# Get the results for a specific comparison with adjusted p-value threshold of 0.05
res_specific <- results(dds, contrast = c("Condition", "Co-culture with PEO4", "Monoculture"), alpha = 0.05)

# Map Ensembl IDs to gene names for specific comparison results
ensembl_ids_specific <- rownames(res_specific)
gene_names_specific <- mapIds(org.Hs.eg.db, keys=ensembl_ids_specific, column="SYMBOL", keytype="ENSEMBL", multiVals="first")

# Add gene names to the specific comparison results
res_specific$gene_name <- gene_names_specific

# Save the specific comparison results to a CSV file
write.csv(as.data.frame(res_specific), "results/DESeq2_specific_results_with_gene_names.csv")

# Perform variance stabilizing transformation (VST)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# Plot PCA using VST data
pcaData_vsd <- plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
percentVar_vsd <- round(100 * attr(pcaData_vsd, "percentVar"))
p <- ggplot(pcaData_vsd, aes(PC1, PC2, color = Condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar_vsd[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar_vsd[2], "% variance")) +
    coord_fixed()
ggsave("figures/pca_plot_vsd.png", plot = p)

# Perform regularized log transformation (rlog) - Optional
rld <- rlogTransformation(dds, blind = FALSE)

# Plot PCA using rlog data - Optional
pcaData_rld <- plotPCA(rld, intgroup = "Condition", returnData = TRUE)
percentVar_rld <- round(100 * attr(pcaData_rld, "percentVar"))
p <- ggplot(pcaData_rld, aes(PC1, PC2, color = Condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar_rld[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar_rld[2], "% variance")) +
    coord_fixed()
ggsave("figures/pca_plot_rld.png", plot = p)

# Create MA plot
ma_data <- as.data.frame(res)
ma_data <- ma_data %>% mutate(significance = ifelse(padj < 0.05, "Significant", "Not Significant"))
p <- ggplot(ma_data, aes(x=baseMean, y=log2FoldChange, color=significance)) +
    geom_point(alpha=0.5) +
    scale_x_log10() +
    scale_color_manual(values=c("red", "black")) +
    theme_minimal() +
    labs(title="MA Plot", x="Mean of normalized counts", y="Log2 Fold Change")
ggsave("figures/ma_plot.png", plot = p)

# Save the top 10 best genes based on adjusted p-value
res_df <- as.data.frame(res)
top_genes <- res_df %>% arrange(padj) %>% head(10)
write_csv(top_genes, "results/Top_10_best_genes.csv")

# Create Volcano plot with gene labels
volcano_data <- res_df %>%
    mutate(significance = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "Not Significant"))

# Label the top 10 significant genes
top_genes <- volcano_data %>% filter(significance == "Significant") %>% arrange(padj) %>% head(10)

p <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("red", "black")) +
    theme_minimal() +
    labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
    geom_text_repel(data = top_genes, aes(label = gene_name), max.overlaps = 10)
ggsave("figures/volcano_plot.png", plot = p)

# Identify the top 30 genes based on adjusted p-value
top_30_genes <- res_df %>%
  arrange(padj) %>%
  head(30)

# Extract the gene IDs, ensuring they are not NA and exist in the normalized counts matrix
top_30_gene_ids <- top_30_genes$gene_name
top_30_gene_ids <- top_30_gene_ids[!is.na(top_30_gene_ids) & top_30_gene_ids %in% rownames(counts(dds, normalized=TRUE))]

# Extract the normalized expression data for these genes
normalized_counts <- counts(dds, normalized = TRUE)
top_30_normalized_counts <- normalized_counts[top_30_gene_ids, ]

# Generate the heatmap
Heatmap(t(scale(t(top_30_normalized_counts))), 
        name = "Z-score",
        cluster_rows = TRUE, 
        cluster_columns = TRUE, 
        show_row_names = TRUE, 
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 10),  # Adjust row names size if necessary
        column_names_gp = gpar(fontsize = 10)  # Adjust column names size if necessary
)
ggsave("figures/heatmap.png")

