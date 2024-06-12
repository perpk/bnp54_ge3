# Generate a synthetic gene expression dataset
set.seed(42)
gene_expression <- matrix(rnorm(1000), nrow = 50, ncol = 20)
rownames(gene_expression) <- paste0("Gene", 1:50)
colnames(gene_expression) <- paste0("Sample", 1:20)

# Perform variance thresholding to select features (genes) with the highest variance
variances <- apply(gene_expression, 1, var)
selected_genes <- names(sort(variances, decreasing = TRUE)[1:30])  # Select top 30 genes with highest variance
gene_expression_selected <- gene_expression[selected_genes,]

# Transpose the data to have samples as rows and selected genes as columns
gene_expression_selected_t <- t(gene_expression_selected)

# Install and load necessary packages
install.packages("umap")
install.packages("ggplot2")
library(umap)
library(ggplot2)

# Perform UMAP
umap_result <- umap(gene_expression_selected_t)

# Convert the result to a data frame for visualization
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Sample <- rownames(gene_expression_selected_t)

# Plot UMAP results
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, label = Sample)) +
  geom_point() +
  geom_text(aes(label = Sample), hjust = 1.5, vjust = 1.5, size = 3) +
  theme_minimal() +
  ggtitle("UMAP of Selected Gene Expression Data")
