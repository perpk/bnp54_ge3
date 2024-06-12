# Generate a synthetic gene expression dataset
set.seed(42)
gene_expression <- matrix(rnorm(1000), nrow = 50, ncol = 20)
rownames(gene_expression) <- paste0("Gene", 1:50)
colnames(gene_expression) <- paste0("Sample", 1:20)

# View the first few rows of the dataset
head(gene_expression)

# Install necessary packages
install.packages("umap")
install.packages("ggplot2")
install.packages("cluster")

# Load the libraries
library(umap)
library(ggplot2)
library(cluster)

# Transpose the data to have samples as rows and genes as columns
gene_expression_t <- t(gene_expression)

# Perform UMAP
umap_result <- umap(gene_expression_t)

# Convert the result to a data frame for visualization
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Sample <- rownames(gene_expression_t)

# Compute the distance matrix on UMAP data
dist_matrix_umap <- dist(umap_df[, 1:2])

# Perform hierarchical clustering using complete linkage
hclust_umap <- hclust(dist_matrix_umap, method = "complete")

# Plot the dendrogram
plot(hclust_umap, main = "Dendrogram of UMAP-Reduced Gene Expression Data", xlab = "Samples", sub = "", cex = 0.9)

# Cut the dendrogram to form clusters (e.g., 4 clusters)
clusters <- cutree(hclust_umap, k = 4)

# Add cluster labels to the UMAP data frame
umap_df$Cluster <- as.factor(clusters)

# Plot UMAP results with clusters
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster, label = Sample)) +
  geom_point(size = 3) +
  geom_text(aes(label = Sample), hjust = 1.5, vjust = 1.5, size = 3) +
  theme_minimal() +
  ggtitle("UMAP of Gene Expression Data with Clusters")

