# Generate a synthetic gene expression dataset
set.seed(42)
gene_expression <- matrix(rnorm(1000), nrow = 50, ncol = 20)
rownames(gene_expression) <- paste0("Gene", 1:50)
colnames(gene_expression) <- paste0("Sample", 1:20)

# View the first few rows of the dataset
## genes are rows, samples columns
head(gene_expression)

# Compute the distance matrix
## after transpose genes are columns, samples rows
dist_matrix <- dist(t(gene_expression))  # Transpose to cluster samples, not genes

# Perform hierarchical clustering using complete linkage
hclust_complete <- hclust(dist_matrix, method = "complete")

# Cut the dendrogram to form 4 clusters
clusters <- cutree(hclust_complete, k = 4)

# Add the cluster labels to the samples
sample_clusters <- data.frame(Sample = colnames(gene_expression), Cluster = clusters)
print(sample_clusters)

# Install and load the cluster package for silhouette analysis
install.packages("cluster")
library(cluster)

# Calculate silhouette widths
silhouette_widths <- silhouette(clusters, dist_matrix)

# Plot the silhouette plot
plot(silhouette_widths, main = "Silhouette Plot for Hierarchical Clustering")

## Silhouette width measures how similar an object is to its own cluster
### compared to other clusters

## Si values close to 1 indicate that an object is well-clustered
## Si values close to 0 indicate that an object is on the boundary between 2 clusters
## Si values close to -1 indicate that an object might be assigned to the wrong cluster
