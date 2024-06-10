data(iris)
iris_data <- iris[, -5]
dist_matrix <- dist(iris_data)
# Perform hierarchical clustering using complete linkage
hclust_complete <- hclust(dist_matrix, method = "complete")
# Plot the dendrogram
plot(hclust_complete, labels = iris$Species, main = "Hierarchical Clustering Dendrogram (Complete Linkage)", xlab = "", sub = "", ylab = "Height")
# Cut the dendrogram to form 3 clusters
clusters <- cutree(hclust_complete, k = 3)

# Add the cluster assignment to the original data
iris$Cluster <- as.factor(clusters)

# Print the first few rows to see the cluster assignments
head(iris)
# Install and load ggplot2 if not already installed
install.packages("ggplot2")
library(ggplot2)
library("uwot")
library(factoextra)


# Plot the clusters
ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Cluster)) +
  geom_point(size = 2) +
  labs(title = "Hierarchical Clustering of Iris Data", x = "Sepal Length", y = "Sepal Width") +
  theme_minimal()

# Set a random seed for reproducibility
set.seed(42)

# Perform UMAP
umap_results <- umap(iris_data, n_neighbors = 15, min_dist = 0.1, n_components = 2)

# Convert UMAP results to a data frame
umap_df <- as.data.frame(umap_results)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Species <- iris$Species

# Plot the UMAP results
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Species)) +
  geom_point(size = 2) +
  labs(title = "UMAP Projection of the Iris Dataset",
       x = "UMAP1", y = "UMAP2") +
  theme_minimal()

kmeans_result <- kmeans(iris_data, centers = 3, nstart = 25)

# Add the cluster assignment to the original data
iris$Cluster <- as.factor(kmeans_result$cluster)

# Compute silhouette information
silhouette_info <- silhouette(kmeans_result$cluster, dist(iris_data))
fviz_silhouette(silhouette_info)
