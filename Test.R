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

# Install and load required packages
install.packages("dendextend")
library(dendextend)

# Load the iris dataset
data(iris)
iris_data <- iris[, -5]  # Remove the Species column for clustering

# Compute the distance matrix
dist_matrix <- dist(iris_data)

# Perform hierarchical clustering using complete linkage
hclust_complete <- hclust(dist_matrix, method = "complete")

# Convert hclust object to a dendrogram object
dend <- as.dendrogram(hclust_complete)

# Color branches by cluster
dend <- color_branches(dend, k = 3)

# Plot the colored dendrogram
plot(dend, main = "Colored Dendrogram", xlab = "", sub = "", ylab = "Height")

# Add rectangles around clusters
rect.dendrogram(dend, k = 3, border = "red")

#########

# Load the iris dataset
data(iris)

# Load necessary libraries
library(caret)
library(rpart)
library(rpart.plot)

# Set seed for reproducibility
set.seed(42)

# Split the data into training (70%) and test (30%) sets
train_index <- createDataPartition(iris$Species, p = 0.7, list = FALSE)
train_data <- iris[train_index, ]
test_data <- iris[-train_index, ]

# Train a decision tree model
model <- rpart(Species ~ ., data = train_data, method = "class")

# Print the model summary
print(summary(model))

# Plot the decision tree
rpart.plot(model, type = 2, extra = 104, fallen.leaves = TRUE, main = "Decision Tree for Iris Dataset")

# Make predictions on the test data
predictions <- predict(model, test_data, type = "class")

# Create the confusion matrix
conf_matrix <- confusionMatrix(predictions, test_data$Species)
print(conf_matrix)

# Example of making a single prediction
new_data <- data.frame(Sepal.Length = 5.1, Sepal.Width = 3.5, Petal.Length = 1.4, Petal.Width = 0.2)
prediction <- predict(model, new_data, type = "class")
print(prediction)

