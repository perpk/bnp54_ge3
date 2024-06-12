# Install and load required packages
install.packages("dendextend")
library(dendextend)

# Generate a synthetic gene expression dataset
set.seed(42)
gene_expression <- matrix(rnorm(1000), nrow = 50, ncol = 20)
rownames(gene_expression) <- paste0("Gene", 1:50)
colnames(gene_expression) <- paste0("Sample", 1:20)

# View the first few rows of the dataset
head(gene_expression)

# Compute the distance matrix
dist_matrix <- dist(t(gene_expression))  # Transpose to cluster samples, not genes

# Perform hierarchical clustering using complete linkage
hclust_complete <- hclust(dist_matrix, method = "complete")

# Convert hclust object to a dendrogram object
dend <- as.dendrogram(hclust_complete)

# Color branches by cluster
dend <- color_branches(dend, k = 4)

dend %>% set("labels_cex", 0.3) %>% plot
# Plot the colored dendrogram
#plot(dend, main = "Colored Dendrogram", xlab = "", sub = "", ylab = "Height")

# Add rectangles around clusters
rect.dendrogram(dend, k = 4, border = "red")
