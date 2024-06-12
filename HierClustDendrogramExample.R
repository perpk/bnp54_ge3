# Generate a synthetic gene expression dataset
set.seed(131)
gene_expression <- matrix(rnorm(1000), nrow = 131, ncol = 600)
rownames(gene_expression) <- paste0("Gene", 1:131)
colnames(gene_expression) <- paste0("Sample", 1:600)

# View the first few rows of the dataset
head(gene_expression)

# Compute the distance matrix
dist_matrix <- dist(t(gene_expression))  # Transpose to cluster samples

# Perform hierarchical clustering using complete linkage
hclust_complete <- hclust(dist_matrix, method = "complete")

# Plot the dendrogram
plot(hclust_complete, labels=FALSE, main = "Dendrogram of Gene Expression Data", xlab = "Samples", sub = "", cex = 0.9)
