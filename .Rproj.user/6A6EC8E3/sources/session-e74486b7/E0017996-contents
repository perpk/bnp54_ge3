install.packages("umap")
install.packages("plotly")
install.packages("ggdendro")
install.packages("factoextra")
install.packages("uwot")
library("dplyr")
library("ggplot2")
library("stats")
library("plotly")
library("uwot")
library("factoextra")
library("stats")
library("ggdendro")
library("cluster")

# Topic 1.
# 1. Read data into dataframes
expression.data.df <- as.data.frame(read.csv("expression_data.csv"))
expression.tag.df <- as.data.frame(read.csv("expression_data_tag.csv"))

temp <- expression.data.df[,-1]
rownames(temp) <- expression.data.df[,1]
expression.data.df <- as.data.frame(temp)
expression.data.df <- t(expression.data.df)

temp <- as.data.frame(expression.tag.df[,-1])
rownames(temp) <- expression.tag.df[,1]
expression.tag.df <- as.data.frame(temp)
colnames(expression.tag.df) <- c('Status')

# 2. Extract and delete all NA values from the imported dataset
na.genes <- apply(expression.data.df, 2, function(c) any(is.na(c)))
clean.expression.data <- expression.data.df[, !na.genes]
na.gene.list <- colnames(expression.data.df[, na.genes])
print(na.gene.list)

# 3. Dimensionality Reduction via UMAP
expressions.umap <- umap(clean.expression.data)
expressions.umap.df <- as.data.frame(expressions.umap)
expressions.umap.df$Status <- expression.tag.df$Status

ggplot(expressions.umap.df, aes(x = V1, y = V2, colour = Status)) + geom_point()

expressions.3component.umap <- umap(clean.expression.data, n_components=3)
layout3d <- as.data.frame(expressions.3component.umap)
layout3d$Status <- expression.tag.df$Status
umap.plot.3d.df <- data.frame(x = layout3d[,1], y = layout3d[,2], z = layout3d[,3], Status = expression.tag.df)
plot_ly(umap.plot.3d.df, x = ~x, y = ~y, z = ~z, color = umap.plot.3d.df$Status)

# 4. Feature Selection via ANOVA
p_values <- c()
for (g in 1:2490) {
  aov_results <- aov(clean.expression.data[,g] ~ expression.tag.df[,1])
  p_values <- c(p_values, summary(aov_results)[[1]][["Pr(>F)"]][[1]])
}
df.pvalues <- as.data.frame(p_values)
genes <- colnames(clean.expression.data)
rownames(df.pvalues) <- genes

stat.expression.data <- data.frame(t(clean.expression.data), df.pvalues, stringsAsFactors = FALSE)
filtered.expression.data <- subset(stat.expression.data, stat.expression.data$p_values < 0.05)
filtered.genes <- rownames(filtered.expression.data)
length(filtered.genes)
cat(paste(shQuote(filtered.genes, type="cmd"), collapse = ", "))

# 5. Plot data after feature selection and filtering
filtered.expression.data <- filtered.expression.data[,1:600]
filtered.expression.data.t <- t(filtered.expression.data)

filtered.expressions.umap <- umap(filtered.expression.data.t)
layout.filtered <- as.data.frame(filtered.expressions.umap)

umap.filtered.plot.df <- data.frame(x = layout.filtered[,1], y = layout.filtered[,2], Status = expression.tag.df)
ggplot(umap.filtered.plot.df, aes(x, y, colour = Status)) + geom_point()

# Topic 2.
# 1. Apply hierarchical clustering on dataset after feature selection
dist.expressions <- dist(filtered.expression.data.t, method = 'euclidean')
hclust.expressions <- hclust(dist.expressions, method = "complete")
plot(hclust.expressions)

#2. Hierarchical clustering and visualisation with complete, average and ward linkage
hclust.complete <- hclust(dist.expressions, method = "complete")
hclust.average <- hclust(dist.expressions, method = "average")
hclust.ward <- hclust(dist.expressions, method = "ward.D")

#p.complete <- ggdendrogram(hclust.complete, rotate = FALSE)
#p.average <- ggdendrogram(hclust.average, rotate = FALSE)
#p.ward <- ggdendrogram(hclust.ward, rotate = FALSE)

plot(hclust.complete)
plot(p.average)
plot(p.ward)

#3. 
fviz_silhouette(silhouette(cutree(hclust.expressions, 2), dist.expressions))
fviz_silhouette(silhouette(cutree(hclust.expressions, 3), dist.expressions))
fviz_silhouette(silhouette(cutree(hclust.expressions, 4), dist.expressions))
fviz_silhouette(silhouette(cutree(hclust.expressions, 5), dist.expressions))
fviz_silhouette(silhouette(cutree(hclust.expressions, 6), dist.expressions))
fviz_silhouette(silhouette(cutree(hclust.expressions, 7), dist.expressions))
fviz_silhouette(silhouette(cutree(hclust.expressions, 8), dist.expressions))
fviz_silhouette(silhouette(cutree(hclust.expressions, 9), dist.expressions))
fviz_silhouette(silhouette(cutree(hclust.expressions, 10), dist.expressions))

