install.packages("umap")
install.packages("plotly")
install.packages("ggdendro")
install.packages("factoextra")
install.packages("uwot")
install.packages("caret")
install.packages("rpart.plot")
library("rpart")
library("rpart.plot")
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
expressions.umap.df <- as.data.frame(expressions.umap$layout)
expressions.umap.df$Status <- expression.tag.df$Status

ggplot(expressions.umap.df, aes(x = V1, y = V2, colour = Status)) + geom_point()

expressions.3component.umap <- umap(clean.expression.data, n_components=3)
layout3d <- as.data.frame(expressions.3component.umap$layout)
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
layout.filtered <- as.data.frame(filtered.expressions.umap$layout)

umap.filtered.plot.df <- data.frame(x = layout.filtered[,1], y = layout.filtered[,2], Status = expression.tag.df)
ggplot(umap.filtered.plot.df, aes(x = x, y = y, colour = Status)) + geom_point()

# Topic 2.
# 1. Apply hierarchical clustering on dataset after feature selection
dist.expressions <- dist(umap.filtered.plot.df[,1:2], method = 'euclidean')
hclust.expressions <- hclust(dist.expressions, method = "complete")

install.packages("dendextend")
library(dendextend)

dend <- as.dendrogram(hclust.expressions)
dend <- color_branches(dend)
dend %>% 
  set("labels", umap.filtered.plot.df$Status) %>%
  set("labels_colors", as.numeric(as.factor(umap.filtered.plot.df$Status)), order_value=TRUE) %>%  
  set("labels_cex", 0.2) %>% 
  plot(main = "Euclidean Distance & Complete Linkage")

#2. Hierarchical clustering and visualisation with complete, average and ward linkage
hclust.complete <- hclust(dist.expressions, method = "complete")
hclust.average <- hclust(dist.expressions, method = "average")
hclust.ward <- hclust(dist.expressions, method = "ward.D")

dend <- as.dendrogram(hclust.expressions)
dend <- color_branches(dend)
dend %>% 
  set("labels", umap.filtered.plot.df$Status) %>%
  set("labels_colors", as.numeric(as.factor(umap.filtered.plot.df$Status)), order_value=TRUE) %>%  
  set("labels_cex", 0.2) %>% 
  plot(main = "Euclidean Distance & Complete Linkage")

dend <- as.dendrogram(hclust.average)
dend <- color_branches(dend)
dend %>% 
  set("labels", umap.filtered.plot.df$Status) %>%
  set("labels_colors", as.numeric(as.factor(umap.filtered.plot.df$Status)), order_value=TRUE) %>%  
  set("labels_cex", 0.2) %>% 
  plot(main = "Average Linkage")

dend <- as.dendrogram(hclust.ward)
dend <- color_branches(dend)
dend %>% 
  set("labels", umap.filtered.plot.df$Status) %>%
  set("labels_colors", as.numeric(as.factor(umap.filtered.plot.df$Status)), order_value=TRUE) %>%  
  set("labels_cex", 0.2) %>% 
  plot(main = "Ward Linkage")

#3. 
library("cluster")
fviz_silhouette(silhouette(cutree(hclust.expressions, 2), dist.expressions))
fviz_silhouette(silhouette(cutree(hclust.expressions, 3), dist.expressions))
fviz_silhouette(silhouette(cutree(hclust.expressions, 4), dist.expressions))
fviz_silhouette(silhouette(cutree(hclust.expressions, 5), dist.expressions))
fviz_silhouette(silhouette(cutree(hclust.expressions, 6), dist.expressions))
fviz_silhouette(silhouette(cutree(hclust.expressions, 7), dist.expressions))
fviz_silhouette(silhouette(cutree(hclust.expressions, 8), dist.expressions))
fviz_silhouette(silhouette(cutree(hclust.expressions, 9), dist.expressions))
fviz_silhouette(silhouette(cutree(hclust.expressions, 10), dist.expressions))

clusters <- cutree(hclust.expressions, 2)
umap.filtered.plot.df$Cluster <- as.factor(clusters)
ggplot(umap.filtered.plot.df, aes(x = x, y = y, color = Cluster)) + geom_point()

# Topic 3.
# 1.
knn.data <- cbind(filtered.expression.data.t, expression.tag.df)
training <- knn.data[sample(nrow(knn.data), 420),]
training.classes <- training$Status
training <- training[,-132]

test <- knn.data[!(row.names(knn.data) %in% row.names(training)),]
test.classes <- test$Status
test <- test[,-132]

# 2.
library(caret)
library(class)
knn.results <- knn(train = training, test = test, cl = training.classes, k = 5)
conf.matrix <- table(knn.results, test.classes)
confusionMatrix(conf.matrix)

# 3.
install.packages("rpart.plot")
library(rpart)
library(rpart.plot)

model <- rpart(training.classes ~ ., data = training)
print(summary(model))
rpart.plot(model)
predictions <- predict(model, test, type="class")
conf.matrix.dt <- confusionMatrix(predictions, factor(test.classes))
print(conf.matrix.dt)

# 4.
importance <- varImp(model, scale=FALSE)
importance <- importance[order(-importance$Overall),, drop=FALSE]
topten <- importance[1:10,,drop=FALSE]
print(topten)
