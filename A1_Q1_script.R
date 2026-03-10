####Initialize libraries####
library(ggplot2)

####Load data####
# Load dataset
load("geneexpression2.rda")

# Inspect object names
ls()

# Rename object to "a1data"
a1data <- dat

####Visualization####
# View structure
str(a1data)

# Dimensions
dim(a1data)

# Check first rows
head(a1data)

# Check row names
head(rownames(a1data))

# Check for missing values
sum(is.na(a1data))

# Summary statistics
summary(a1data)

#####Plot based visualizations#####
# Boxplot of dataset, looking for outliers and large shifts in expression levels
boxplot(a1data,
        outline = FALSE,
        las = 2,
        main = "Distribution of Gene Expression Across Samples",
        ylab = "Expression Level")

# Heatmap of gene expression, looking for potential groupings
heatmap(as.matrix(a1data),
        scale = "row",
        col = colorRampPalette(c("blue","white","red"))(100))

# Sample correlation heatmap, showing whether subjects affect expression patterns
sample_cor <- cor(t(a1data))

heatmap(sample_cor,
        col = colorRampPalette(c("blue","white","red"))(100),
        main = "Sample Correlation Heatmap")

####Extract cell type and subject status####
# Extract row names
sample_names <- rownames(a1data)

# Split row names
split_names <- strsplit(sample_names, "_")

# Extract components
subject_status <- sapply(split_names, "[", 1)
cell_type <- sapply(split_names, "[", 2)

# Create metadata dataframe
metadata <- data.frame(
  Subject = subject_status,
  CellType = cell_type
)

head(metadata)

####PCA####
# Perform PCA
pca_result <- prcomp(a1data, scale. = TRUE)

# Summary of PCA
summary(pca_result)

####Variance explanation through PCs####
# Variance explained
pca_var <- pca_result$sdev^2

# Proportion of variance explained
var_explained <- pca_var / sum(pca_var)

# View first PCs
var_explained[1:10]

# Cumulative variance
cumsum(var_explained)[1:10]

#####Scree Plot#####
plot(var_explained,
     type = "b",
     xlab = "Principal Component",
     ylab = "Variance Explained",
     main = "Scree Plot")

####PCA Plots####
#####Prepare PCA Data#####
pca_result <- prcomp(a1data, scale. = TRUE)

sample_pca <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  CellType = cell_type,
  Subject = subject_status
)

gene_pca <- data.frame(
  PC1 = pca_result$rotation[,1],
  PC2 = pca_result$rotation[,2],
  Gene = rownames(pca_result$rotation)
)

#####PCA Plot of Samples#####
ggplot(sample_pca, aes(PC1, PC2, color = CellType)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Samples",
       x = "PC1",
       y = "PC2")

#####PCA Plot of Genes#####
ggplot(gene_pca, aes(PC1, PC2)) +
  geom_point(color = "darkred") +
  theme_minimal() +
  labs(title = "PCA of Genes (Loadings)",
       x = "PC1",
       y = "PC2")
