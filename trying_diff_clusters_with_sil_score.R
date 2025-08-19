library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)

seurat_cl <- readRDS("seurat_normalised_cellClustered_foundAllGenes.rds")


# Cell clustering 

## determining how many PCA to choose 
ElbowPlot(seurat_cl) # i'll go with 1:12 PCA


# Test multiple resolutions (adjust based on your expected cell type diversity)

# Find clusters at different resolutions
seurat_cl <- FindNeighbors(seurat_cl, dims = 1:20, reduction = "harmony")  

# Leiden clustering (algorithm 4)
# Try and play with the resolution value 
# resolution 0.4
seurat_Eqc <- FindClusters(
  seurat_cl,
  resolution = 0.6,  
  algorithm = 4,     # 4 = Leiden
  method = "igraph", # Uses `igraph` backend
  random.seed = 42   # For reproducibility
)

# validate clusters 
# Visualize
# DimPlot(seurat_Eqc, 
#         reduction = "umap",  # or "tsne"
#         group.by = "seurat_clusters", 
#         label = TRUE,        # Show cluster IDs
#         repel = TRUE,        # Avoid label overlap
#         pt.size = 0.5) +     # Adjust point size
#   ggtitle("Leiden Clustering (Resolution = 0.2)") +
#   theme_minimal()
#############################################################################
# cal. sil. score
library(cluster)  # For silhouette()

# Get cluster labels
clusters <- seurat_Eqc@meta.data$seurat_clusters

# Get the reduced dimensions (e.g., Harmony/PCA)
embeddings <- Embeddings(seurat_Eqc, reduction = "harmony")[, 1:20]  # Use same dims as FindNeighbors

# Compute Euclidean distance matrix
dist_matrix <- dist(embeddings, method = "euclidean")

# step 2. 
# Calculate silhouette scores
sil_scores <- silhouette(as.numeric(clusters), dist_matrix)

# Summary statistics
summary(sil_scores)

# Average silhouette width per cluster
avg_sil_width <- aggregate(sil_scores[, "sil_width"], 
                           by = list(cluster = sil_scores[, "cluster"]), 
                           FUN = mean)
print(avg_sil_width)

# Plot silhouette scores
plot(sil_scores, col = 1:length(unique(clusters)), main = "Silhouette Plot")

## interpreting scores 
FindMarkers(seurat_Eqc, ident.1 = 5, min.pct = 0.25) 

DimPlot(seurat_Eqc, 
        reduction = "umap", 
        cells.highlight = WhichCells(seurat_Eqc, idents = 5),
        sizes.highlight = 0.5) +
  ggtitle("Cluster 5 (T Cells)")

