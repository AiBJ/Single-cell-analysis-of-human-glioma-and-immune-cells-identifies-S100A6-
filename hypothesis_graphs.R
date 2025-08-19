library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(CellChat)
library(Seurat)
library(patchwork)

library(Seurat)
library(ggplot2)
library(patchwork) # For combining plots
library(ComplexHeatmap) # For CellChat heatmap
library(CellChat)
setwd("~/Desktop/raw_data_S100A4/")

seurat_Eqc <- readRDS(file = "seurat_annotated_cluster.rds")

cellchat <- readRDS(file = "cellchat_up2date_annot_doc.rds") 

# S100A6 gene expression 
# Create the violin plot
s100a6_vln <- VlnPlot(seurat_Eqc,
                      features = "S100A6",
                      pt.size = 0.1, # Remove the dots to see distribution clearly
                      group.by = "annotated_clusters", # Use your metadata column name for cell identities
                      cols = c("B_Cells" = "skyblue", # Color code for clarity
                               "Glioma_Cells" = "red3",
                               "Myeloid_Cells" = "orange2",
                               "T_Cells" = "forestgreen",
                               "Neurons" = "grey80",
                               "Oligodendrocytes" = "grey70",
                               "Pericytes" = "grey60")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-labels
  ggtitle("S100A6 Expression by Cell Type") +
  ylab("Expression Level") + xlab("Cell Type")

# Display the plot
print(s100a6_vln)

# Save it
ggsave("new_cluster_cellChat_graph/s100a6_violin_plot.png", plot = s100a6_vln, width = 8, height = 6, dpi = 300)



# UMAP S100A6 feature 
# Create the feature plot
s100a6_umap <- FeaturePlot(seurat_Eqc,
                           features = "S100A6",
                           order = TRUE, # Plot positive cells on top
                           label = TRUE,
                           repel = TRUE,
                           cols = c("lightgrey", "darkred")) + # Gradient from low to high
  ggtitle("Spatial Localization of S100A6 Expression")

# Display the plot
print(s100a6_umap)

# Save it
ggsave("new_cluster_cellChat_graph/s100a6_umap_plot.png", plot = s100a6_umap, width = 8, height = 6, dpi = 300)



# Define your checkpoint genes
checkpoint_genes <- c("PDCD1", "CTLA4", "LAG3", "TIGIT")

# Create the dot plot
checkpoint_dotplot <- DotPlot(seurat_Eqc,
                              features = checkpoint_genes,
                              group.by = "annotated_clusters", # Use your metadata column name
                              col.min = 0,
                              dot.scale = 5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") + # Custom color scale
  labs(x = "Gene", y = "Cell Type", title = "T Cell Exhaustion Checkpoint Expression") +
  guides(size = guide_legend(title = "Percent Expressed"),
         color = guide_colorbar(title = "Average Expression"))

# Display the plot
print(checkpoint_dotplot)


# 1. Compute the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# 2. Identify the top pathways based on overall signaling weight
pathways.show <- cellchat@netP$pathways
# Print to console to choose from, e.g., "MIF", "GALECTIN", "SPP1"
print(pathways.show)

# 3. Create a heatmap for a specific pathway (Replace "MIF" with your chosen pathway)
cellchat_heatmap <- netVisual_heatmap(cellchat,
                                      signaling = "MIF",
                                      measure = "weight", # Use 'weight' for communication strength
                                      color.heatmap = "Reds",
                                      title.name = "MIF Signaling Pathway Strength",
                                      font.size = 8,
                                      font.size.title = 10)

# The heatmap is a ComplexHeatmap object. Draw it to display.
draw(cellchat_heatmap, heatmap_legend_side = "left")
png("new_cluster_cellChat_graph/cellchat_mif_heatmap.png", width = 1000, height = 800, res = 150)
draw(cellchat_heatmap, heatmap_legend_side = "left")
dev.off()


# Combine the plots into a multi-panel figure
combined_plot <- (s100a6_vln | s100a6_umap) /
  (checkpoint_dotplot)

# Display the combined plot
combined_plot
##################################################################################################
mat <- cellchat@net$weight

# Find the index of the glioma cell type
glioma_index <- which(rownames(mat) == "Glioma_Cells")

# Get group sizes from the CellChat object
groupSize <- as.numeric(table(cellchat@idents))

# Check if glioma was found before proceeding
if (length(glioma_index) > 0) {
  # Create an empty matrix
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  # Populate only the row for the glioma cell type
  mat2[glioma_index, ] <- mat[glioma_index, ]
  # Plot the interactions for glioma cells only
  netVisual_circle(mat2,
                   vertex.weight = groupSize,
                   weight.scale = T,
                   edge.weight.max = max(mat),
                   title.name = rownames(mat)[glioma_index])
} else {
  print("Cell type 'Glioma_Cells' not found. Check rownames(mat).")
}

groupSize <- as.numeric(table(cellchat@idents))

df.net <- subsetCommunication(cellchat)

# Or, visualize the significance in a heatmap (This is what you NEED for your poster)
print(
  netVisual_heatmap(cellchat,
                    measure = "count",
                    color.heatmap = "Reds",
                    title.name = "Significance of Communication (p-value)")
)


vertex.receiver = seq(1,4)
netVisual_aggregate(cellchat, signaling = "MIF",  vertex.receiver = vertex.receiver, layout = "hierarchy")
