

library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(CellChat)
library(Seurat)
library(patchwork)

seurat_Eqc <- readRDS(file = "seurat_annotated_cluster.rds")

cellchat <- readRDS(file = "cellchat_up2date_annot_doc.rds") 
#######################################################################

## plot annotated clusters 
library(ggplot2)
library(RColorBrewer)

# Create a color palette with enough distinct colors
n_types <- length(unique(seurat_Eqc$annotated_clusters))
cluster_colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_types)

# Generate the plot
p <- DimPlot(seurat_Eqc, 
             group.by = "annotated_clusters",
             label = TRUE,
             repel = TRUE,
             pt.size = 0.5,
             cols = cluster_colors) +
  ggtitle("Annotated Cell Clusters") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right")

print(p)
############################################################################
# cellChat heatMap plot 
pathways.show <- c("BAFF") 
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# Compute the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# Identify signaling roles (e.g., major sources/targets)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# bird's-eye view.
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction strength/weight")

# Rank pathways based on overall information flow
pathways.show <- cellchat@netP$pathways
# View the list, the top ones are globally important


# Create the heatmap for the pathway
library(NMF)
library(ComplexHeatmap)

# This heatmap shows the complete communication pattern for the MIF pathway
netVisual_heatmap(cellchat,
                  signaling = "MIF", 
                  measure = "weight", # Use 'weight' for strength, 'count' for number of interactions
                  color.heatmap = "Reds",
                  title.name = "MIF Signaling Pathway Strength",
                  font.size = 10,
                  font.size.title = 12)





#######################################################################


# VinPlot 
VlnPlot(seurat_Eqc, features = c("S100A4", "S100A8", "S100A9", "S100A6"), 
        group.by = "annotated_clusters", pt.size = 0.1)

# proteins related to S100A6
VlnPlot(seurat_Eqc, features = c("CACYBP", "ANXA11", "AGER", "S100A6"), 
        group.by = "annotated_clusters", pt.size = 0.1)

VlnPlot(seurat_Eqc, features = c("S100A4", "S100A6"), 
        group.by = "annotated_clusters", pt.size = 0.1)

# Functional Signatures (e.g., Hypoxia, EMT)
# 
# Plot gene sets associated with tumor microenvironments:

hypoxia_genes <- c("HIF1A", "VEGFA", "LDHA", "SLC2A1")
emt_genes <- c("VIM", "ZEB1", "CDH2", "FN1")

VlnPlot(seurat_Eqc, features = c(hypoxia_genes, emt_genes), 
        group.by = "annotated_clusters",
        ncol = 4)

# Immune checkpoints 
checkpoint_genes <- c("PDCD1", "CTLA4", "LAG3", "TIGIT")

VlnPlot(seurat_Eqc, features = checkpoint_genes,
        group.by = "annotated_clusters", # Subset clusters
        pt.size = 0.1)

VlnPlot(seurat_Eqc, features = c("HIF1A", "VIM"),
        group.by = "annotated_clusters",
        split.by = "Grade_group")  
########################################################################
# Dot plot 

# marker genes 
markers <- FindAllMarkers(
  seurat_Eqc,
  only.pos = TRUE,
  min.pct = 0.25,    # Gene detected in at least 25% of cells in a cluster
  logfc.threshold = 0.5,  # Minimum log2 fold-change
  test.use = "wilcox"     # Wilcoxon rank-sum test (default)
)

# Define the desired order of cell types
celltype_order <- c(
  "Myeloid_Cells", "Tumor_Associated_Macrophages", "Microglia",
  "Glioma_Cells", "Progenitor_Cells",
  "T_Cells", "B_Cells", "NK_Cells",
  "Oligodendrocytes", "Astrocytes", "Neurons",
  "Endothelial_Cells", "Pericytes"
)

# Update Idents with the correct order
Idents(seurat_Eqc) <- factor(
  seurat_Eqc$annotated_clusters, 
  levels = celltype_order
)

# # Dot plot for Genes
# # Get unique genes (removes duplicates)
# genes_to_plot <- unique(unlist(markers, use.names = FALSE))
# 
# # Proceed with DotPlot
# DotPlot(seurat_Eqc, 
#         features = genes_to_plot,
#         cols = c("blue", "red"),
#         dot.scale = 6,
#         cluster.idents = FALSE
# ) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(title = "Glioma Marker Genes (Unique Genes Only)")


## DotPlot for highly expressed genes 
# Get expressed genes (non-zero in at least 1 cell)
expressed_genes <- rownames(seurat_Eqc)[Matrix::rowSums(GetAssayData(seurat_Eqc, slot = "counts")) > 0]

# Filter your gene list to keep only expressed genes
genes_to_plot <- unique(unlist(markers, use.names = FALSE))
genes_to_plot <- genes_to_plot[genes_to_plot %in% expressed_genes]  # Keep only expressed

DotPlot(seurat_Eqc, 
        features = genes_to_plot,
        cols = c("blue", "red"),
        dot.scale = 6,
        cluster.idents = FALSE
) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "Glioma Marker Genes (Expressed Only)")

## Exploring GLUL gene 
# Plot GLUL expression across clusters
VlnPlot(seurat_Eqc, features = "GLUL", group.by = "annotated_clusters")

# Check co-expression with myeloid/glioma markers
FeaturePlot(seurat_Eqc, features = c("GLUL", "CD68"), blend = TRUE) # , "SOX2"
##############################################################################
# GLUL Questions 
#Spatial Localization of GLUL+ Cells

netVisual_aggregate(cellchat, signaling = "glutamate", sources.use = c("Myeloid_Cells"), targets.use = c("Glioma_Cells"))

# Violin plot to compare expression
VlnPlot(seurat_Eqc, features = c("GLUL", "CD68"), 
        group.by = "annotated_clusters", pt.size = 0)
###

# We don't have spatial data 
# Step 1: Identify GLUL+CD68+ cells
seurat_Eqc$glul_cd68 <- ifelse(
  seurat_Eqc@assays$RNA$counts["GLUL", ] > 0 & 
    seurat_Eqc@assays$RNA$counts["CD68", ] > 0, 
  "Double+", "Other"
)

# Step 2: Compare their transcriptomes
markers <- FindMarkers(seurat_Eqc, 
                       ident.1 = "Double+", 
                       group.by = "glul_cd68")
head(markers)  # Check enriched pathways (e.g., glycolysis, mTOR)

# Step 3: Correlate with outcomes (if clinical data exists)
VlnPlot(seurat_Eqc, features = "GLUL",  # e.g., "Responder" vs "Non-responder"
        group.by = "glul_cd68")

# 
# # Install packages (if not already installed)
# install.packages(c("survival", "survminer"))
# 
# # Load libraries
# library(survival)    # For survfit()
# library(survminer)   # For ggsurvplot() (visualization)
# survfit(Surv(time, status) ~ GLUL_CD68_overlap, data = clinical_df)
############################################################################################
# Immune Checkpoints PDCD1

# Define your checkpoint genes
checkpoint_genes <- c("PDCD1", "CTLA4", "LAG3", "TIGIT")

# 1. VIOLIN PLOTS
# This shows the distribution of expression for each gene across all cell types.
violin_plots <- VlnPlot(seurat_Eqc,
                        features = checkpoint_genes,
                        pt.size = 0.1, # Size of the dots (set to 0 to remove them)
                        group.by = "annotated_clusters", # Specify your metadata column for cell types
                        ncol = 2) # Arrange plots in 2 columns
print(violin_plots)

# 2. FEATURE PLOTS (on UMAP)
feature_plot_list <- list()
for (gene in checkpoint_genes) {
  feature_plot_list[[gene]] <- FeaturePlot(seurat_Eqc,
                                           features = gene,
                                           order = TRUE, # Plot positive cells on top
                                           label = TRUE,
                                           repel = TRUE) +
    ggtitle(gene) +
    theme(legend.position = "right")
}
# Combine all feature plots into one image
combined_feature_plots <- wrap_plots(feature_plot_list, ncol = 2)
print(combined_feature_plots) # they are all there 


# 3. DOT PLOT (Most informative summary)
dot_plot <- DotPlot(seurat_Eqc,
                    features = checkpoint_genes,
                    group.by = "annotated_clusters", # Specify your metadata column for cell types
                    col.min = 0, # Minimum expression level to color
                    dot.scale = 6 # Scale the size of the dots
) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) # Rotate x-axis labels
print(dot_plot)


# Exhausted T cells are a type of immune cell that have lost their ability to function properly 
# due to prolonged stimulation, often seen in chronic infections and cancer
t_cell_subset <- subset(seurat_Eqc, subset = annotated_clusters == "T_Cells" & Number_Checkpoints >= 2)
# Now we can visualize just these potentially exhausted T cells
DimPlot(t_cell_subset, label = TRUE) + ggtitle("T Cells Expressing >=2 Checkpoints")


# How to Interpret the Results:
# 
#     Violin Plots: The "shape" of the plot shows the distribution. A tall, thin plot means only a small subset of cells in that group express the gene. A short, wide plot means many cells express it at varying levels.
# 
#     Feature Plots: Look for genes that are expressed in specific clusters. You expect PDCD1, LAG3, and TIGIT to be in T cell clusters. CTLA4 might be more specific to a regulatory T cell subset.
# 
#     Dot Plot: This is your best friend. Large, dark dots represent cell types where a high percentage of cells express a high level of the gene. This is a strong signal of biological importance. For example, you might see that TIGIT is highly expressed in most NK cells, while PDCD1 is expressed in a smaller subset of T cells.
# 
#     Co-expression Plot: Cells that are red in the Number_Checkpoints plot are highly likely to be "exhausted" T cells. The bar plot gives you a precise quantitative breakdown of this population.












