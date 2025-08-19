library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)

# # joining meta data with seurat object 
# meta <- read.csv("Meta_GBM.csv")
# 
# # The first row contains the TYPE information
# type_row <- meta[1, ]
# 
# # Get the current column names
# current_colnames <- colnames(meta)
# 
# # Combine column names with TYPE information
# new_colnames <- paste(current_colnames, type_row, sep = "_")
# 
# # Assign the new column names
# colnames(meta) <- new_colnames
# 
# # Remove the TYPE row (now redundant as it's incorporated in column names)
# meta <- meta[-1, ]
# 
# # Write the reformatted data
# write.csv(meta, "reformatted_metadata.csv", row.names = FALSE)

meta <- read.csv("reformatted_metadata.csv")

# object from early QC process 
seurat_Eqc <- readRDS("seurat_subset_removed_mt_&_low_nFeatures.rds")

# 1. First ensure your metadata is properly formatted with cell barcodes as row names
rownames(meta) <- meta$NAME_TYPE

# 2. Select only the columns you want to add
columns_to_add <- c("biosample_id_group", "donor_id_group", "Grade_group", 
                    "sex_group", "Phase_group", "GSMID_group")  # Changed 'phase' to match your metadata
meta_subset <- meta[, columns_to_add]

# 3. Match cells between metadata and Seurat object
common_cells <- intersect(rownames(meta_subset), colnames(seurat_Eqc))
meta_subset <- meta_subset[common_cells, ]
seurat_Eqc <- seurat_Eqc[, common_cells]

# 4. Add the metadata to your Seurat object
seurat_Eqc <- AddMetaData(seurat_Eqc, metadata = meta_subset)

# 5. Verify the metadata was added correctly
head(seurat_Eqc@meta.data)

# # 6. Example analyses using GSMID grouping:
# 
# # Calculate average expression by GSMID
# avg_expression <- AverageExpression(seurat_Eqc,
#                                     assays = "RNA")
# 
# # Visualization
# DimPlot(seurat_Eqc, group.by = "orig.ident") +
#   ggtitle("Cells by GSMID (orig.ident)")
# 
# # Find markers for each GSMID group
# markers <- FindAllMarkers(seurat_Eqc, assay = "RNA")
# 
# # You can also create new groupings combining multiple factors
# seurat_Eqc$group_combo <- paste(seurat_Eqc$orig.ident, 
#                                 seurat_Eqc$Grade_group, 
#                                 sep = "_")

##############################################################################
# Normalizing the data
seurat_Eqc <- NormalizeData(seurat_Eqc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection) 
seurat_Eqc <- FindVariableFeatures(seurat_Eqc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_Eqc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_Eqc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling the data 
all.genes <- rownames(seurat_Eqc)
seurat_Eqc <- ScaleData(seurat_Eqc, features = all.genes)

# Perform linear dimensional reduction 

seurat_Eqc <- RunPCA(seurat_Eqc, features = VariableFeatures(object = seurat_Eqc))

# Visualize PCA results colored by potential batch variables
p1 <- DimPlot(seurat_Eqc, reduction = "pca", group.by = "orig.ident") + 
  ggtitle("PCA by Sample (orig.ident)")
p2 <- DimPlot(seurat_Eqc, reduction = "pca", group.by = "donor_id_group") + 
  ggtitle("PCA by Donor")
p3 <- DimPlot(seurat_Eqc, reduction = "pca", group.by = "Phase_group") + 
  ggtitle("PCA by Cell Cycle Phase")

# Combine plots
library(patchwork)
(p1 + p2) / p3

# Examine PCA loadings to see which features drive variation
VizDimLoadings(seurat_Eqc, dims = 1:2, reduction = "pca")

#Heatmap of first PCA components
DimHeatmap(seurat_Eqc, dims = 1:6, cells = 500, balanced = TRUE)

# Run UMAP for non-linear visualization
seurat_Eqc <- RunUMAP(seurat_Eqc, dims = 1:20)

# Plot UMAP colored by batch variables
u1 <- DimPlot(seurat_Eqc, reduction = "umap", group.by = "orig.ident") + 
  ggtitle("UMAP by Sample")
u2 <- DimPlot(seurat_Eqc, reduction = "umap", group.by = "donor_id_group") + 
  ggtitle("UMAP by Donor")
u1 + u2

################################################################################
# removing Batch effect 
library(harmony)

# Run Harmony integration
# Now run Harmony - using the correct function name
seurat_Eqc <- RunHarmony(
  object = seurat_Eqc,
  group.by.vars = "donor_id_group",
  reduction = "pca",
  assay = "RNA",  # Note: 'assay' instead of 'assay.use'
  reduction.save = "harmony"  # Explicitly save the reduction
)

# Create UMAP using Harmony embeddings
seurat_Eqc <- RunUMAP(seurat_Eqc, 
                      reduction = "harmony", 
                      dims = 1:30)

# Visualize corrected data
DimPlot(seurat_Eqc, 
        reduction = "umap", 
        group.by = "donor_id_group") + 
  ggtitle("Post-Harmony Integration")

#################################################################################################
# Cell clustering 

## determining how many PCA to choose 
ElbowPlot(seurat_Eqc) # i'll go with 1:12 PCA


# Test multiple resolutions (adjust based on your expected cell type diversity)

# Find clusters at different resolutions
seurat_Eqc <- FindNeighbors(seurat_Eqc, dims = 1:20, reduction = "harmony")  

# Leiden clustering (algorithm 4)
# Try and play with the resolution value 

seurat_Eqc <- FindClusters(
  seurat_Eqc,
  resolution = 0.6,  # Start here, adjust later
  algorithm = 4,     # 4 = Leiden
  method = "igraph", # Uses `igraph` backend
  random.seed = 42   # For reproducibility
)

# validate clusters 
# Visualize
DimPlot(seurat_Eqc, group.by = "RNA_snn_res.0.6", label = TRUE)  

# Check marker genes per cluster
markers <- FindAllMarkers(seurat_Eqc, only.pos = TRUE, min.pct = 0.25)
head(markers)

# Top 5 markers per cluster (adjust 'n' as needed)
top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

View(top_markers)

# Plot
DotPlot(seurat_Eqc, 
        features = unique(top_markers$gene),
        cols = c("lightgrey", "red"),  # Customize colors
        dot.scale = 6) + 
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster2.markers <- FindMarkers(seurat_Eqc, ident.1 = 2, ident.2 = c(1, 3))
head(cluster2.markers, n = 5)


# find markers for every cluster compared to all remaining cells, report only the positive
# ones

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


# features plot 
FeaturePlot(seurat_Eqc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))
############################################################################
# Features plot 

## Anti-Tumor Myeloid Subsets (Better Survival)
# Inflammatory Microglia (i-Mic)
FeaturePlot(seurat_Eqc, features = c("CCL3", "CCL4", "TNF", "IL1B", "NFKBIZ"), 
            ncol = 3, order = TRUE, cols = c("lightgrey", "red"))

# Homeostatic Microglia (h-Mic)
FeaturePlot(seurat_Eqc, features = c("CST3", "P2RY12", "TMEM119"), 
            ncol = 3, order = TRUE, cols = c("lightgrey", "blue"))

# Activated Microglia (a-Mic)
FeaturePlot(seurat_Eqc, features = c("SPRY1", "P2RY13"), 
            ncol = 2, order = TRUE, cols = c("lightgrey", "purple"))

## Pro-Tumor Myeloid Subsets (Worse Survival)
# M2-like Macrophages (s-Mac1)
FeaturePlot(seurat_Eqc, features = c("CD14", "CD163", "MSR1"), 
            ncol = 3, order = TRUE, cols = c("lightgrey", "darkred"))

# MDSCs
FeaturePlot(seurat_Eqc, features = c("MIF", "ITGA4"), 
            ncol = 2, order = TRUE, cols = c("lightgrey", "orange"))

# Immunosuppressive Macrophages (s-Mac2)
FeaturePlot(seurat_Eqc, features = c("S100A4", "VEGFA", "TGFB1", "IL10"), 
            ncol = 2, order = TRUE, cols = c("lightgrey", "brown"))

# Antigen-Presenting Microglia (AP-Mic)
FeaturePlot(seurat_Eqc, features = c("CX3CR1", "CD86", "IFNGR1"), 
            ncol = 3, order = TRUE, cols = c("lightgrey", "cyan"))

## glioma cell status 
# Mesenchymal (MES-like)
FeaturePlot(seurat_Eqc, features = c("EGFR", "CHI3L1", "CD44"), 
            ncol = 3, order = TRUE, cols = c("lightgrey", "darkgreen"))

# Neural Progenitor-like (NPC-like)
FeaturePlot(seurat_Eqc, features = c("SOX2", "OLIG2"), 
            ncol = 2, order = TRUE, cols = c("lightgrey", "darkblue"))

# Astrocyte-like (AC-like)
FeaturePlot(seurat_Eqc, features = c("GFAP", "S100B"), 
            ncol = 2, order = TRUE, cols = c("lightgrey", "pink"))

# Oligodendrocyte-like (OPC-like)
FeaturePlot(seurat_Eqc, features = c("PDGFRA", "OLIG1"), 
            ncol = 2, order = TRUE, cols = c("lightgrey", "gold"))


saveRDS(seurat_Eqc, "seurat_normalised_cellClustered_foundAllGenes.rds")
