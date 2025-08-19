install.packages("remotes")  # To install from GitHub
install.packages("NMF")      # Required for CellChat
install.packages("ggalluvial")

library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)

if (!requireNamespace("CellChat", quietly = TRUE)) {
  remotes::install_github("sqjin/CellChat")
}
library(CellChat)
library(Seurat)
library(patchwork)

seurat_Eqc <- readRDS("seurat_normalised_cellClustered_foundAllGenes.rds")

all_markers <- list(
  `Microglia` = c("CCL3", "CCL4", "TNF", "IL1B", "NFKBIZ"), #  (i-Mic)
  `Microglia` = c("CST3", "P2RY12", "TMEM119"), # (h-Mic)
  `Macrophages` = c("CD14", "CD163", "MSR1"), # (M2-like)
  `MDSCs` = c("MIF", "ITGA4"),
  `Glioma` = c("EGFR", "CHI3L1", "CD44"), # (MES)
  `Tregs` = c("FOXP3", "IL2RA", "TNFRSF4")
)
####################################################################
# plotting cluster Diagram 
DimPlot(seurat_Eqc, reduction = "umap", label = TRUE, group.by = "RNA_snn_res.0.6")

# Assigning Identity to each cluster 
markers <- FindAllMarkers(seurat_Eqc, only.pos = TRUE, min.pct = 0.25)
# check top markers per cluster 
# View top 5 markers per cluster
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# Print as a formatted table
library(knitr)
kable(top_markers, caption = "Top Markers by Cluster")

# Check specific markers via FeaturePlot
FeaturePlot(seurat_Eqc, features = c("GFAP", "AQP4", "SPOCK3", "OLIG2"), order = TRUE)

# Compare with known cell-type signatures
celltype_markers <- list(
  "Microglia" = c("P2RY12", "TMEM119"),
  "GBM_MES" = c("EGFR", "CD44"),
  "T_Cells" = c("CD3D", "CD8A")
)
DotPlot(seurat_Eqc, features = unlist(celltype_markers), group.by = "RNA_snn_res.0.6")


new_ids <- c(
  "1" = "Resident_Macrophages",
  "2" = "GBM_MES",
  "3" = "Activated_Microglia",
  "4" = "MDSCs",
  "5" = "Astrocytes",
  "6" = "Hypoxic_GBM",
  "7" = "Cytotoxic_T_NK",
  "8" = "Proliferating_GBM",
  "9" = "OPCs",
  "10" = "Neurons",
  "11" = "Fibroblasts",
  "12" = "Glutamatergic_Neurons",
  "13" = "CD4_T_Cells",
  "14" = "B_Plasma_Cells"
)

#Extract MKI67 expression values for Cluster 9
mk167_expr <- GetAssayData(seurat_Eqc, assay = "RNA", slot = "data")["MKI67", ]
cluster9_cells <- WhichCells(seurat_Eqc, idents = 9)  # Get cells in Cluster 9
mk167_mean <- mean(mk167_expr[cluster9_cells])

# Update annotation based on MKI67 expression
if (mk167_mean > 0.5) {
  new_ids["9"] <- "Proliferating_OPCs"
} else {
  new_ids["9"] <- "Quiescent_OPCs"
}

# Apply to Seurat object
seurat_Eqc$celltype <- plyr::mapvalues(
  x = seurat_Eqc$RNA_snn_res.0.6,
  from = names(new_ids),
  to = new_ids
)

# visual confirmation 
VlnPlot(seurat_Eqc, features = "MKI67", group.by = "RNA_snn_res.0.6", idents = 9)

# visulize after assigning the identities 
DimPlot(seurat_Eqc, group.by = "celltype", label = TRUE)
#####################################################################

# # Cell chat 
# 
# ## prepare the data 
# # Ensure your Seurat object has the celltype annotations
# data.input <- GetAssayData(seurat_Eqc, assay = "RNA", slot = "scale.data")  # Normalized data
# metadata <- seurat_Eqc@meta.data
# celltypes <- metadata$celltype  # Use your annotated cell types
# 
# ## create cellChat object 
# cellchat <- createCellChat(
#   object = data.input,
#   meta = metadata,
#   group.by = "celltype"  # Column name in metadata
# )
# 
# ## set the ligned-receptor dataset 
# CellChatDB <- CellChatDB.human  # Use 'CellChatDB.mouse' for mouse data
# cellchat@DB <- CellChatDB
# 
# ## preprocess data
# cellchat <- subsetData(cellchat)  # Optional: subset for speed
# cellchat <- identifyOverExpressedGenes(cellchat)
# cellchat <- identifyOverExpressedInteractions(cellchat)
# 
# # complete communication probabilities 
# cellchat <- computeCommunProb(cellchat)
# cellchat <- filterCommunication(cellchat, min.cells = 10)  # Filter low-confidence interactions
# 
# ## aggregate and visulize network 
# cellchat <- aggregateNet(cellchat)
# groupSize <- as.numeric(table(cellchat@idents))
# netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE)
# 
# ## plot specific signaling pathwayas 
# 
# # check components are expressed 
# # Check expression of TGFβ pathway genes
# tgf_components <- c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3")
# DotPlot(seurat_Eqc, features = tgf_components, group.by = "celltype") + # change group by to display by patients or rank
#   theme(axis.text.x = element_text(angle = 45, hjust=1))
# 
# 
# # Example: Plot TGFb signaling (Immune cells)
# pathways.show <- c("TGFb")
# netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# 
# all_tgf <- subsetCommunication(cellchat, signaling = "TGFb", thresh = 0.05) # Lower threshold
# View(all_tgf) # Look for ANY interactions involving your cell types
# 
# 
# # cancer cells 
# # Filter for TGFβ pairs where tumor clusters are senders
# mif_pairs <- subsetCommunication(
#   cellchat,
#   signaling = "MIF",
#   sources.use = c("Hypoxic_GBM", "GBM_MES", "Proliferating_GBM"),
#   targets.use = c("MDSCs", "Activated_Microglia", "Fibroblasts", "Astrocytes"),
#   thresh = 0.7  # Lower the probability threshold
# )
# 
# nrow(mif_pairs)
# 
# # Plot top interactions (adjust thickness.scale for clarity)
# # 1. First confirm your MIF interactions
# View(mif_pairs)  # Should show MIF-related L-R pairs
# 
# # 2. Plot with corrected parameters
# netVisual_individual(
#   cellchat,
#   signaling = "MIF",  # Must match your filtered pairs!
#   pairLR.use = mif_pairs[1:12, ],  
#   layout = "circle",
#   sources.use = c("Hypoxic_GBM", "GBM_MES", "Proliferating_GBM"),
#   targets.use = c("MDSCs", "Fibroblasts", "Activated_Microglia"),
#   signaling.name = "MIF in GBM-TME Crosstalk"  # Update title
# )
# ################################################
# #
# # heatmap for   interaction strength 
# netVisual_heatmap(cellchat, signaling = "MIF", color.heatmap = "Reds")
# 
# # save and export results 
# saveRDS(cellchat, file = "cellchat_results.rds")  # Save for later
# write.csv(cellchat@net$weight, file = "cellchat_interaction_weights.csv")
# 
# 
# 
# 
# 
# ###########################################
# # Dot plots and specific genes cluster plots 
# 
# # Flatten into a single vector
# marker_genes <- unlist(all_markers)
# 
# # Create DotPlot
# DotPlot(seurat_Eqc, 
#         features = marker_genes,
#         group.by = "RNA_snn_res.0.6",  # Replace with your cluster column
#         cols = c("lightgrey", "red"),
#         dot.scale = 6) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(title = "GBM Marker Expression by Cluster")
# 
# # Feature Plot 
# marker_genes <- c("C1QC", "CSF1R", "TYROBP", "CD3D", "EGFR")
# 
# # Plot on UMAP
# FeaturePlot(seurat_Eqc, 
#             features = marker_genes,
#             reduction = "umap",
#             order = TRUE,          # Highlight expressing cells
#             blend = FALSE,         # Set TRUE to compare 2 genes
#             ncol = 3)    
