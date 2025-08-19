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

# Find markers with stricter thresholds (adjust based on your data)
markers <- FindAllMarkers(
  seurat_Eqc,
  only.pos = TRUE,
  min.pct = 0.25,    # Gene detected in at least 25% of cells in a cluster
  logfc.threshold = 0.5,  # Minimum log2 fold-change
  test.use = "wilcox"     # Wilcoxon rank-sum test (default)
)
# check top markers per cluster 
# View top 5 markers per cluster
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

########################################################################################### 
# Annotating clusters 

# updated gene markers for glioma markers 
glioma_markers <- list(
  "Myeloid_Cells" = c("CD14", "CD68", "AIF1", "CSF1R", "TYROBP", "C1QA", "C1QB", "C1QC", "LYZ", "S100A8", "S100A9"),
  "Glioma_Cells" = c("SOX2", "OLIG2", "GFAP", "EGFR", "PDGFRA", "VIM", "NES", "S100B", "IDH1", "TNC", "CHI3L1"),
  "T_Cells" = c("CD3D", "CD3E", "CD8A", "CD4", "IL7R", "FOXP3", "CD2", "TRAC", "TRBC1", "GZMA", "GZMB"),
  "Oligodendrocytes" = c("MBP", "PLP1", "MOBP", "MOG", "CNP", "MAG", "OLIG1", "OLIG2", "ERMN"),
  "Endothelial_Cells" = c("PECAM1", "CDH5", "CLDN5", "VWF", "FLT1", "ESAM", "CD34", "PLVAP", "SLC2A1"),
  "Pericytes" = c("PDGFRB", "RGS5", "CSPG4", "ACTA2", "DES", "MCAM", "NOTCH3", "KCNJ8"),
  "B_Cells" = c("CD19", "MS4A1", "CD79A", "PAX5", "CD20", "CD22", "CD79B", "BANK1"),
  "Neurons" = c("SLC17A7", "SLC17A6", "SYT1", "RBFOX3", "MAP2", "DCX", "NEFL", "GAD1", "GAD2"),
  "Astrocytes" = c("AQP4", "ALDH1L1", "SLC1A2", "SLC1A3", "GJA1", "GLUL", "FGFR3"),
  "Microglia" = c("TMEM119", "P2RY12", "CX3CR1", "SALL1", "GPR34", "HEXB"),
  "Progenitor_Cells" = c("PROM1", "CD44", "SOX9", "NESTIN", "HES1", "ASCL1"),
  "Tumor_Associated_Macrophages" = c("CD163", "MRC1", "MSR1", "TREM2", "APOE", "LGALS3")
)

# convert to a searchable dataframe 
marker_df <- stack(glioma_markers) %>%
  setNames(c("gene", "cell_type")) %>%
  mutate(cell_type = as.character(cell_type))

# annotate clusters 
# Find matches with known markers
annotations <- markers %>%
  inner_join(marker_df, by = "gene") %>%
  group_by(cluster, cell_type) %>%
  summarize(
    n_markers = n(),
    avg_logfc = mean(avg_log2FC),
    .groups = "drop"
  ) %>%
  arrange(cluster, desc(n_markers))

# View top matches per cluster
annotations %>% group_by(cluster) %>% slice_max(n_markers, n = 3)
#################
## add annotation to seurat object 
# First, extract the dominant cell type for each cluster
dominant_annotations <- annotations %>%
  group_by(cluster) %>%
  slice_max(n_markers, n = 1, with_ties = FALSE) %>%
  # Handle cases where multiple cell types tie for top markers
  group_by(cluster) %>%
  mutate(cell_type = ifelse(n() > 1, 
                            paste(cell_type, collapse = "/"), 
                            cell_type)) %>%
  distinct(cluster, cell_type)

# Convert cluster numbers to match Seurat's format
dominant_annotations$cluster <- as.character(dominant_annotations$cluster)

# Add annotations to metadata
seurat_Eqc[["annotated_clusters"]] <- dominant_annotations$cell_type[match(
  as.character(seurat_Eqc$seurat_clusters),
  dominant_annotations$cluster
)]

# Verify the annotation was added
head(seurat_Eqc[[]])

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

## Split by cell type graphs
DimPlot(seurat_Eqc,
        group.by = "annotated_clusters",
        split.by = "annotated_clusters",
        ncol = 4)  # Adjust columns as needed



##############################################################################################

# ## validate manually in T-cells
# # Example: Check T-cell markers in Cluster 3
# cluster5_markers <- FindMarkers(seurat_Eqc, ident.1 = 3, only.pos = TRUE)
# tcell_markers <- human_cell_markers[grepl("T cell", human_cell_markers$cell_name), ]
# 
# # Check overlap
# intersect(cluster5_markers$gene, tcell_markers$gene)


# # Load built-in CellMarker data (human)
# data(human_cell_markers)
# 
# # Search for T-cell markers
# tcell_markers <- human_cell_markers %>%
#   filter(str_detect(cell_name, "T cell")) %>%
#   pull(marker_gene)
# 
# # Check overlap with your Cluster 5 markers
# cluster5_markers <- FindMarkers(seurat_Eqc, ident.1 = 3, only.pos = TRUE)
# intersect(cluster5_markers$gene, tcell_markers)  # Should return CD3E, CD8A, etc.




# Clustering based on the paper's gene markers 

# new_ids <- c(
#   "1" = "Myeloid_Cells",        # C1 (PTPRC/CD45+, ITGAM/CD11B+, CD68+)
#   "2" = "Glioma_Cells",         # C2 (SOX2+, OLIG1+, GFAP+, S100B+)
#   "3" = "T_Cells",             # C3 (PTPRC/CD45+, CD3E+, CD4/CD8+)
#   "4" = "Myeloid_Cells",       # C4 
#   "5" = "Oligodendrocytes",    # C5 (OLIG2+, MBP+)
#   "6" = "Glioma_Cells",        # C6
#   "7" = "Myeloid_Cells",       # C7
#   "8" = "Pericytes",           # C8 (ACTA2+, PDGFRB+)
#   "9" = "Glioma_Cells",        # C9
#   "10" = "Endothelial_Cells",  # C10 (PECAM+)
#   "11" = "B_Cells",            # C11 (CD79A+, CD19+)
#   "12" = "Glutamatergic_Neurons",
#   "13" = "CD4_T_Cells",
#   "14" = "B_Plasma_Cells"
# )
# 
# # Apply new IDs
# seurat_Eqc$celltype <- plyr::mapvalues(
#   x = seurat_Eqc$RNA_snn_res.0.6,
#   from = names(new_ids),
#   to = new_ids
# )
# 
# DimPlot(seurat_Eqc, group.by = "celltype", label = TRUE)


# validate with key markers 
# Key marker visualization
markers_to_check <- list(
  Myeloid = c("PTPRC", "ITGAM", "CD68"),
  Glioma = c("SOX2", "OLIG1", "GFAP", "S100B"),
  T_Cells = c("CD3E", "CD4", "CD8A"),
  B_Cells = c("CD79A", "CD19"),
  Pericytes = c("ACTA2", "PDGFRB"),
  Endothelial = c("PECAM1"),
  Oligodendrocytes = c("OLIG2", "MBP")
)

DotPlot(seurat_Eqc, 
        features = unlist(markers_to_check),
        group.by = "celltype",
        cols = c("lightgrey", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#############################################################################################
# Cell chat 

data.input <- GetAssayData(seurat_Eqc, assay = "RNA", slot = "data")  # Normalized data but not scaled data 
metadata <- seurat_Eqc@meta.data
celltypes <- metadata$celltype  # Use your annotated cell types

# CellChatDB <- CellChatDB.human  # Explicitly use human database
# cellchat@DB <- CellChatDB
# 
# # Verify database loaded correctly
# print(dim(cellchat@DB$interaction))
# 
# ## create cellChat object 
# # Recompute with ultra-sensitive settings
# cellchat <- computeCommunProb(
#   cellchat,
#   type = "truncatedMean",
#   trim = 0.2,           # Increased from 0.1 (include more cells)
#   raw.use = TRUE,       # Use non-scaled data if scaling removed signals
#   population.size = TRUE, # Account for group size differences
#   thresh = 0.01         # Lower probability threshold
# )
# 
# # Analyze glioma-immune interactions
# glioma_targets <- c("Myeloid_Cells", "T_Cells", "B_Cells", "Pericytes")
# pathways_of_interest <- c("TGF-beta", "GALECTIN", "ANNEXIN", "PD-L1")  # "TGF-beta" instead of "TGFb"
# 
# for(pathway in pathways_of_interest){
#   netVisual_aggregate(
#     cellchat, 
#     signaling = pathway,
#     sources.use = c("Glioma_Cells"),
#     targets.use = glioma_targets,
#     layout = "circle"
#   )
# }
# 
# #saveRDS(seurat_Eqc, "seurat_paper_cluster_annotation.rds")


###########################################################################
# Doc playaround 

cellchat <- createCellChat(object = data.input, group.by = "annotated_clusters", assay = "RNA", meta = metadata)

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 4)

# Part II: Inference of cell-cell communication network

## Compute the communication probability and infer cellular communication network
cellchat@meta <- data.frame(
  row.names = colnames(cellchat@data),  # assumes cell barcodes are column names
  labels = metadata$annotated_clusters  # assign celltype labels
)

options(future.globals.maxSize = 900 * 1024^2)  # Set to 900 MiB (adjust as needed)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

options(future.globals.maxSize = 1.5 * 1024^3)  # Set to 900 MiB (adjust as needed)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# displaying for all clusters/cell types
mat <- cellchat@net$weight
for (i in 1:nrow(mat)) {  
  mat2 <- matrix(0, nrow = nrow(mat), 
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))  
  mat2[i, ] <- mat[i, ]  
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])}



saveRDS(seurat_Eqc, "seurat_annotated_cluster.rds")
saveRDS(cellchat, file = "cellchat_up2date_annot_doc.rds") 