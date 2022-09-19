# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----        Distinct cerebrospinal fluid immune perturbations           -----
# -----          in healthy brain aging and cognitive impairment           -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 04-18-2022
# Written by: Natalie Piehl
# Summary: Run supervised clustering with CITEseq reference atlas
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("SeuratDisk")
})

# Initialize paths
seurat_object <- "path/to/seurat_object/"
ref_atlas_path <- "path/to/ref/atlas"
output_dir <- "path/to/export/results"

# Source helper functions
source("code/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#------------------------------------------------------------------------------
# Format data

# Load in seurat object
load(seurat_object)
s

# Load in reference
reference <- LoadH5Seurat(ref_atlast_path)

# Run standard normalization
DefaultAssay(object = s) <- "SCT"

# Find transfer anchors
anchors <- FindTransferAnchors(
  reference = reference,
  query = s,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

# Map to reference atlas
s <- MapQuery(
  anchorset = anchors,
  query = s,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)
s

# Save updated object
save(s, file = paste0(output_dir, "s_sup_clustering"))

# Visualize results
p1 = DimPlot(s, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE)
p2 = DimPlot(s, reduction = "umap", group.by = "predicted.celltype.l2",
             label = TRUE, label.size = 3 ,repel = TRUE)
set_panel_size(p1, file = paste0(output_dir, "og_umap_L1.pdf"), width = unit(5, "in"), height = unit(4, "in"))
set_panel_size(p2, file = paste0(output_dir, "og_umap_L2.pdf"), width = unit(5, "in"), height = unit(4, "in"))

p1 = DimPlot(s, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE)
p2 = DimPlot(s, reduction = "ref.umap", group.by = "predicted.celltype.l2",
             label = TRUE, label.size = 3 ,repel = TRUE)
set_panel_size(p1, file = paste0(output_dir, "umap_L1.pdf"), width = unit(5, "in"), height = unit(4, "in"))
set_panel_size(p2, file = paste0(output_dir, "umap_L2.pdf"), width = unit(5, "in"), height = unit(4, "in"))

# Visualize with clonality
s_clonal <- subset(s, clonal %in% c("NC", "C"))
p2 = DimPlot(s_clonal, reduction = "ref.umap", group.by = "clonal",
             label = TRUE, label.size = 3 ,repel = TRUE,
             cols = c("aquamarine2", "deeppink"))
p2
set_panel_size(p2, file = paste0(output_dir, "umap_L2_clonal.pdf"), width = unit(5, "in"), height = unit(4, "in"))
