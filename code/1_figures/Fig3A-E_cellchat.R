# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----     Single cell transcriptomics reveals CD8 T cell recruitment     -----
# -----       to the cerebrospinal fluid during cognitive impairment       -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 06-14-2022
# Written by: Natalie Piehl
# Summary: Run cellchat
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("CellChat")
})

# Initialize paths
seurat_object <- "path/to/seurat/object"
output_dir <- "path/to/export/results"

# Source helper functions
source("../0_preprocessing/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Data formatting

# Load Seurat object
load(seurat_object)

# Run standard normalization
DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE)

# Subset for celltypes of interest
s <- subset(s, cluster_ident %!in% c("CD4+/CD8+ T Cells", "Undetermined"))

# Subset for disease or age bin
s <- subset(s, Diagnosis != "HC")

#-------------------------------------------------------------------------------
# Run CellChat

# Create cellchat object
cellchat <- createCellChat(object = s,
                           group.by = "cluster_ident")

# Define ligand-receptor pair database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

# Set the used database in the object
cellchat@DB <- CellChatDB

# Preprocessing
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Calculate communication probability and network
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Save object
save(cellchat, file = paste0(output_dir, "cellchat"))

#-------------------------------------------------------------------------------
# Generate interaction strength scatter plot (Fig 3B)

# Define celltype colors
colors_use <- c(
  "darkviolet", "gold1", "red", 
  "darkorange", "darkturquoise", "lawngreen", "lightpink", "hotpink",
  "violet", "navy"
)

# Visualize outgoing vs incoming signal
pdf(paste0(output_dir, "celltype_strength_scatter.pdf"))
netAnalysis_signalingRole_scatter(cellchat,
                                  color.use = colors_use) +
  theme_Publication_blank()
dev.off()

#-------------------------------------------------------------------------------
# Rerun on HC and CI separately

# Subset for disease
s_ci <- subset(s, Diagnosis != "HC")
s_hc <- subset(s, Diagnosis == "HC")

# Run cellchat
cellchat_ci <- createCellChat(object = s_ci,
                              group.by = "cluster_ident")
cellchat_hc <- createCellChat(object = s_hc,
                              group.by = "cluster_ident")

# Set the used database in the object
cellchat_ci@DB <- CellChatDB
cellchat_hc@DB <- CellChatDB

# Preprocessing
cellchat_hc <- subsetData(cellchat_hc)
cellchat_ci <- subsetData(cellchat_ci)
cellchat_hc <- identifyOverExpressedGenes(cellchat_hc)
cellchat_hc <- identifyOverExpressedInteractions(cellchat_hc)
cellchat_ci <- identifyOverExpressedGenes(cellchat_ci)
cellchat_ci <- identifyOverExpressedInteractions(cellchat_ci)

# Calculate communication probability and network
cellchat_hc <- computeCommunProb(cellchat_hc)
cellchat_hc <- computeCommunProbPathway(cellchat_hc)
cellchat_hc <- aggregateNet(cellchat_hc)
cellchat_hc <- netAnalysis_computeCentrality(cellchat_hc, slot.name = "netP")
cellchat_ci <- computeCommunProb(cellchat_ci)
cellchat_ci <- computeCommunProbPathway(cellchat_ci)
cellchat_ci <- aggregateNet(cellchat_ci)
cellchat_ci <- netAnalysis_computeCentrality(cellchat_ci, slot.name = "netP")

# Save object
save(cellchat_hc, file = paste0(output_dir, "cellchat_hc"))
save(cellchat_ci, file = paste0(output_dir, "cellchat_ci"))

# Create cellchat object
cellchat <- mergeCellChat(list(cellchat_hc, cellchat_ci),
                          add.names = c("HC", "CI"))
cell_chat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Save merged
save(cellchat, file = paste0(output_dir, "cellchat"))

#-------------------------------------------------------------------------------
# Generate dot plot of CD68hi -> T cell interactions (Fig 3C)

# Generate bubble of CD68hi signaling to t cells
pdf(paste0(output_dir, "bubble_cd68hi_to_t_cells.pdf"))
netVisual_bubble(cellchat_ci,
                 sources.use = 3,
                 targets.use = c(5,6,10),
                 remove.isolate = FALSE)
dev.off()

#-------------------------------------------------------------------------------
# Generate interaction plots of HC vs CI (Fig 3A)

# Make interaction circle plot for HC only
pdf(paste0(output_dir, "hc_number_of_interactions.pdf"))
netVisual_circle(cellchat_hc@net$count, weight.scale = T,
                 label.edge= F, title.name = "HC Cell-Cell Interactions",
                 color.use = colors_use)
dev.off()

# Make interaction circle plot for CI only
pdf(paste0(output_dir, "ci_number_of_interactions.pdf"))
netVisual_circle(cellchat_hc@net$count, weight.scale = T,
                 label.edge= F, title.name = "CI Cell-Cell Interactions",
                 color.use = colors_use)
dev.off()

#-------------------------------------------------------------------------------
# Generate interaction signaling in CI plot (Fig 3D)

# Generate bubble of CD68hi signaling to t cells increased in CI
pdf(paste0(output_dir, "bubble_cd68hi_to_t_cells_ci_increased.pdf"))
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(6),
                 comparison = c(1, 2), max.dataset = 2,
                 title.name = "Increased signaling in CI",
                 angle.x = 45, remove.isolate = T,)
dev.off()

#-------------------------------------------------------------------------------
# Generate plot of CXCL signaling (Fig 3E)

# Make CI CXCL signaling circle plot
pdf(paste0(output_dir, "celltype_signalling_cxcl.pdf"))
netVisual_aggregate(cellchat_ci,
                    signaling = c("CXCL"),
                    color.use = colors_use)
dev.off()