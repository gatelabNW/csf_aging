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
# Date: 10-11-2021
# Written by: Natalie Piehl
# Summary: Visualize contamination removal with SoupX on representative sample
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("gridExtra")
  library("Matrix")
  library("SoupX")
  library("Seurat")
  library("DropletUtils")
  library("Hmisc")
})

# Initialize paths
gex_dir <- "path/to/cellranger_count/output"
output_dir <- "path/to/export/results"

# Source helper functions
source("../0_preprocessing/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Generate tSNEs visualizing contamination removal (Fig S1B)

# Create lists of directories to load as SoupX objects and Seurat objects
sample_dirs <- list.dirs(gex_dir, recursive = FALSE)

# Isolate sample A2
dir <- sample_dirs[2]
sample <- "A2"

# Initialize gene list to estimate contamination fraction
# Selecting for monocyte/dendritic markers which are highly specific and abundant
M.genes <- c("CD14", "CD68", "MS4A7", "FCGR3A")

# Load in SoupX and Seurat data
soupx <- load10X(paste0(dir, "/outs"))
seurat <- Read10X(paste0(dir, "/outs/filtered_feature_bc_matrix")) %>%
  CreateSeuratObject()
seurat

# Normalize and run PCA on Seurat object
seurat <- SCTransform(seurat, variable.features.n = 2000, verbose = TRUE) %>%
  RunPCA()

# Generate clusters
seurat <- RunTSNE(seurat, dims = 1:10) %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.3)

# Add cluster info to SoupX object
soupx <- setClusters(soupx, setNames(seurat$seurat_clusters, colnames(seurat)))

# Estimate contamination fraction
useToEst <- estimateNonExpressingCells(soupx, nonExpressedGeneList = list(M.genes = M.genes))
soupx <- calculateContaminationFraction(soupx,
                                        list(M.genes = M.genes),
                                        useToEst = useToEst)

# Create adjusted counts
adj_counts <- adjustCounts(soupx)

# Generate corrected Seurat object
seurat_adj <- CreateSeuratObject(adj_counts) %>%
  SCTransform(variable.features.n = 2000, verbose = TRUE) %>%
  RunPCA() %>%
  RunTSNE(dims = 1:10) %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.3)

# Isolate non zero CD14 cells
cd14 <- FetchData(seurat, vars = "CD14")
cd14 <- cd14[ which(cd14$CD14 > 0), ,drop = FALSE]
cd14_highlight <- data.frame(Embeddings(seurat[["tsne"]])[ rownames(cd14), ]) %>%
  merge(cd14, by = 0)
cd14_outline <- cd14_highlight[ which(cd14_highlight$tSNE_1 < 25), ]

# Plot uncorrected CD14 expression
p1 <- FeaturePlot(seurat, features = "CD14", reduction = "tsne", pt.size = 1) +
  theme_Publication_blank() +
  geom_point(data=cd14_highlight,
             aes(x=tSNE_1, y=tSNE_2, color=CD14),
             size=3,
             alpha=1) +
  geom_point(data=cd14_outline,
             aes(x=tSNE_1, y=tSNE_2),
             shape=1,
             color='black',
             size=3,
             alpha=1) +
  labs(title = "CD14 Expression pre-Correction",
       color = "CD14 Expression")

# Plot soup fraction of CD14
p <- plotChangeMap(soupx, adj_counts, geneSet = "CD14", DR = Embeddings(seurat[["tsne"]]))
soupfrac_highlight <- p$data[ cd14_highlight$Row.names, ]
soupfrac_outline <- soupfrac_highlight[ which(soupfrac_highlight$RD1 < 25), ]

p2 <- plotChangeMap(soupx, adj_counts, geneSet = "CD14", DR = Embeddings(seurat[["tsne"]])) +
  theme_Publication_blank() +
  geom_point(data=soupfrac_highlight,
             aes(x=RD1, y=RD2, color=relChange),
             size=3,
             alpha=1) +
  geom_point(data=soupfrac_outline,
             aes(x=RD1, y=RD2),
             shape=1,
             color='black',
             size=3,
             alpha=1) +
  scale_colour_gradient(low = "blue", high = "red", limits = c(0,1)) +
  labs(title = "Soup Fraction of CD14 Expression")

# Isolate non zero CD14 cells
cd14_adj <- FetchData(seurat_adj, vars = "CD14")
cd14_adj <- cd14_adj[ which(cd14_adj$CD14 > 0), ,drop = FALSE]
cd14_adj_highlight <- data.frame(Embeddings(seurat[["tsne"]])[ rownames(cd14_adj), ]) %>%
  merge(cd14_adj, by = 0)

# Plot corrected CD14 expression
p3 <- FeaturePlot(seurat_adj, features = "CD14", reduction = "tsne", pt.size = 1) +
  theme_Publication_blank() +
  geom_point(data=cd14_adj_highlight,
             aes(x=tSNE_1, y=tSNE_2, color=CD14),
             size=3,
             alpha=1) +
  labs(title = "CD14 Expression post-Correction")

# Export plots
set_panel_size(p1, file=paste0(output_dir, "cd14_precorrection.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))

set_panel_size(p2, file=paste0(output_dir, "cd14_soupfrac.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))

set_panel_size(p3, file=paste0(output_dir, "cd14_postcorrection.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))