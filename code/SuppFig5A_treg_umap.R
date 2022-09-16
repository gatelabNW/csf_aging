# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----         Activated monocytes recruit CD8 T cells to the             -----
# -----         cerebrospinal fluid during cognitive impairment            -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 03-04-2022
# Written by: Natalie Piehl
# Summary: Generate UMAPs of Tregs
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("ggrepel")
  library("ggthemes")
  library("grid")
  library("Seurat")
  library("scales")
  library("immunarch")
})

# Initialize paths
seurat_object <- "path/to/seurat_object/"
output_dir <- "path/to/export/results"

# Source helper functions
source("code/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#------------------------------------------------------------------------------
# Generate Treg UMAPs (Fig S5A)

# Subset normalized frequency for HC Tregs
s@meta.data$normalized_frequency[
  which(s@meta.data$cluster_ident != "T Regulatory Cells") ] <- NA

# Generate clonality UMAP
p <- FeaturePlot(s, features = "normalized_frequency", pt.size = 1,
                 cols = c("grey88", "deeppink"), reduction = "umap",
                 order = TRUE, max.cutoff = 8)
p

# Export plot
set_panel_size(p,
               file = paste0(output_dir, "umap_treg_norm_freq.pdf"),
               width = unit(5, "in"), height = unit(4, "in"))


# Generate FOXP3 UMAP
DefaultAssay(s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE)
p <- FeaturePlot(s, features = "FOXP3", pt.size = 1,
                 cols = c("grey88", "#1000FF"), reduction = "umap",
                 order = TRUE)
p

# Export plot
set_panel_size(p,
               file = paste0(output_dir, "umap_treg_foxp3.pdf"),
               width = unit(5, "in"), height = unit(4, "in"))

# FOXP3 of Tregs only
s@assays$RNA@data['FOXP3',][which(s@meta.data$cluster_ident != "T Regulatory Cells")] <- 0
p <- FeaturePlot(s, features = "FOXP3", pt.size = 1,
                 cols = c("grey88", "#1000FF"), reduction = "umap",
                 order = TRUE)
p

# Export plot
set_panel_size(p,
               file = paste0(output_dir, "umap_treg_foxp3_tregsonly.pdf"),
               width = unit(5, "in"), height = unit(4, "in"))