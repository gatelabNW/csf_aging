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
# Date: 01-18-2022
# Written by: Natalie Piehl
# Summary: Generate UMAP of normalized clone frequency
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("ggrepel")
  library("ggthemes")
  library("ggpubr")
})

# Initialize paths
seurat_object <- "path/to/seurat_object/"
output_dir <- "path/to/export/results"

# Source helper functions
source("code/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Generate clonal frequency UMAP (Fig 4A)

# Load in seurat object
load(seurat_object)

# Set Ident to normalized frequency
s <- SetIdent(s, value = "normalized_frequency")

# Generate UMAP
myplot <- FeaturePlot(s,
                      features = c("normalized_frequency"), pt.size = 1,
                      cols = c("grey88", "deeppink"), reduction = "umap",
                      order = TRUE
)
set_panel_size(myplot,
               file = paste0(output_dir, "normalized_frequency.pdf"),
               width = unit(5, "in"), height = unit(4, "in")
)