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
# Summary: Generate bubbleplot of enrichment analysis
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
library_vector <- c("plyr", "tidyverse", "ggpubr", "ggrepel", "ggthemes",
                    "stringr", "xlsx", "DEswan", "scales")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = library_vector)

# Initialize paths
seurat_object <- "path/to/seurat_object/"
output_dir <- "path/to/export/results"

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