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
# Date: 04-11-2022
# Written by: Natalie Piehl
# Summary: Generate UMAP of C and NC cells
#
# - Coexpression UMAP of CXCR6 and CXCR4 generated with the ShinyCell app
#   at https://gatelabnu.shinyapps.io/csf_aging/
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

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Generate clonal frequency UMAP (Fig 4A)

# Load in seurat object
load(seurat_object)

# Set Ident to normalized frequency
s <- SetIdent(s, value = "clonal")

# Remove NA cells
s <- subset(s, clonal %in% c("C", "NC"))

# Generate UMAP
myplot <- FeaturePlot(s,
                      features = c("clonal"), pt.size = 1,
                      cols = c("aquamarine2", "deeppink"), reduction = "umap",
                      order = TRUE
)
set_panel_size(myplot,
               file = paste0(output_dir, "clonality.pdf"),
               width = unit(5, "in"), height = unit(4, "in")
)