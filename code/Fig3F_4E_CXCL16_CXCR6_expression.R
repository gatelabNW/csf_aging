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
# Date: 04-05-2022
# Written by: Natalie Piehl
# Summary: Generate violin/box plots of CXCR4, CXCR6, and CXCL16 expression
#
# - Coexpression UMAP of CXCR6 and CXCL16 generated with the ShinyCell app
#   at https://gatelabnu.shinyapps.io/csf_aging/
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
})

# Initialize paths
seurat_object <- "path/to/seurat_object/"
output_dir <- "path/to/export/results"

# Source helper functions
source("code/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#------------------------------------------------------------------------------
# Generate CXCR6 violin/box plot (Fig 4E)

# Load in seurat object
load(seurat_object)

# Normalize with standard processing (log transformation)
DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst",
                       verbose = FALSE)

# Subset for C and NC in CD8+ T Cells
s_sub <- subset(s, clonal %in% c("C", "NC") & cluster_ident == "CD8+ T Cells")
s_sub$clonal <- factor(s_sub$clonal, levels = c("NC", "C"))

# Generate violin plot
p <- VlnPlot(s_sub, features = "CXCR6", split.by = c("clonal"), pt.size = 0,
             cols = c("#E5007F", "#7CC8AB"), log = TRUE) +
  geom_boxplot(width=0.5, outlier.shape = NA)
set_panel_size(p,
               file = paste0(output_dir, "CD8+_T_Cells_CXCR6.pdf"),
               width = unit(4, "in"), height = unit(6, "in"))

#------------------------------------------------------------------------------
# Generate CXCL16 violin/box plot (Fig 3F)

# Set ident to cell type
Idents(s) <- "cluster_ident"
levels(s) <- c(
  "CD4+ T Cells", "CD8+ T Cells", "CD4+/CD8+ T Cells", "T Regulatory Cells",
  "NK Cells", "DC",
  "CD14+/CD16-/CD68lo Monocytes", "CD14+/CD16+/CD68mid Monocytes", "CD14+/CD16+/CD68hi Monoctyes",
  "B Cells", "Plasma Cells", "Undetermined"
)

# Subset for clusters of interest
s <- subset(s, cluster_ident %!in% c("Undetermined", "CD4+/CD8+ T Cells"))

# Generate violin plot
p <- VlnPlot(s, features = "CXCL16", pt.size = 0,
             log = TRUE, cols = c(
               "darkturquoise", "lawngreen", "navy",
               "hotpink", "lightpink", 
               "gold1", "darkorange", "red",
               "darkviolet", "violet"
             )) +
  geom_boxplot(width=0.5, outlier.shape = NA)
p
set_panel_size(p,
               file = paste0(output_dir, "cxcl16_violin.pdf"),
               width = unit(4, "in"), height = unit(6, "in"))