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
# Date: 04-05-2022
# Written by: Natalie Piehl
# Summary: Generate violin plots of CXCR4 and CXCR6 expression
#
# - Fig S4D created using pseudobulked data with Prism 9.2.0
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
  library("ggpubr")
  library("scales")
  library("immunarch")
})

# Initialize paths
seurat_object <- "path/to/seurat_object/"
bulk_data <- "path/to/pseudobulked/data/"
output_dir <- "path/to/export/results"

# Source helper functions
source("../0_preprocessing/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Generate violin/box plots of CXCR4 and CXCR6 in CD4+ T Cells (Fig S5C)

# Load in seurat object
load(seurat_object)

# Normalize with standard processing (log transformation)
DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE)

# Subset for C and NC in CD4+ T Cells
s_sub <- subset(s, clonal %in% c("C", "NC") & cluster_ident == "CD4+ T Cells")
s_sub$clonal <- factor(s_sub$clonal, levels = c("NC", "C"))

# Generate violin plot
Idents(s_sub) <- "Diagnosis"
p <- VlnPlot(s_sub, features = "CXCR4", split.by = c("clonal"), pt.size = 0,
             log = TRUE, cols = c("#E5007F", "#7CC8AB")) +
  geom_boxplot(width=0.5, outlier.shape = NA)
set_panel_size(p,
               file = paste0(output_dir, "CD4+_T_Cells_CXCR4.pdf"),
               width = unit(4, "in"), height = unit(6, "in"))

p <- VlnPlot(s_sub, features = "CXCR6", split.by = c("clonal"), pt.size = 0,
             cols = c("#E5007F", "#7CC8AB"), log = TRUE) +
  geom_boxplot(width=0.5, outlier.shape = NA)
set_panel_size(p,
               file = paste0(output_dir, "CD4+_T_Cells_CXCR6.pdf"),
               width = unit(4, "in"), height = unit(6, "in"))