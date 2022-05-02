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
# Date: 10-11-2021
# Written by: Natalie Piehl
# Summary: Visualize QC on diagnosis and sample ID as well as celltype 
#          composition
#
# - Fig 2SF visual created with Prism 9.2.0
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
  library("gridExtra")
  library("Matrix")
  library("DropletUtils")
  library("Hmisc")
})

# Initialize paths
seurat_object <- "path/to/seurat_object"
output_dir <- "path/to/export/results"

# Source helper functions
source("code/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Generate violin plots across Diagnosis (Fig S2D)

# Load in seurat object
load(seurat_object)

# nFeatures per Diagnosis
Idents(s) <- 'Diagnosis'
myplot <- VlnPlot(s, features = c("nFeature_RNA"), pt.size = 0, cols = c("black", "red")) +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  labs(x = "Diagnosis",
       y = "Number of Features",
       title = "Number of Features per Diganosis")
print(myplot$data)
set_panel_size(myplot, file=paste0(output_dir, "nFeature_Diagnosis.pdf"),
               width=unit(3, "in"), height=unit(4, "in"))

# nCounts per Diagnosis
myplot <- VlnPlot(s, features = c("nCount_RNA"), pt.size = 0, cols = c("black", "red")) +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  labs(x = "Diagnosis",
       y = "Counts",
       title = "Number of Counts per Diganosis")
set_panel_size(myplot, file=paste0(output_dir, "nCount_Diagnosis.pdf"),
               width=unit(3, "in"), height=unit(4, "in"))

# percent.mt per Diagnosis
myplot <- VlnPlot(s, features = c("percent.mt"), pt.size = 0, cols = c("black", "red")) +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  labs(x = "Diagnosis",
       y = "Percent Mitochondrial",
       title = "Percent Mitochondrial per Diagnosis")
set_panel_size(myplot, file=paste0(output_dir, "percent.mt_Diagnosis.pdf"),
               width=unit(3, "in"), height=unit(4, "in"))

#-------------------------------------------------------------------------------
# Generate violin plots across IDs (Fig S2E)

# nFeatures per ID
Idents(s) <- 'ID'
myplot <- VlnPlot(s, features = c("nFeature_RNA"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8)) +
  labs(x = "Sample ID",
       y = "Number of Features",
       title = "Number of Features per Sample ID")
set_panel_size(myplot, file=paste0(output_dir, "nFeature_ID.pdf"),
               width=unit(8, "in"), height=unit(4, "in"))

# nCounts per ID
myplot <- VlnPlot(s, features = c("nCount_RNA"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8)) +
  labs(x = "Sample ID",
       y = "Number of Counts",
       title = "Number of Counts per Sample ID")
set_panel_size(myplot, file=paste0(output_dir, "nCount_ID.pdf"),
               width=unit(8, "in"), height=unit(4, "in"))

#-------------------------------------------------------------------------------
# Calculate cell type composition across age (Fig S2F)

# Isolate age and cell type
age_celltype <- s[[c("age", "cluster_ident")]]

# Generate table of percent cell type per age
table <- table(age_celltype)
table <- table / rowSums(table)
table

# Export table
write.csv(table, paste0(output_dir, "celltypes_over_age.csv"))

#-------------------------------------------------------------------------------
# Generate split UMAPs on sort day (Fig S2G)

# View diagnosis split
s <- SetIdent(s, value = "cluster_ident")
myplot <- DimPlot(s, reduction = "umap", label = TRUE, pt.size = 1,
                  cols = c(
                    "darkturquoise", "lawngreen", "dodgerblue", "red",
                    "lightpink", "hotpink", "gold1", "navy", "violet",
                    "darkorange", "darkviolet", "gray"
                  ),
                  shuffle = TRUE, split.by = "sort_day")
set_panel_size(myplot,
               file = paste0(output_dir, "sort_day_split.pdf"),
               width = unit(5, "in"), height = unit(4, "in")
)