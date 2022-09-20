# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----   Single cell transcriptomics reveals cerebrospinal fluid immune   -----
# -----  dysregulation during healthy brain aging and cognitive impairment -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 09-08-2022
# Written by: Natalie Piehl
# Summary: Look at CXCR6 in CD8TEMs C vs NC, CI vs HC
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("car")
  library("FSA")
})

# Initialize paths
seurat_object <- "path/to/seurat_object/"
output_dir <- "path/to/export/results"

# Source helper functions
source("../0_preprocessing/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#------------------------------------------------------------------------------
# Format data

# Load in seurat object
load(seurat_object)

# Run standard normalization
DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE)

# Isolate data
data <- GetAssayData(s, slot = "data")
data <- as.matrix(data[c("CXCR6", "CXCL16"),]) %>% t()
data <- merge(data, s[[c("Diagnosis", "predicted.celltype.l2", "clonal")]], by = 0)
data <- data[which(data$clonal == "C"), ]
data <- data[which(data$predicted.celltype.l2 == "CD8 TEM"), ]

# Run Wilcox Rank Sum
wilcox_test <- function(celltype, gene) {
  data_sub <- data[which(data$predicted.celltype.l2 == celltype),]
  pval <- tryCatch({
    wilcox <- wilcox.test(data_sub[which(data_sub$Diagnosis == "HC"), gene],
                          data_sub[which(data_sub$Diagnosis != "HC"), gene])
    return(wilcox$p.value)
  }, error = function(e) {
    return(1)
  })
  return(pval)
}

wilcox_cxcr6 <- sapply(unique(data$predicted.celltype.l2),
                       wilcox_test, gene = "CXCR6")

#------------------------------------------------------------------------------
# Visualize expression on only clonal CD8 TEM cells (Fig 4H)

# Violin Plot CXCR6
Idents(s) <- "predicted.celltype.l2"
s <- subset(s, predicted.celltype.l2 == "CD8 TEM")
s <- subset(s, clonal == "C")

p <- VlnPlot(s, features = "CXCR6",
             split.by = 'Diagnosis',
             cols = c("gray", "red"),
             pt.size = 0.5,
             log = FALSE)
  geom_boxplot(width=0.5, alpha=0.5, outlier.shape = NA, position = position_dodge(width = 0.9))

# Export
set_panel_size(p,
               file = paste0(output_dir, "cxcr6_violin_cd8tem_clonal.pdf"),
               width = unit(6, "in"), height = unit(6, "in"))
