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
# Date: 06-14-2022
# Written by: Natalie Piehl
# Summary: DE on diagnosis in advanced age only
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("ggpubr")
  library("ggrepel")
  library("ggthemes")
  library("scales")
  library("grid")
  library("UpSetR")
  library("doMC")
})

# Initialize input parameters
seurat_object <- "path/to/seurat/object"
output_dir <- "path/to/export/results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Source helper functions
source("code/00_helper_functions.R")

# Set thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

# Set core number for parallel model fitting
registerDoMC(cores = 6)

# Load Seurat object
load(seurat_object)

#-------------------------------------------------------------------------------
# Run DE on diagnosis in advanced age and generate volcano plot (Fig 3B)

# Run standard normalization
DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE)

# Subset for advanced only
s <- subset(s, age_bin == "advanced")
print(table(s[["age_bin"]]))

# Set Ident to Diagnosis
s <- SetIdent(s, value = "Diagnosis")
print(table(s[["Diagnosis"]]))

run_de <- function(cell_type) {
  print(cell_type)
  cell_type_label <- gsub("/", "", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  
  # Find DEGs b/w MCI/AD and HC
  degs <-FindMarkers(object = subset(s, cluster_ident == cell_type),
                     ident.1 = "MCI/AD",
                     ident.2 = "HC",
                     latent.vars = c("sex"),
                     test.use = "MAST",
                     logfc.threshold = -Inf,
                     min.pct = 0.1,
                     assay = "RNA"
  )
  
  # Remove ribosomal, mitochondrial, and HLA genes
  degs <- degs[-grep(pattern = "^RPS|^RPL|^MT-|^HLA-", x = rownames(degs)),]
  
  # Run Benjamini-Hochberg adjustment
  degs$BH <- p.adjust(degs$p_val, method = "BH")
  print(head(degs))
  
  # Write out results
  write.csv(degs, paste0(output_dir, cell_type_label, "_advancedonly_MCIAD_vs_HC_degs.csv"))
  
  # Create volcano plot
  volcano_plot(degs, title = paste0("MCI/AD vs HC (Advanced only)\nin ", cell_type),
               file = paste0(output_dir, cell_type_label, "_advancedonly_MCIAD_vs_HC_volcano.pdf"),
               padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
}

# Run DE on all celltypes
cell_types <- unique(s[["cluster_ident"]])[,1]
cell_types <- cell_types[ cell_types %!in% c("CD4+/CD8+ T Cells", "Undetermined")]
mclapply(cell_types, run_de, mc.cores = 6)