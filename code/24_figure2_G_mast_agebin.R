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
# Summary: Run DE with MAST on age bins
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
library_vector <- c("plyr", "tidyverse", "ggrepel", "ggthemes", "Seurat", "grid",
                    "scales", "doMC", "UpSetR")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = library_vector)

# Initialize paths
seurat_object <- "path/to/seurat_object/"
output_dir <- "path/to/export/results"

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set core number for parallel model fitting
registerDoMC(cores = 12)

#-------------------------------------------------------------------------------
# Run DE (Fig 2G)

# Load Seurat object
load(seurat_object)

# Subset for HC samples only
s <- subset(s, Diagnosis == "HC")

# Run standard normalization
DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE)

# Set Ident to age bin
s <- SetIdent(s, value = "age_bin")

run_de <- function(cell_type) {
  # Define cell type label
  cell_type_label <- gsub("/", "", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  
  # Find DEGs b/w MCI/AD and middle
  degs <-FindMarkers(object = subset(s, cluster_ident == cell_type),
                     ident.1 = "advanced",
                     ident.2 = "middle",
                     test.use = "MAST",
                     logfc.threshold = -Inf,
                     latent.vars = c("sex"),
                     min.pct = 0.1,
                     assay = "RNA"
  )
  
  # Remove ribosomal, mitochondrial, and HLA genes
  degs <- degs[-grep(pattern = "^RPS|^RPL|^MT-|^HLA-", x = rownames(degs)),]
  
  # Run Benjamini-Hochberg adjustment
  degs$BH <- p.adjust(degs$p_val, method = "BH")
  
  # Write out results
  write.csv(degs, paste0(output_dir, cell_type_label, "_advanced_vs_middle_degs.csv"))
  
  # Create volcano plot
  volcano_plot(degs, title = paste0("advanced vs middle in ", cell_type),
               file = paste0(output_dir, cell_type_label, "_advanced_vs_middle_volcano.pdf"),
               padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
}

# Run DE on all celltypes
cell_types <- unique(s[["cluster_ident"]])[,1]
cell_types <- cell_types[ cell_types %!in% c("CD4+/CD8+ T Cells", "Undetermined")]
mclapply(cell_types, run_de, mc.cores = 12)