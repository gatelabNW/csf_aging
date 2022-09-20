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
# Date: 03-02-2022
# Written by: Natalie Piehl
# Summary: Run DE with MAST on CI vs HC in clonal and nonclonal T cells
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
  library("Seurat")
  library("scales")
  library("doMC")
  library("UpSetR")
  library("colorspace")
})

# Initialize paths
seurat_object <- "path/to/seurat_object/"
output_dir <- "path/to/export/results"

# Source helper functions
source("../0_preprocessing/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set core number for parallel model fitting
registerDoMC(cores = 3)

# Specify thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

#------------------------------------------------------------------------------
# Run DE (Fig 4D)

# Load Seurat object
load(seurat_object)

# Run standard normalization
DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE)

# Subset clonal and nonclonal cells
s_clonal <- subset(s, clonal == "C")
s_nonclonal <- subset(s, clonal == "NC")

run_de <- function(cell_type, s) {
  # Generate labels for naming
  cell_type_label <- gsub("/", "", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  clonality <- gsub("s_", "", deparse(substitute(s)))
  
  # Find DEGs b/w CI and NC
  s <- SetIdent(s, value = "Diagnosis")
  degs <-FindMarkers(object = subset(s, cluster_ident == cell_type),
                     ident.1 = "CI",
                     ident.2 = "HC",
                     latent.vars = c("sex", "age"),
                     test.use = "MAST",
                     logfc.threshold = -Inf,
                     min.pct = 0.1,
                     assay = "RNA"
  )
  
  # Remove ribosomal, mitochondrial, and HLA genes
  degs <- degs[-grep(pattern = "^RPS|^RPL|^MT-|^HLA-", x = rownames(degs)),]
  
  # Run Benjamini-Hochberg adjustment
  degs$BH <- p.adjust(degs$p_val, method = "BH")
  
  # Write out results
  write.csv(degs, paste0(output_dir, cell_type_label, "_", clonality, "_CI_vs_HC_degs.csv"))
  
  # Create volcano plot
  volcano_plot(degs, title = paste0(clonality, " CI vs HC in ", cell_type),
               file = paste0(output_dir, cell_type_label, "_", clonality, "_CI_vs_HC_volcano.pdf"),
               padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
}

# Run DE on all celltypes
cell_types <- c("CD4+ T Cells", "CD8+ T Cells", "T Regulatory Cells")
mclapply(cell_types, run_de, s = s_clonal, mc.cores = 3)
mclapply(cell_types, run_de, s = s_nonclonal, mc.cores = 3)

#------------------------------------------------------------------------------
# Create Upset plot (Fig 4C)

# Initialize sig gene lists
sig_genes_clonal_ls <- list()

# Create list with sig genes for each cell type
for (cell_type in cell_types) {
  # Define cell type label
  cell_type_label <- gsub("/", "", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  
  # Load in degs
  degs_c <- read.csv(paste0(output_dir, cell_type_label, "_clonal_CI_vs_HC_degs.csv"))

  # Identify sig genes
  sig_genes_c <- degs_c[which(degs_c$BH <= padj.thresh & abs(degs_c$avg_log2FC) >= lfc.thresh),]

  # Add sig genes to list
  sig_genes_clonal_ls[[cell_type]] <- sig_genes_c$X
}

# Find number of genes in each set
num_degs_c <- sapply(sig_genes_clonal_ls, length)
sig_genes_clonal_ls <- sig_genes_clonal_ls[order(-num_degs_c)]
sig_genes_clonal_ls <- sig_genes_clonal_ls[ lapply(sig_genes_clonal_ls, length) > 0 ]

# Define colors
color_vector_c <- mapvalues(
  names(sig_genes_clonal_ls),
  from = c("CD4+ T Cells", "CD8+ T Cells", "T Regulatory Cells"),
  to = c("darkturquoise", "lawngreen", "navy")
)

# Create upset plot and export
pdf(file = paste0(output_dir, "/clonal_CI_vs_HC_upset_celltype.pdf"))
upset(
  fromList(sig_genes_clonal_ls),
  nsets = length(sig_genes_clonal_ls),
  order.by = "freq",
  sets.bar.color = color_vector_c,
  text.scale = 2
)
dev.off()
