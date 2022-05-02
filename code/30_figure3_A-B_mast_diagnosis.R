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
# Summary: Run DE with MAST on diagnosis
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
  library("scales")
  library("doMC")
  library("UpSetR")
})

# Initialize paths
seurat_object <- "path/to/seurat_object/"
output_dir <- "path/to/export/results"

# Source helper functions
source("code/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set core number for parallel model fitting
registerDoMC(cores = 12)

# Specify thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

#-------------------------------------------------------------------------------
# Run DE (Fig 3B)

# Load Seurat object
load(seurat_object)

# Run standard normalization
DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE)

# Set Ident to Diagnosis
s <- SetIdent(s, value = "Diagnosis")

run_de <- function(cell_type) {
  # Generate cell typ label
  cell_type_label <- gsub("/", "", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  
  # Find DEGs b/w CI and HC
  degs <-FindMarkers(object = subset(s, cluster_ident == cell_type),
                     ident.1 = "CI",
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

  # Write out results
  write.csv(degs, paste0(output_dir, cell_type_label, "_CI_vs_HC_degs.csv"))
  
  # Create volcano plot
  volcano_plot(degs, title = paste0("CI vs HC in ", cell_type),
               file = paste0(output_dir, cell_type_label,
                             "_CI_vs_HC_volcano.pdf"),
               padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
}

# Run DE on all celltypes
cell_types <- unique(s[["cluster_ident"]])[,1]
cell_types <- cell_types[
  cell_types %!in% c("CD4+/CD8+ T Cells", "Undetermined")]
mclapply(cell_types, run_de, mc.cores = 12)

#------------------------------------------------------------------------------
# Generate Upset plot (Fig 3A)

# Initialize sig gene list
sig_genes_ls <- list()

# Create list with sig genes for each cell type
for (cell_type in cell_types) {
  # Generate cell type label
  cell_type_label <- gsub("/", "", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  
  # Load in degs
  degs <- read.csv(paste0(output_dir, cell_type_label, "_CI_vs_HC_degs.csv"))
  
  # Identify sig genes
  sig_genes <- degs[
    which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
  
  # Add sig genes to list
  sig_genes_ls[[cell_type]] <- sig_genes$X
}

# Find number of genes in each set
num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
for (key in names(sig_genes_ls)) {
  num_degs <- rbind(num_degs, data.frame(key, length(sig_genes_ls[[key]])))
}

# Get order of sets for coloring
num_degs <- num_degs[order(-num_degs[, 2]),]

# Define colors
color_vector <- num_degs[, 1]
color_vector <- mapvalues(
  color_vector,
  from = c(
    "CD4+ T Cells", "CD8+ T Cells", "CD4+/CD8+ T Cells",
    "CD14+/CD16+/CD68hi Monoctyes", "DC", "NK Cells",
    "CD14+/CD16-/CD68lo Monocytes", "T Regulatory Cells", "Plasma Cells",
    "CD14+/CD16+/CD68mid Monocytes", "B Cells", "Undetermined"
  ),
  to = c(
    "darkturquoise", "lawngreen", "dodgerblue", "red",
    "lightpink", "hotpink", "gold1", "navy", "violet",
    "darkorange", "darkviolet", "gray"
  )
)

# Create upset plot and export
pdf(file = paste0(output_dir, "/upset_celltype.pdf"))
upset(
  fromList(sig_genes_ls),
  nsets = length(sig_genes_ls),
  order.by = "freq",
  sets.bar.color = color_vector
)
dev.off()