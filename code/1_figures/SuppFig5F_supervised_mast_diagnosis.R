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
# Date: 04-18-2022
# Written by: Natalie Piehl
# Summary: Run DE on diagnosis across supervised annotated celltypes
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("doMC")
  library("UpSetR")
})

# Initialize paths
seurat_object <- "path/to/seurat_object/"
output_dir <- "path/to/export/results"

# Source helper functions
source("../0_preprocessing/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

# Set core number for parallel model fitting
registerDoMC(cores = 12)

#-------------------------------------------------------------------------------
# Run DE on HC vs MCI/AD

# Load Seurat object
load(seurat_object)

# Run standard normalization
DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE)

# Set Ident to Diagnosis
s <- SetIdent(s, value = "Diagnosis")
print(table(s[["Diagnosis"]]))

run_de <- function(cell_type) {
  print(cell_type)
  cell_type_label <- gsub("/", "", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)

  # Find DEGs b/w MCI/AD and HC
  degs <-FindMarkers(object = subset(s, predicted.celltype.l2 == cell_type),
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
  write.csv(degs, paste0(output_dir, cell_type_label, "_MCIAD_vs_HC_degs.csv"))

  # Create volcano plot
  volcano_plot(degs, title = paste0("MC/AD vs HC in ", cell_type),
               file = paste0(output_dir, cell_type_label, "_MCIAD_vs_HC_volcano.pdf"),
               padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
}

# # Run DE on all celltypes
cell_types <- unique(s[["predicted.celltype.l2"]])[,1]
mclapply(cell_types, run_de, mc.cores = 12)

#------------------------------------------------------------------------------
# Create Upset plot (SuppFig 5F)

# Initialize empty list
sig_genes_ls <- list()

# Create list with sig genes for each cell type
for (cell_type in cell_types) {
  print(cell_type)
  cell_type_label <- gsub("/", "", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  
  # Load in degs
  tryCatch({
    degs <- read.csv(paste0(output_dir, cell_type_label, "_MCIAD_vs_HC_degs.csv"))
    
    # Identify sig genes
    sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
    print(head(sig_genes))
    
    # Add sig genes to list
    sig_genes_ls[[cell_type]] <- sig_genes$X
  }, error = function(e) {
    NULL
  })
  
}

sig_genes_ls
sig_genes_ls <- sig_genes_ls[ lapply(sig_genes_ls, length) > 0 ]

# Find number of genes in each set
num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
for (key in names(sig_genes_ls)) {
  num_degs <- rbind(num_degs, data.frame(key, length(sig_genes_ls[[key]])))
}

# Create upset plot and export
pdf(file = paste0(output_dir, "/upset_celltype.pdf"))
upset(
  fromList(sig_genes_ls),
  nsets = length(sig_genes_ls),
  order.by = "freq"
  )
dev.off()
