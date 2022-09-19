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
# Date: 02-28-2022
# Written by: Natalie Piehl
# Summary: Generate pseudobulk values
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

#-------------------------------------------------------------------------------
# Generate pseuodbulk values

# Load in seurat object
load(seurat_object)

bulk_by_id <- function(s, clonal = FALSE) {
  # Extract counts
  counts <- GetAssayData(object = s, assay = "RNA", slot = "counts")
  
  # Transpose counts and convert to matrix
  counts <- t(as.matrix(counts))
  
  # Extract age and ID per barcode
  meta_data <- s[[c("age", "ID", "Diagnosis", "sex", "clonal", "age_bin")]]
  
  # Merge rna data and age over barcode
  counts <- merge(counts, meta_data, by = 0)
  counts$age <- as.numeric(as.character(counts$age))
  
  # Sum counts across samples
  gene_names <- names(counts)[names(counts) %!in% c("Row.names", "age", "ID", "Diagnosis", "sex", "clonal", "age_bin")]
  
  if (clonal == FALSE) {
    count_sum <- counts %>%
      dplyr::group_by(ID) %>%
      dplyr::summarize(across(all_of(gene_names), sum))
    
    col <- 2
  } else {
    count_sum <- counts %>%
      dplyr::group_by(ID, clonal) %>%
      dplyr::summarize(across(all_of(gene_names), sum))
    count_sum <- count_sum[count_sum$clonal != "NA",]
    
    col <- 3
  }
  
  # LogNormalize counts
  data <- count_sum
  data[,col:ncol(data)] <- log1p(data[,col:ncol(data)]  / rowSums(data[,col:ncol(data)]) * 10000)
  
  # Add metadata back on
  data <- merge(data, unique(counts[,c("ID", "age", "Diagnosis", "sex", "age_bin")]), all.x = TRUE, all.y = FALSE)
  
  return(data)
}

#------------------------------------------------------------------------------
# Run bulking

# Run function on all cells
data_all <- bulk_by_id(s)
data_all_clonal <- bulk_by_id(s, clonal = TRUE)
data_all$cell_type <- rep("all", nrow(data_all))
data_all_clonal$cell_type <- rep("all", nrow(data_all_clonal))

# Isolate cell types
cell_types <- unique(s[["cluster_ident"]])[,1]
cell_types <- cell_types[ cell_types %!in% c("CD4+/CD8+ T Cells", "Undetermined")]

# Reserve s of all cell_types
s_all <- s

# Run in cell_type specific manner
for (cell_type in cell_types) {
  print(cell_type)
  
  # Subset cell_type of interest
  s <- subset(s_all, cluster_ident == cell_type)
  
  # Run bulking
  data <- bulk_by_id(s)
  data_clonal <- bulk_by_id(s, clonal = TRUE)
  data$cell_type <- rep(cell_type, nrow(data))
  data_clonal$cell_type <- rep(cell_type, nrow(data_clonal))
  
  # Append to running total
  data_all <- merge(data_all, data, all = TRUE)
  data_all_clonal <- merge(data_all_clonal, data_clonal, all = TRUE)
}

# Add clonal column to data all and merge
data_all$clonal <- rep(NA, nrow(data_all))
data <- merge(data_all, data_all_clonal, all = TRUE) %>%
  tidyr::drop_na(ID)

# Export data
save(data, file = paste0(output_dir, "bulk_expression"))

#-------------------------------------------------------------------------------
# Generate lists of genes expressed in at least 10% of cells per cell type

# Extract counts
data <- GetAssayData(object = s, assay = "RNA", slot = "counts") %>% as.matrix()

# Subset over 10%
data <- data[ which((rowSums(data) / ncol(data)) > .1), ]
over10_genes <- rownames(data)

# Remove junk
over10_genes <- over10_genes[-grep(pattern = "^RPS|^RPL|^MT-|^HLA-", x = over10_genes)]

# Export genes
write.csv(over10_genes, file.path(output_dir, "all_over10_genes.csv"))
save(over10_genes, file = file.path(output_dir, "all_over10_genes"))

# Get over 10% for each cell type
for (cell_type in unique(s[["cluster_ident"]])[,1]) {
  # Generate cell_type label for naming
  cell_type_label <- gsub("/","", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  
  # Subset cell type of interest
  s_sub <- subset(s, cluster_ident == cell_type)
  
  # Extract counts
  data <- GetAssayData(object = s_sub, assay = "RNA", slot = "counts") %>% as.matrix()
  
  # Subset over 10%
  data <- data[ which((rowSums(data) / ncol(data)) > .1), ]
  over10_genes <- rownames(data)
  
  # Remove junk
  over10_genes <- over10_genes[-grep(pattern = "^RPS|^RPL|^MT-|^HLA-", x = over10_genes)]
  
  # Export genes
  write.csv(over10_genes, paste0(output_dir, cell_type_label, "_over10_genes.csv"))
  save(over10_genes, file = paste0(output_dir, cell_type_label, "_over10_genes"))
}
