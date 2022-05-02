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
# Summary: Calculate Spearman correlation between HC and CI expression over age
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
  library("scales")
  library("pheatmap")
})

# Initialize paths
data_path <- "path/to/pseudobulk_data"
clust_path <- "path/to/loess/clustering"
output_dir <- "path/to/export/results"

# Source helper functions
source("code/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Run correlation (Fig 3C)

# Define non-gene column names
non_gene_cols <- c("ID", "age", "Diagnosis", "sex",
                   "age_bin", "clonal", "cell_type")

# Define type of correlation
corr_type <- "spearman"

# Define correlation function
run_corr <- function(data, clust, i) {
  # Identify cluster genes
  genes <- rownames(clust[ which(clust[,1] == i), ,drop = FALSE])
  data <- data[,c(non_gene_cols, genes)]
  
  # Subset on HC and MCIAD
  hc <- data[which(data$Diagnosis == "HC"),]
  mciad <- data[which(data$Diagnosis == "CI"),]
  
  # Average samples of same age
  hc <- hc %>%
    group_by(age) %>%
    summarise_at(names(hc)[names(hc) %!in% non_gene_cols], mean) %>%
    arrange(age)
  
  mciad <- mciad %>%
    group_by(age) %>%
    summarise_at(names(mciad)[names(mciad) %!in% non_gene_cols], mean) %>%
    arrange(age)
  
  # Convert to vectors
  corr <- c()
  for (col in colnames(hc)[colnames(hc) %!in% non_gene_cols]) {
    hc_vec <- as.vector(unlist(hc[,col]))
    mciad_vec <- as.vector(unlist(mciad[,col]))
    corr_tmp <- cor(hc_vec, mciad_vec, method = corr_type)
    corr <- c(corr, corr_tmp)
  }
  return(mean(corr, na.rm = TRUE))
}

# Define function to run correlation on each celltype
corr_on_celltype <- function(cell_type) {
  # Generate cell_type label for naming
  cell_type_label <- gsub("/","", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  
  # Load data
  data <- read.csv(paste0(data_path, cell_type_label, "/",
                          cell_type_label, "_data_scaled.csv"))
  clust <- read.csv(paste0(clust_path, cell_type_label, "/",
                           cell_type_label, "_hclust_cut.csv"), row.names = 1)
  clust <- clust[,1, drop = FALSE]
  
  # Subset on shared ages b/w HC and MCIAD
  shared_ages <- base::intersect(data[which(data$Diagnosis == "HC"),"age"],
                                 data[which(data$Diagnosis == "CI"),"age"])
  data <- data[which(data$age %in% shared_ages),]
  
  # Isolate 12 clusters
  corr_vec <- sapply(1:12, run_corr, data = data, clust = clust)
  
  return(corr_vec)
}

# Run correlation
cell_types <- c("CD4+ T Cells", "CD8+ T Cells", "T Regulatory Cells",
                "NK Cells", "DC", "CD14+/CD16-/CD68lo Monocytes",
                "CD14+/CD16+/CD68mid Monocytes", "CD14+/CD16+/CD68hi Monoctyes",
                "B Cells", "Plasma Cells")
corr_df <- data.frame(lapply(cell_types, corr_on_celltype))
names(corr_df) <- cell_types
rownames(corr_df) <- 1:12

# Set plotting zlim
lim <- 0.15

# Plot heatmap
pheatmap(t(corr_df), cluster_rows = FALSE, cluster_cols = FALSE,
         breaks = seq(-lim,lim,length.out = 30),
         main = paste0("    ", corr_type, " corr. b/w HC & CI"),
         cellwidth = 20, cellheight = 20,
         angle_col = 0, fontsize = 10,
         fontsize_row = 10, fontsize_col = 10,
         color = colorRampPalette(c("blue", "white", "red"))(30),
         show_rownames = TRUE, show_colnames = TRUE,
         filename = paste0(output_dir, "/HC_CI_cluster_",
                           corr_type, "_correlation_heatmap.pdf"))