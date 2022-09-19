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
# Date: 02-14-2022
# Written by: Natalie Piehl
# Summary: Compare DE-SWAN to linear model results
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("ggpubr")
  library("ggthemes")
  library("ggrepel")
  library("Seurat")
  library("UpSetR")
  library("DEswan")
  library("scales")
})

# Initialize paths
lm_deg_path <- "path/to/linear_modeling/results/lm_age_stats.csv"
deswan_dir <- "path/to/deswan/store/"
output_dir <- "path/to/export/results"

# Source helper functions
source("../0_preprocessing/00_helper_functions.R")

# Generate output/store directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Specify thresholds
padj.thresh <- 0.01
beta.thresh <- 0.005
padj_deswan.thresh <- 1e-4
beta_deswan.thresh <- 1e-4

#------------------------------------------------------------------------------
# Load data

# Define cell types
cell_types <- c("CD4+ T Cells", "CD8+ T Cells", "T Regulatory Cells",
                "NK Cells", "DC", "CD14+/CD16-/CD68lo Monocytes",
                "CD14+/CD16+/CD68mid Monocytes", "CD14+/CD16+/CD68hi Monoctyes")


extract_data <- function(cell_type) {
  # Generate cell_type label for naming
  cell_type_label <- gsub("/","", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  
  # Load data
  load(paste0(deswan_dir, cell_type_label,
              "_DEswan_res_singlecell_10centers_10bootstrap20cells_bucket4"))
  
  # Convert to wide and calculate BH
  res_wide_b <- reshape.DEswan(res, parameter = 2, factor = "qt")
  res_wide_b$cell_type <- rep(cell_type, nrow(res_wide_b))
  res_wide_p <- reshape.DEswan(res, parameter = 1, factor = "qt")
  res_wide_q <- q.DEswan(res_wide_p, method = "BH")
  res_wide_q$cell_type <- rep(cell_type, nrow(res_wide_q))
  return(list(res_wide_b, res_wide_q))
}

# Extract data and format
res <- lapply(cell_types, extract_data)
b <- lapply(res, `[[`, 1)
q <- lapply(res, `[[`, 2)

# Load lm data
lm_data <- read.csv(lm_deg_path, stringsAsFactors = FALSE)

#------------------------------------------------------------------------------
# Run comparison and generate upset plot (Fig 2G)

# Iterate through up, down, and any direction
for (direction in c("up", "down", "any")) {
  # Initialize gene list
  gene_list <- list()
  
  # Iterate through cell types
  for (i in seq(length(cell_types))) {
    # Generate label for filenaming
    cell_type <- cell_types[i]
    cell_type_label <- gsub("/","", cell_type)
    cell_type_label <- gsub(" ", "_", cell_type_label)
    print(paste0("Analyzing ", cell_type))
    
    # Subset data for celltype of interest
    lm_sub <- lm_data[ which(lm_data$cell_type == cell_type), ]
    q_sub <- q[[i]]
    b_sub <- b[[i]]
    
    # Subset for sig genes
    if (direction == "up") {
      lm_sig <- lm_sub[which(lm_sub$BH <= padj.thresh & lm_sub$avg_log2FC >= beta.thresh),"gene"]
      q_66 <- q_sub[which(q_sub[,3] <= padj_deswan.thresh),c(1,3)]
      b_66 <- b_sub[which(b_sub[,3] >= beta_deswan.thresh),c(1,3)]
      deswan_66 <- intersect(q_66$variable, b_66$variable)
      q_72 <- q_sub[which(q_sub[,6] <= padj_deswan.thresh),c(1,6)]
      b_72 <- b_sub[which(b_sub[,6] >= beta_deswan.thresh),c(1,6)]
      deswan_72 <- intersect(q_72$variable, b_72$variable)
      q_78 <- q_sub[which(q_sub[,9] <= padj_deswan.thresh),c(1,9)]
      b_78 <- b_sub[which(b_sub[,9] >= beta_deswan.thresh),c(1,9)]
      deswan_78 <- intersect(q_78$variable, b_78$variable)
    } else if (direction == "down") {
      lm_sig <- lm_sub[which(lm_sub$BH <= padj.thresh & lm_sub$avg_log2FC <= -beta.thresh),"gene"]
      q_66 <- q_sub[which(q_sub[,3] <= padj_deswan.thresh),c(1,3)]
      b_66 <- b_sub[which(b_sub[,3] <= -beta_deswan.thresh),c(1,3)]
      deswan_66 <- intersect(q_66$variable, b_66$variable)
      q_72 <- q_sub[which(q_sub[,6] <= padj_deswan.thresh),c(1,6)]
      b_72 <- b_sub[which(b_sub[,6] <= -beta_deswan.thresh),c(1,6)]
      deswan_72 <- intersect(q_72$variable, b_72$variable)
      q_78 <- q_sub[which(q_sub[,9] <= padj_deswan.thresh),c(1,9)]
      b_78 <- b_sub[which(b_sub[,9] <= -beta_deswan.thresh),c(1,9)]
      deswan_78 <- intersect(q_78$variable, b_78$variable)
    } else {
      lm_sig <- lm_sub[which(lm_sub$BH <= padj.thresh & abs(lm_sub$avg_log2FC) >= beta.thresh),"gene"]
      q_66 <- q_sub[which(q_sub[,3] <= padj_deswan.thresh),c(1,3)]
      b_66 <- b_sub[which(abs(b_sub[,3]) >= beta_deswan.thresh),c(1,3)]
      deswan_66 <- intersect(q_66$variable, b_66$variable)
      q_72 <- q_sub[which(q_sub[,6] <= padj_deswan.thresh),c(1,6)]
      b_72 <- b_sub[which(abs(b_sub[,6]) >= beta_deswan.thresh),c(1,6)]
      deswan_72 <- intersect(q_72$variable, b_72$variable)
      q_78 <- q_sub[which(q_sub[,9] <= padj_deswan.thresh),c(1,9)]
      b_78 <- b_sub[which(abs(b_sub[,9]) >= beta_deswan.thresh),c(1,9)]
      deswan_78 <- intersect(q_78$variable, b_78$variable)
    }
    
    # Define list of sig genes and number of sig genes for each cell type
    sets <- list(LM = lm_sig, DESWAN_66 = deswan_66, DESWAN_72 = deswan_72, DESWAN_78 = deswan_78)
    gene_list[[cell_type]] <- sets
    sets <- sets[ lapply(sets, length) > 0 ]
    
    # Get order of sets for coloring
    num_degs <- sapply(sets, length)
    sets <- sets[order(-num_degs)]
    
    # Define colors
    color_vector <- mapvalues(
      names(sets),
      from = c("LM", "DESWAN_66", "DESWAN_72", "DESWAN_78"),
      to = c("skyblue1", "mistyrose", "hotpink", "brown1")
    )
    
    # Create upset plot and export
    pdf(file = paste0(output_dir, direction, "_", cell_type_label, "_upset.pdf"))
    upset(
      fromList(sets),
      nsets = length(sets),
      order.by = "freq",
      sets.bar.color = color_vector,
      text.scale = 2
    )
    dev.off()
  }
  
  # Write out gene lists
  sink(paste0(output_dir, "/deg_lists_", direction, ".txt"))
  print(gene_list)
  sink()
  saveRDS(gene_list, paste0(output_dir, "/deg_lists_", direction, ".rds"))
}