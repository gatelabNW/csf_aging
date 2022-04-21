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
# Date: 02-07-2022
# Written by: Natalie Piehl
# Summary: Generate LOESS curves of LM + MAST shared DEGs
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
library_vector <- c("plyr", "tidyverse", "ggrepel", "ggthemes", "Seurat", "grid",
                    "scales", "doMC", "UpSetR")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = library_vector)

# Initialize paths
mast_path <- "path/to/mast/age_bin/results/"
lm_path <- "path/to/lm/age/results"
output_dir <- "path/to/export/results"

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Specify thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25
beta.thresh <- 0.005

#------------------------------------------------------------------------------
# Data preparation

# List deg csvs
mast_files <- list.files(path = mast_deg_path, pattern = "*_advanced_vs_middle_degs.csv")

# Load lm data
lm_data <- read.csv(lm_deg_path, stringsAsFactors = FALSE)
lm_data$cell_type <- sapply(lm_data$cell_type, function(x) gsub("/", "", x))
lm_data$cell_type <- sapply(lm_data$cell_type, function(x) gsub(" ", "_", x))

# Initialize de_data
mast_data <- setNames(data.frame(matrix(ncol = 8, nrow = 0)),
                      c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "BH", "cell_type"))

for (file in mast_files) {
  # Read in data
  data <- read.csv(file.path(mast_deg_path, file), stringsAsFactors = FALSE) %>% dplyr::rename(gene = X)
  
  # Add cell-type column
  cell_type_label <- sub("_advanced_vs_middle_degs.csv", "", file)
  data$cell_type <- cell_type_label
  
  # Merge data together
  mast_data <<- merge(mast_data, data, all = TRUE)
}

#------------------------------------------------------------------------------
# Identify shared MAST + LM genes for each cell type

# Iterate through up and down
for (direction in c("up", "down")) {
  
  # Initialize gene list
  gene_list <- list()
  
  # Iterate through cell types
  for (cell_type in unique(lm_data$cell_type)) {
    print(paste0("Analyzing ", cell_type))
    
    # Subset data for celltype of interest
    lm_sub <- lm_data[ which(lm_data$cell_type == cell_type), ]
    mast_sub <- mast_data[ which(mast_data$cell_type == cell_type), ]
    
    # Subset for sig genes
    if (direction == "up") {
      lm_sig <- lm_sub[which(lm_sub$BH <= padj.thresh & lm_sub$avg_log2FC > beta.thresh),]
      mast_sig <- mast_sub[which(mast_sub$BH <= padj.thresh & mast_sub$avg_log2FC > lfc.thresh),]
    } else {
      lm_sig <- lm_sub[which(lm_sub$BH <= padj.thresh & lm_sub$avg_log2FC < -beta.thresh),]
      mast_sig <- mast_sub[which(mast_sub$BH <= padj.thresh & mast_sub$avg_log2FC < -lfc.thresh),]
    }
    
    # Make gene sets
    lm_all <- lm_sig$gene
    mast_all <- mast_sig$gene
    lm_unique <- setdiff(lm_all, mast_all)
    mast_unique <- setdiff(mast_all, lm_all)
    shared <- intersect(lm_all, mast_all)
    
    # Append genes to gene_list
    gene_list[[cell_type]] <- list(lm_unique = lm_unique,
                                   mast_unique = mast_unique,
                                   shared = shared)
  }
  
  # Write out gene lists
  sink(paste0(output_dir, "/deg_lists_", direction, ".txt"))
  print(gene_list)
  sink()
  saveRDS(gene_list, paste0(output_dir, "/deg_lists_", direction, ".rds"))
  
  if (direction == "up") {
    gene_list_up <<- gene_list
  } else {
    gene_list_down <<- gene_list
  }
}

#------------------------------------------------------------------------------
# Generate LOESS curves (Fig S3D)

# Specify desired age range
age_min <- 62
age_max <- 82

# Define non-gene column names
non_gene_cols <- c("ID", "age", "Diagnosis", "sex", "age_bin", "cell_type")

# Load in data
down <- readRDS(paste0(output_dir, "/deg_lists_down.rds"))
up <- readRDS(paste0(output_dir, "/deg_lists_up.rds"))

plot_degs <- function(cell_type, dir = "up") {
  # Define cell type label
  cell_type_label <- gsub("/", "", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  
  # Select for DEGs
  if (dir == "up") {
    degs <- up[[cell_type_label]]$shared
  } else {
    degs <- down[[cell_type_label]]$shared
  }
  
  # Skip if no degs
  if (length(degs) == 0) {
    return(NULL)
  }
  
  # Subset data
  data <- data[ which(data$cell_type == cell_type), ]
  data <- data[ ,names(data) %in% c(non_gene_cols, degs) ] %>%
    dplyr::relocate(all_of(non_gene_cols))
  data <- data[ which(data$Diagnosis == "HC" & data$age >= age_min & data$age <= age_max), ]
  
  # Scale data
  data_scaled <- data
  data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)] <- scale(data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)])
  data_scaled$age <- as.numeric(as.character(data_scaled$age))
  
  # Convert to long
  data_scaled <- data_scaled[ ,names(data) %in% c("age", degs) ]
  data_long <- gather(data_scaled, gene, expression, -age)

  # Generate plot
  plot <- ggplot(data_long, aes(x = age, y = expression, group = gene,
                                color = gene)) +
    theme_Publication_blank() +
    geom_line(stat="smooth", method = "loess", span = 0.75, se = FALSE) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    theme(legend.position = "right",
          aspect.ratio = 1)
  
  # Export plot
  set_panel_size(plot, file = paste0(output_dir, cell_type_label, "_", dir, "_MAST+LMdegs_loess.pdf"),
                 width = unit(4, "inch"), height = unit(4, "inch"))
}

# Define celltypes
celltypes <- c(
  "CD4+ T Cells", "CD8+ T Cells", "T Regulatory Cells",
  "NK Cells", "DC",
  "CD14+/CD16-/CD68lo Monocytes", "CD14+/CD16+/CD68mid Monocytes", "CD14+/CD16+/CD68hi Monoctyes"
)

# Plot DEGs for each celltype
sapply(celltypes, plot_degs)