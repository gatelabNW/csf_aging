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
# Date: 02-14-2022
# Written by: Natalie Piehl
# Summary: Run DE-SWAN on healthy cells
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
  library("DEswan")
  library("car")
  library("lsr")
})

# Initialize paths
seurat_object <- "path/to/seurat_object"
output_dir <- "path/to/export/results"
store_dir <- "path/to/store/objects"

# Source helper functions
source("code/00_helper_functions.R")

# Generate output/store directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(store_dir, showWarnings = FALSE, recursive = TRUE)

# Specify thresholds
padj_deswan.thresh <- 1e-4
beta_deswan.thresh <- 1e-4

# Set seed
set.seed(123)

#-------------------------------------------------------------------------------
# Data normalizing and set up

# Load Seurat object
load(seurat_object)

# Run standard normalization
DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE)

# Subset healthy and aged 62-82
s <- subset(s, Diagnosis == "HC")
s <- subset(s, age %in% seq(62, 82))

#------------------------------------------------------------------------------
# Define bootstrap and DE-SWAN functions

bootstrap_cells <- function(s) {
  # Extract counts
  counts <- GetAssayData(object = s, assay = "RNA", slot = "counts")
  
  # Transpose counts and convert to matrix
  counts <- t(as.matrix(counts))
  meta_data <- s[[c("age", "ID", "Diagnosis", "sex", "clonal", "age_bin")]]
  
  # Merge rna data and age over barcode
  counts <- merge(counts, meta_data, by = 0)
  counts$age <- as.numeric(as.character(counts$age))
  names(counts)[names(counts) == "Row.names"] <- "barcode"
  
  # Find minimum cell number per age
  barcode_age <- counts[,c("barcode", "age")]
  barcode_age_sampled <- barcode_age %>%
    group_by(age) %>%
    slice_sample(n = 200, replace = TRUE)
  
  # Merge data together
  counts <- dplyr::left_join(barcode_age_sampled, counts)
  counts$sampling <- paste(counts$age, rep(seq(10), each = 20, times = 20), sep = "_")
  
  # Sum and normalize bootstrap cells
  gene_names <- names(counts)[names(counts) %!in% c(non_gene_cols, "sampling", "barcode")]
  counts_grouped <- counts %>%
    group_by(sampling) %>%
    summarize_at(all_of(gene_names), sum)
  
  # Normalize counts with standard Seurat method
  data <- counts_grouped[,-c(1)] / rowSums(counts_grouped[,-c(1)]) * 10000 %>% log1p()
  
  # Merge data back with metadata
  data$sampling <- counts_grouped$sampling
  data$age <- unlist(sapply(counts_grouped$sampling, function(x) {sub("_.*", "", x)}))
  data$age <- as.numeric(as.character(data$age))
  return(data)
}

run_deswan <- function(cell_type) {
  # Generate cell_type label for naming
  cell_type_label <- gsub("/","", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  print(paste0("Analyzing ", cell_type))
  
  # Load over10_genes
  load(paste0("results/pseudobulk/main/out/", cell_type_label,"_over10_genes"))
  
  # Subset s
  if (cell_type != "all") {
    s_sub <- subset(s, cluster_ident == cell_type)
  } else {
    s_sub <- s
  }
  
  # Create bootstrapped data
  data <- bootstrap_cells(s_sub)
  
  # Subset for genes expressed in >10% of cells
  data <- data[, c("sampling", "age", over10_genes)] %>% dplyr::relocate(all_of(c("sampling", "age")))
  data <- data.frame(data)
  
  # Run DE-SWAN
  res <- DEswan(data.df = data[ ,(names(data) %!in% non_gene_cols)],
                qt = data[,2],
                buckets.size = 4,
                window.center = seq(62,82,2)) # covariates = data[,4]
  
  # Save results
  save(res, file = paste0(store_dir, cell_type_label, "_DEswan_res_singlecell_10centers_10bootstrap20cells_bucket4"))
}

# Run DE-SWAN on each cell type
cell_types <- unique(s[["cluster_ident"]])[,1]
cell_types <- cell_types[cell_types %!in% c("CD4+/CD8+ T Cells", "Undetermined", "B Cells", "Plasma Cells")]
lapply(cell_types, run_deswan)

#------------------------------------------------------------------------------
# Generate line plot of DEG number per cell type (Fig 2B)

# Initialize list of sig genes
sig_gene_list <- list()

for (cell_type in cell_types) {
  # Generate cell_type label for naming
  cell_type_label <- gsub("/","", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  
  # Load data
  load(paste0(store_dir,
              cell_type_label,
              "_DEswan_res_singlecell_10centers_10bootstrap20cells_bucket4"))
  
  # Convert to wide and calculate BH
  res_wide_p <- reshape.DEswan(res, parameter = 1, factor = "qt")
  res_wide_q <- q.DEswan(res_wide_p, method = "BH")
  
  # Find beta
  res_wide_b <- reshape.DEswan(res, parameter = 2, factor = "qt")
  
  # Subset BH and beta
  q <- lapply(res_wide_q, function(x) {
    res_wide_q$variable[which(x <= padj_deswan.thresh)]
    })
  b <- lapply(res_wide_b, function(x) {
    res_wide_b$variable[which((x >= beta_deswan.thresh) | (x <= -beta_deswan.thresh))]
    })
  q[["variable"]] <- NULL
  b[["variable"]] <- NULL
  
  # Subset sig genes
  sig_genes <- list()
  for (age in names(b)) {
    sig_genes[[age]] <- intersect(q[[age]], b[[age]])
  }
  sig_genes <- sapply(sig_genes, length)
  
  # Add to list
  sig_gene_list[[cell_type_label]] <- sig_genes
}

# Generate dataframe of number of sig genes
sig_gene_df <- data.frame(sig_gene_list)
names(sig_gene_df) <- cell_types
sig_gene_df$age <- gsub("X", "", rownames(sig_gene_df))
sig_gene_df_long <- gather(sig_gene_df, cell_type, num_genes, -age)
sig_gene_df_long$cell_type <- factor(sig_gene_df_long$cell_type, levels = cell_types)

# Generate color vector
color_vector <- mapvalues(
  cell_types,
  from = c(
    "CD4+ T Cells", "CD8+ T Cells", "CD4+/CD8+ T Cells", "CD14+/CD16+/CD68hi Monoctyes",
    "DC", "NK Cells", "CD14+/CD16-/CD68lo Monocytes", "T Regulatory Cells", "Plasma Cells",
    "CD14+/CD16+/CD68mid Monocytes", "B Cells", "Undetermined"
  ),
  to = c(
    "darkturquoise", "lawngreen", "dodgerblue", "red",
    "lightpink", "hotpink", "gold1", "navy", "violet",
    "darkorange", "darkviolet", "gray"
  )
)

# Generate plot
p <- ggplot(sig_gene_df_long, aes(x = age, y = num_genes, color = cell_type, group = cell_type)) +
  geom_line(size = 2) +
  scale_color_manual(values = color_vector) +
  theme_Publication_blank() +
  theme(legend.position = "right") +
  labs(title = "Num degs with q < 1e-04")

set_panel_size(p, 
               file = paste0(output_dir,
                                "allcelltypes_q1e-04_b1e-04_numgenes_10centers_10bootstrap20cells_bucket4.pdf"))

#------------------------------------------------------------------------------
# Generate manhattan plot of adjusted p-values at age 78 (Fig 2C)

# Function to load data for each cell type
extract_data <- function(cell_type) {
  # Generate cell_type label for naming
  cell_type_label <- gsub("/","", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  
  # Load data
  load(paste0(store_dir, cell_type_label, "_DEswan_res_singlecell_10centers_10bootstrap20cells_bucket4"))
  
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
q <- lapply(res, `[[`, 2)

# Extract age 78 data
q_78 <- lapply(q, `[`, c(1,9,11))
q_78_df <- rbind.fill(q_78)
names(q_78_df) <- c("gene", "padj", "cell_type")
q_78_df$cell_type <- factor(q_78_df$cell_type, cell_types)
q_78_df <- q_78_df[which(q_78_df$gene != "MALAT1"),]

# Add label for top three genes in each cell type
q_78_df <- q_78_df[order(q_78_df$cell_type, q_78_df$padj),]
indeces <- c()
for (cell_type in cell_types) {
  df <- q_78_df[which(q_78_df$cell_type == cell_type),]
  indeces <- c(indeces, which(q_78_df$cell_type == cell_type)[1:3])
}
q_78_df$label <- rep(FALSE, nrow(q_78_df))
q_78_df$label[indeces] <- TRUE
q_78_df$padj[q_78_df$padj < 10^-40] <- 10^-40

# Generate plot
p <- ggplot(q_78_df, aes(x = cell_type, y = -log10(padj), color = cell_type, size = label)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 41)) +
  geom_jitter(alpha = 0.8) +
  scale_size_manual(values = c(1,2)) +
  theme_Publication_blank() +
  labs(title = "DE-SWAN at age 78") +
  scale_color_manual(values = color_vector) +
  geom_text_repel(data = q_78_df[which(q_78_df$label == TRUE), ],
                  label = q_78_df[which(q_78_df$label == TRUE), "gene"],
                  inherit.aes = T,
                  color = 'black', size = 3, force = 3) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank())
p

set_panel_size(p, file = paste0(output_dir, "manhattan_padj_78_celltypes.pdf"),
               height = unit(3, "inch"), width = unit(5, "inch"))