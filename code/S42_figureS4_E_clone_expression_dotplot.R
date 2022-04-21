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
# Date: 03-04-2022
# Written by: Natalie Piehl
# Summary: Generate dotplot of genes of interest across HC and CI clones
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
library_vector <- c("plyr", "tidyverse", "ggrepel", "ggthemes", "Seurat",
                    "grid", "ggpubr", "data.table")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = library_vector)

# Initialize paths
seurat_object <- "path/to/seurat_object/"
output_dir <- "path/to/export/results"

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set seed
set.seed(123)

#------------------------------------------------------------------------------
# Data preparation

# Load data
load(seurat_object)

# Normalize with standard processing (log transformation)
DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst",
                       verbose = FALSE)
s

# Define cell type
cell_type <- "CD8+ T Cells"
cell_type_label <- gsub("/","", cell_type)
cell_type_label <- gsub(" ", "_", cell_type_label)


# Isolate biggest clones from HC and MCI/AD
extract_clones <- function(diagnosis, cell_type) {
  # Extract clonotypes
  clonotypes <- s[[c("ID", "trb_cdr3s", "frequency_filtered", "age", "Diagnosis", "cluster_ident")]] %>%
    dplyr::filter(Diagnosis == diagnosis & cluster_ident == cell_type) %>%
    tidyr::drop_na(trb_cdr3s) %>%
    dplyr::add_count(trb_cdr3s, sort = TRUE) %>%
    dplyr::group_by(trb_cdr3s) %>%
    dplyr::mutate(IDs = paste(ID, sep = "|")) %>%
    dplyr::mutate(ages = paste(age, sep = "|")) %>%
    dplyr::select(trb_cdr3s, n, IDs, ages) %>%
    dplyr::distinct()
  clonotypes <- clonotypes[1:5,]
  
  # Isolate singletons
  singletons <- s[[c("ID", "trb_cdr3s", "frequency_filtered", "Diagnosis", "cluster_ident")]] %>%
    dplyr::filter(Diagnosis == diagnosis & cluster_ident == cell_type & frequency_filtered == 1) %>%
    dplyr::slice_sample(n = sample_num) %>%
    dplyr::rename(n = frequency_filtered) %>%
    dplyr::select(trb_cdr3s, n)
  
  # Define sampling column
  singletons$sampling <- rep(1:5, each = 10)
  
  # Group singletons
  singletons <- singletons %>%
    dplyr::group_by(sampling) %>%
    dplyr::summarize(trb_cdr3s = paste(trb_cdr3s, collapse = "|"))
  
  # Merge clones and singletons
  clonotypes <- dplyr::full_join(clonotypes, singletons)
  
  return(clonotypes)
}

# Run function
clones <- lapply(c("HC", "MCI/AD"), extract_clones, cell_type = cell_type)

#------------------------------------------------------------------------------
# Find avg gene expression of each clone

# Extract normalized data
data <- GetAssayData(s, assay = "RNA", slot = "data")

# Define genes of interest
genes <- c("CXCR4", "CXCR3", "CXCR6", "IL32", "CD69", "NFKBIA")

extract_expression <- function(cdr3b, diagnosis, cell_type) {
  # Find average expression per clonotype
  cdr3b <- unlist(strsplit(cdr3b, split = "|", fixed = TRUE))
  barcodes <- s[[c("trb_cdr3s", "cluster_ident", "Diagnosis")]] %>%
    dplyr::filter(trb_cdr3s %in% cdr3b & Diagnosis == diagnosis & cluster_ident == cell_type)
  barcodes <- rownames(barcodes)
  
  # Avg data and find percentage of expressing cells
  if (length(barcodes) > 1) {
    clone_data <- data[genes, barcodes]
    clone_avg <- rowMeans(clone_data)
    pct_cells <- rowSums(clone_data != 0) / ncol(clone_data)
  } else {
    clone_data <- data[genes, barcodes, drop = FALSE]
    clone_avg <- clone_data[,1]
    pct_cells <- rep(1, length(genes))
    pct_cells[which(clone_avg == 0)] <- 0
    names(pct_cells) <- names(clone_avg)
  }
  
  # Add meta data
  clone_avg$type <- "expression"
  pct_cells$type <- "pct_cells"
  data <- rbind(data.frame(clone_avg), data.frame(pct_cells))
  if (length(cdr3b) == 1) {
    data$trb_cdr3s <- cdr3b
  } else {
    data$trb_cdr3s <- "NC Bootstrap"
  }
  data$diagnosis <- diagnosis
  return(data)
}

# Generate data
clone_expression_hc <- lapply(clones[[1]]$trb_cdr3s, extract_expression,
                              cell_type = cell_type, diagnosis = "HC")
clone_expression_mciad <- lapply(clones[[2]]$trb_cdr3s, extract_expression,
                                 cell_type = cell_type, diagnosis = "MCI/AD")
hc_df <- rbindlist(clone_expression_hc)
hc_df <- dplyr::left_join(hc_df, clones[[1]], by = "trb_cdr3s")
mciad_df <- rbindlist(clone_expression_mciad)
mciad_df <- dplyr::left_join(mciad_df, clones[[2]], by = "trb_cdr3s")
df <- dplyr::full_join(hc_df, mciad_df) %>% data.frame()

# Scale gene expression across clonotypes
df[which(df$type == "expression"), 1:length(genes)] <- scale(df[which(df$type == "expression"), 1:length(genes)])
df$trb_cdr3s[ which(df$trb_cdr3s == "NC Bootstrap")] <- paste0(df$trb_cdr3s[ which(df$trb_cdr3s == "NC Bootstrap")],
                                                               " ", rep(1:5, each = 2, times = 2))
df$sampling <- NULL

# Convert to long for plotting
df_long <- gather(df, gene, expression, -type, -trb_cdr3s, -diagnosis, -n, -IDs, -ages)
df_long <- dplyr::inner_join(df_long[which(df_long$type == "expression"),],
                             df_long[which(df_long$type == "pct_cells"),],
                             by = c("trb_cdr3s", "diagnosis", "gene", "n", "IDs", "ages"))
df_long[,c("type.x", "type.y")] <- NULL
names(df_long)[ names(df_long) %in% c("expression.x", "expression.y")] <- c("expression", "pct_cells")

#------------------------------------------------------------------------------
# Generate exports (Fig S4E)

# Export data
write.csv(df_long, paste0(output_dir, cell_type_label, "_topclones_and_bootstrap_singletons_HCvsMCIAD_avgexpression.csv"))

# Order and label CDR3bs
df_long$n[ which(is.na(df_long$IDs)) ] <- 10
df_long$trb_cdr3s <- factor(df_long$trb_cdr3s, levels = unique(df$trb_cdr3s))
df_long$trb_cdr3s_label <- paste0(df_long$diagnosis, ": ", df_long$trb_cdr3s, " (", df_long$n, ")")
df_long$trb_cdr3s_label <- factor(df_long$trb_cdr3s_label, levels = rev(unique(df_long$trb_cdr3s_label)))
df_long$gene <- factor(df_long$gene, levels = genes)

expression_lim <- 1
df_long$expression[ df_long$expression > expression_lim] <- expression_lim
df_long$expression[ df_long$expression < -expression_lim] <- -expression_lim

# Generate plot
p <- ggplot(df_long, aes(x = gene, y = trb_cdr3s_label, color = expression, size = pct_cells)) +
  geom_point() +
  theme_Publication_blank() +
  geom_hline(yintercept = 10.5, color = "gray") +
  geom_hline(yintercept = 5.5, color = "gray") +
  geom_hline(yintercept = 15.5, color = "gray") +
  scale_x_discrete(position = 'top') +
  scale_y_discrete(position = 'right') +
  scale_color_gradient(low = "darkorchid4", high = "yellow",
                       limits = c(-expression_lim, expression_lim)) +
  scale_size_continuous(range = c(2, 8)) +
  labs(title = paste0(cell_type, " most expanded clones")) +
  theme(legend.position = "left",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank())

set_panel_size(p, file = paste0(output_dir, cell_type_label, "_topclones_and_bootstrap_singletons_gex_HCvsMCIAD_scatterplot.pdf"),
               width = unit(3.5, "in"), height = unit(6, "in"))