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
# Date: 11-05-2021
# Written by: Natalie Piehl
# Summary: Generate linear models of expression ~ age + sex per cell type
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
  library("UpSetR")
  library("car")
  library("lsr")
})

# Initialize paths
seurat_object <- "path/to/seurat_object"
output_dir <- "path/to/export/results"

# Source helper functions
source("code/00_helper_functions.R")

# Generate output directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Specify thresholds
threshold.padj <- 0.01
threshold.beta <- 0.005

# Set seed
set.seed(123)

#------------------------------------------------------------------------------
# Load and normalize data

# Load Seurat object
load(seurat_object)

# Select for HC samples only
s <- subset(s, subset = Diagnosis == "HC")

# Normalize with standard processing (log transformation)
DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst",
                       verbose = FALSE)
s

#------------------------------------------------------------------------------
# Generate lm of expression ~ age + sex For each cell type

# List cell types
cell_types <- unique(s[["cluster_ident"]])[,1]
cell_types <- cell_types[ cell_types %!in% c("CD4+/CD8+ T Cells", "Undetermined")]

# Reserve s object of all cell types
s_all <- s

# Initialize stats output for each LM
stats <- data.frame(matrix(ncol = 4, nrow = 0))

for (cell_type in cell_types) {
  # Format cell type label
  cell_type_label <- gsub("/", "", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  
  # Subset for cells of interest
  s <- subset(s_all, cluster_ident == cell_type)
  
  # Extract normalized RNA data and demographics
  rna_data <- GetAssayData(object = s, assay = "RNA", slot = "data")
  demographics <- s[[c("age", "sex", "ID")]]
  
  # Remove ribosomal, mitochondrial, and HLA genes
  rna_data <-
    rna_data[-grep(pattern = "^RPS|^RPL|^MT-|^HLA-", x = rownames(rna_data)), ]
  
  # Transpose rna data and convert to matrix
  rna_data <- t(as.matrix(rna_data))
  
  # Remove genes not expressed in at least 10% of cells
  rna_data <- rna_data[ ,which(colSums(rna_data != 0) / nrow(rna_data) > .1) ]
  
  # Merge RNA data and demographics together
  rna_data <- merge(rna_data, demographics, by = 0)
  rna_data$age <- as.numeric(as.character(rna_data$age))
  
  # Initialize beta coefficient + p-value df
  age_stats_df <- data.frame(matrix(ncol = 3, nrow = 0))
  
  for (gene in colnames(rna_data)) {
    # Move to next iteration if column is age
    if (gene == "age" |
        gene == "ID" |
        gene == "sex" |
        gene == "Row.names") {
      next
    }
    
    # Subset gene of interest
    tmp <- rna_data[, c(gene, "age", "sex", "ID")]
    
    # Generate linear model
    lm <- lm(reformulate(c("age", "sex"), gene), data = tmp)
    
    # Run ANOVA
    anova <- Anova(lm, type = "II")
    
    # Isolate age coefficient and pval
    pval <- anova["age", "Pr(>F)"]
    beta <- coef(summary(lm))["age", "Estimate"]
    df <- data.frame(gene, beta, pval)
    
    # Add stats to df
    age_stats_df <- rbind(age_stats_df, df)
  }
  # Rename columns
  colnames(age_stats_df) <- c("gene", "avg_log2FC", "pval")
  rownames(age_stats_df) <- age_stats_df$gene
  
  # BH correct p values
  age_stats_df$BH <- p.adjust(age_stats_df$pval, method = "BH")
  
  #------------------------------------------------------------------------------
  # Plot lm results (Fig 1G)
  
  # Generate volcano plot
  volcano_plot(age_stats_df, 
               title = paste0("LM Age in ", cell_type),
               file = paste0(output_dir, cell_type_label, "_lm_age_volcano.pdf"),
               x_title = "beta",
               padj.thresh = threshold.padj,
               lfc.thresh = threshold.beta)
  
  # Save stats to cumulative df
  age_stats_df <- age_stats_df[order(age_stats_df$BH),]
  age_stats_df$cell_type <- rep(cell_type, nrow(age_stats_df))
  stats <- rbind(stats, age_stats_df)
}

# Write out results
write.csv(stats, paste0(output_dir, "/lm_age_stats.csv"))

#------------------------------------------------------------------------------
# Create Upset plot (Fig 1F)

# Isolate sig genes
sig_genes <- stats[which(stats$BH <= threshold.padj & abs(stats$avg_log2FC) >= threshold.beta),]

# Define list of sig genes and number of sig genes for each cell type
sets <- list()
num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
for (cell_type in unique(sig_genes$cell_type)) {
  sets[[cell_type]] <-
    sig_genes[which(sig_genes$cell_type == cell_type), "gene"]
  num_degs <-
    rbind(num_degs, data.frame(cell_type, length(sets[[cell_type]])))
}

# Get order of sets for coloring
num_degs <- num_degs[order(-num_degs[, 2]),]

# Define colors
color_vector <- num_degs[, 1]
color_vector <- mapvalues(
  color_vector,
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

# Create upset plot and export
pdf(file = paste0(output_dir, "/upset_celltype.pdf"))
upset(
  fromList(sets),
  nsets = length(sets),
  order.by = "freq",
  sets.bar.color = color_vector
)
dev.off()
