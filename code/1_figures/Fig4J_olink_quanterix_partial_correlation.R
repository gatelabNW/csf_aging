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
# Date: 08-22-2022
# Written by: Natalie Piehl
# Summary: Generate partial correlations b/w Olink CXCL16 and Quanterix data
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("ppcor")
  library("xlsx")
  library("ggnewscale")
  library("corrplot")
})

# Initialize paths
protein_path <- "path/to/olink-quanterix/protein/data/"
output_dir <- "path/to/export/results"

# Source helper functions
source("../0_preprocessing/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Run correlation between genes and markers

# Load in data
data <- read.xlsx(protein_path, sheetIndex = 2)

# Convert Gender to binary
data$Gender <- plyr::mapvalues(data$Gender,
                         from = c("F", "M"),
                         to = c(0,1))
data$Gender <- as.numeric(data$Gender)

# Select markers and genes
markers <- c( "NFL",
             "GFAP", "UCHL1",
             "Ab40", "Ab42", "pTau181")
genes <- c("CXCL16", "PLTP")

# Define diagnoses
diagnoses <- list("HC", "MCI-AD", "AD", c("MCI-AD", "AD"))
diagnoses_labels <- c("HC", "MCI-AD", "AD", "CI")

# Function to run for each marker
corr_marker <- function(marker, gene, data, diagnosis) {
  
  # Subset data for diagnosis and marker
  data <- data[which(data$Diagnosis %in% diagnosis),
               c("Age", "Gender", gene, marker)]
  
  # Remove any incomplete rows
  data <- data[complete.cases(data),]
  
  # Run partial correlation
  results <- pcor(data, method = "spearman")
}

corr_diagnosis <- function(diagnosis, data, gene) {
  
  # Generate result for each marker
  results <- lapply(markers,  corr_marker,
                    diagnosis = diagnosis, 
                    data = data, 
                    gene = gene)
  
  return(results)
  
}

# Function to run for each gene
corr_gene <- function(data, gene) {
  
  # Generate result for each gene
  results <- lapply(diagnoses, corr_diagnosis, 
                    data = data, 
                    gene = gene)
  
  return(results)
}

# Run partial correlation
results <- lapply(genes, corr_gene, data = data)
pvals <- lapply(results, `[[`, 3)

#-------------------------------------------------------------------------------
# Format for plotting

# Extract R and pval for each comparison
extract_rpval <- function(element) {
  return(list(R = element$estimate[3,4],
              pval = element$p.value[3,4]))
}

# Run extraction
rpval <- list()
for (i in seq(2)) {
  rpval[[genes[i]]] <- sapply(results[[i]], function(list) {sapply(list, extract_rpval)}) %>% as.data.frame
}

# Rbind and add gene names
rpval <- bind_rows(rpval)

# Add type column and rename columns and add genes
names(rpval) <- diagnoses_labels
rpval$gene <- rep(genes, each = length(markers)*2)
rpval$type <- rep(c("R","pval"), length(markers))
rpval$marker <- rep(markers, each = 2)

# Reshape data
rpval_shaped <- gather(rpval, Diagnosis, val, HC:CI)
rpval_shaped <- spread(rpval_shaped, type, val)

# Add levels
rpval_shaped$Diagnosis <- factor(rpval_shaped$Diagnosis,
                                 levels = diagnoses_labels)
rpval_shaped$marker <- factor(rpval_shaped$marker,
                                 levels = markers)

# Unlist
rpval_shaped$pval <- unlist(rpval_shaped$pval)
rpval_shaped$R <- unlist(rpval_shaped$R)

# Make gene_diagnosis column
group_order <- as.vector(sapply(genes, function(gene) {
                sapply(diagnoses_labels, function(diag) {paste(gene, diag)} ) }))
rpval_shaped$gene_diagnosis <- paste(rpval_shaped$gene, rpval_shaped$Diagnosis)
rpval_shaped$gene_diagnosis <- factor(rpval_shaped$gene_diagnosis,
                                      levels = group_order)

#-------------------------------------------------------------------------------
# Make correlation heatmap (Fig 4J)

# Shape into matrix
rpval_shaped <- rpval_shaped[ order(rpval_shaped$marker), ]

r_mat <- matrix(rpval_shaped$R, nrow = length(unique(rpval_shaped$gene))*4)
colnames(r_mat) <- unique(rpval_shaped$marker)
rownames(r_mat) <- unique(rpval_shaped$gene_diagnosis)

p_mat <- matrix(rpval_shaped$pval, nrow = length(unique(rpval_shaped$gene))*4)
colnames(p_mat) <- unique(rpval_shaped$marker)
rownames(p_mat) <- unique(rpval_shaped$gene_diagnosis)

# Reorder
r_mat <- r_mat[c(3,4,1,2,7,8,5,6),]
p_mat <- p_mat[c(3,4,1,2,7,8,5,6),]

# Generate plot
pdf(paste0(output_dir, "Olink_pcor_markers_corrplot.pdf"),
    width = 4, height = 8)
corrplot(r_mat, tl.col = "black",
         p.mat = p_mat,
         insig = 'label_sig',
         pch.cex = 1,
         sig.level = c(0.001, 0.01, 0.05),
         method = "square",
         col=colorRampPalette(c("blue","white","red"))(200))
dev.off()
