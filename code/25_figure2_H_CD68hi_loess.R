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
# Date: 03-22-2022
# Written by: Natalie Piehl
# Summary: Generate LOESS curves of MAST+LM DEGs in CD68hi monocytes
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
})

# Initialize paths
input_path <- "path/to/pseudbulk/data"
input_clustering <- "path/to/cd68hi_monocyte/loess/clustering"
output_dir <- "path/to/export/results"

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Load and format data

# desired age range
age_min <- 62
age_max <- 82

# Define non-gene column names
non_gene_cols <- c("ID", "age", "Diagnosis", "sex", "age_bin", "cell_type")

# Load in data
load(input_path)
cut <- read.csv(input_clustering)

# Identify yellow cluster genes
clust_genes <- cut[which(cut[,2] == 2), "X"]
genes <- c("APOE", "APOC1", "PLTP")

# Subset data
data <- data[ which(data$cell_type == "CD14+/CD16+/CD68hi Monoctyes"), ]
data <- data[ ,names(data) %in% c(non_gene_cols, clust_genes, genes) ] %>%
  dplyr::relocate(all_of(non_gene_cols))
data <- data[ which(data$Diagnosis == "HC" & data$age >= age_min & data$age <= age_max), ]

# Scale data
data_scaled <- data
data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)] <- scale(data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)])
data_scaled$age <- as.numeric(as.character(data_scaled$age))
data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)] <- data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)][ , colSums(is.na(data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)])) == 0]

# Center samples
data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)] <- apply(data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)], 1,
                                                                 function(x) {scale(x, scale = FALSE)})

# Convert to long
data_scaled <- data_scaled[ ,names(data) %in% c("age", clust_genes) ]
data_long <- gather(data_scaled, gene, expression, -age)

#------------------------------------------------------------------------------
# Generate LOESS plot (Fig 2H)

# Add color column
colors <- c("#ff6d60", "#00fab7", "#ffb922", "grey30")
data_long$color <- rep("Remaining\ncluster genes", nrow(data_long))
data_long$color[ which(data_long$gene == "APOE") ] <- "APOE"
data_long$color[ which(data_long$gene == "APOC1") ] <- "APOC1"
data_long$color[ which(data_long$gene == "PLTP") ] <- "PLTP"
data_long$color <- factor(data_long$color, 
                          levels = c("APOE", "APOC1", "PLTP", "Remaining\ncluster genes"))

# Add size column
sizes <- c(0.3, 3)
data_long$size <- rep(sizes[1], nrow(data_long))
data_long$size[ which(data_long$gene %in% genes) ] <- sizes[2]

# Add alpha column
alphas <- c(0.05, 0.8)
data_long$alpha <- rep(alphas[1], nrow(data_long))
data_long$alpha[ which(data_long$gene %in% genes) ] <- alphas[2]

# Order genes of interests first
data_long$gene <- factor(data_long$gene, levels = rev(union(genes, unique(data_long$gene))))

# Generate plot
plot <- ggplot(data_long, aes(x = age, y = expression, group = gene, 
                              color = color, size = as.factor(size),
                              alpha = as.factor(alpha))) +
  theme_Publication_blank() +
  scale_color_manual(values = colors) +
  scale_size_manual(values = sizes) +
  scale_alpha_manual(values = alphas) +
  geom_line(stat="smooth", method = "loess", span = 0.75, se = FALSE) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  guides(size = "none", alpha = "none") +
  theme(legend.position = "right",
        aspect.ratio = 1)
plot

# Export plot
set_panel_size(plot, file = paste0(output_dir, "all_cluster_genes_loess_cd68hi_mono.pdf"),
               width = unit(4, "inch"), height = unit(4, "inch"))
