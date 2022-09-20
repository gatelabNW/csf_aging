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
# Date: 06-14-2022
# Written by: Natalie Piehl
# Summary: Look at APOE, APOC1, and PLTP in HC and CI
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("ggpubr")
  library("ggrepel")
  library("ggthemes")
  library("scales")
})

# Organize inputs
input_path <- "path/to/pseudbulk/data"
output_dir <- "path/to/export/results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Source helper functions
source("../0_preprocessing/00_helper_functions.R")

# Desired age range
age_min <- 62
age_max <- 82

# Define non-gene column names
non_gene_cols <- c("ID", "age", "Diagnosis", "sex", "age_bin", "cell_type")

#------------------------------------------------------------------------------
# Format data

# Load in data
load(input_path)

# Specify genes
genes <- c("APOE", "APOC1", "PLTP")

# Subset data
data <- data[ which(data$cell_type == "CD14+/CD16+/CD68hi Monoctyes"), ]
data <- data[ ,names(data) %in% c(non_gene_cols, genes) ] %>%
  dplyr::relocate(all_of(non_gene_cols))
data <- data[ which(data$age >= age_min & data$age <= age_max), ]

# Scale data
data_scaled <- data
data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)] <- scale(data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)])
data_scaled$age <- as.numeric(as.character(data_scaled$age))
data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)] <- data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)][ , colSums(is.na(data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)])) == 0]

# Convert to long
data_scaled <- data_scaled[ ,names(data) %in% c("age", "Diagnosis", genes) ]
data_long <- gather(data_scaled, gene, expression, -age, -Diagnosis)

#------------------------------------------------------------------------------
# Plot APOC1, APOE, and PLTP in CI and HC (Fig 2I)

# Add color column
colors <- c("#ff6d60", "#00fab7", "#ffb922")

# Add gene diagnosis column
data_long$gene <- factor(data_long$gene, levels = c("APOE", "APOC1", "PLTP"))
data_long$gene_diagnosis <- paste0(data_long$gene, "_", data_long$Diagnosis)

# Generate plot
plot <- ggplot(data_long, aes(x = age, y = expression, color = gene)) +
  theme_Publication_blank() +
  geom_line(aes(linetype=Diagnosis), stat="smooth", size = 2,
            method = "loess", span = 0.75, se = FALSE) +
  scale_color_manual(values = colors) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  guides(size = "none", alpha = "none") +
  ggtitle("CD14+/CD16+/CD68hi\nMonocytes") +
  theme(legend.position = "right",
        text = element_text(size=24),
        legend.key.size =  unit(0.5, "in"))

# Export plot
set_panel_size(plot, file = paste0(output_dir,
                                   "apoe_apoc1_pltp_loess_cd68hi_mono.pdf"),
               width = unit(5, "inch"), height = unit(3, "inch"))
