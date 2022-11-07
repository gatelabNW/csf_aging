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
# Date: 09-19-2022
# Written by: Natalie Piehl
# Summary: Look at sig of AD risk genes
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("corrplot")
})

# Initialize paths
supervised_degs_dir <- "path/to/supervised_celltypes/MAST_DE_results"
manual_degs_dir_ <- "path/to/manual_celltypes/MAST_DE_results"
output_dir <- "path/to/export/results"

# Source helper functions
source("../0_preprocessing/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#------------------------------------------------------------------------------
# Format data

# Read in degs
deg_reader <- function(path) {
  # Read in degs
  degs <- read.csv(paste0(supervised_degs_dir, path))
  
  # Isolate celltype
  celltype <- gsub("_MCIAD_vs_HC_degs.csv", "", path)
  
  # Add celltype to degs
  degs$celltype <- rep(celltype, nrow(degs))
  
  return(degs)
}

# Define paths to degs
degs_paths <- list.files(path = supervised_degs_dir, pattern = "*_degs.csv")

# Generate merged degs
degs <- lapply(degs_paths, deg_reader)
degs <- do.call("rbind", degs)

# Define risk genes in order
risk_genes <- c("APOE", "BIN1", "MS4A6A", "PICALM", "CR1", "CLU", "TREM2",
                "ABCA7", "NYAP1", "PTK2B", "PLCG2", "SPI1", "SORL1", "HLA-DRB1",
                "CD2AP", "SLC24A4", "RIN3", "ADAMTS1", "ADAMTS4", "CASS4",
                "ADAM10", "FERMT2", "HAVCR2", "SCIMP", "CLNK", "ECHDC3",
                "TNIP1", "ABCA1", "CNTNAP2", "USP6NL", "INPP5D", "CD33", "ACE",
                "IQCK", "WWOX", "ABI3", "HESX1", "FHL2", "APH1B", "HS3ST1",
                "CHRNE", "CCDC6", "AGRN", "KAT8", "IL34")

# Select for risk genes
degs <- degs[which(degs$X %in% risk_genes),]
degs$X <- factor(degs$X, levels = risk_genes)

# Make sig col
degs <- dplyr::mutate(degs, sig = ifelse(BH <= .01, TRUE, FALSE))

# Isolate innate and adaptive
innate_types <- c("ASDC", "CD14_Mono", "CD16_Mono", "cDC1", "cDC2", "ILC",
                  "NK_CD56bright", "NK", "NK_Proliferating", "pDC")
adaptive_types <- c("B_intermediate", "B_memory",
                    "CD4_CTL", "CD4_Naive", "CD4_Proliferating", "CD4_TCM", "CD4_TEM",
                    "CD8_Naive", "CD8_Proliferating", "CD8_TCM", "CD8_TEM",
                    "dnT", "gdT", "MAIT", "Treg")
innate <- degs[which(degs$celltype %in% innate_types),]
adaptive <- degs[which(degs$celltype %in% adaptive_types),]
adaptive <- merge(adaptive, unique(innate[,"X",drop=FALSE]), all.y = TRUE)
adaptive[which(is.na(adaptive$celltype)), "avg_log2FC"] <- 0
adaptive[which(is.na(adaptive$celltype)), "celltype"] <- "Treg"

#------------------------------------------------------------------------------
# Generate plots on supervised celltypes (Supp Fig 7A)

# Generate plot
p <- ggplot(innate, aes(x = celltype, y = X,
                      fill = avg_log2FC)) +
  geom_tile() +
  geom_point(aes(size=ifelse(sig, "dot", "no_dot")), shape = 8) +
  scale_size_manual(values=c(dot=1, no_dot=NA), guide="none") +
  scale_fill_gradient2(low = "blue", high = "red", limits = c(-1.2,.85)) +
  theme_Publication_blank() +
  scale_x_discrete(position = 'top') +
  scale_y_discrete(limits = rev) +
  theme(legend.position = "bottom",
        legend.direction="vertical",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0, hjust=0))
p

# Export plot
set_panel_size(p, file = paste0(output_dir, "innate_ad_risk_genes.pdf"),
               width = unit(4, "in"), height = unit(7, "in"))

# Generate plot
p <- ggplot(adaptive, aes(x = celltype, y = X,
                        fill = avg_log2FC)) +
  geom_tile() +
  geom_point(aes(size=ifelse(sig, "dot", "no_dot")), shape = 8) +
  scale_size_manual(values=c(dot=1, no_dot=NA), guide="none") +
  scale_fill_gradient2(low = "blue", high = "red", limits = c(-1.2,.85)) +
  theme_Publication_blank() +
  scale_x_discrete(position = 'top') +
  scale_y_discrete(limits = rev) +
  theme(legend.position = "bottom",
        legend.direction="vertical",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0, hjust=0))
p

# Export plot
set_panel_size(p, file = paste0(output_dir, "adaptive_ad_risk_genes.pdf"),
               width = unit(6, "in"), height = unit(7, "in"))

#------------------------------------------------------------------------------
# Generate plots on manual celltypes (Fig 5A)

# Read in degs
deg_reader <- function(path) {
  # Read in degs
  degs <- read.csv(paste0(manual_degs_dir, path))
  
  # Isolate celltype
  celltype <- gsub("_MCIAD_vs_HC_degs.csv", "", path)
  
  # Add celltype to degs
  degs$celltype <- rep(celltype, nrow(degs))
  
  return(degs)
}

# Define paths to degs
degs_paths <- list.files(path = manual_degs_dir, pattern = "*_degs.csv")

# Generate merged degs
degs <- lapply(degs_paths, deg_reader)
degs <- do.call("rbind", degs)

# Select for risk genes
degs <- degs[which(degs$X %in% risk_genes),]
degs$X <- factor(degs$X, levels = risk_genes)

# Make sig col
degs <- dplyr::mutate(degs, sig = ifelse(BH <= .01, TRUE, FALSE))

# Generate plot
p <- ggplot(degs, aes(x = celltype, y = X,
                        fill = avg_log2FC)) +
  geom_tile() +
  geom_point(aes(size=ifelse(sig, "dot", "no_dot")), shape = 8) +
  scale_size_manual(values=c(dot=1, no_dot=NA), guide="none") +
  scale_fill_gradient2(low = "blue", high = "red") +
  theme_Publication_blank() +
  scale_x_discrete(position = 'top') +
  scale_y_discrete(limits = rev) +
  theme(legend.position = "bottom",
        legend.direction="vertical",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0, hjust=0))
p

# Export plot
set_panel_size(p, file = paste0(output_dir, "manual_celltypes_ad_risk_genes.pdf"),
               width = unit(4, "in"), height = unit(7, "in"))