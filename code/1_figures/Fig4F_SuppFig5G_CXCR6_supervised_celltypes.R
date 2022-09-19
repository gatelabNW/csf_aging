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
# Date: 04-18-2022
# Written by: Natalie Piehl
# Summary: Compare supervised annotation celltypes
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("car")
  library("FSA")
})

# Initialize paths
seurat_object <- "path/to/seurat_object/"
output_dir <- "path/to/export/results"

# Source helper functions
source("code/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#------------------------------------------------------------------------------
# Run stats on CXCL16 and CXCR6 across cell types

# Load in seurat object
load(seurat_object)

# Subset for T cells
s <- subset(s, cluster_ident %in% c("CD4+ T Cells",
                                    "CD8+ T Cells",
                                    "T Regulatory Cells"))

# Isolate data
data <- GetAssayData(s, slot = "data")
data <- as.matrix(data[c("CXCR6", "CXCL16"),]) %>% t()
data <- merge(data, s[[c("Diagnosis", "predicted.celltype.l2")]], by = 0)

# Generate test function
wilcox_test <- function(celltype, gene) {
  data_sub <- data[which(data$predicted.celltype.l2 == celltype),]
  pval <- tryCatch({
    wilcox <- wilcox.test(data_sub[which(data_sub$Diagnosis == "HC"), gene],
                          data_sub[which(data_sub$Diagnosis != "HC"), gene])
    return(wilcox$p.value)
  }, error = function(e) {
    return(1)
  })
  return(pval)
}

# Run Wilcox Rank Sum
wilcox_cxcr6 <- sapply(unique(data$predicted.celltype.l2), wilcox_test, gene = "CXCR6")
wilcox_cxcl16 <- sapply(unique(data$predicted.celltype.l2), wilcox_test, gene = "CXCL16")
wilcox_cxcr6_adj <- p.adjust(wilcox_cxcr6, method = "BH")
wilcox_cxcl16_adj <- p.adjust(wilcox_cxcl16, method = "BH")
results <- data.frame(celltype = unique(data$predicted.celltype.l2),
                      cxcr6 = wilcox_cxcr6_adj,
                      cxcl16 = wilcox_cxcl16_adj
)
write.csv(results, paste0(output_dir, "wilcox_ranksum_withBH_CXCR6-CXCL16.csv"))

#------------------------------------------------------------------------------
# Visualize expression on only CD8 T cells (Fig 4F)

# Violin Plot CXCR6
Idents(s) <- "predicted.celltype.l2"
cd8_subsets <- c("CD8 TEM", "CD8 TCM", 
                 "CD8 Naive", "CD8 Proliferating")
s_cd8 <- subset(s, predicted.celltype.l2 %in% cd8_subsets)

# Generate p values
pvals <- compare_means(CXCR6 ~ Diagnosis,
                       group.by = "predicted.celltype.l2",
                       method = "wilcox.test",
                       p.adjust.method = "BH",
                       data = data)
pvals_cd8 <- pvals[which(pvals$predicted.celltype.l2 %in% cd8_subsets),]
pvals_cd8$y.position <- rep(max(data$CXCR6), length(cd8_subsets))
pvals_cd8$group1 <- cd8_subsets
pvals_cd8$group2 <- cd8_subsets

p <- VlnPlot(s_cd8, features = "CXCR6",
             split.by = 'Diagnosis',
             cols = c("gray", "red"),
             pt.size = 0.5,
             log = FALSE) +
  scale_y_continuous(limits = c(NA,4)) +
  geom_boxplot(width=0.5, alpha=0.5, outlier.shape = NA,
               position = position_dodge(width = 0.9)) +
  stat_pvalue_manual(pvals_cd8,
                     inherit.aes = FALSE,
                     label = "p = {p.adj}",
                     remove.bracket = TRUE)
p

# Export
set_panel_size(p,
               file = paste0(output_dir, "cxcr6_violin_cd8only.pdf"),
               width = unit(8, "in"), height = unit(6, "in"))

#------------------------------------------------------------------------------
# Visualize expression on only CD4 T cells (SuppFig 5G)

# Violin Plot CXCR6
Idents(s) <- "predicted.celltype.l2"
cd4_subsets <- c("CD4 TEM", "CD4 TCM",  "Treg", "CD4 CTL", "CD4 Naive",
                 "CD4 Proliferating")
s_cd4 <- subset(s, predicted.celltype.l2 %in% cd4_subsets)
s_cd4@meta.data$predicted.celltype.l2 <- factor(s_cd4@meta.data$predicted.celltype.l2,
                                                cd4_subsets)

# Generate p values
pvals <- compare_means(CXCR6 ~ Diagnosis,
                       group.by = "predicted.celltype.l2",
                       method = "wilcox.test",
                       p.adjust.method = "BH",
                       data = data)
pvals_cd4 <- pvals[which(pvals$predicted.celltype.l2 %in% cd4_subsets),]
pvals_cd4[5,] <- list("CD4 Naive", "", "", "", NA, NA, NA, NA, "Wilcoxon")
pvals_cd4[6,] <- list("CD4 Proliferating", "", "", "", NA, NA, NA, NA, "Wilcoxon")
pvals_cd4$y.position <- rep(max(data$CXCR6), length(cd4_subsets))
pvals_cd4$group1 <- cd4_subsets
pvals_cd4$group2 <- cd4_subsets

p <- VlnPlot(subset(s, predicted.celltype.l2 %in% cd4_subsets),
             features = "CXCR6",
             split.by = 'Diagnosis',
             cols = c("gray", "red"),
             pt.size = 0.5,
             log = FALSE) +
  scale_y_continuous(limits = c(NA,3.6)) +
  geom_boxplot(width=0.5, alpha=0.5, outlier.shape = NA,
               position = position_dodge(width = 0.9)) +
  stat_pvalue_manual(pvals_cd4,
                     inherit.aes = FALSE,
                     label = "p = {p.adj}",
                     remove.bracket = TRUE)
p

# Export
set_panel_size(p,
               file = paste0(output_dir, "cxcr6_violin_cd4only.pdf"),
               width = unit(8, "in"), height = unit(6, "in"))