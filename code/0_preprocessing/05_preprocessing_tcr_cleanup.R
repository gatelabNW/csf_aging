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
# Date: 02-10-2022
# Written by: Natalie Piehl
# Summary: Remove TCR data for non T-cell annotated cells and update frequency
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
})

# Initialize paths
seurat_object <- "path/to/seurat_object/"
output_dir <- "path/to/export/results"

# Source helper functions
source("00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# General formatting fixes

# Load in seurat object
load(seurat_object)

# Remove old pANN columns
s@meta.data <- s@meta.data %>% select(-contains("pANN"))

# Add age bin metadata
s@meta.data$age_bin <- rep("middle", nrow(s@meta.data))
s@meta.data$age_bin[
  which(s@meta.data$age %in% seq(70,82,1))] <- "advanced"

#------------------------------------------------------------------------------
# Filter TCRs

# Remove TCR annotations from Non T cells
t_cell_clusters <- c("CD4+ T Cells", "CD4+/CD8+ T Cells", "CD8+ T Cells", "T Regulatory Cells")

print(table(s[[c("cluster_ident", "clonal")]]))

s@meta.data$clonal[
  which(s@meta.data$cluster_ident %!in% t_cell_clusters)] <- NA
s@meta.data$clonotype_id[
  which(s@meta.data$cluster_ident %!in% t_cell_clusters)] <- NA
s@meta.data$trb_cdr3s[
  which(s@meta.data$cluster_ident %!in% t_cell_clusters)] <- NA
s@meta.data$tra_cdr3s[
  which(s@meta.data$cluster_ident %!in% t_cell_clusters)] <- NA
s@meta.data$frequency[
  which(s@meta.data$cluster_ident %!in% t_cell_clusters)] <- NA
s@meta.data$disease.clonal[
  which(s@meta.data$cluster_ident %!in% t_cell_clusters)] <- NA
s@meta.data$Genderclonal[
  which(s@meta.data$cluster_ident %!in% t_cell_clusters)] <- NA

print(table(s[[c("cluster_ident", "clonal")]]))

# Calculate filtered frequency
filtered_freq <- data.frame(table(s[["clonotype_id"]]))
names(filtered_freq) <- c("clonotype_id", "frequency_filtered")
filtered_freq_merged <- plyr::join(s[["clonotype_id"]], filtered_freq)

# If in same order, merge
if (all.equal(filtered_freq_merged$clonotype_id, s[["clonotype_id"]]$clonotype_id)) {
  s@meta.data$frequency_filtered <- filtered_freq_merged$frequency_filtered
} else {
  stop("Unable to merge.")
}

# Reassign NC
s@meta.data$clonal[
  which(s@meta.data$frequency_filtered == 1)] <- "NC"
print("After NC reassignmnet")
print(table(s[[c("cluster_ident", "clonal")]]))

#------------------------------------------------------------------------------
# Generate normalized frequency

# Calculate TCR counts per sample
tcr_count <- s[[c("orig.ident", "clonotype_id", "age", "frequency_filtered")]] %>%
  dplyr::filter(!is.na(frequency_filtered)) %>%
  dplyr::distinct(clonotype_id, .keep_all = TRUE) %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarize(tcr_count = sum(frequency_filtered), across(age)) %>%
  dplyr::distinct()
head(tcr_count)

# Calculate median TCR counts
median_tcr_count <- median(tcr_count$tcr_count)
median_tcr_count

# Initialize column
s@meta.data$normalized_frequency <- rep(NA, nrow(s@meta.data))

# Normalize frequency for each sample
for (sample in tcr_count$orig.ident) {
  sample_tcr_count <- tcr_count %>% dplyr::filter(orig.ident == sample) %>% dplyr::pull(tcr_count)
  s@meta.data$normalized_frequency[
    which(s@meta.data$orig.ident == sample) ] <- s@meta.data$frequency_filtered[
      which(s@meta.data$orig.ident == sample) ] * (median_tcr_count / sample_tcr_count)
}

# Check output
head(s[[c("frequency", "frequency_filtered", "normalized_frequency", "age_bin")]])

# Save modified seurat object
save(s, file = paste0(output_dir, "s_tcrclean_new"))
