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
# Date: 10-11-2021
# Written by: Natalie Piehl, Emma Tapp
# Summary: Integrate GEX data and generate cell type annotated clusters
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
library_vector <- c("tidyverse", "ggrepel", "ggthemes", "grid",
                    "Seurat", "scales")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = library_vector)

# Initialize paths
seurat_object <- "path/to/seurat_object"
output_dir <- "path/to/export/results"

# Generate output directories
pca_plot_dir <- paste0(output_dir, "/pca_plots/")
umap_plot_dir <- paste0(output_dir, "/umap_plots/")
dir.create(pca_plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(umap_plot_dir, showWarnings = FALSE, recursive = TRUE)

#------------------------------------------------------------------------------
# Integrate on SCT assay

# Load Seurat object
load(seurat_object)

# Split by sample
s_split <- SplitObject(s, split.by = "orig.ident")

# Normalize each sample separately
for (i in names(s_split)) {
  s_split[[i]] <- SCTransform(s_split[[i]],
                              variable.features.n = 1000,
                              ncells = 10000,
                              vars.to.regress = c(
                                "nCount_RNA",
                                "nFeature_RNA",
                                "percent.mt"
                              )
  )
  DefaultAssay(s_split[[i]]) <- "SCT"
}


# Prep for integration
integration_features <- SelectIntegrationFeatures(object.list = s_split, assay = rep("SCT", length(s_split)))
s_split <- PrepSCTIntegration(object.list = s_split, anchor.features = integration_features)
anchors <- FindIntegrationAnchors(object.list = s_split,
                                  normalization.method = "SCT",
                                  anchor.features = integration_features,
                                  reference = c(3,5)) # sample A4(HC) and A6(MCI) as reference for memory saving, too many anchors leads to errors

# Integrate data
s <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Checkpoint save/load
save(s, file = paste0(output_dir, "s_integrated"))

#------------------------------------------------------------------------------
# PCA Visualization (Fig S1C)

# Set assay to integrated
DefaultAssay(object = s) <- "integrated"

# Generate PCA
s <- RunPCA(object = s)

# Initialize identities to look at on PCA
idents <- c("ID", "age", "sort_day", "Diagnosis")

# Plot PCA for each identity
for (ident in idents) {
  s <- SetIdent(s, value = ident)
  myplot <- DimPlot(s, reduction = "pca", pt.size = 2, shuffle = TRUE)
  set_panel_size(myplot,
                 file = paste0(pca_plot_dir, "PCA_", ident, ".pdf"),
                 width = unit(5, "in"), height = unit(4, "in")
  )
}

# Visualize PC Dimensionality of Dataset
myplot <- ElbowPlot(s)
set_panel_size(myplot,
               file = paste0(output_dir, "/Elbowplot.pdf"),
               width = unit(5, "in"), height = unit(4, "in")
)

#------------------------------------------------------------------------------
# Clustering and UMAP visualization

# Cluster Cells on SCT integrated assay
s <- FindNeighbors(s, dims = 1:15)
s <- FindClusters(s, resolution = 0.3)

# Run UMAP on SCT integrated assay
s <- RunUMAP(s, dims = 1:15, n.neighbors = 30)

# Clustering done on integrated assay, but visualize expression with SCT
DefaultAssay(object = s) <- "SCT"

# List cell type markers
markers <- c(
  "CD3D", "CD4", "CD8A", "CD8B", "FOXP3",
  "NKG7", "GNLY", "FCER1A", "CD68",
  "CD14", "FCGR3A",
  "CD19", "JCHAIN", "CD44"
)

# Plot markers on UMAP
for (marker in markers) {
  myplot <- FeaturePlot(s,
                        features = c(marker), pt.size = 1,
                        cols = c("grey88", "deeppink"), reduction = "umap"
  )
  set_panel_size(myplot,
                 file = paste0(umap_plot_dir, marker, ".pdf"),
                 width = unit(5, "in"), height = unit(4, "in")
  )
}

# Add seurat clusters and clonality to idents
idents <- c(idents, "seurat_clusters", "clonal")

# Visualize each identity on UMAP
for (ident in idents) {
  s <- SetIdent(s, value = ident)
  myplot <- DimPlot(s, reduction = "umap", label = TRUE, pt.size = 1, shuffle = TRUE)
  set_panel_size(myplot,
                 file = paste0(umap_plot_dir, ident, ".pdf"),
                 width = unit(5, "in"), height = unit(4, "in")
  )
}

# View diagnosis split
s <- SetIdent(s, value = "Diagnosis")
myplot <- DimPlot(s, reduction = "umap", label = TRUE, pt.size = 1,
                  shuffle = TRUE, split.by = "Diagnosis")
set_panel_size(myplot,
               file = paste0(umap_plot_dir, "Diagnosis_split.pdf"),
               width = unit(5, "in"), height = unit(4, "in")
)

# Plot clonality with both C and NC and NA removed
s <- SetIdent(s, value = "clonal")
clonality <- subset(s, idents = c("C", "NC"))
myplot <- DimPlot(clonality,
                  reduction = "umap", pt.size = 1, shuffle = TRUE,
                  cols = c("aquamarine2", "deeppink")
)
set_panel_size(myplot,
               file = paste0(umap_plot_dir, "clonality_NA_removed.pdf"),
               width = unit(5, "in"), height = unit(4, "in")
)

# Identify marker genes of seurat clusters
s <- SetIdent(s, value = "seurat_clusters")
cluster_markers <- FindAllMarkers(s, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster_markers, paste0(output_dir, "seurat_allcells_cluster_markers.csv"))

# Annotate non-subclustered cell types
s@meta.data$cluster_ident <- "markers"
s@meta.data$cluster_ident[
  which(s@meta.data$seurat_clusters %in% c("10"))
] <- "Plasma Cells"
s@meta.data$cluster_ident[
  which(s@meta.data$seurat_clusters %in% c("7"))
] <- "NK Cells"
s@meta.data$cluster_ident[
  which(s@meta.data$seurat_clusters %in% c("0", "1", "2", "4", "6", "9"))
] <- "T Cells"
s@meta.data$cluster_ident[
  which(s@meta.data$seurat_clusters %in% c("8", "3"))
] <- "Monocytes"
s@meta.data$cluster_ident[
  which(s@meta.data$seurat_clusters %in% c("5"))
] <- "DC"
s@meta.data$cluster_ident[
  which(s@meta.data$seurat_clusters %in% c("11"))
] <- "B Cells"
s@meta.data$cluster_ident[
  which(s@meta.data$seurat_clusters %in% c("12"))
] <- "Undetermined"

# Save seurat clusters for all cells
s@meta.data$seurat_clusters_allcells <- s@meta.data$seurat_clusters
head(s@meta.data)
table(s[["cluster_ident"]])

#------------------------------------------------------------------------------
# Recluster T cells and monocytes separately

# Define celltypes to recluster
celltypes <- c("T Cells", "Monocytes")

for (celltype in celltypes) {
  # Define celltype label
  cell_type_label <- gsub("/", "", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  
  # Initialize output directory
  celltype_dir <- paste0(output_dir, "/", celltype_label, "_reclustering/")
  dir.create(celltype_dir, showWarnings = FALSE)
  
  # Subset celltype
  s_celltype <- subset(s, cluster_ident == celltype)
  print(s_celltype)
  
  #------------------------------------------------------------------------------
  # Visualize PCA on celltype only
  
  # Set assay to integrated
  DefaultAssay(object = s_celltype) <- "integrated"
  
  # Re-calculate variable features
  s_celltype <- FindVariableFeatures(s_celltype, assay = "integrated")
  
  # Generate PCA
  s_celltype <- RunPCA(object = s_celltype)
  
  # Initialize identities to look at on PCA
  idents <- c("ID", "age", "sort_day", "sort_round", "Diagnosis")
  
  # Plot PCA for each identity
  for (ident in idents) {
    s_celltype <- SetIdent(s_celltype, value = ident)
    
    myplot <- DimPlot(s_celltype, reduction = "pca", pt.size = 2, shuffle = TRUE)
    
    set_panel_size(myplot,
                   file = paste0(celltype_dir, "PCA_", ident, ".pdf"),
                   width = unit(5, "in"), height = unit(4, "in")
    )
  }
  
  # Visualize PC Dimensionality of Dataset
  myplot <- ElbowPlot(s_celltype)
  set_panel_size(myplot,
                 file = paste0(celltype_dir, "/Elbowplot.pdf"),
                 width = unit(5, "in"), height = unit(4, "in")
  )
  
  #------------------------------------------------------------------------------
  # Clustering and UMAP visualization on celltype only
  
  # Cluster Cells on SCT integrated assay
  s_celltype <- FindNeighbors(s_celltype, dims = 1:15)
  s_celltype <- FindClusters(s_celltype, resolution = 0.3)
  
  # Run UMAP on SCT integrated assay
  s_celltype <- RunUMAP(s_celltype, dims = 1:15, n.neighbors = 30)
  
  # Clustering done on integrated assay, but visualize expression with SCT
  DefaultAssay(object = s_celltype) <- "SCT"
  
  # List cell type markers
  markers <- c("CD3D", "CD4", "CD8A", "CD8B", "FOXP3",
               "CD14", "FCGR3A", "CD68")
  
  # Plot markers
  for (marker in markers) {
    if (marker %!in% rownames(GetAssayData(s_celltype))) {
      print(paste0(marker, " is not in the seurat object."))
      next
    }
    myplot <- FeaturePlot(s_celltype,
                          features = c(marker), pt.size = 1,
                          cols = c("grey88", "deeppink"), reduction = "umap"
    )
    set_panel_size(myplot,
                   file = paste0(celltype_dir, marker, ".pdf"),
                   width = unit(5, "in"), height = unit(4, "in")
    )
  }
  
  # Add to identities
  idents <- c(idents, "seurat_clusters", "clonal")
  
  # Visualize each identity on UMAP
  for (ident in idents) {
    s_celltype <- SetIdent(s_celltype, value = ident)
    myplot <- DimPlot(s_celltype, reduction = "umap", label = TRUE, pt.size = 1, shuffle = TRUE)
    set_panel_size(myplot,
                   file = paste0(celltype_dir, ident, ".pdf"),
                   width = unit(5, "in"), height = unit(4, "in")
    )
  }
  
  # Set Ident to seurat clusters
  s_celltype <- SetIdent(s_celltype, value = "seurat_clusters")
  
  if (celltype == "T Cells") {
    # Save seurat clusters for celltype
    # Changed cluster nums to +1 because of weird transformation
    s@meta.data$seurat_clusters_tcells <- "NA"
    s@meta.data$seurat_clusters_tcells[
      which(s@meta.data$cluster_ident == celltype)
    ] <- s_celltype@meta.data$seurat_clusters
    print(table(s@meta.data$seurat_clusters_tcells))
    print(table(s_celltype@meta.data$seurat_clusters))
    
    # Assign cell types
    s@meta.data$cluster_ident[
      which(s@meta.data$seurat_clusters_tcells %in% c("3"))
    ] <- "CD8+ T Cells"
    s@meta.data$cluster_ident[
      which(s@meta.data$seurat_clusters_tcells %in% c("1", "2", "5"))
    ] <- "CD4+ T Cells"
    s@meta.data$cluster_ident[
      which(s@meta.data$seurat_clusters_tcells %in% c("4", "7"))
    ] <- "CD4+/CD8+ T Cells"
    s@meta.data$cluster_ident[
      which(s@meta.data$seurat_clusters_tcells %in% c("6"))
    ] <- "T Regulatory Cells"
    
  } else if (celltype == "Monocytes") {
    # Save seurat clusters for celltype
    s@meta.data$seurat_clusters_monocytes <- "NA"
    s@meta.data$seurat_clusters_monocytes[
      which(s@meta.data$cluster_ident == celltype)
    ] <- s_celltype@meta.data$seurat_clusters
    print(table(s@meta.data$seurat_clusters_monocytes))
    print(table(s_celltype@meta.data$seurat_clusters))
    
    # Assign cell types
    # Changed cluster nums to +1 because of weird transformation
    s@meta.data$cluster_ident[
      which(s@meta.data$seurat_clusters_monocytes %in% c("1", "3"))
    ] <- "CD14+/CD16+/CD68hi Monoctyes"
    s@meta.data$cluster_ident[
      which(s@meta.data$seurat_clusters_monocytes %in% c("2"))
    ] <- "CD14+/CD16-/CD68lo Monocytes"
    s@meta.data$cluster_ident[
      which(s@meta.data$seurat_clusters_monocytes %in% c("4", "5", "6"))
    ] <- "CD14+/CD16+/CD68mid Monocytes"
    
  }
  
  # Print cluster markers
  cluster_markers <- FindAllMarkers(s_celltype, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(cluster_markers, paste0(celltype_dir, "seurat_cluster_markers.csv"))
}

# Print table of cells per cell type
print(table(s[["cluster_ident"]]))

# Save final seurat object with celltypes labelled
save(s, file = "data/seurat_object_withcelltype_allsamples")

#------------------------------------------------------------------------------
# Visualize cell types on UMAP (Fig 1C)

# Assign levels to cell type
s <- SetIdent(s, value = "cluster_ident")
levels(s) <- c(
  "CD4+ T Cells", "CD8+ T Cells", "CD4+/CD8+ T Cells", "T Regulatory Cells",
  "NK Cells", "DC", "CD14+/CD16-/CD68lo Monocytes", "CD14+/CD16+/CD68mid Monocytes",
  "CD14+/CD16+/CD68hi Monoctyes", "B Cells", "Plasma Cells", "Undetermined"
)

# Generate umap of cell types
myplot <- DimPlot(s,
                  cols = color_vector, shuffle = TRUE,
                  reduction = "umap", pt.size = 1, label = FALSE
)
set_panel_size(myplot,
               file = paste0(umap_plot_dir, "clusterid.pdf"),
               width = unit(10, "in"), height = unit(8, "in")
)


#------------------------------------------------------------------------------
# Visualize cell types on cell marker heatmap (Fig 1D)

# Sample maximum cells from each celltype
barcode_celltype <- s[["cluster_ident"]]
barcode_celltype$barcode <- rownames(barcode_celltype)
barcode_sampled <- barcode_celltype %>%
  group_by(cluster_ident) %>%
  slice_sample(n = 161, replace = FALSE)
s_sub <- subset(s, cells = barcode_sampled$barcode)

# Assign levels to cell type
s_sub <- SetIdent(s_sub, value = "cluster_ident")
levels(s_sub) <- c(
  "CD4+ T Cells", "CD8+ T Cells", "CD4+/CD8+ T Cells", "T Regulatory Cells",
  "NK Cells", "DC", "CD14+/CD16-/CD68lo Monocytes", "CD14+/CD16+/CD68mid Monocytes",
  "CD14+/CD16+/CD68hi Monoctyes", "B Cells", "Plasma Cells", "Undetermined"
)

# Define cell type colors
color_vector <- c(
  "darkturquoise", "lawngreen", "dodgerblue", "red",
  "lightpink", "hotpink", "gold1", "navy", "violet",
  "darkorange", "darkviolet", "gray"
)

# Visualize heatmap of cell type markers
p <- DoHeatmap(s_sub,
               group.colors = color_vector,
               features = as.vector(c(
                 "CD3D", "CD4", "CD8A", "CD8B", "FOXP3",
                 "NKG7", "GNLY", "FCER1A", "CD68",
                 "CD14", "FCGR3A",
                 "CD19", "JCHAIN", "CD44"
               )),
               slot = "data",
               assay = "RNA") + 
  scale_fill_gradientn(colors = c("white", "darkorchid2", "navy")) +
  theme(plot.margin = unit(c(3,2,1,1), "cm"),
        legend.position = "bottom")
p

ggsave(p, filename = paste0(output_dir, "heatmap_celltype_markers.pdf"),
       width = unit(10, "inch"))

#------------------------------------------------------------------------------
# Look at cell type composition (Fig 1E)

# Generate table of number of cells per type
celltype_composition <- table(s[["cluster_ident"]])

# Generate percentages
celltype_composition_pct <- celltype_composition / sum(celltype_composition_pct)
celltype_composition_pct

#------------------------------------------------------------------------------
# General formatting fixes

# Remove old pANN columns
s@meta.data <- s@meta.data %>% select(-contains("pANN"))

# Add age bin metadata
s@meta.data$age_bin <- rep("middle", nrow(s@meta.data))
s@meta.data$age_bin[
  which(s@meta.data$age %in% seq(70,82,1))] <- "advanced"

#------------------------------------------------------------------------------
# Filter TCRs based on cell type annotations

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
save(s, file = paste0(output_dir, "s_tcrclean"))