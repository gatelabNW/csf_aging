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
# Summary: Remove doublets and perform initial QC on scRNA-TCRseq data
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
  library("DoubletFinder")
})

# Organize inputs
soupx_dir <- "/path/to/soupx_corrected/gex_matrices"
contigs_merged_path <- "/path/to/merged_and_formatted/contigs_matrix"
paired_clonotypes_path <- "/path/to/merged_and_formatted/clonotypes_matrix"
qc_plot_dir <- "/path/to/export/qc_plots/to"
seurat_obj_dir <- "/path/to/export/seurat_object/to"

# Create output directories
dir.create(qc_plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(seurat_obj_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Load in SoupX corrected counts

# List SoupX corrected counts sample dirs
soupx_sample_dirs <- list.dirs(path = soupx_dir, full.names = TRUE, recursive = FALSE)
soupx_sample_dirs

# Create function to load seurat objects
load_seurat <- function(dir) {
  counts <- Read10X(dir)
  project <- tail(unlist(strsplit(dir, "/")), 1)
  return(CreateSeuratObject(counts = counts, project = project,
                            min.cells = 3, min.features = 200))
}

# Apply Seurat loading function to all samples
seurat_object_list <- sapply(soupx_sample_dirs, load_seurat)

#-------------------------------------------------------------------------------
# Identify doublets with DoubletFinder

# Define doublet formation rate
doublet_formation_rate <- 0.010
print(paste0("Using doublet formation rate of ", doublet_formation_rate))

run_doubletfinder <- function(s) {
  # Process normally
  s <- s %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
  
  # Run TSNE clustering
  s <- RunPCA(s, features = VariableFeatures(object = s))
  s <- FindNeighbors(s, dims = 1:12)
  s <- FindClusters(s, resolution = 0.3)
  s <- RunTSNE(s, dims = 1:12)
  
  # pK Identification (no ground-truth) 
  sweep.res.list <- paramSweep_v3(s, PCs = 1:12, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  max_index <- which.max(bcmvn$BCmetric)
  optimal_pK <- as.numeric(as.character(bcmvn[max_index, "pK"]))
  
  # Homotypic Doublet Proportion Estimate
  annotations <- s@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(doublet_formation_rate*nrow(s@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  print(s)
  
  # Run DoubletFinder
  s <- doubletFinder_v3(s, PCs = 1:12, pN = 0.25, pK = optimal_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  # Rename column name for consistency
  colnames(s@meta.data)[ grep("DF.classifications*", colnames(s@meta.data)) ] <- "DF.classifications"
  print(head(s@meta.data))
  
  return(s)
}

# Apply DoubletFinder function to all samples
seurat_object_list <- sapply(seurat_object_list, run_doubletfinder)

#-------------------------------------------------------------------------------
# Merge samples togther
s <- merge(seurat_object_list[[1]], unlist(seurat_object_list, use.names = FALSE)[2:length(seurat_object_list)],
           add.cell.ids = c('A1', 'A2', 'A4', 'A5', 'A6', 'A7', 'A8',
                            'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8',
                            'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C8',
                            'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D8',
                            'E1', 'E2', 'E3', 'E4', 'E5', 'E7', 'E8',
                            'F1', 'F2', 'F3', 'F4', 'F6', 'F7', 'F8',
                            'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8',
                            'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8'),
           project = "CSF_aging")
s

#-------------------------------------------------------------------------------
# Add TCR Data

# Add cell ID to metadata
ID <- sapply(strsplit(rownames(s@meta.data), split="_"), "[[", 1)
s <- AddMetaData(object=s,
                 metadata=data.frame(ID=ID,
                                     row.names=rownames(s@meta.data)))

# Load in TCR information
contigs_merged <- read.csv(contigs_merged_path)
TCRs_paired <- read.csv(paired_clonotypes_path)

# Remove duplicate rows per cell
tcr <- contigs_merged[!duplicated(contigs_merged$barcode), ]

# Only keep the barcode and clonotype columns
# We'll get additional clonotype info from the clonotype table.
tcr <- tcr[,c("barcode", "raw_clonotype_id")]
names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"

# Slap the AA sequences onto our original table by clonotype_id
tcr <- merge(tcr, TCRs_paired[, c("clonotype_id", "trb_cdr3s", "tra_cdr3s", "frequency")])

# Reorder so barcodes are first column and set them as rownames.
tcr <- tcr[, c(2,1,3,4, 5)]
rownames(tcr) <- tcr[,1]
tcr[,1] <- NULL
print("Head of TCR")
head(tcr)

# Add to the Seurat object's metadata.
s <- AddMetaData(object=s,
                 metadata=tcr)

#-------------------------------------------------------------------------------
# Add Metadata

# Define Diagnosis Groups
AD_samples <- c("E1", "B4", "B6", "G6", "B7", "C5")
MCI_samples <- c("A6", "B1", "H1", "E4", "E7", "E5", "D2", "D3")

# Add diagnosis to metadata
s@meta.data$Diagnosis<- rep("HC", nrow(s@meta.data))
s@meta.data$Diagnosis[
  which(s@meta.data$ID %in% c(AD_samples, MCI_samples))] <- "CI"

# Check out number of cells per diagnosis
table(s@meta.data$Diagnosis)

# Calculate Percent Mitochondrial Genes
s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")

# Add clonal info to metadata
s@meta.data$clonal <- "NA"
s@meta.data$clonal[
  s@meta.data$frequency > 1] <- "C"
s@meta.data$clonal[
  s@meta.data$frequency == 1] <- "NC"
table(s@meta.data$clonal)

# Define Sex Groups
f <- c('A1', 'A7', 'A8', 'B2', 'B6', 'B7', 'B8',
       'C2', 'C3', 'C8', 'D1', 'D5', 'D6', 'D8',
       'E1', 'E2', 'E3', 'E4', 'E5', 'E7', 'F1', 'F2', 'F4', 'F8',
       'G2', 'G3', 'G4', 'G6', 'G7', 'G8', 'H2', 'H3', 'H6', 'H7')
m <- c('A2', 'A4', 'A5', 'A6', 'B1', 'B3', 'B4', 'B5',
       'C1', 'C4', 'C5', 'C6', 'D2', 'D3', 'D4',
       'E8', 'F3', 'F6', 'F7',
       'G1', 'G5', 'H1', 'H4', 'H5', 'H8')

# Add sex to metadata
s@meta.data$sex <- "NA"
s@meta.data$sex[
  which(s@meta.data$ID %in% m)] <- "m"
s@meta.data$sex[
  which(s@meta.data$ID %in% f)] <- "f"
table(s@meta.data$sex)

# Define Age Groups
x48 <- c('D3')
x54 <- c('H7', 'C4')
x55 <- c('H2')
x56 <- c('E5')
x62 <- c('G3', 'H6')
x63 <- c('A7', 'C3', 'E1', 'G1')
x64 <- c('B7', 'D8', 'G2', 'G6')
x65 <- c('D1', 'E3')
x66 <- c('C2', 'E2', 'G7', 'G5')
x67 <- c('A1', 'B2', 'F8', 'A4', 'F7')
x69 <- c('F4', 'H3', 'H8')
x70 <- c('D5', 'B1', 'B5', 'C1')
x71 <- c('A8', 'C8', 'B3')
x72 <- c('F1', 'G4', 'A5', 'D2')
x73 <- c('F2', 'G8', 'B4')
x74 <- c('E4', 'E8')
x75 <- c('H5')
x76 <- c('B6', 'B8', 'D6')
x77 <- c('A2', 'C5')
x78 <- c('A6', 'F6')
x79 <- c('E7', 'D4', 'H1')
x80 <- c('H4')
x81 <- c('F3')
x82 <- c('C6')

# Add age to metadata
s@meta.data$age <- "NA"
s@meta.data$age[
  which(s@meta.data$ID %in% x48)] <- 48
s@meta.data$age[
  which(s@meta.data$ID %in% x54)] <- 54
s@meta.data$age[
  which(s@meta.data$ID %in% x55)] <- 55
s@meta.data$age[
  which(s@meta.data$ID %in% x56)] <- 56
s@meta.data$age[
  which(s@meta.data$ID %in% x57)] <- 57
s@meta.data$age[
  which(s@meta.data$ID %in% x58)] <- 58
s@meta.data$age[
  which(s@meta.data$ID %in% x62)] <- 62
s@meta.data$age[
  which(s@meta.data$ID %in% x63)] <- 63
s@meta.data$age[
  which(s@meta.data$ID %in% x64)] <- 64
s@meta.data$age[
  which(s@meta.data$ID %in% x65)] <- 65
s@meta.data$age[
  which(s@meta.data$ID %in% x66)] <- 66
s@meta.data$age[
  which(s@meta.data$ID %in% x67)] <- 67
s@meta.data$age[
  which(s@meta.data$ID %in% x69)] <- 69
s@meta.data$age[
  which(s@meta.data$ID %in% x70)] <- 70
s@meta.data$age[
  which(s@meta.data$ID %in% x71)] <- 71
s@meta.data$age[
  which(s@meta.data$ID %in% x72)] <- 72
s@meta.data$age[
  which(s@meta.data$ID %in% x73)] <- 73
s@meta.data$age[
  which(s@meta.data$ID %in% x74)] <- 74
s@meta.data$age[
  which(s@meta.data$ID %in% x75)] <- 75
s@meta.data$age[
  which(s@meta.data$ID %in% x76)] <- 76
s@meta.data$age[
  which(s@meta.data$ID %in% x77)] <- 77
s@meta.data$age[
  which(s@meta.data$ID %in% x78)] <- 78
s@meta.data$age[
  which(s@meta.data$ID %in% x79)] <- 79
s@meta.data$age[
  which(s@meta.data$ID %in% x80)] <- 80
s@meta.data$age[
  which(s@meta.data$ID %in% x81)] <- 81
s@meta.data$age[
  which(s@meta.data$ID %in% x82)] <- 82
table(s@meta.data$age)

# Define sort day groups
sort_day1 <- c('A1', 'A2', 'A4',
               'B1', 'B2', 'B3', 'B4',
               'C1', 'C2', 'C3', 'C4',
               'D1', 'D2', 'D3', 'D4',
               'E1', 'E2', 'E3', 'E4',
               'F1', 'F2', 'F3', 'F4',
               'G1', 'G2', 'G3', 'G4',
               'H1', 'H2', 'H3', 'H4')
sort_day2 <- c('A5', 'A6', 'A7', 'A8',
               'B5', 'B6', 'B7', 'B8',
               'C5', 'C6', 'C8',
               'D5', 'D6', 'D8',
               'E5', 'E7', 'E8',
               'F6', 'F7', 'F8',
               'G5', 'G6', 'G7', 'G8',
               'H5', 'H6', 'H7', 'H8')

# Add sort day to metadata
s@meta.data$sort_day <- "NA"
s@meta.data$sort_day[
  which(s@meta.data$ID %in% sort_day1)] <- "sort_day1"
s@meta.data$sort_day[
  which(s@meta.data$ID %in% sort_day2)] <- "sort_day2"
table(s@meta.data$sort_day)

#-------------------------------------------------------------------------------
# Visualize initial QC

# nFeatures per ID
Idents(s) <- 'ID'
myplot <- VlnPlot(s, features = c("nFeature_RNA"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8)) +
  labs(x = "Sample ID",
       y = "Number of Features",
       title = "Number of Features per Sample ID")
set_panel_size(myplot, file=paste0(qc_plot_dir, "nFeature_ID.pdf"),
               width=unit(8, "in"), height=unit(4, "in"))

# nCounts per ID
myplot <- VlnPlot(s, features = c("nCount_RNA"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8)) +
  labs(x = "Sample ID",
       y = "Number of Counts",
       title = "Number of Counts per Sample ID")
set_panel_size(myplot, file=paste0(qc_plot_dir, "nCount_ID.pdf"),
               width=unit(8, "in"), height=unit(4, "in"))

# percent.mt per ID
myplot <- VlnPlot(s, features = c("percent.mt"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8)) +
  labs(x = "Sample ID",
       y = "Percent Mitochondrial",
       title = "Percent Mitochondrial per Sample ID")
set_panel_size(myplot, file=paste0(qc_plot_dir, "percent.mt_ID.pdf"),
               width=unit(8, "in"), height=unit(4, "in"))

# Add levels to age so it plots in order
s@meta.data$age <- factor(s@meta.data$age, levels = c(48, 54, 55, 56, 57, 58,
                                                      62, 63, 64, 65, 66, 67, 69,
                                                      70, 71, 72, 73, 74, 75, 76,
                                                      77, 78, 79, 80, 81, 82))
# nFeatures per age
Idents(s) <- 'age'
myplot <- VlnPlot(s, features = c("nFeature_RNA"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  labs(x = "Age",
       y = "Number of Features (Genes)",
       title = "nFeatures per Age")
set_panel_size(myplot, file=paste0(qc_plot_dir, "nFeature_Age.pdf"),
               width=unit(8, "in"), height=unit(4, "in"))

# nCounts per age
myplot <- VlnPlot(s, features = c("nCount_RNA"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  labs(x = "Age",
       y = "Number of Counts",
       title = "nCounts per Age")
set_panel_size(myplot, file=paste0(qc_plot_dir, "nCount_Age.pdf"),
               width=unit(8, "in"), height=unit(4, "in"))

# percent.mt per age
myplot <- VlnPlot(s, features = c("percent.mt"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  labs(x = "Age",
       y = "Percent Mitochondrial",
       title = "Percent MT per Age")
set_panel_size(myplot, file=paste0(qc_plot_dir, "percent.mt_Age.pdf"),
               width=unit(8, "in"), height=unit(4, "in"))

# nFeatures per Diagnosis
Idents(s) <- 'Diagnosis'
myplot <- VlnPlot(s, features = c("nFeature_RNA"), pt.size = 0, cols = c("black", "red")) +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  labs(x = "Diagnosis",
       y = "Number of Features",
       title = "Number of Features per Diganosis")
print(myplot$data)
set_panel_size(myplot, file=paste0(qc_plot_dir, "nFeature_Diagnosis.pdf"),
               width=unit(3, "in"), height=unit(4, "in"))

# nCounts per Diagnosis
myplot <- VlnPlot(s, features = c("nCount_RNA"), pt.size = 0, cols = c("black", "red")) +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  labs(x = "Diagnosis",
       y = "Counts",
       title = "Number of Counts per Diganosis")
set_panel_size(myplot, file=paste0(qc_plot_dir, "nCount_Diagnosis.pdf"),
               width=unit(3, "in"), height=unit(4, "in"))

# percent.mt per Diagnosis
myplot <- VlnPlot(s, features = c("percent.mt"), pt.size = 0, cols = c("black", "red")) +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  labs(x = "Diagnosis",
       y = "Percent Mitochondrial",
       title = "Percent Mitochondrial per Diagnosis")
set_panel_size(myplot, file=paste0(qc_plot_dir, "percent.mt_Diagnosis.pdf"),
               width=unit(3, "in"), height=unit(4, "in"))

# nFeatures per sort day
Idents(s) <- 'sort_day'
myplot <- VlnPlot(s, features = c("nFeature_RNA"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  labs(x = "Sort Day",
       y = "Number of Features (Genes)",
       title = "nFeatures per Sort Day")
set_panel_size(myplot, file=paste0(qc_plot_dir, "nFeature_sort_day.pdf"),
               width=unit(8, "in"), height=unit(4, "in"))

# nCounts per sort day
myplot <- VlnPlot(s, features = c("nCount_RNA"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  labs(x = "Sort Day",
       y = "Number of Counts",
       title = "nCounts per Sort Day")
set_panel_size(myplot, file=paste0(qc_plot_dir, "nCount_sort_day.pdf"),
               width=unit(8, "in"), height=unit(4, "in"))

# percent.mt per sort day
myplot <- VlnPlot(s, features = c("percent.mt"), pt.size = 0) +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  labs(x = "Sort Day",
       y = "Percent Mitochondrial",
       title = "Percent MT per Sort Day")
set_panel_size(myplot, file=paste0(qc_plot_dir, "percent.mt_sort_day.pdf"),
               width=unit(8, "in"), height=unit(4, "in"))

# percent.mt vs nCount per ID
Idents(s) <- 'ID'
myplot <- FeatureScatter(s, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 10, linetype = 2, color = "gray") +
  labs(x = "Number of Counts",
       y = "Percent Mitochondrial",
       title = "Percent MT by nCounts")
set_panel_size(myplot, file=paste0(qc_plot_dir, "percentmt_by_nCount.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))

# nFeature vs nCount per ID
myplot <- FeatureScatter(s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  theme_Publication_blank() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 20000, linetype = 2, color = "gray") +
  geom_hline(yintercept = 4000, linetype = 2, color = "gray") +
  labs(x = "Number of Counts",
       y = "Number of Features (Genes)",
       title = "nFeatures per nCounts")
set_panel_size(myplot, file=paste0(qc_plot_dir, "nFeature_by_nCount.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))

# Generate PCA
s <- SCTransform(s,
                 variable.features.n = 1000,
                 vars.to.regress = c("nCount_RNA",
                                     "nFeature_RNA"))
s <- RunPCA(object = s)

# PCA labelled by sort day
s <- SetIdent(s, value = 'sort_day')
myplot <- DimPlot(s, reduction = "pca", pt.size = 0)
set_panel_size(myplot, file=paste0(qc_plot_dir, "PCA_sort_day.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))

# PCA labelled by diagnosis
s <- SetIdent(s, value = 'Diagnosis')
myplot <- DimPlot(s, reduction = "pca", pt.size = 0, cols = c("black", "red"))
set_panel_size(myplot, file=paste0(qc_plot_dir, "PCA_diagnosis.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))

#-------------------------------------------------------------------------------
# Filter and export results

# Subset for singlets
s <- subset(s, subset = DF.classifications == "Singlet")

# Remove cells with above threshold mitochondrial reads
s <- subset(s, subset = percent.mt < 10)

# Export results
save(s, file = paste0(seurat_obj_dir, "/seurat_object"))