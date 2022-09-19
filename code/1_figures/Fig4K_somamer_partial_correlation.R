# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----     Single cell transcriptomics reveals CD8 T cell recruitment     -----
# -----       to the cerebrospinal fluid during cognitive impairment       -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 08-23-2022
# Written by: Natalie Piehl
# Summary: Generate partial correlations b/w genes oi using Somalogic data
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
protein_dir <- "path/to/somamer/protein/data/"
output_dir <- "path/to/export/results"

# Source helper functions
source("../0_preprocessing/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Format data

# Load in data
data <- read.xlsx(paste0(protein_dir, "ADRC_cleaned_CSF_Wagner.xlsx"), sheetIndex = 1)
cdr <- read.xlsx(paste0(protein_dir, "ADRC_cleaned_CDR.xlsx"), sheetIndex = 1)
moca <- read.xlsx(paste0(protein_dir, "ADRC_cleaned_MOCA.xlsx"), sheetIndex = 1)

# Only complete rows
cdr <- cdr[complete.cases(cdr),]
moca <- moca[complete.cases(moca),]

# Clean up cdr and moca
cdr <- cdr[,c("ADRC_ID", "Age_at_dementia_scoring", "Clinical.Dementia.Rating..CDR..CDRSUM")]
names(cdr) <- c("ADRC_ID", "Age_CDR", "CDR")
moca <- moca[,c("ADRC_ID", "Age_at_MOCA_scoring", "MOCATOTS")]
names(moca) <- c("ADRC_ID", "Age_MOCA", "MoCA")

# Merge onto data
data <- merge(data, cdr, all.x = TRUE, by = "ADRC_ID")
data <- merge(data, moca, all.x = TRUE, by = "ADRC_ID")

# Remove columns not of interest
data <- data[,-c(1:6,10:11)]

# Rename columns
names(data) <- c("Age_GENE", "Sex", "Diagnosis",
                 "CXCL16", "ApoE", "ApoE3", "ApoE4", "ApoE2",
                 "PLTP", "APOC1", "NEFL", "GFAP", "TREM2", "UCHL1",
                 "Age_CDR", "CDR", "Age_MOCA", "MoCA")

# Convert Sex to binary
data$Sex <- plyr::mapvalues(data$Sex,
                         from = c("F", "M"),
                         to = c(0,1))
data$Sex <- as.numeric(data$Sex)

# Fix dumb "MCI " value
data[which(data$Diagnosis == "MCI "), "Diagnosis"] <- "MCI"

# Select markers and genes
markers <- c("NEFL")
genes <- c("CXCL16")

# Define diagnoses
diagnoses <- list("HC", "MCI", "AD", c("MCI", "AD"))
diagnoses_labels <- c("HC", "MCI", "AD", "CI")

#-------------------------------------------------------------------------------
# Run correlation between genes and markers

# Function to run for each marker
corr_marker <- function(marker, gene, data, diagnosis) {
  
  # Select for appropriate age variable
  if (marker == "MoCA") {
    age_var <- "Age_MOCA"
  } else if (marker == "CDR") {
    age_var <- "Age_CDR"
  } else {
    age_var <- "Age_GENE"
  }
  
  # Subset data for diagnosis and marker
  data <- data[which(data$Diagnosis %in% diagnosis),
               c(age_var, "Sex", gene, marker)]
  print(dim(data))
  
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

#-------------------------------------------------------------------------------
# Format for plotting

# Extract R and pval for each comparison
extract_rpval <- function(element) {
  return(list(R = element$estimate[3,4],
              pval = element$p.value[3,4]))
}

# Run extraction
rpval <- list()
for (i in seq(length(genes))) {
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

# Unlist
rpval_shaped$pval <- unlist(rpval_shaped$pval)
rpval_shaped$R <- unlist(rpval_shaped$R)
rpval_shaped$R_sig <- rpval_shaped$R

#-------------------------------------------------------------------------------
# Make corr scatter plot (Fig 4K)

# Make diagnosis in right order
data$Diagnosis <- factor(data$Diagnosis, levels = c("HC", "MCI", "AD"))

# Subset for non na
data <- data[which(!is.na(data$Diagnosis)),]

# Generate plot
p <- ggplot(data, aes(x = NEFL, y = CXCL16,
                      fill = Diagnosis, color = Diagnosis)) +
  theme_Publication_blank() +
  geom_point(size = 2, shape = 21, color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("gray", "lightblue", "red")) +
  scale_color_manual(values = c("gray", "lightblue", "red")) +
  geom_smooth(method='lm', se = FALSE) +
  theme(aspect.ratio = 1,
        text = element_text(size=20))
p

# Export plot
set_panel_size(p, file = paste0(output_dir, "Somalogic_CXCL16_NFL_scatterplot.pdf"),
               width = unit(3.5, "in"), height = unit(3, "in"))

# Write out corr values
write.csv(rpval_shaped, file = paste0(output_dir, "Somalogic_pcor_stats.csv"))