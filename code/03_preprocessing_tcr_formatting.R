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
# Summary: Format TCR data to match Seurat and merge samples into one file
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("data.table")
  library("tidyverse")
})

# Organize inputs
tcr_dir <- "/path/to/cellranger_vdj/results"
output_dir <- "/path/to/output/directory"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Merge and format contig matrices

# Define function to grab nth element of a directory
list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}

# Generate vector of "outs" directories
out_dirs <- unlist(lapply(list.dirs.depth.n(tcr_dir, 2),
                          function(x) grep("*/outs", x, value = TRUE)))

# Create list of paths to filtered_contig_annotations.csv files
filtered_contig_paths <- list.files(path = out_dirs,
                                    pattern = "filtered_contig_annotations.csv",
                                    full.names = TRUE, recursive = TRUE)

# Initialize list to contain contig dataframes
filtered_contig_dfs <- list()

# Generate list of dataframes of filtered_contig_annotations.csv files
for (path in filtered_contig_paths) {
  # Read in csv
  tmp <- read.csv(path)
  
  # Isolate sample ID from path
  id <- unlist(strsplit(path, "/")) %>%
    tail(3) %>%
    pluck(1)
  
  # Create ID  column
  tmp$id <- rep(id, dim(tmp)[1])
  
  # Add dataframe to list
  filtered_contig_dfs[[id]] <- tmp
}

# Merge all dataframes into one
filtered_contig_merged <- rbindlist(filtered_contig_dfs)
head(filtered_contig_merged)

# Add sample ID to barcode and raw_clonotype_id
filtered_contig_merged$barcode <- paste0(filtered_contig_merged$id,"_",
                                         filtered_contig_merged$barcode)
filtered_contig_merged$raw_clonotype_id <- paste0(filtered_contig_merged$id,"_",
                                                  filtered_contig_merged$raw_clonotype_id)

# Remove cells with no consensus clonotype
filtered_contig_merged <- filtered_contig_merged[
  which(filtered_contig_merged$raw_consensus_id != 'None'),]

# Write out merged contigs
write.csv(filtered_contig_merged,
          file = paste0(output_dir, "/contigs_merged.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------
# Merge and format clonotype matrices

# Create list of paths to Clonotype<ID>.csv files
clonotype_paths <- list.files(path = out_dirs,
                              pattern = "clonotypes.csv",
                              full.names = TRUE, recursive = TRUE)

# Initialize list to cointain contig dataframes
clonotype_dfs <- list()

# Generate list of dataframes of filtered_contig_annotations.csv files
for (path in clonotype_paths) {
  # Read in csv
  tmp <- read.csv(path)
  
  # Isolate sample name from path
  id <- unlist(strsplit(path, "/")) %>%
    tail(3) %>%
    pluck(1)
  
  # Create id column
  tmp$id <- id
  
  # Add dataframe to list
  clonotype_dfs[[id]] <- tmp
}

# Merge all dataframes into one
clonotype_merged <- rbindlist(clonotype_dfs)
head(clonotype_merged)

# Break apart CDR3 column into separate consensus sequences
clonotype_merged <- separate(data = clonotype_merged, col = cdr3s_aa,
                             into = c("seq1", "seq2", "seq3", "seq4"), sep = ";")
clonotype_merged$tra_cdr3s <- rep(NA, dim(clonotype_merged)[1])
clonotype_merged$trb_cdr3s <- rep(NA, dim(clonotype_merged)[1])

# Remove unneccesary columns
clonotype_merged$cdr3s_nt <- NULL
clonotype_merged$inkt_evidence <- NULL
clonotype_merged$mait_evidence <- NULL

# ITERATE THROUGH EACH ROW
for (row_num in seq(nrow(clonotype_merged))) {
  # Create TRA column containing <tra_consensus_1>;<tra_consensus_2>
  row <- clonotype_merged[row_num,]
  tra_seqs <- grep("^TRA:", row, value = TRUE) %>%
    str_remove_all("TRA:") %>%
    paste(collapse = ";")
  
  # clonotype_merged[ row_num, "tra_cdr3s"] <- tra_seqs
  clonotype_merged$tra_cdr3s[ row_num ] <- tra_seqs
  
  # Create TRB column containing <trb_consensus_1>;<trb_consensus_2>
  trb_seqs <- grep("^TRB:", row, value = TRUE) %>%
    str_remove_all("TRB:") %>%
    paste(collapse = ";")
  
  clonotype_merged$trb_cdr3s[ row_num ] <- trb_seqs
}

# Remove seq columns
clonotype_merged <- subset(clonotype_merged, select=-c(seq1, seq2, seq3, seq4))
head(clonotype_merged)

# Isolate paired clonotypes
paired_clonotypes_merged <- clonotype_merged[
  which(clonotype_merged$trb_cdr3s != "" & clonotype_merged$tra_cdr3s != "" & clonotype_merged$tra_cdr3s != 0 ),]

# Add sample ID to clonotype_id
paired_clonotypes_merged$clonotype_id <- paste0(paired_clonotypes_merged$id,"_",
                                                paired_clonotypes_merged$clonotype_id)

# Write out merged paired clonotypes
write.csv(paired_clonotypes_merged, paste0(output_dir, "/paired_clonotypes_merged.csv"), row.names = FALSE)