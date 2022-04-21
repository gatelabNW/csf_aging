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
# Date: 03-21-2022
# Written by: Natalie Piehl
# Summary: Identify number of >0.9 Lsim connections between HC and CI samples
#
# - Fig 4C visual created with Prism 9.2.0
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("ggrepel")
  library("ggthemes")
  library("RecordLinkage")
  library("Seurat")
})

# Initialize paths
seurat_object <- "path/to/seurat_object/"
output_dir <- "path/to/export/results"

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#------------------------------------------------------------------------------
# Data preparation

# Load in 
load(seurat_object)

# Extract clonotype data
df <- s[[c("cluster_ident", "clonotype_id", "Diagnosis", "ID", "age",
           "frequency_filtered", "tra_cdr3s", "trb_cdr3s")]] %>%
  dplyr::distinct(clonotype_id, .keep_all = TRUE)

# Filter unpaired TCRs
df <- df[ which(df$tra_cdr3s!="" & df$trb_cdr3s!=""), ]

# Filter ambiguous TCRs
df <- df[ -grep(";",df$tra_cdr3s), ]
df <- df[ -grep(";",df$trb_cdr3s), ]

# Filter nonclonal TCRs
df <- df[ which(df$frequency_filtered > 1), ]

# create clonotype_id.long as unique id for each clonotype
df$clonotype_id_long <- paste(df$tra_cdr3s, df$trb_cdr3s, sep="_")

# Add bins
age_bins <- c("54_64", "65_67", "68_72", "73_82")
samples_54_64 <- unique((subset(df, age >= 54 & age <= 64))$ID)
samples_65_67 <- unique((subset(df, age >= 65 & age <= 67))$ID)
samples_68_72 <- unique((subset(df, age >= 68 & age <= 72))$ID)
samples_73_82 <- unique((subset(df, age >= 73 & age <= 82))$ID)
df <- mutate(df, bin =
               case_when(ID %in% samples_54_64 ~ age_bins[1],
                         ID %in% samples_65_67 ~ age_bins[2],
                         ID %in% samples_68_72 ~ age_bins[3],
                         ID %in% samples_73_82 ~ age_bins[4]))
# Add diagnostic values
df$bin[ which(df$Diagnosis == "CI") ] <- "CI"
table(df$bin)

#------------------------------------------------------------------------------
# Calculate number of CI and HC connections (Fig 4C)

# Define function to find number of connections
find_num_connections <- function(age_bin, cell_type = NULL) {
  # Subset disease and age bins
  df_age <- df[ which(df$bin == age_bin), ]
  df_disease <- df[ which(df$bin == "CI"), ]
  
  # Subset for celltype (if applicable)
  if (!is.null(cell_type)) {
    df_age <- df_age[ which(df_age$cluster_ident == cell_type), ]
    df_disease <- df_disease[ which(df_disease$cluster_ident == cell_type), ]
  }
  
  # Calculate Lsim
  lsim_df <- sapply(df_age$clonotype_id_long, function(seq_1) {
    sapply(df_disease$clonotype_id_long, function(seq_2) {
      levenshteinSim(seq_1, seq_2)
      })
    })
  connections <- which(lsim_df >= 0.9, arr.ind = T)
  
  # Generate lists of pairs meeting threshold
  disease_tcrs <- df_disease[connections[,"row"],"clonotype_id_long"]
  age_tcrs <- df_age[connections[,"col"],"clonotype_id_long"]
  lsim <- lsim_df[ which(lsim_df >= 0.9) ]
  
  return(list(disease_tcrs, age_tcrs, lsim, rep(age_bin, length(lsim))))
}

# Generate results
res_54_64 <- find_num_connections("54_64") %>% data.frame
res_65_67 <- find_num_connections("65_67") %>% data.frame
res_68_72 <- find_num_connections("68_72") %>% data.frame
res_73_82 <- find_num_connections("73_82") %>% data.frame

# Merge results together
connections <- data.frame("bin" = age_bins,
                          "connections" = c(nrow(res_54_64), nrow(res_65_67),
                                            nrow(res_68_72), nrow(res_73_82)))
connections <- dplyr::mutate(connections,
                             norm_connections = connections / table(df$bin[which(df$bin != "CI")]))

# Export csv
write.csv(connections, paste0(output_dir, "ci_agebin_connections.csv"))