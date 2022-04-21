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
# Date: 01-18-2022
# Written by: Natalie Piehl
# Summary: Generate bubbleplot of enrichment analysis
#
# - Enrichment was performed with https://metascape.org
# - Gene lists from clusters of interests were copied from "hclust_cut.csv"
#   output of LOESS analysis
# - Enrichment was run using GO MF, GO BP, GO CC, Reactome, and KEGG databases
# - All other default Metascape parameters were used
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
  library("xlsx")
  library("DEswan")
  library("scales")
})

# Initialize paths
input_dir <- "path/to/metascape/results"
output_dir <- "path/to/export/results"

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#------------------------------------------------------------------------------
# Generate bubble plot of enrichment (Fig 3E)

# Identify directories containing metadata results
dirs <- list.dirs(input_dir, full.names = FALSE, recursive = FALSE)

for (dir in dirs) {
  # Read in xlsx
  enrichment <- read.xlsx(file.path(input_dir, dir, "metascape_result.xlsx"),
                          sheetIndex = 2)
  enrichment <- enrichment[grep("Summary", enrichment$GroupID),]
  
  # Sort and subset ORA
  enrichment <- enrichment[order(enrichment$LogP),]
  enrichment <- enrichment[1:5,]
  enrichment$term_description <- paste0(enrichment$Term, ":\n", enrichment$Description)
  
  # Format GeneRatio to be numeric
  enrichment["InTerm_InList"] <- sapply(sapply(enrichment["InTerm_InList"],
                                               function(y) strsplit(y, split = "/")),
                                        function(x) as.numeric(x[1])/as.numeric(x[2]))
  enrichment <- enrichment[order(-enrichment$LogP),]
  enrichment$term_description <- str_wrap(enrichment$term_description,
                                          width = 25)
  enrichment$term_description <- factor(enrichment$term_description,
                                        levels = enrichment$term_description)
  
  # Generate dotplot
  p <- ggplot(data = enrichment, aes(x = -LogP, y = term_description,
                                     color = -LogP, size = InTerm_InList)) +
    theme_Publication_blank() +
    theme(aspect.ratio = 2,
          panel.grid.major.y = element_line(color = "grey80", size = 0.3),
          legend.direction = "vertical",
          legend.box = "vertical",
          legend.position = "right",
          text = element_text(size=20)) +
    xlim(-max(enrichment$LogP) - 2, NA) +
    scale_size_continuous(range = c(4, 12)) +
    scale_color_gradient(low = "blue", high = "red") +
    geom_point()
  
  # Export plot
  set_panel_size(p, file=paste0(output_dir, "ORA_", dir, ".pdf"),
                 width=unit(3, "in"), height=unit(4, "in"))
}