# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----         Activated monocytes recruit CD8 T cells to the             -----
# -----         cerebrospinal fluid during cognitive impairment            -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 04-20-2022
# Written by: Emma Tapp, Natalie Piehl
# Summary: Define helper functions used across the project
#
#-------------------------------------------------------------------------------
# Define functions

# Load in libraries
suppressMessages({
  library("ggrepel")
  library("tidyverse")
  library("ggpubr")
  library("ggthemes")
  library("grid")
})

# Negation of %in%
'%!in%' <- Negate('%in%')

# Define publication theme plotting function
theme_Publication_blank <- function(base_size=12, base_family="") {
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(fill = "transparent",colour = NA),
           plot.background = element_rect(fill = "transparent",colour = NA),
           panel.border = element_rect(colour = NA, fill = "transparent"),
           axis.title = element_text(size = rel(1)),
           axis.title.y = element_text(angle=90,margin=margin(0,10,0,0)),
           axis.title.x = element_text(margin=margin(10,0,0,0)),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(size = 0.3),
           axis.line.x = element_line(size = 0.3, linetype = "solid", colour = "black"),
           axis.line.y = element_line(size = 0.3, linetype = "solid", colour = "black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA, fill="transparent"),
           legend.position = "bottom",
           legend.margin = margin(t = 10, unit='pt'),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#d8d8d8",fill="#d8d8d8")
   ))
} 

# Define plot exporting helper function
set_panel_size <- function(p=NULL, g=ggplotGrob(p), file=NULL, 
                           margin = unit(1,"mm"),
                           width=unit(7, "inch"), 
                           height=unit(5, "inch")){
  panels <- grep("panel", g$layout$name)
  panel_index_w<- unique(g$layout$l[panels])
  panel_index_h<- unique(g$layout$t[panels])
  nw <- length(panel_index_w)
  nh <- length(panel_index_h)
  
  g$widths[panel_index_w] <-  rep(width,  nw)
  g$heights[panel_index_h] <- rep(height, nh)
  
  if(!is.null(file)) {
    ggplot2::ggsave(file, g,
                    width = convertWidth(sum(g$widths) + margin,
                                         unitTo = "in", valueOnly = TRUE),
                    height = convertHeight(sum(g$heights) + margin,
                                           unitTo = "in", valueOnly = TRUE), useDingbats=F,
                    dpi=300)
    invisible(g)
  }
}

# Define volcano plotting function
volcano_plot <- function(data, file = NULL, title = NULL,
                         padj.lim = NULL, lfc.lim = NULL, lfc.asymmetric = NULL,
                         padj.thresh = 0.01, lfc.thresh = 0.25, x_title = "avg_log2FC",
                         width=unit(4, "inch"), height=unit(4, "inch")) {
  # Create gene name column
  data$gene <- rownames(data)
  
  # Generate PFC scores
  data$PFC <- -log10(data$BH) * abs(data$avg_log2FC)
  PFC <- unique(data$PFC)
  data$PFC[data$PFC == Inf] <- sort(PFC, partial=length(PFC)-1)[length(PFC)-1]
  
  #Generate log padj column
  data$log_padj <- -log10(data$BH)
  
  # Define limits if not provided
  if (is.null(padj.lim)) {
    log.padj.lim <- unique(data$log_padj)[order(-unique(data$log_padj))][2]
  } else {
    log.padj.lim <- -log10(padj.lim)
  }
  if (is.null(lfc.lim)) {
    lfc.lim <- abs(data[order(-abs(data$avg_log2FC)),"avg_log2FC"][1])
  }
  
  # Generate color column
  data$color <- rep("black", nrow(data))
  data[which(data$BH <= padj.thresh &
               data$avg_log2FC > lfc.thresh), 'color'] <- "red"
  data[which(data$BH <= padj.thresh &
               data$avg_log2FC < -lfc.thresh), 'color'] <- "blue"
  
  # Scale down genes outside of bounds
  data[which(data$log_padj > log.padj.lim), 'log_padj'] <- log.padj.lim
  data[which(data$avg_log2FC > lfc.lim), 'avg_log2FC'] <- lfc.lim
  data[which(data$avg_log2FC < -lfc.lim), 'avg_log2FC'] <- -lfc.lim
  
  # Generate asymmetric lfc limits if necessary
  if (is.null(lfc.asymmetric)) {
    lfc_lims <- c(-lfc.lim, lfc.lim)
  } else {
    lfc_lims <- lfc.asymmetric
  }
  
  # Plot data
  p <-
    ggplot(data,
           aes(
             x = avg_log2FC,
             y = log_padj,
             color = color,
             label = gene,
             size = PFC
           )) +
    theme_Publication_blank() +
    geom_hline(
      yintercept = -log10(padj.thresh),
      linetype = 2,
      color = "gray"
    ) +
    geom_vline(xintercept = lfc.thresh,
               linetype = 2,
               color = "gray") +
    geom_vline(
      xintercept = -lfc.thresh,
      linetype = 2,
      color = "gray"
    ) +
    geom_point(aes(size = PFC), alpha = 0.5) +
    scale_color_manual(values = c("red" = "red",
                                  "black" = "black",
                                  "blue" = "blue")) +
    geom_text_repel(
      data = data[which(data$color != 'black'), ],
      inherit.aes = T,
      color = 'black',
      size = 4,
      force = 3
    ) +
    theme(legend.position = "none") +
    labs(title = title,
         x = x_title) +
    scale_x_continuous(limits = lfc_lims) +
    scale_y_continuous(limits = c(0, log.padj.lim)) 
  theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12))
  
  # Export plot
  if (is.null(file)) {
    return(p)
  } else {
    set_panel_size(
      p,
      file = file,
      width = width,
      height = height
    )
    return(NULL)
  }
}
