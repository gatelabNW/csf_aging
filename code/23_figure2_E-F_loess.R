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
# Date: 02-11-2022
# Written by: Natalie Piehl
# Summary: Generate LOESS trajectories
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
library_vector <- c("ggplot2", "ggpubr", "ggrepel", "ggthemes", "Seurat",
                    "tidyverse", "plyr", "dplyr", "factoextra", "scales",
                    "patchwork", "clusterProfiler", "pheatmap", "colorspace")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = library_vector)

# Initialize paths
data_path <- "path/to/pseudobulk/data"
output_parent_dir <- "path/to/export/results"

# Generate output directory
dir.create(output_parent_dir, showWarnings = FALSE, recursive = TRUE)

# Organize parameters
# clonal: "C", "NC", "compare", "all"
clonal <- "all"
# Diagnosis: "HC", "MCI/AD", "compare", "all"
Diagnosis <- "compare"
# sex: "f", "m", "compare", "all"
sex <- "all"
# cell_type: "all", celltype, "allandcelltype"
cell_type <- "allandcelltype"

# Specific reference clustering folder if doing a comparison
ref_clustering <- "HC"

# Define number of clusters to generate
cluster_nums = 6:12

# desired age range
age_min <- 62
age_max <- 82

# Define non-gene column names
non_gene_cols <- c("ID", "age", "Diagnosis", "sex", "age_bin", "clonal", "cell_type")

#------------------------------------------------------------------------------
# Define functions

subset_parameters <- function(data, param) {
  # Select for param settings
  group <- deparse(substitute(param))
  if (group == "clonal") {
    if (param == "compare") {
      data <- data[ which(data[,group] %in% c("C", "NC")), ]
    } else if (param == "all") {
      data <- data[ which(data[,group] %!in% c("C", "NC")), ]
    } else {
      data <- data[ which(data[,group] == param), ]
    }
  } else if (group == "cell_type") {
    if (param != "allandcelltype") {
      data <- data[ which(data[,group] == param), ]
    }
  } else {
    if (param %!in% c("all", "compare")) {
      data <- data[ which(data[,group] == param), ]
    }
  }
  return(data)
}

subset_data <- function(data, clonal, Diagnosis, sex, cell_type, genes = NULL) {
  # Subset based on features
  if (!is.null(genes)) {
    data <- data[ ,c(genes, non_gene_cols)]
    if ((ncol(data) - length(non_gene_cols)) == length(genes)) {
      print(paste0("Successfully selected for ", length(genes), " feature(s)"))
    }
  }
  
  # Subset based on parameters
  data <- subset_parameters(data, clonal)
  data <- subset_parameters(data, Diagnosis)
  data <- subset_parameters(data, sex)
  data <- subset_parameters(data, cell_type)
  
  # Subset based on age range
  data <- data[which(data$age %in% seq(age_min, age_max)),]
  
  # Return data subset
  return(data)
}

run_loess <- function(gene, data) {
  # Extract gene specific data
  data <- data[, c(gene, non_gene_cols)]
  
  # Generate LOESS model and save to lo_data
  lo <- loess(get(gene) ~ age, data, span = 0.75)
  
  # Generate predicted values
  lo_predict <- predict(lo, data.frame(age = seq(age_min, age_max, 1))) %>%
    # lo_predict <- predict(lo, data.frame(age = seq(age_min, age_max, 0.0833))) %>%
    as.data.frame()
  colnames(lo_predict) <- gene
  
  return(lo_predict)
}

heatmap_of_loess <- function(lo_predict, output_dir, cell_type_label) {
  # Convert wide to long
  lo_predict_long <- gather(lo_predict, gene, expression, 2:ncol(lo_predict))
  
  # Order genes by hclustering
  clust_order <- hclust(dist(t(lo_predict[,2:ncol(lo_predict)])))$order
  lo_predict_long$gene <- factor(lo_predict_long$gene, levels = unique(lo_predict_long$gene)[clust_order])
  
  # Add limits for heatmap plotting
  zscore_limit <- 0.5
  lo_predict_long_plotting <- lo_predict_long
  lo_predict_long_plotting[ which(lo_predict_long_plotting$expression < -zscore_limit), "expression"] <- -zscore_limit
  lo_predict_long_plotting[ which(lo_predict_long_plotting$expression > zscore_limit), "expression"] <- zscore_limit
  
  # Generate heatmap
  p <- ggplot(lo_predict_long_plotting , aes(x = age, y = gene, fill = expression)) +
    theme_Publication_blank() +
    geom_raster(interpolate=TRUE) +
    theme(axis.text.y=element_blank()) +
    scale_fill_gradient2(low = "cyan", mid = "black", high = "yellow",
                         breaks = c(-zscore_limit, zscore_limit),
                         limits = c(-zscore_limit, zscore_limit))
  
  # Export heatmap
  ggplot2::ggsave(paste0(output_dir, "/", cell_type_label, "_loess_predicted_heatmap.pdf"), p, width = 5, height = 4)
}

line_plot_of_loess <- function(data, color, ylim = NULL, data_compare = NULL,
                               alpha = 0.05, size = 0.3, plot_ind_genes = TRUE) {
  # Find average expression over age of all genes in cluster
  avg <- data %>%
    group_by(age) %>%
    summarize_each(funs(mean, sd, se=sd(.)/sqrt(n())), expression)
  
  # Generate LOESS model and save to lo_data
  lo <- loess(mean ~ age, avg, span = 0.75)
  lo_predict <- predict(lo, data.frame(age = avg$age))
  avg <- add_column(avg, lo = lo_predict)
  
  # Define gene number annotation df
  annotations <- data.frame(
    xpos = c(-Inf),
    ypos =  c(Inf),
    annotateText = c(paste0("n=", length(unique(data$gene)))),
    hjustvar = c(-0.5) ,
    vjustvar = c(1))
  
  if (is.null(data_compare)) {
    # Generate plot
    plot <- ggplot(data, aes(x = age, y = expression)) +
      theme_Publication_blank() +
      geom_line(stat="smooth", method = "loess", span = 0.75, se = FALSE,
                aes(group = gene),
                alpha = .05, color = color, size = size) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      geom_line(data = avg,
                aes(x = age, y = avg_exp),
                stat="smooth", method = "loess", span = 0.75, se = FALSE,
                color = darken(color, amount = 0.3), size = size * 3) +
      theme(legend.position = "none",
            aspect.ratio = 1,
            text = element_text(size = 8),
            panel.margin = unit(c(0, 0, 0, 0), "null")) +
      geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) +
      {if(!is.null(ylim))ylim(ylim)}
  } else {
    # Find average expression over age of all genes in cluster
    avg_compare <- data_compare %>%
      group_by(age) %>%
      summarize_each(funs(mean, sd, se=sd(.)/sqrt(n())), expression)
    
    # Generate LOESS model and save to lo_data
    lo <- loess(mean ~ age, avg_compare, span = 0.75)
    lo_predict <- predict(lo, data.frame(age = avg_compare$age))
    avg_compare <- add_column(avg_compare, lo = lo_predict)
    
    # Find shared age range
    min_age_compare <- max(min(data$age), min(data_compare$age))
    max_age_compare <- min(max(data$age), max(data_compare$age))
    
    # Extract regression values
    p_tmp <- ggplot(data, aes(x = age, y = expression)) +
      geom_line(data = avg, aes(x = age, y = mean),
                stat="smooth", method = "loess", span = 0.75)
    stats <- ggplot_build(p_tmp)$data[[1]][seq(1,80,by=8),]
    
    p_tmp <- ggplot(data_compare, aes(x = age, y = expression)) +
      geom_line(data = avg_compare, aes(x = age, y = mean),
                stat="smooth", method = "loess", span = 0.75)
    compare_stats <- ggplot_build(p_tmp)$data[[1]][seq(1,80,by=8),]
    
    # Generate plot
    plot <- ggplot(data, aes(x = age, y = expression)) +
      theme_Publication_blank() +
      {if(plot_ind_genes)
        geom_line(stat="smooth", method = "loess", span = 0.75, se = FALSE,
                  aes(group = gene),
                  alpha = alpha, color = "gray", size = size)} +
      {if(plot_ind_genes)
        geom_line(data = data_compare,
                  aes(x = age, y = expression, group = gene),
                  stat="smooth", method = "loess", span = 0.75, se = FALSE,
                  alpha = alpha, color = color, size = size)} +
      geom_line(data = avg,
                aes(x = age, y = mean),
                stat="smooth", method = "loess", span = 0.75, se = FALSE,
                color = "black", size = size * 3) + #gray40
      geom_line(data = avg_compare,
                aes(x = age, y = mean),
                stat="smooth", method = "loess", span = 0.75, se = FALSE,
                color = color, size = size * 3) +
      geom_errorbar(inherit.aes = FALSE, data = stats, color = "black",
                    mapping = aes(x = x, ymin = y-se, ymax=y+se), width = 0.5, size = 0.3) +
      geom_errorbar(inherit.aes = FALSE, data = compare_stats, color = color,
                    mapping = aes(x = x, ymin = y-se, ymax=y+se), width = 0.5, size = 0.3) +
      theme(legend.position = "none",
            aspect.ratio = 1,
            text = element_text(size = 8),
            panel.margin = unit(c(0, 0, 0, 0), "null")) +
      coord_cartesian(xlim=c(min_age_compare, max_age_compare)) +
      geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) +
      {if(!is.null(ylim))ylim(ylim)}
  }
  return(plot)
}

cluster_loess <- function(data_scaled, output_dir, cell_type_label,
                          cluster_nums = 1:5, lo_predict = NULL, hclust_cut_merged = NULL) {
  # If not performing comparison
  if (!is.null(lo_predict)) {
    # Perform distance calculation on all trajectories (removing age column)
    dist_mat <- dist(t(lo_predict[,2:ncol(lo_predict)]))
    
    # Plot hierarchal clusters
    pdf(paste0(output_dir, "/", cell_type_label, "_dendogram_loess_predicted.pdf"))
    clust <- hclust(dist_mat, method = "complete")
    plot(clust, labels = FALSE)
    dev.off()
    
    # Initialize hclust_cut_list
    hclust_cut_list <- list()
  }
  
  # Generate clustering with k = cluster_nums
  for (k in cluster_nums) {
    # Check if we have enough genes
    if (!is.null(lo_predict)) {
      if (k > (ncol(lo_predict) - 1)) {
        print(paste0("Insufficient gene number for ", k, " clusters"))
        next
      }
      
      # Generate hierarchical clusters
      hclust_cut <- data.frame(cutree(clust, k = k))
    } else {
      # Isolate hierarchical clusters
      hclust_cut <- hclust_cut_merged[,match(k, cluster_nums),drop = FALSE]
    }
    
    # Generate color vector
    color_vector <- hue_pal()(k)
    
    # Initialize plot list
    hclust_plot_list <- list()
    
    for (i in seq(k)) {
      # Subset genes of interest
      cols <- c(rownames(hclust_cut[ which(hclust_cut[,1] == i), ,drop = FALSE]), "age")
      data_scaled_subset <- data_scaled[ ,cols ]
      
      if (!is.null(lo_predict)) {
        data_long <- gather(data_scaled_subset, gene, expression, -age)
        
        # Run plotting
        plot <- line_plot_of_loess(data_long, color = color_vector[i])
      } else {
        # Extract comparison parameters of interest
        group <- c("clonal", "Diagnosis", "sex")[grep("compare", c(clonal, Diagnosis, sex))]
        if (group == "clonal") {
          group1 <- "NC"
          group2 <- "C"
        } else if (group == "Diagnosis") {
          group1 <- "HC"
          group2 <- "MCI/AD"
        } else if (group == "sex") {
          group1 <- "m"
          group2 <- "f"
        }
        
        # Subset two groups
        data_tmp <- data_scaled_subset[which(data_scaled[,group] == group1),]
        data_tmp_compare <- data_scaled_subset[which(data_scaled[,group] == group2),]
        data_long <- gather(data_tmp, gene, expression, -age)
        data_long_compare <- gather(data_tmp_compare, gene, expression, -age)
        
        # Run plotting
        plot <- line_plot_of_loess(data_long, color = color_vector[i],
                                   data_compare = data_long_compare)
      }
      
      # Append plot to list
      hclust_plot_list[[i]] <- plot
    }
    
    # Generate figure of clustered trajectories
    if (k == 1) {
      pdf(paste0(output_dir, "/", cell_type_label, "_hclust_k", k, ".pdf"))
      print(hclust_plot_list[[1]])
      dev.off()
    } else {
      patchwork <- wrap_plots(hclust_plot_list, guides = "collect")
      if (is.null(lo_predict)) {
        patchwork <- patchwork + plot_annotation(subtitle = paste0(cell_type_label, ": gray = ", group1, ", color = ", group2))
      } else {
        patchwork <- patchwork + plot_annotation(subtitle = cell_type_label)
      }
      pdf(paste0(output_dir, "/", cell_type_label, "_hclust_k", k, ".pdf"))
      print(patchwork)
      dev.off()
    }
    
    # Append cut to list
    if (!is.null(lo_predict)) {hclust_cut_list[[k]] <- hclust_cut}
  }
  # Write out cluster cut
  if (!is.null(lo_predict)) {
    hclust_cut_merged <- as.data.frame(do.call(cbind, hclust_cut_list[cluster_nums]))
    write.csv(hclust_cut_merged, file = paste0(output_dir, "/", cell_type_label, "_hclust_cut.csv"))
  }
}

execute_on_celltype <- function(data, clonal, Diagnosis, sex, cell_type, output_parent_dir) {
  # Generate cell_type label for naming
  cell_type_label <- gsub("/","", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  
  # Define cell_type output dir
  output_dir <- file.path(output_parent_dir, cell_type_label)
  dir.create(output_dir, showWarnings = FALSE)
  
  # Choose features of interest
  load(paste0("/path/to/pseudobulk/results/", cell_type_label, "_over10_genes"))
  genes <- over10_genes
  names(data) <- gsub("-", ".", names(data))
  genes <- gsub("-", ".", genes)
  
  # Subset based on clonal, diagnosis, and sex
  data_subset <- subset_data(data,
                             clonal, Diagnosis, sex, cell_type,
                             genes = genes)
  
  # Scale genes
  data_scaled <- data_subset
  data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)] <- scale(data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)])
  data_scaled$age <- as.numeric(as.character(data_scaled$age))
  data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)] <- data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)][ , colSums(is.na(data_scaled[ ,(names(data_scaled) %!in% non_gene_cols)])) == 0]

  # Reassign genes
  genes <- names(data_scaled)
  genes <- genes[ genes %!in% non_gene_cols ]
  
  # If not running comparison
  if (length(grep("compare", c(clonal, Diagnosis, sex))) == 0) {
    # Run LOESS on all genes
    lo_predict_list <- lapply(genes, run_loess, data = data_scaled)
    
    # Merge lo_predict
    lo_predict <- as.data.frame(lo_predict_list)
    lo_predict$age <- seq(age_min, age_max, 1)
    lo_predict <- dplyr::relocate(lo_predict, age)
    
    # Generate heatmap of lo_predict, if not comparing
    heatmap_of_loess(lo_predict, output_dir, cell_type_label)
    
    # Generate clusters of lo_predict and plot data_scaled of clusters with regression
    cluster_loess(lo_predict = lo_predict, data_scaled = data_scaled,
                  output_dir = output_dir, cell_type_label = cell_type_label,
                  cluster_nums = cluster_nums)
    
    # Export lo_predict
    write.csv(lo_predict, paste0(output_dir, "/", cell_type_label, "_lo_predict.csv"))
  } else {
    # Load HC data for gene clusters
    hclust_cut_merged <- read.csv(paste0("path/to/loess/results/", ref_clustering,
                                         "/", cell_type_label, "/", cell_type_label, "_hclust_cut.csv"),
                                  row.names = 1)
    
    # Generate clusters of data_scaled using HC clusters
    cluster_loess(hclust_cut_merged = hclust_cut_merged, data_scaled = data_scaled,
                  output_dir = output_dir, cell_type_label = cell_type_label,
                  cluster_nums = cluster_nums)
  }
  # Export data scaled
  write.csv(data_scaled, paste0(output_dir, "/", cell_type_label, "_data_scaled.csv"))
}

#------------------------------------------------------------------------------
# Run Analysis

main <- function() {
  # Load data
  load(data_path)
  print(head(data[,1:10]))
  
  # Run LOESS
  if (length(cell_type) != 1) {
    for (type in cell_type) {
      execute_on_celltype(data, clonal, Diagnosis, sex, type, output_parent_dir)
    }
  } else if (cell_type != "allandcelltype") {
    execute_on_celltype(data, clonal, Diagnosis, sex, cell_type, output_parent_dir)
  } else {
    if (clonal == "all") {
      for (type in unique(data$cell_type)) {
        print(paste0("Analyzing ", type))
        execute_on_celltype(data, clonal, Diagnosis, sex, type, output_parent_dir)
      }
    } else {
      for (type in c("T Regulatory Cells", "CD8+ T Cells", "CD4+ T Cells", "all")) {
        print(paste0("Analyzing ", type))
        execute_on_celltype(data, clonal, Diagnosis, sex, type, output_parent_dir)
      }
    }
  }
}

main()