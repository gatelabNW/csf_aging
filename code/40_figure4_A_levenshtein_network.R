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
# Date: 03-02-2022
# Written by: Benoit Lehallier, Hamilton Oh, Natalie Piehl
# Summary: Generate TCR network of Levenshtein Similarities
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
  library("Seurat")
  library("RecordLinkage")
  library("qgraph")
})

# Initialize paths
seurat_object <- "path/to/seurat_object/"
output_dir <- "path/to/export/results"

# Source helper functions
source("code/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set levenshtein similarity cutoff 
lev_cutoff <- 0.9

#------------------------------------------------------------------------------
# Data preparation

# Load in seurat object
load(seurat_object)

# Load TCR data from seurat object
clonotypes <- s[[c("frequency_filtered", "tra_cdr3s", "trb_cdr3s", "clonotype_id", "ID", "age", "Diagnosis")]]
clonotypes <- clonotypes[!duplicated(clonotypes[,"clonotype_id"]),]

# Add age bins
samples_54_64 <- unique((subset(clonotypes, age >= 54 & age <= 64))$ID)
samples_65_67 <- unique((subset(clonotypes, age >= 65 & age <= 67))$ID)
samples_68_72 <- unique((subset(clonotypes, age >= 68 & age <= 72))$ID)
samples_73_82 <- unique((subset(clonotypes, age >= 73 & age <= 82))$ID)
clonotypes <- mutate(clonotypes, tmp =
                       case_when(ID %in% samples_54_64 ~ "54_64",
                                 ID %in% samples_65_67 ~ "65_67",
                                 ID %in% samples_68_72 ~ "68_72",
                                 ID %in% samples_73_82 ~ "73_82"))

# Add diagnostic values
clonotypes$tmp[ which(clonotypes$Diagnosis == "CI") ] <- "CI"

# Subset and rename columns of interest
clonotypes <- clonotypes[,c("frequency_filtered", "tra_cdr3s", "trb_cdr3s", "ID", "tmp")]
names(clonotypes) <- c("Frequency","CDR3a","CDR3b","id","Diagnosis")

# Filter
clonotypes <- clonotypes[ which(clonotypes$Frequency > 1), ] # Select for clonal
clonotypes <- clonotypes[ which(clonotypes$CDR3a!="" & clonotypes$CDR3b!=""), ]  # Select for paired TCRs
clonotypes <- clonotypes[ -grep(";",clonotypes$CDR3a), ] # Select for unambiguous CDR3a
clonotypes <- clonotypes[ -grep(";",clonotypes$CDR3b), ] # Select for unambiguous CDR3b

# Create clonotype_id.long as unique id for each clonotype
clonotypes$clonotype_id.long <- paste(clonotypes$CDR3a,
                                      clonotypes$CDR3b, sep = "_")

#------------------------------------------------------------------------------
# Network construction (Fig 4A)

# Load homemade function (from Hamilton/Benoit)
mat2edge<-function(x,cutoff){
  if(missing(cutoff)==T){cutoff<-0}
  res=NULL
  for(i in 1:nrow(x)){
    for(k in 1:ncol(x)){
      if(x[i,k]>cutoff){
        res=rbind(res,data.frame(X=rownames(x)[i],Y=colnames(x)[k],weight=x[i,k],stringsAsFactors = F))
      }
    }
    print(i)
  }
  res
}

# Generate table of unique IDs with corresponding Diagnosis
ForNetwork1 <- data.frame(Diagnosis=clonotypes$Diagnosis, id=clonotypes$id, stringsAsFactors = F)
ForNetwork1 <- ForNetwork1[which(duplicated(paste(ForNetwork1[,1], ForNetwork1[,2], sep="_"))==F),]
colnames(ForNetwork1) <- c("x", "y")

# Calculate levenshtein similarity between all possible pairs of TCRab
lsim_df <- NULL
print(length(1:length(clonotypes$clonotype_id)))
for(i in 1:length(clonotypes$clonotype_id)){
  lsim_df <- rbind(lsim_df,
                   levenshteinSim(clonotypes$clonotype_id.long[i],
                                  clonotypes$clonotype_id.long))
  if (i %% 100 == 0) {
    print(i)
  }
}

# rename rows/cols
colnames(lsim_df) <- clonotypes$clonotype_id.long
rownames(lsim_df) <- clonotypes$clonotype_id.long
write.csv(lsim_df, paste0(output_dir, "lsim_matrix.csv"))

# Make heatmap of TCR similarity
lsim_df.2 <- lsim_df
lsim_df.2[lsim_df.2 < lev_cutoff] <- 0  #only visualize pairs with levenshteinSim > threshold
diag(lsim_df.2) <- 0
lsim_df.2 <- lsim_df.2[which(apply(lsim_df.2,1,function(x) sum(x>0,na.rm=T))>0),
                       which(apply(lsim_df.2,1,function(x) sum(x>0,na.rm=T))>0)]
pairs.breaks <- seq(0, 1, by=0.01)
mycol <- colorpanel(n=length(pairs.breaks)-1,low="darkslateblue",mid="red",high="yellow")
pdf(paste0(output_dir, "Lsim_heatmap_over", lev_cutoff, ".pdf"))
heatmap.2(lsim_df.2,
          cexRow=.2,cexCol=.2,
          trace="none",
          dendrogram="both",
          breaks=pairs.breaks,
          col=mycol,
          Rowv=T,key=F,
          Colv=T,
          lhei=c(0.2,4),
          lwid=c(.2,3))
dev.off()

# Generate network
{
  # Generate table of ID, clonotype_id.long, and frequency
  ForNetwork2 = data.frame(clonotypes$id,
                           clonotypes$clonotype_id.long,
                           clonotypes$Frequency,
                           stringsAsFactors = F)
  ForNetwork2 = ForNetwork2[which(duplicated(paste(
    ForNetwork2[, 1], ForNetwork2[, 2], ForNetwork2[, 3], sep = "_"
  )) == F), ]
  colnames(ForNetwork2) = c("x", "y", "weight")
  
  # Format lsim table
  ForNetwork3 = lsim_df
  ForNetwork3[upper.tri(ForNetwork3)] <- NA
  diag(ForNetwork3) <- NA
  ForNetwork3[is.na(ForNetwork3)] <- 0
  
  # Generate table of clonotype_id.long pairs and lsims
  ForNetwork3 = mat2edge(x = ForNetwork3, cutoff = lev_cutoff)
  colnames(ForNetwork3) = c("x", "y", "weight")
  
  {
    # Combine all three
    ForNetwork1$weight = 5
    toQgraph = rbind(ForNetwork1, ForNetwork2)
    toQgraph = rbind(toQgraph, ForNetwork3)
    
    # Customize qgraph
    {
      # Set color of disease groups
      l=unique(c(toQgraph[,1],toQgraph[,2]))
      
      col.tmp=rep(adjustcolor("grey85",.3),length(l))
      col.tmp[which(l %in% c("54_64"))]<-adjustcolor("yellow",.8)
      col.tmp[which(l %in% c("65_67"))]<-adjustcolor("orange",.8)
      col.tmp[which(l %in% c("68_72"))]<-adjustcolor("red",.8)
      col.tmp[which(l %in% c("73_82"))]<-adjustcolor("purple",.8)
      col.tmp[which(l %in% c("CI"))]<-adjustcolor("green",.8)
      
      # Set color of subjects nodes
      col.tmp[which(l %in% unique(clonotypes$id[which(clonotypes$Diagnosis=="54_64")]))]<-adjustcolor("yellow",.3)
      col.tmp[which(l %in% unique(clonotypes$id[which(clonotypes$Diagnosis=="65_67")]))]<-adjustcolor("orange",.3)
      col.tmp[which(l %in% unique(clonotypes$id[which(clonotypes$Diagnosis=="68_72")]))]<-adjustcolor("red",.3)
      col.tmp[which(l %in% unique(clonotypes$id[which(clonotypes$Diagnosis=="73_82")]))]<-adjustcolor("purple",.3)
      col.tmp[which(l %in% unique(clonotypes$id[which(clonotypes$Diagnosis=="CI")]))]<-adjustcolor("green",.3)
      
      # Set variable size of nodes
      scaling.factor <- 5
      size.tmp = rep(1, length(l))
      size.tmp[grep("1|2|3|4|5|6|7|8", l)] <- 7
      size.tmp[which(l %in% c("54_64", "65_67", "68_72", "73_82", "CI"))] <- 12
      for (i in seq(length(l))) {
        # Extract node value
        node = l[i]
        
        # Average frequencies for each TCR seq
        if (node %in% ForNetwork2[, 2]) {
          node_freq <- ForNetwork2[which(ForNetwork2[, 2] == node), "weight"]
          mean_freq <- mean(node_freq)
          # Set frequency limit of 20
          if (mean_freq > 20) {
            mean_freq <- 20
          }
        } else {
          next
        }
        
        # Update size.tmp with mean frequency divided by scaling factor
        size.tmp[i] <- mean_freq / scaling.factor + 0.5
      }
      
      # Set color of edges
      line.tmp=rep("black",nrow(toQgraph))
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("54_64"))]<-"yellow"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("65_67"))]<-"orange"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("68_72"))]<-"red"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("73_82"))]<-"purple"
      line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("CI"))]<-"green"
      
      # Essentially remove labels
      labels.cex.tmp=l
      labels.cex.tmp=rep(0.00000000001,length(l))
      labels.cex.tmp[grep("1|2|3|4|5|6|7|8",l)]<-0.00000000001
      labels.cex.tmp[which(l %in% c("54_64", "65_67", "68_72", "73_82", "CI"))]<-0.00000000001
    }
  }
  
  # Set weights
  toQgraph$weight <- 1
  toQgraph$weight[ 
    which(toQgraph$x %in% c("54_64", "65_67", "68_72", "73_82", "CI")) ] <- 1.3
  
  # Generate plot
  pdf(paste0(output_dir, "Lsim_network_over", lev_cutoff, ".pdf"))
  qgraph(toQgraph, 
         color=col.tmp,
         vsize=size.tmp/2,
         edge.color=line.tmp, 
         labels=TRUE, 
         label.cex=labels.cex.tmp, 
         directed=F,
         layout="spring")
  dev.off()
}