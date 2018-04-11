#!/usr/bin/env Rscript

# This script takes a loom object, builds a model to predict separability, and
# outputs the separability predictions for 1k-50k cells (by 1k)

################################ Arguments #####################################
# args[1]    Path to loom file
# args[2]    Name of cluster 1 
# args[3]    Fraction of cluster 1 in the full dataset
# args[4]    Name of cluster 2 
# args[5]    Fraction of cluster 2 in the full dataset
# args[6]    Path for output tsv file
################################################################################

suppressMessages(library(Seurat))
#suppressMessages(library(loomR))
suppressMessages(library(rhdf5))
suppressMessages(library(RANN))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
################################################################################

# Parameters
B0 = -4.8
B_logit_score_sil = 2.17
B_frequency_cluster_tot = 0.51
B_log_ncell = 0.71

args <- commandArgs(trailingOnly = TRUE)
loom_file <- args[1]
freq1 <- as.numeric(args[3])
freq2 <- as.numeric(args[5])
cluster_1 <- args[2]
cluster_2 <- args[4]
outfile <- args[6]

#lf <- connect(filename = args[1])
data <- h5read(loom_file,name = "/")
data_matrix <- t(data$matrix)
dimnames(data_matrix) <- list(data$row_attrs$gene_names, data$col_attrs$cell_names)
clusters <- data$col_attrs$cluster; 
names(clusters) <- colnames(data_matrix)
rm(data)
umi.matrix <- Matrix(data_matrix,sparse=T)

# Get UMI matrix and cluster information from the loom object
#umi.matrix <- t(lf[["matrix"]][,])
#rownames(umi.matrix) <- lf[["row_attrs/gene_names"]][]
#colnames(umi.matrix) <- lf[["col_attrs/cell_names"]][]
#clusters <- lf[["col_attrs/cluster"]][]

seurat.object <- CreateSeuratObject(raw.data = umi.matrix)
seurat.object <- SetIdent(object = seurat.object,
                          cells.use = seurat.object@cell.names,
                          ident.use = clusters)
seurat.object <- NormalizeData(seurat.object,
                               display.progress = F,
                               scale.factor = 1e4,
                               normalization.method = "LogNormalize")

# Silhouette Score Metric
SilScore <- function(seurat.object, cluster_1, cluster_2, clusters) {
  
  # Compute Cluster Averages
  avg_cluster <- AverageExpression(seurat.object,return.seurat = T,show.progress = F)
  
  # Identify member cells
  is_compared <- clusters %in% c(cluster_1,cluster_2)
  
  # Projective Average Sil Width
  x <- t(t(avg_cluster@data[,cluster_1,drop = FALSE] -
           avg_cluster@data[,cluster_2,drop = FALSE]) %*% seurat.object@data[,is_compared])      
  d <- dist(x)
  g <- as.numeric(factor(clusters[is_compared]))
  
  return(mean(cluster::silhouette(dist = d,x=g)[,"sil_width"]))
}

score_sil <- SilScore(seurat.object = seurat.object,
                      cluster_1 = cluster_1,
                      cluster_2 = cluster_2,
                      clusters = clusters)
logit_score_sil <- boot::logit(score_sil/2 + 1/2)

# Total Cluster Frequency Metric
frequency_cluster_tot <- log(freq1) + log(freq2)

# dependence on number of cells
num.cells <- seq(1000, 100000, 1000)
predicted.separability <- boot::inv.logit(B0 +
                                          B_logit_score_sil*logit_score_sil +
                                          B_frequency_cluster_tot*frequency_cluster_tot +
                                          B_log_ncell*log(num.cells))
  
output <- cbind(num.cells, predicted.separability)
colnames(output) <- c("Number of Cells", "Predicted Separability")
write.table(x = output, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
