#!/usr/bin/env Rscript

# This script takes a loom object, builds a model to predict separability, and
# outputs the separability predictions for 1k-50k cells (by 1k)

################################ Arguments #####################################
# args[1]    Path to loom file
# args[2]    Fraction of cluster 1 in the full dataset
# args[3]    Fraction of cluster 2 in the full dataset
# args[4]    Path for output tsv file
################################################################################

suppressMessages(library(Seurat))
#suppressMessages(library(loomR))
suppressMessages(library(rhdf5))
suppressMessages(library(RANN))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
################################################################################

# Parameters
B_freq_tot = 1
B_freq_rel = 1
B_sil_curve = 1
B_sil_asym = 1
Int_asym = 1
B_ncells = 1
B_cov_rel = 1

args <- commandArgs(trailingOnly = TRUE)

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
seurat.object <- NormalizeData(seurat.object,display.progress = F,scale.factor = 1e4,normalization.method = "LogNormalize")
cluster.avgs <- AverageExpression(object = seurat.object, display.progress = FALSE)

#Silhouette Score Function
score_sil <- SilScore(seurat.object, cluster.avgs,
                      cluster_comparison = paste(cluster_1, ":", cluster_2),
                      cluster_ids = seurat.object@ident)

# saturation point
sat_point <- Int_asym + B_sil_asym * score_sil


# dependence on number of cells
num.cells <- seq(1000, 100000, 1000)
predicted.separability <- #SOMETHING
  
out <- cbind(num.cells, predicted.separability)
write.table(out, args[4])

# Simple prediction based soley on DE genes
num.cells <- seq(1000, 100000, 1000)
predicted.separability <- 1 / (1 + exp(-num.cells * num.de.genes / 1e7))
output <- data.frame(num.cells, predicted.separability)
colnames(output) <- c("Number of Cells", "Predicted Separability")
write.table(x = output, file = args[4], sep = "\t", quote = FALSE, row.names = FALSE)
