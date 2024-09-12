#Code for Wellington et al. 2024 Figure 2, 10X 3' sc-RNA-seq of CD34+, CD43+ or 
#CD45+ iPSC-derived cells

library(monocle3)
library(dplyr)
library(ggplot2)
library(Matrix)

#Note that there are areas in this code that contain underlying variability
#that may lead to slight variations in the distribution of cells in UMAP space
mat_path <- "~/cellranger/iPS_10X/matrix.mtx.gz"
feature_anno_path <- "~/cellranger/iPS_10X/features.tsv.gz"
cell_anno_path <- "~/cellranger/iPS_10X/barcodes.tsv.gz"

cds <- load_mm_data(mat_path = mat_path, 
                    feature_anno_path = feature_anno_path, 
                    cell_anno_path = cell_anno_path)

#Add a column to the metadata containing the number of UMIs for each cell
colData(cds)$n.umi <- Matrix::colSums(counts(cds))

#Look at the distribution of UMIs
#Note that these plots were originally generated with qplot, which has been
#deprecated
ggplot(data = colData(cds), aes(x=n.umi)) + geom_density()

#Set UMI Cutoffs
cds <- cds[,colData(cds)$n.umi > 1000 & colData(cds)$n.umi < 30000]

#Preprocess data
cds <- preprocess_cds(cds, num_dim = 5)

#Reduce dimensions of the data using UMAP
cds <- reduce_dimension(cds)

#Cluster cells
cds_clust <- cluster_cells(cds, resolution = 3.5e-3)

#Look at cells in UMAP
#Note that there is a random element in one of the monocle3 functions that may
#lead to some differences
plot_cells(cds_clust)