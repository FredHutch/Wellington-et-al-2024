#Conversion of the integrated Seurat object to Monocle3 for further analysis

library(Seurat)
library(monocle3)

#Load Seurat object
#Note: don't forget to add the seurat cluster info as a metadata column
int_obj <- readRDS("~/RDS_objects/integrated_data_Seurat.RDS")

#Pull counts, metadata, UMAP and PCA information from Seurat object
counts <- int_obj@assays$RNA@counts
metadata <- int_obj@meta.data
umap_mat <- int_obj@reductions$umap@cell.embeddings
pca_mat <- int_obj@reductions$pca@cell.embeddings

#Create gene name data.frame (to be used for gene_metadata argument)
genes <- rownames(counts)
genes_df <- data.frame(genes)
colnames(genes_df) <- c("gene_short_name")
rownames(genes_df) <- genes_df$gene_short_name

#Create Monocle3 object
cds <- new_cell_data_set(counts,
                         cell_metadata = metadata, 
                         gene_metadata = genes_df)

#Add UMAP and PCA to monocle object
reducedDim(x = cds, type = "UMAP") <- umap_mat
reducedDim(x = cds, type = "PCA") <- pca_mat

#Calculate new clusters via Monocle for later analysis
cds <- cluster_cells(cds)