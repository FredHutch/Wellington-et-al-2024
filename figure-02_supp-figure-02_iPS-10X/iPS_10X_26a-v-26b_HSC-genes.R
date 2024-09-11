#Look at the 6 HSC-related genes identified by Calvanese et al. 2022 in clusters
#26a vs. 26b

cds <- readRDS("~/iPS_10X_preprocessed_marker_scores.RDS")

#Select only cluster 26
clust26 <- cds[,colData(cds)$new_cluster %in% c("26a", "26b")]

#Create a vector of the markers associated with HSCs from Calvanese et al. 2022
hsc_genes <- c("HLF","HOXA9","MECOM", "MLLT3", "RUNX1", "SPINK2")

#Plot the expression of each gene in each subcluster for comparison
hsc_exp <- plot_genes_by_group(ips_clust26, 
                               hsc_genes, 
                               group_cells_by = "new_cluster",
                               ordering_type = "none",
                               axis_order = "group_marker",
                               max.size = 10) +
  theme(text = element_text(size = 20)) +
  xlab("Cluster")