#Marker scoring for sci-RNA-seq clusters

library(monocle3)
library(VGAM) 
library(viridis)
library(stringr)
library(tibble) 
library(reticulate)
library(pheatmap)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(sf)
library(Rcpp)
library(scales)
library(ComplexHeatmap)
library(tidyverse)
library(pals)
library(svglite)
library(cowplot)

#Load the preprocessed RDS file from github
#cds <- getURL("https://raw.github.com/...")
#cds <- readRDS(text = cds)
cds <- r"(C:\Users\rachw\Documents\GitHub\code_wellington_etal2024_EBprofiling\objects\sci_cds_BBI_preprocessed.RDS)"
cds <- gsub("\\\\", "/", cds)
cds <- readRDS(cds)

#Add cluster numbers as a metdata column
colData(cds)$cluster <- clusters(cds)

#Markers are from supplemental table S1A from Xu, Y., et al. 2023. Nat Cell Biol.
#"A Single Cell Transcriptome Atlas Profiles Early Organogenesis in Human Embyros"
neural_progenitor <- c("SOX2")
neuron <- c("ELAVL4", "TUBB3")
epidermis <- c("PERP", "KRT18")
sensory_neuron <- c("SIX1", "ISL1", "TLX3")
schwann <- c("FOXD3", "MPZ", "PLP1", "S100B")
craniofacial <- c("ALX1", "ALX3", "DLX1", "DLX2")
head_mesoderm <- c("SNAI2", "TWIST1", "FOXC1", "FOXC2")
somite <- c("PAX1", "PAX9", "SOX9", "PAX3", "FOXC1", "FOXC2", "MYF5")
im <- c("PAX8", "WT1", "PAX2")
somatic_lpm <- c("IRX3")
limb <- c("TBX5", "PITX1")
splanchnic_lpm <- c("FOXF1")
endothelium <- c("PLVAP", "ECSCR")
blood <- c("CD53", "PLAC8", "CD68", "HBA1", "CORO1A")
endoderm <- c("FOXA2", "AFP", "CLDN6")
pgc <- c("POU5F1", "DDX4", "DPPA3", "DAZL", "NANOS3", "DND1", "ALPL", "SALL4")
epithelium <- c("KRT18", "KRT19", "CLDN4")
fibroblast <- c("COL1A2", "COL3A1")

gene_groups <- list(neural_progenitor, neuron, epidermis, sensory_neuron, schwann,
                    craniofacial, head_mesoderm, somite, im, somatic_lpm, limb,
                    splanchnic_lpm, endothelium, blood, endoderm, pgc, epithelium,
                    fibroblast)
group_names <- c("neural_progenitor", "neuron", "epidermis", "sensory_neuron", 
                 "schwann", "craniofacial", "head_mesoderm", "somite", "im",
                 "somatic_lpm", "limb", "splanchnic_lpm", "endothelium", "blood",
                 "endoderm", "pgc", "epithelium", "fibroblast")

#marker set scoring function; required for multi_group_scoring
gene_group_scoring <- function(cds, gene_group, column_name) {
  cds_gene_group <- cds[fData(cds)$gene_short_name %in% gene_group,]
  aggregate_expression <- exprs(cds_gene_group)
  aggregate_expression <- Matrix::t(Matrix::t(aggregate_expression) / pData(cds_gene_group)$Size_Factor)
  aggregate_expression <- Matrix::colSums(aggregate_expression)
  pData(cds)[,column_name] <- log(aggregate_expression+1)
  return(cds)
}

#multiple marker set scoring function
multi_group_scoring <- function(cds, 
                                gene_groups, 
                                group_names, 
                                cds_save_directory = "NONE") {
  x = 1
  cds_copy = cds
  for (gene_group in gene_groups) {
    column_name <- paste(sub("_", " ", group_names[x]), " Score", sep = "")
    cds_copy <- gene_group_scoring(cds_copy, 
                                   gene_group, 
                                   column_name = column_name)
    x = x + 1
  }
  if (cds_save_directory != "NONE") {
    saveRDS(cds_copy, file = paste(cds_save_directory, "cds_with_gene-set_scores.RDS", 
                                   sep = ""))
  }
  return(cds_copy)
}

cds <- multi_group_scoring(cds,
                          gene_groups = gene_groups,
                          group_names = group_names)

paletteLength <- 100
mycol <- colorRampPalette(c("navy", "blue", "lightblue", "white", "yellow", "red", "red4"))(paletteLength)
heatmap_save_dir <- "~/save_dir/" #unique for user

#Heatmap fucnction
avg_score_matrix <- function(cds, cell_group, scores){
  avg_score_list<-lapply(levels(factor(cds[[cell_group]])), function(x){
    sub<-cds@colData[cds[[cell_group]] == x, scores] %>% as.data.frame()
    avg<-colMeans(sub)
  })
  names(avg_score_list)<-levels(factor(cds[[cell_group]]))
  mat<-do.call(rbind, avg_score_list)
  mat<-t(scale(mat))
  mat
}

#Create heatmap of scores
score_names <- colnames(colData(cds))[-c(1:17,36)]
mat_sci_scores <- avg_score_matrix(cds = cds,
                                   cell_group = "cluster",
                                   scores = score_names)
mat_rownames <- rownames(mat_sci_scores)
mat_rownames <- gsub("\\.", " ", mat_rownames)
mat_rownames <- gsub(" Score", "", mat_rownames)
rownames(mat_sci_scores) <- mat_rownames

#Set colors for heatmap
mybreaks <- c(c(seq(min(floor(mat_sci_scores)), 0, length.out=ceiling(paletteLength/2)+1), 
                seq(max(ceiling(mat_sci_scores))/paletteLength, max(mat_sci_scores), length.out=floor(paletteLength/2))))

#Create svg file containing heatmap of Xu et al S1A scores per cluster
svglite(filename = paste0(heatmap_save_dir, 
                          "Xu_etal_2023_Heatmap_FigS1A.svg"),
        height = 5,
        width = 10)
my_pheatmap <- pheatmap(mat_sci_scores[row.names(mat_sci_scores)
                                       ,c(1:max(as.numeric(colData(cds)$cluster)))],
                        legend=T, 
                        show_rownames = T, 
                        show_colnames = T, 
                        cluster_cols=F, 
                        cluster_rows = F,
                        cellheight = 10, 
                        cellwidth = 10,
                        color = mycol,
                        breaks = mybreaks) 
draw(my_pheatmap)
dev.off()

#Create png file containing heatmap of Xu et al S1A scores per cluster
png(filename = paste0(heatmap_save_dir, 
                      "Xu_etal_2023_Heatmap_FigS1A.png"),
    height = 5,
    width = 10,
    units = "in",
    res = 1200)
my_pheatmap <- pheatmap(mat_sci_scores[row.names(mat_sci_scores)
                                       ,c(1:max(as.numeric(colData(cds)$cluster)))],
                        legend=T, 
                        show_rownames = T, 
                        show_colnames = T, 
                        cluster_cols=F, 
                        cluster_rows = F,
                        cellheight = 10, 
                        cellwidth = 10,
                        color = mycol,
                        breaks = mybreaks) 
draw(my_pheatmap)
dev.off()

#annotate Xu et al S1A maximal cell types in sci data as cell_type in metadata
df <- data.frame(mat_sci_scores)
cols <- colnames(df)
max_cell_types <- c()

for (col in cols) {
  max_cell_type <- rownames(df[df[col] == max(df[col]),])
  max_cell_types <- c(max_cell_types, max_cell_type)
}

clusters <- as.numeric(unname(colData(cds)$cluster))
cell_types <- c()
for (clust in clusters) {
  cell_type <- max_cell_types[clust]
  cell_types <- c(cell_types, cell_type)
}

colData(cds)$s1a_cell_type <- cell_types

#Plot Xu et al S1A cluster max scored cell type in UMAP
simple_theme <-  theme(axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       panel.grid = element_blank(),
                       axis.title = element_blank(),
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       axis.line.x = element_blank(),
                       axis.line.y = element_blank(), 
                       panel.border = element_rect(color = "black", 
                                                   linewidth = 2, 
                                                   fill = NA)) 

plt <- plot_cells(cds, 
                  color_cells_by = "cell_type",
                  show_trajectory_graph = FALSE,
                  group_label_size = 4,
                  rasterize = TRUE) + 
  simple_theme +
  scale_color_manual(values = as.vector(stepped(length(unique(colData(cds)$cell_type)))))

save_dir <- "~/save_dir/" #different for each user
ggsave(plt, filename = paste0(save_dir, "XuFigS1A_cell-types.png"),
       width = 8.5,
       height = 8,
       dpi = 600)

#Markers are from supplemental table S1B from Xu, Y., et al. 2023. Nat Cell Biol.
#"A Single Cell Transcriptome Atlas Profiles Early Organogenesis in Human Embyros"

s1b_labels <- read.csv("~/path_to_csv.csv", header = TRUE)
s1b_labels_df <- data.frame(s1b_labels)
s1b_labels_df <- s1b_labels_df %>% select(annotation, marker_genes_in_literature)

for (i in 1:dim(s1b_labels_df)[1]) {
  marker_genes <- s1b_labels_df$marker_genes_in_literature[i]
  marker_genes <- strsplit(s1b_labels_df$marker_genes_in_literature[i], split = "\\+")[[1]]
  assign(s1b_labels_df$annotation[i], marker_genes)
}

gene_groups <- s1b_labels_df$annotation

#marker set scoring function; required for multi_group_scoring;
#modified for calling gene_group string as variable
gene_group_scoring <- function(cds, gene_group, column_name) {
  cds_gene_group <- cds[fData(cds)$gene_short_name %in% get(gene_group),]
  aggregate_expression <- exprs(cds_gene_group)
  aggregate_expression <- Matrix::t(Matrix::t(aggregate_expression) / pData(cds_gene_group)$Size_Factor)
  aggregate_expression <- Matrix::colSums(aggregate_expression)
  pData(cds)[,column_name] <- log(aggregate_expression+1)
  return(cds)
}

save_dir <- "~/save_dir/" #different for each user
heatmap_save_dir <- "~/save_dir/" #different for each user

#calculate marker set score for each cluster for each potential cell type
cds <- multi_group_scoring(cds = cds,
                          gene_groups = gene_groups,
                          group_names = gene_groups)

#Pull score names as set of columns
score_names <- colnames(colData(cds))[-c(1:17,158)] #may need to be modified
mat_sci_scores <- avg_score_matrix(cds = cds,
                                   cell_group = "cluster",
                                   scores = score_names)
mat_rownames <- rownames(mat_sci_scores)
mat_rownames <- gsub("\\.", " ", mat_rownames)
mat_rownames <- gsub(" Score", "", mat_rownames)
rownames(mat_sci_scores) <- mat_rownames

#Set breaks for coloring in heatmap
mybreaks <- c(c(seq(min(floor(mat_sci_scores)), 0, length.out=ceiling(paletteLength/2)+1), 
                seq(max(ceiling(mat_sci_scores))/paletteLength, max(mat_sci_scores), length.out=floor(paletteLength/2))))

#Create svg file containing heatmap of Xu et al S1B scores per cluster
svglite(filename = paste0(heatmap_save_dir, 
                          "Xu_etal_2023_Heatmap_FigS1B_clustered.svg"),
        height = 20,
        width = 10)
my_pheatmap <- pheatmap(mat_sci_scores[row.names(mat_sci_scores)
                                       ,c(1:max(as.numeric(colData(cds)$cluster)))],
                        legend=T, 
                        show_rownames = T, 
                        show_colnames = T, 
                        cluster_cols=F, 
                        cluster_rows = T,
                        cellheight = 10, 
                        cellwidth = 10,
                        color = mycol,
                        breaks = mybreaks) 
draw(my_pheatmap)
dev.off()

#Create png file containing heatmap of Xu et al S1B scores per cluster
png(filename = paste0(heatmap_save_dir, 
                      "Xu_etal_2023_Heatmap_FigS1B.png"),
    height = 20,
    width = 10,
    units = "in",
    res = 1200)
my_pheatmap <- pheatmap(mat_sci_scores[row.names(mat_sci_scores)
                                       ,c(1:max(as.numeric(colData(cds)$cluster)))],
                        legend=T, 
                        show_rownames = T, 
                        show_colnames = T, 
                        cluster_cols=F, 
                        cluster_rows = F,
                        cellheight = 10, 
                        cellwidth = 10,
                        color = mycol,
                        breaks = mybreaks) 
draw(my_pheatmap)
dev.off()

#determine the maximal cell type score for each cluster and store in dataframe
df <- data.frame(mat_sci_scores)
cols <- colnames(df)
max_cell_types <- c()

for (col in cols) {
  max_cell_type <- rownames(df[df[col] == max(df[col]),])
  max_cell_types <- c(max_cell_types, max_cell_type)
}

clusters <- as.numeric(unname(colData(cds)$cluster))
cell_types <- c()
for (clust in clusters) {
  cell_type <- max_cell_types[clust]
  cell_types <- c(cell_types, cell_type)
}

#add maximal cell type as cell_type metadata column
colData(cds)$cell_type <- cell_types

#Plot cluster cell types in UMAP
simple_theme <-  theme(axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       panel.grid = element_blank(),
                       axis.title = element_blank(),
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       axis.line.x = element_blank(),
                       axis.line.y = element_blank(), 
                       panel.border = element_rect(color = "black", 
                                                   linewidth = 2, 
                                                   fill = NA)) 

plt <- plot_cells(cds, 
                  color_cells_by = "cell_type",
                  show_trajectory_graph = FALSE,
                  group_label_size = 4,
                  rasterize = TRUE) + 
  simple_theme +
  scale_color_manual(values = as.vector(ocean.phase(length(unique(colData(cds)$cell_type)))))

save_dir <- "~/save_dir/" #different for each user
ggsave(plt, filename = paste0(save_dir, "XuFigS1B_cell-types.png"),
       width = 8.5,
       height = 8,
       dpi = 600)
