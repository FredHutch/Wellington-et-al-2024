#Label transfer of previously identified cell types for iPS Day 8, Crosse et al 
#and Calvanese et al (more complete cell types)
#(Zeng labels were added pre-integration, see iPS_Human-Embryo_Preprocessing_Integration.RDS)

library(monocle3)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(pals)
library(cowplot)

#Load Monocle integrated object; note that the Google Drive version already has 
#these labels transferred

cds <- readRDS("~/RDS_objects/integrated_data_Monocle.RDS")

#####Add Crosse et al. labels (acquired from Edie Crosse) to the objects so they will####
#be included in the final integrated object
#Note that this only includes labels for Aov embryos (no AoMid) integrated
crosse_ids <- read.csv("~/csv_files/Crosse-et-al_Cluster_IDs.csv")
#store ids associated with ventral 1 and ventral 2 in two variabless
labels_ventral1 <- crosse_ids[str_detect(crosse_ids$index, 
                                         "Human10XVentral1"),]
labels_ventral2 <- crosse_ids[str_detect(crosse_ids$index, 
                                         "Human10XVentral2"),]
#create corrected ids (without all the extra text and matching barcode formatting)
#Note that this was done separately for the two included embryos (ventral 1 and
#ventral 2) to avoid issues with potential shared barcodes
new_cells_v1 <- sub("[^:]*:([^x]*)x", 
                    "\\1", 
                    labels_ventral1$index)
new_cells_v2 <- sub("[^:]*:([^x]*)x", 
                    "\\1", 
                    labels_ventral2$index)
new_clusters_v1 <- sub("/[^/]*", 
                       "", 
                       labels_ventral1$clusters)
new_clusters_v2 <- sub("/[^/]*", 
                       "", 
                       labels_ventral2$clusters)
v1.mod1 <- labels_ventral1["clusters"]
v1.mod2 <- data.frame(new_clusters_v1)
v1.mod2 <- rename(v1.mod2, 
                  "cluster"="new_clusters_v1")
v1.mod2$barcode <- new_cells_v1
v2.mod1 <- labels_ventral2["clusters"]
v2.mod2 <- data.frame(new_clusters_v2)
v2.mod2 <- rename(v2.mod2, 
                  "cluster" = "new_clusters_v2")
v2.mod2$barcode <- new_cells_v2
labels_all <- rbind(v1.mod2, 
                    v2.mod2)
labels_all_nodup <- labels_all[!duplicated(labels_all$barcode),]
labels_all_nodup$barcode <- paste0(labels_all_nodup$barcode, 
                                   "_5") #_5 was assigned by Seurat
row.names(labels_all_nodup) <- c(labels_all_nodup$barcode)
labels_all_nodup <- subset(labels_all_nodup, 
                           select=-(barcode))
labels_crosse <- labels_all_nodup
#recode clusters to cell type; cell types were expanded based on annotation in 
#the text of Crosse et al. 2020
labels_crosse$cluster <- dplyr::recode(labels_crosse$cluster,
                                       "1" = "Arterial Endothelium",
                                       "2" = "Endothelium",
                                       "3" = "Endothelium",
                                       "4" = "Endothelium",
                                       "5" = "Hemogenic Endothelium",
                                       "6" = "Endothelium",
                                       "7" = "Endothelium",
                                       "8" = "Endothelium",
                                       "9" = "Endothelium",
                                       "10" = "HSPC",
                                       "11" = "Hematopoietic",
                                       "12" = "Mesenchymal",
                                       "13" = "Mesenchymal",
                                       "14" = "Pericytes",
                                       "15" = "Mesenchymal",
                                       "16" = "SNS",
                                       "17" = "EMX2+ Epithelium",
                                       "18" = "Unknown PERPhi",
                                       "19" = "Unknown PTH+",
                                       "20" = "PGC")
colnames(labels_crosse) <- c("crosse_cell_types_complete")


######Add iPS Day 8 labels######
cds_ips <- readRDS("~/RDS_objects/iPS_10X_preprocessed_marker_scores.RDS")
#extract cell types assigned in Figure 02
labels_ips <- colData(cds_ips)["cell_type"]
colnames(labels_ips) <- "45-iPS Cell Type"
labels_ips_rownames <- paste0(rownames(labels_ips), "_4") #_4 was assigned by Seurat
rownames(labels_ips) <- labels_ips_rownames

#####Add more complete form of Calvanese et al. labels####
#Clusters identified by Calvanese et al. pulled from Fig1_17/seurat_metadata.csv
#at https://ucla.app.box.com/s/koqumtfbja7ll66ucffd291eexppro8c

#open saved file from above (re-saved for easy access)
calv_labels <- read.csv("~/csv_files/Calvanese-et-al_CS14-15_cellidents.csv", 
                        row.names = 1)

#add new metadata column with cell identities from Calvanese et al. Figure 1 and
#Extended Data Figure 1
calv_labels$calvanese_cell_types_complete <- as.character(calv_labels$seurat_clusters)
calv_labels$calvanese_cell_types_complete <- dplyr::recode(calv_labels$calvanese_cell_types_complete,
                                                 "0" = "Endothelium",
                                                 "1" = "Endothelium",
                                                 "2" = "Other",
                                                 "3" = "Endothelium",
                                                 "4" = "Stroma",
                                                 "5" = "Stroma",
                                                 "6" = "Stroma",
                                                 "7" = "Stroma",
                                                 "8" = "Stroma",
                                                 "9" = "Other",
                                                 "10" = "Stroma",
                                                 "11" = "Stroma",
                                                 "12" = "HSC",
                                                 "13" = "Endothelium",
                                                 "14" = "Stroma",
                                                 "15" = "Epithelium",
                                                 "16" = "Epithelium",
                                                 "17" = "Stroma",
                                                 "18" = "Erythroid Progenitor",
                                                 "19" = "Stroma")
#suffix _7 was added to all calvanese cell barcodes along with
#prefixes wk4_ for -3 cells, wk5_1_ for -1 cells, wk5_2_ for -2 cells
#Convert rownames of calv_labels to match integrated barcode form
new_names <- c()
for (row_name in rownames(calv_labels)) {
  if (grepl("(*)-1", row_name)) {
    new_name <- paste0("wk5_1_", row_name, "_7")
    new_names <- c(new_names, new_name)
  } else if (grepl("(*)-2", row_name)) {
    new_name <- paste0("wk5_2_", row_name, "_7")
    new_names <- c(new_names, new_name)
  } else if (grepl("(*)-3", row_name)) {
    new_name <- paste0("wk4_", row_name, "_7")
    new_names <- c(new_names, new_name)
  }
}

#convert to new rownames
rownames(calv_labels) <- new_names

#remove all rows except cell type
calv_labels <- calv_labels[dim(calv_labels)[2]]

####Adding cell type labels to cds####

#Extract integrated metadata as dataframe
cds_metadata <- colData(cds)
#Add column containing a number corresponding to current row order
cds_metadata$id <- 1:length(rownames(cds_metadata))
#merge the cell type labels into the metadata
cds_metadata <- merge(cds_metadata, 
                      labels_ips,
                      by = "row.names",
                      all.x = TRUE)
rownames(cds_metadata) <- cds_metadata$Row.names
cds_metadata <- subset(cds_metadata, select = -c(Row.names))

cds_metadata <- merge(cds_metadata, 
                      labels_crosse, 
                      by = "row.names",
                      all.x = TRUE)
rownames(cds_metadata) <- cds_metadata$Row.names
cds_metadata <- subset(cds_metadata, select = -c(Row.names))

cds_metadata <- merge(cds_metadata, 
                      calv_labels, 
                      by = "row.names",
                      all.x = TRUE)
rownames(cds_metadata) <- cds_metadata$Row.names
cds_metadata <- subset(cds_metadata, select = -c(Row.names))

#Ensure that the order of cells is correct before adding metadata back to cds
cds_metadata <- cds_metadata[order(cds_metadata$id),]

#Create copy of cds
cds_copy <- cds

#Add new metadata columns to cds_copy
colData(cds_copy)$"crosse_cell_types_complete" <- cds_metadata$"crosse_cell_types_complete"
colData(cds_copy)$"ips_cell_type" <- cds_metadata$"45-iPS Cell Type"
colData(cds_copy)$"calvanese_cell_types_complete" <- cds_metadata$"calvanese_cell_type_complete"

#Subset cds for each label transfer plot
cds_crosselabels <- cds_copy[,!is.na(colData(cds_copy)$"crosse_cell_types_complete")]
cds_ipslabels <- cds_copy[,!is.na(colData(cds_copy)$"ips_cell_type")]
cds_zenglabels <- cds_copy[,!is.na(colData(cds_copy)$"Zeng_Label")]
cds_calvlabels <- cds_copy[,!is.na(colData(cds_copy)$"calvanese_cell_types_complete")]

#For iPS labels
plt_ipslabels <- plot_cells(cds_ipslabels, 
                            color_cells_by = "ips_cell_type",
                            cell_size = 0.8,
                            label_cell_groups = FALSE,
                            rasterize = TRUE,
                            show_trajectory_graph = FALSE) + 
  simple_theme +
  scale_color_manual(values = as.vector(glasbey(11)))
ggsave(plot = plt_ipslabels + theme(legend.position = "none"), 
       filename = paste0(save_dir, 
                         "int_ips_labels.svg"), 
       width = 7, 
       height = 7)

plt_ipslegend <- cowplot::get_legend(plt_ipslabels)
svg(filename = paste0(save_dir, "int_ips_legend.svg"),
    width = 4,
    height = 4)
grid.newpage()
grid.draw(plt_ipslegend)
dev.off()

#For Zeng lables
plt_zenglabels <- plot_cells(cds_zenglabels, 
                             color_cells_by = "Zeng_Label",
                             cell_size = 1,
                             label_cell_groups = FALSE,
                             rasterize = TRUE,
                             show_trajectory_graph = FALSE) + 
  simple_theme +
  scale_color_manual(values = as.vector(glasbey(11)))
ggsave(plot = plt_zenglabels + theme(legend.position = "none"), 
       filename = paste0(save_dir, 
                         "int_zeng_labels.svg"), 
       width = 7, 
       height = 7)

plt_zenglegend <- cowplot::get_legend(plt_zenglabels)
svg(filename = paste0(save_dir, "int_zeng_legend.svg"),
    width = 4,
    height = 4)
grid.newpage()
grid.draw(plt_zenglegend)
dev.off()

#For Calvanese labels (from Calvanese et al. Figure 1)
plt_calvlabels <- plot_cells(cds_calvlabels, 
                             color_cells_by = "calvanese_cell_types_complete",
                             cell_size = 1,
                             label_cell_groups = FALSE,
                             rasterize = TRUE,
                             show_trajectory_graph = FALSE) + 
  simple_theme +
  scale_color_manual(values = as.vector(glasbey(6)))
ggsave(plot = plt_calvlabels + theme(legend.position = "none"), 
       filename = paste0(save_dir, 
                         "int_calv_labels.svg"), 
       width = 7, 
       height = 7)

plt_calvlegend <- cowplot::get_legend(plt_calvlabels)
svg(filename = paste0(save_dir, "int_calv_legend.svg"),
    width = 4,
    height = 4)
grid.newpage()
grid.draw(plt_calvlegend)
dev.off()

#For Crosse labels (there were only labels for AoV, no AoMid)
plt_crosselabels <- plot_cells(cds_crosselabels, 
                               color_cells_by = "crosse_cell_types_complete",
                               cell_size = 1,
                               label_cell_groups = FALSE,
                               rasterize = TRUE,
                               show_trajectory_graph = FALSE) + 
  simple_theme +
  scale_color_manual(values = as.vector(glasbey(20)))
ggsave(plot = plt_crosselabels + theme(legend.position = "none"), 
       filename = paste0(save_dir, 
                         "int_crosse_labels.svg"), 
       width = 7, 
       height = 7)

plt_crosselegend <- cowplot::get_legend(plt_crosselabels)
svg(filename = paste0(save_dir, "int_crosse_legend.svg"),
    width = 4,
    height = 4)
grid.newpage()
grid.draw(plt_crosselegend)
dev.off()

#####To label single label transfer populations#####

#Extract UMAP coordinates for each cell
colData(cds_copy)$umap_1 = reducedDims(cds_copy)[["UMAP"]][,1]
colData(cds_copy)$umap_2 = reducedDims(cds_copy)[["UMAP"]][,2]

#create vectors of columns containing labels and the abbreviation to be used
cols <- c("ips_cell_type", "Zeng_Label", "crosse_cell_types_complete", "calvanese_cell_types_complete")
abbr <- c("ips_", "zeng_", "crosse_", "calv_")

#Remove all strange characters (+, /, -, " ") from ips and crosse labels (to 
#avoid regex issues in grepl)
crosse_updated <- str_replace_all(colData(cds_copy)$"crosse_cell_types_complete", 
                                  fixed("+"), 
                                  "")
crosse_updated <- str_replace_all(crosse_updated, 
                                  fixed(" "), 
                                  "_")
ips_updated <- str_replace_all(colData(cds_copy)$ips_cell_type, 
                               fixed(" "), 
                               "_")
ips_updated  <- str_replace_all(ips_updated, 
                                fixed("("), 
                                "")
ips_updated  <- str_replace_all(ips_updated, 
                                fixed(")"), 
                                "")
ips_updated  <- str_replace_all(ips_updated, 
                                fixed("-"), 
                                "_")
ips_updated  <- str_replace_all(ips_updated, 
                                fixed("/"), 
                                "_")

#Add updated cell types to metadata
colData(cds_copy)$crosse_cell_types_complete <- crosse_updated
colData(cds_copy)$ips_cell_type <- ips_updated

#Create new metadata columns containing calls for single populations (i.e. HE (Calvanese))
for (i in 1:length(cols)){
  cell_types <- unique(colData(cds_copy)[,cols[i]])
  cell_types <- cell_types[!is.na(cell_types) & cell_types != "NA"]
  for (cell_type in cell_types) {
    colData(cds_copy)[,paste0(abbr[i], cell_type)] <- ifelse(grepl(cell_type, 
                                                                   colData(cds_copy)[,cols[i]]), 
                                                             "TRUE", 
                                                             "FALSE")
  }
}

#Create dataframe containing metadata
cds_colData <- colData(cds_copy) %>% as.data.frame()
col_names <- gsub("[.]","_",colnames(cds_colData))
col_names <- gsub("[..]","_",col_names)
colnames(cds_colData) <- col_names

#Pull column names for columns containing data for individual labels
cell_type_cols <- colnames(cds_colData)[c(19:length(colnames(cds_colData)))]

#Set colors for FALSE and TRUE calls
colors <- c("FALSE" = "grey",
            "TRUE" = "#F8766D")

#Function for creating a clean plot
simple_theme <-  theme(axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       panel.grid = element_blank(),
                       axis.title = element_blank(),
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       axis.line.x = element_blank(),
                       axis.line.y = element_blank(), 
                       panel.border = element_blank(),
                       plot.margin = grid::unit(c(0,0,0,0), "mm"),
                       legend.position = "none",
                       plot.background = element_rect(fill = "white", 
                                                      color = "white"))

#Set save directory destination for plots
save_dir <- "~/save_dir/"

#Create plots for each individual label that was transferred to view where the
#population with that label ends up in the integrated UMAP space
#This plots TRUE cells on top of the FALSE cells so they are more visible
for(cell_type_col in cell_type_cols) {
  plt <- ggplot() +
    geom_point(data = cds_colData[cds_colData[,cell_type_col] == "FALSE",],
               aes(x = umap_1,
                   y = umap_2),
               color = "black",
               stroke = 0,
               size = 2) +
    geom_point(data = cds_colData[cds_colData[,cell_type_col] == "FALSE",],
               aes(x = umap_1,
                   y = umap_2),
               color = "grey",
               stroke = 0,
               size = 1.9) +
    geom_point(data = cds_colData[cds_colData[,cell_type_col] == "TRUE",],
               aes(x = umap_1,
                   y = umap_2),
               color = "black",
               stroke = 0,
               size = 2) +
    geom_point(data = cds_colData[cds_colData[,cell_type_col] == "TRUE",],
               aes(x = umap_1,
                   y = umap_2,
                   color = get(cell_type_col)),
               stroke = 0,
               size = 1.9) +
    theme_void() +
    scale_color_manual(values = colors) +
    simple_theme
  ggsave(plot = plt, 
         filename = paste0(save_dir, cell_type_col, ".png"),
         width = 5,
         height = 5)
}