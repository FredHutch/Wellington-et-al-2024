#Label Endothelial and Hematopoietic populations

library(viridis)
library(viridisLite)
library(monocle3)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

#Load the preprocessed RDS file from Google Drive
cds <- readRDS("~/RDS_objects/sci_cds_BBI_preprocessed.RDS")
save_dir <- "~/save_dir/" #set by individual user

#set up simple_theme for plotting
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

#function for aggregate gene scoring for a marker set
gene_group_scoring <- function(cds, gene_group, column_name) {
  cds_gene_group <- cds[fData(cds)$gene_short_name %in% gene_group,]
  aggregate_expression <- exprs(cds_gene_group)
  aggregate_expression <- Matrix::t(Matrix::t(aggregate_expression) / pData(cds_gene_group)$Size_Factor)
  aggregate_expression <- Matrix::colSums(aggregate_expression)
  pData(cds)[,column_name] <- log(aggregate_expression+1)
  return(cds)
}

#Hematopoietic signature (RUNX1, CD43, CD45)
hemato_markers <- c("RUNX1","SPN", "PTPRC") #These are from Calvanese Fig 1
cds <- gene_group_scoring(cds, hemato_markers, column_name = "hemato_score")

#Create plot showing hematopoietic score (not shown in paper)
col <- colorRampPalette(c("lightgrey", "blue", "darkblue"))(256)
hemato_score <- plot_cells(cds, 
                           color_cells_by = "hemato_score", 
                           cell_size = 0.1, 
                           show_trajectory_graph = FALSE) + 
  scale_color_gradientn(colors = col, 
                        limits = c(0, 
                                   max(colData(cds)$hemato_score))) +
  labs(color = "Hematopoietic Score") +
  simple_theme 

ggsave(hemato_score, 
       filename = paste0(save_dir, "hemato_score.png"), 
       width = 6, 
       height = 4)

#Endothelial signature (RUNX1, CD34, CDH5)
endo_markers <- c("CDH5", "KDR", "TIE1", "CD34") #These are Calvanese Fig 1 Endo genes + CD34
cds <- gene_group_scoring(cds, endo_markers, column_name = "endo_score")
endo_not_hemato_scores <- c()
for (x in 1:nrow(colData(cds))) {
  if (colData(cds)$hemato_score[x] != 0) {
    endo_not_hemato_scores <- c(endo_not_hemato_scores, 0)
  }
  else {
    endo_not_hemato_scores <- c(endo_not_hemato_scores, 
                                colData(cds)$endo_score[x])
  }
}
endo_not_hemato_scores <- unname(endo_not_hemato_scores)
colData(cds)$endo_not_hemato_score <- endo_not_hemato_scores

#Create plot showing endothelial score (not shown in paper)
col <- colorRampPalette(c("lightgrey", "blue", "darkblue"))(256)
he_not_hemato_score_plt <- plot_cells(cds, 
                                      color_cells_by = "endo_not_hemato_score", 
                                      cell_size = 0.1, 
                                      show_trajectory_graph = FALSE) + 
  scale_color_gradientn(colors = col, 
                        limits = c(0, 
                                   max(colData(cds)$endo_score))) +
  labs(color = "Endothelial Score") +
  simple_theme 

ggsave(endo_not_hemato_score_plt, 
       filename = paste0(save_dir, "endo_not_hemato_score.png"), 
       width = 5.25, 
       height = 4)

#Create endo_or_hemato metadata column with cell type (either Endo, Hemato or NA)
cds_coldata <- data.frame(colData(cds))
cds_coldata <- cds_coldata %>% mutate(endo_or_hemato = case_when(endo_not_hemato_score > 0 ~ "Endo",
                                                                 endo_not_hemato_score == 0 & hemato_score > 0 ~ "Hemato",
                                                                 endo_not_hemato_score == 0 & hemato_score == 0 ~ NA))
colData(cds)$endo_or_hemato <- cds_coldata$endo_or_hemato

#Label cell type calls in UMAP
he_or_hemato_plt <- plot_cells(cds, 
                               color_cells_by = "endo_or_hemato",
                               cell_size = 0.1,
                               show_trajectory_graph = FALSE,
                               label_cell_groups = FALSE) +
  scale_color_manual(values = c("#005eb8", "#990000", "#999999")) +
  simple_theme

#Label cell type calls in individual plots for each day of differentiation
d7 <- cds[,colData(cds)$time.point %in% c("Day 07")]
d8 <- cds[,colData(cds)$time.point %in% c("Day 08")]
d11 <- cds[,colData(cds)$time.point %in% c("Day 11")]
d14 <- cds[,colData(cds)$time.point %in% c("Day 14")]
d18 <- cds[,colData(cds)$time.point %in% c("Day 18")]

#Plot for day 7
d7_plt <- plot_cells(d7, 
                     color_cells_by = "endo_or_hemato",
                     cell_size = 0.5,
                     show_trajectory_graph = FALSE,
                     label_cell_groups = FALSE) +
  scale_color_manual(values = c("#005eb8", "#990000", "#999999")) +
  simple_theme + 
  theme(legend.position = "none", plot.margin=grid::unit(c(0,0,0,0), "mm"))

#Plot for day 8
d8_plt <- plot_cells(d8, 
                     color_cells_by = "endo_or_hemato",
                     cell_size = 0.5,
                     show_trajectory_graph = FALSE,
                     label_cell_groups = FALSE) +
  scale_color_manual(values = c("#005eb8", "#990000", "#999999")) +
  simple_theme + 
  theme(legend.position = "none", plot.margin=grid::unit(c(0,0,0,0), "mm"))

#Plot for day 11
d11_plt <- plot_cells(d11, 
                      color_cells_by = "endo_or_hemato",
                      cell_size = 0.5,
                      show_trajectory_graph = FALSE,
                      label_cell_groups = FALSE) +
  scale_color_manual(values = c("#005eb8", "#990000", "#999999")) +
  simple_theme + 
  theme(legend.position = "none", plot.margin=grid::unit(c(0,0,0,0), "mm"))

#Plot for day 14
d14_plt <- plot_cells(d14, 
                      color_cells_by = "endo_or_hemato",
                      cell_size = 0.5,
                      show_trajectory_graph = FALSE,
                      label_cell_groups = FALSE) +
  scale_color_manual(values = c("#005eb8", "#990000", "#999999")) +
  simple_theme + 
  theme(legend.position = "none", plot.margin=grid::unit(c(0,0,0,0), "mm"))

#Plot for day 18
d18_plt <- plot_cells(d18, 
                      color_cells_by = "endo_or_hemato",
                      cell_size = 0.5,
                      show_trajectory_graph = FALSE,
                      label_cell_groups = FALSE) +
  scale_color_manual(values = c("#005eb8", "#990000", "#999999")) +
  simple_theme + 
  theme(legend.position = "none", plot.margin=grid::unit(c(0,0,0,0), "mm"))

#Create a plot of all days separated into panels
combo_plt <- d7_plt + 
  d8_plt + 
  d11_plt + 
  d14_plt + 
  d18_plt + 
  plot_layout(ncol = 5)

ggsave(plot = combo_plt, 
       filename = paste0(save_dir, "days_separate_colored_by_endo_hemato.png"), 
       width = 25,
       height = 5,
       units = "in")

#For hematopoietic-related cell populations only
cds_hemato_only <- readRDS("~/sci_cds_BBI_preprocessed_hemato-only.RDS")

d7 <- cds[,colData(cds)$time.point %in% c("Day 07") & colData(cds)$cell %in% colData(cds_hemato_only)$cell]
d8 <- cds[,colData(cds)$time.point %in% c("Day 08") & colData(cds)$cell %in% colData(cds_hemato_only)$cell]
d11 <- cds[,colData(cds)$time.point %in% c("Day 11") & colData(cds)$cell %in% colData(cds_hemato_only)$cell]
d14 <- cds[,colData(cds)$time.point %in% c("Day 14") & colData(cds)$cell %in% colData(cds_hemato_only)$cell]
d18 <- cds[,colData(cds)$time.point %in% c("Day 18") & colData(cds)$cell %in% colData(cds_hemato_only)$cell]
d21 <- cds[,colData(cds)$time.point %in% c("Day 21") & colData(cds)$cell %in% colData(cds_hemato_only)$cell]

#Plot for day 7
d7_plt <- plot_cells(d7, 
                     color_cells_by = "endo_or_hemato",
                     cell_size = 0.5,
                     show_trajectory_graph = FALSE,
                     label_cell_groups = FALSE) +
  scale_color_manual(values = c("#005eb8", "#990000", "#999999")) +
  simple_theme + 
  theme(legend.position = "none", plot.margin=grid::unit(c(0,0,0,0), "mm"))

#Plot for day 8
d8_plt <- plot_cells(d8, 
                     color_cells_by = "endo_or_hemato",
                     cell_size = 0.5,
                     show_trajectory_graph = FALSE,
                     label_cell_groups = FALSE) +
  scale_color_manual(values = c("#005eb8", "#990000", "#999999")) +
  simple_theme + 
  theme(legend.position = "none", plot.margin=grid::unit(c(0,0,0,0), "mm"))

#Plot for day 11
d11_plt <- plot_cells(d11, 
                      color_cells_by = "endo_or_hemato",
                      cell_size = 0.5,
                      show_trajectory_graph = FALSE,
                      label_cell_groups = FALSE) +
  scale_color_manual(values = c("#005eb8", "#990000", "#999999")) +
  simple_theme + 
  theme(legend.position = "none", plot.margin=grid::unit(c(0,0,0,0), "mm"))

#Plot for day 14
d14_plt <- plot_cells(d14, 
                      color_cells_by = "endo_or_hemato",
                      cell_size = 0.5,
                      show_trajectory_graph = FALSE,
                      label_cell_groups = FALSE) +
  scale_color_manual(values = c("#005eb8", "#990000", "#999999")) +
  simple_theme + 
  theme(legend.position = "none", plot.margin=grid::unit(c(0,0,0,0), "mm"))

#Plot for day 18
d18_plt <- plot_cells(d18, 
                      color_cells_by = "endo_or_hemato",
                      cell_size = 0.5,
                      show_trajectory_graph = FALSE,
                      label_cell_groups = FALSE) +
  scale_color_manual(values = c("#005eb8", "#990000", "#999999")) +
  simple_theme + 
  theme(legend.position = "none", plot.margin=grid::unit(c(0,0,0,0), "mm"))

combo_plt <- d7_plt + d8_plt + d11_plt + d14_plt + d18_plt + plot_layout(ncol = 5)
ggsave(plot = combo_plt, 
       filename = paste0(save_dir, "hemato-only_days_separate_colored_by_endo_hemato.png"), 
       width = 9,
       height = 2)