#Code for Wellington et al. 2024 Figure 1, sci-RNA-seq of EBs on Days 7-21

library(monocle3)
library(stringr)
library(viridis)
library(viridisLite)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(RCurl)

#load the cell dataset object containing all data as provided by the Brotman
#Baty Institute post-sequencing
#Load the RDS file from Google Drive
cds <- readRDS("~/RDS_objects/sci_cds_BBI.RDS") 

#Create sample metadata column containing sample information (BBI sample codename)
cells <- colnames(cds)
samples <- str_extract(cells, "(?=LIG).*")
samples <- str_extract(samples, "[^_]+$")
colData(cds)$sample <- samples

#Create metadata column containing the iPS line used to generate each sample
colData(cds)$line <- samples
colData(cds)$line <- dplyr::recode(colData(cds)$line, 
                                       "SD34" = "MSC-iPS", 
                                       "SD35" = "MSC-iPS", 
                                       "SD36" = "45-iPS", 
                                       "SD37" = "45-iPS", 
                                       "SD38" = "MSC-iPS", 
                                       "SD39" = "MSC-iPS", 
                                       "SD40" = "45-iPS", 
                                       "SD41" = "45-iPS", 
                                       "SD46" = "MSC-iPS", 
                                       "SD47" = "MSC-iPS", 
                                       "SD48" = "45-iPS", 
                                       "SD49" = "45-iPS", 
                                       "SD50" = "MSC-iPS", 
                                       "SD51" = "MSC-iPS", 
                                       "SD52" = "45-iPS", 
                                       "SD53" = "45-iPS", 
                                       "SD54" = "MSC-iPS", 
                                       "SD55" = "MSC-iPS", 
                                       "SD56" = "45-iPS", 
                                       "SD57" = "45-iPS", 
                                       "SD58" = "MSC-iPS", 
                                       "SD59" = "MSC-iPS", 
                                       "SD60" = "45-iPS", 
                                       "SD61" = "45-iPS") 

#Create metadata column including both the iPS line and day of differentiation
colData(cds)$line_and_time.point <- samples
colData(cds)$line_and_time.point <- dplyr::recode(colData(cds)$line_and_time.point, 
                                                      "SD34" = "MSC-iPS Day 07", 
                                                      "SD35" = "MSC-iPS Day 07", 
                                                      "SD36" = "45-iPS Day 07", 
                                                      "SD37" = "45-iPS Day 07", 
                                                      "SD38" = "MSC-iPS Day 08", 
                                                      "SD39" = "MSC-iPS Day 08", 
                                                      "SD40" = "45-iPS Day 08", 
                                                      "SD41" = "45-iPS Day 08", 
                                                      "SD46" = "MSC-iPS Day 11", 
                                                      "SD47" = "MSC-iPS Day 11", 
                                                      "SD48" = "45-iPS Day 11", 
                                                      "SD49" = "45-iPS Day 11", 
                                                      "SD50" = "MSC-iPS Day 14", 
                                                      "SD51" = "MSC-iPS Day 14", 
                                                      "SD52" = "45-iPS Day 14", 
                                                      "SD53" = "45-iPS Day 14", 
                                                      "SD54" = "MSC-iPS Day 18", 
                                                      "SD55" = "MSC-iPS Day 18", 
                                                      "SD56" = "45-iPS Day 18", 
                                                      "SD57" = "45-iPS Day 18", 
                                                      "SD58" = "MSC-iPS Day 21", 
                                                      "SD59" = "MSC-iPS Day 21", 
                                                      "SD60" = "45-iPS Day 21", 
                                                      "SD61" = "45-iPS Day 21") 

#Create metadata column containing day of differentiation for each cell
colData(cds)$time.point <- samples
colData(cds)$time.point <- dplyr::recode(colData(cds)$time.point, 
                                             "SD34" = "Day 07", 
                                             "SD35" = "Day 07", 
                                             "SD36" = "Day 07", 
                                             "SD37" = "Day 07", 
                                             "SD38" = "Day 08", 
                                             "SD39" = "Day 08", 
                                             "SD40" = "Day 08", 
                                             "SD41" = "Day 08", 
                                             "SD46" = "Day 11", 
                                             "SD47" = "Day 11", 
                                             "SD48" = "Day 11", 
                                             "SD49" = "Day 11", 
                                             "SD50" = "Day 14", 
                                             "SD51" = "Day 14", 
                                             "SD52" = "Day 14", 
                                             "SD53" = "Day 14", 
                                             "SD54" = "Day 18", 
                                             "SD55" = "Day 18", 
                                             "SD56" = "Day 18", 
                                             "SD57" = "Day 18", 
                                             "SD58" = "Day 21", 
                                             "SD59" = "Day 21", 
                                             "SD60" = "Day 21", 
                                             "SD61" = "Day 21") 

#Calculate mitochondrial umi percentage and add metadata column
mito_genes <- c("ENSG00000198727","ENSG00000198695","ENSG00000198786", 
                "ENSG00000212907","ENSG00000198886","ENSG00000198840",
                "ENSG00000198938","ENSG00000198899","ENSG00000228253",
                "ENSG00000198712","ENSG00000198804","ENSG00000198763",
                "ENSG00000198888")
colData(cds)$n.mito <- Matrix::colSums(counts(cds[mito_genes]))
cds$perc_mito_umi <- 100*(cds$n.mito / cds$n.umi)

#minimum n.umi cutoff was already set by the BBI to 100
#maximum n.umi is 62698, which is reasonable so no upper limit of n.umi was set
#only include cells that contain < 10% mitochondrial reads
cds <- cds[, colData(cds)$perc_mito_umi < 10]

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds, resolution = 1e-4)

#Look at the distribution of cells within UMAP space
#Note: one of the functions in Monocle has a random element, so it may look a
#bit different than the publication if you run the above newly
plot_cells(cds)
plot_cells(cds, color_cells_by = "line", label_cell_groups = FALSE)
plot_cells(cds, color_cells_by = "line_and_time.point", label_cell_groups = FALSE)
plot_cells(cds, color_cells_by = "time.point", label_cell_groups = FALSE)