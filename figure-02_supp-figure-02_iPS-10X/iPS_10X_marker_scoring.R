#Marker scoring for iPSC Day 8 10X clusters

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

#colors to be used in heatmap
mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")

#Function for calculating the aggregate gene score for a set of expressed markers
gene_group_scoring <- function(cds, gene_group, column_name) {
  cds_gene_group <- cds[fData(cds)$gene_short_name %in% gene_group,]
  aggregate_expression <- exprs(cds_gene_group)
  aggregate_expression <- Matrix::t(Matrix::t(aggregate_expression) / pData(cds_gene_group)$Size_Factor)
  aggregate_expression <- Matrix::colSums(aggregate_expression)
  pData(cds)[,column_name] <- log(aggregate_expression+1)
  return(cds)
}

#Function for running the gene_group_scoring function on multiple marker sets
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

#Heatmap matrix function
avg_score_matrix<-function(cds, cell_group, scores){
  avg_score_list<-lapply(levels(factor(cds[[cell_group]])), function(x){
    sub<-cds@colData[cds[[cell_group]] == x, scores] %>% as.data.frame()
    avg<-colMeans(sub)
  })
  names(avg_score_list)<-levels(factor(cds[[cell_group]]))
  mat<-do.call(rbind, avg_score_list)
  mat<-t(scale(mat))
  mat
}

###Gene sets are from:
#Calvanese et al.: https://www.nature.com/articles/s41586-022-04571-x
#Zeng et al.: https://www.nature.com/articles/s41422-019-0228-6

#Load the preprocessed iPSC day 8 10X data
cds <- readRDS("~/iPS_10X_preprocessed.RDS")

#Add column to metadata containing original clustering info
colData(cds)$cluster <- clusters(cds)

#Separate top and bottom sections of cluster 26
cds_26 <- cds[,colData(cds)$cluster %in% c("26")] #only cluster 26 cells
cds_26a <- choose_cells(cds_26) #all cells in top section
cds_26b <- choose_cells(cds_26) #all cells in bottom section

#Get cell names for selected cells in 26a and 26b
cells_26a <- rownames(colData(cds_26a))
cells_26b <- rownames(colData(cds_26b))

#Create a dataframe of the cds metadata
cds_metadata <- colData(cds)
cds_metadata <- data.frame(cds_metadata)

#Create new column in metadata dataframe for identifying cells that are 26a/b
cds_metadata$new_cluster <- cds_metadata$cluster
cds_metadata$new_cluster <- ifelse(cds_metadata$barcode %in% cells_26a,
                                  "26a", cds_metadata$new_cluster)
cds_metadata$new_cluster <- ifelse(cds_metadata$barcode %in% cells_26b,
                                  "26b", cds_metadata$new_cluster)

#Create a factor containing cell to cluster info
new_cluster <- factor(cds_metadata$new_cluster, 
                      levels = c(1:25, "26a", "26b", 27:28))

#Create a new cds containing the new cluster info
cds_26_split <- cds
cds_26_split$new_cluster <- new_cluster

#Add scores to cds with cluster 26 split: 26a = top, 26b = bottom
#Commented figures are from the associated paper (either Calvanese or Zeng)
#Calvanese et al. markers
endo <- c("CDH5", "KDR", "TIE1") #fig 1
vec <- c("APLNR", "NRP2", "NR2F2") #composite of fig 1&3
calv_aec <- c("GJA5", "CXCR4") #fig 1
calv_hema <- c("RUNX1", "PTPRC", "SPN") #fig 1
hsc <- c("HLF", "SPINK2", "HOXA9", "RUNX1", "MLLT3", "MECOM") #general paper
mo <- c("C1QA", "CD14", "LYVE1") #fig 1
gr <- c("LYZ", "RNASE2") #fig 1
stroma <- c("COL1A1", "PDGFRA", "CXCL12", "POSTN", "PTN") #fig 1
fibro_sub <- c("PAX1", "SOX9", "HAND1", "NKX2-3", "DCN", "ALDH1A2", 
               "COL14A1", "CRABP1", "LUM", "PAX3", "TBX5", "HAND2") #fig 1
peric <- c("REN", "GATA3") #fig 1
sm <- c("ACTC1", "ACTA2") #fig 1
podoc <- c("NPHS2", "DSC1", "NPHS1") #fig 1
calv_epi <- c("EPCAM", "AFP") #fig 1
liver <- c("FGB", "APOA1") #fig 1
kidney <- c("MAL", "CALB1") #fig 1
ery <- c("HBE1", "HBZ", "GYPA") #fig 1
prolif <- c("MKI67", "TOP2A", "AURKB") #fig 1
pre_he <- c("PALMD", "TMEM100", "EDN1", "LTBP4", "HEY2", "SULF1", "IL33", 
            "CYP26B1", "ADGRG6", "COL23A1", "GATA6", "BMX", "TMCC3", "DKK1", 
            "AGTR2", "FBN2", "ELN", "ALDH1A1", "NKX2-3", "PROCR", "GATA3", 
            "GBP4") #fig 3
calv_he <- c("ADGRG6", "COL23A1", "GATA6", "BMX", "TMCC3", "DKK1", "AGTR2", "FBN2", 
             "ELN", "ALDH1A1", "NKX2-3", "PROCR", "GATA3", "GBP4", "MYCN", "KCNK17", 
             "MYB", "STAT5A", "SMIM24", "RAB27B", "SPINK2") #fig 3

#Zeng et al. markers
zeng_aec <- c("TFPI", "GJA5", "IGFBP4", "RDX", "RAMP2", "HSPG2", "CLEC14A", "TMEM100",
              "PLVAP", "CYR61") #Fig 1F
zeng_he <- c("SYNHG16", "RPL5", "RPL12", "RUNX1", "GAS5", "NPM1", "RPL6", "RPSAP58",
             "EIF3E", "LIDHB") #Fig 1F
zeng_hema <- c("HCST", "RAC2", "FYB", "LCP1", "TYROBP", "HOPX", "NKG7", "SPINK2", 
               "CD52", "PLAC8") #Fig 1F
twist_mes <- c("NRP2", "TWIST2", "CXCL12") #Fig 1C
dlk_mes <- c("CXCL12", "DLK1", "HAND1") #Fig 1C
fbln5_mes <- c("NRP2", "CXCL12", "HAND1", "FBLN5") #Fig 1C
wt1_mes <- c("NRP2", "WT1") #Fig 1C
zeng_epi <- c("NRP2", "CDX2") #Fig 1C

#Add marker set scores to cds
gene_groups <- list(endo, vec, calv_aec, calv_hema, hsc, mo, gr, stroma, fibro_sub,
                    peric, sm, podoc, calv_epi, liver, kidney, ery, prolif, pre_he,
                    calv_he, zeng_aec, zeng_he, zeng_hema, twist_mes, dlk_mes,
                    fbln5_mes, wt1_mes, zeng_epi)
group_names <- c("Endothelium (Calvanese)", "Venous Endothelium (Calvanese)", 
                 "Arterial Endothelium (Calvanese)", 
                 "Hematopoeitic (Calvanese)", "HSC (Calvanese)", 
                 "Monocyte/Macrophage (Calvanese)", "Granulocyte (Calvanese)", 
                 "Stroma (Calvanese)", "Fibroblast (Calavanese)",
                 "Pericyte (Calvanese)", "SM (Calvanese)", 
                 "Podocyte (Calvanese)", "Epithelium (Calvanese)", 
                 "Liver (Calvanese)", "Kidney (Calvanese)", 
                 "Erythrocyte (Calvanese)", "Proliferation (Calvanese)", 
                 "Pre-HE (Calvanese)", "HE (Calvanese)", 
                 "Arterial Endothelium (Zeng)", "HE (Zeng)", 
                 "Hematopoeitic (Zeng)", "TWIST+ Mesoderm (Zeng)", 
                 "DLK+ Mesoderm (Zeng)", "FBLN5+ Mesoderm (Zeng)", 
                 "WT1+ Mesoderm (Zeng)", "Epithelium (Zeng)")

cds_save_directory <- "~/save_dir_for_cds_output/"

cds_scores <- multi_group_scoring(cds = cds_26_split, 
                                  gene_groups = gene_groups,
                                  group_names = group_names,
                                  cds_save_directory = cds_save_directory)

#Now separate the hematopoiteic-related populations from the contaminating ones
cds_hemato <- choose_cells(cds_scores) #upper populations
cds_other <- choose_cells(cds_scores) #lower populations

#Get score column names
score_names <- names(cds_scores@colData)[-c(1:4)]

score_names_hemato <- score_names[c(1:7, 16:22)]
score_names_other <- score_names[c(8:15, 23:27)]

#Create heatmap matrix for hematopoietic-related populations
mat_hemato <- avg_score_matrix(cds_hemato, "new_cluster", score_names_hemato)
mat_hemato <- mat_hemato[sort(rownames(mat_hemato)),]

#Select hematopoietic-related score rownames and reformat them
mat_hemato_rownames <- rownames(mat_hemato)
new_hemato_rownames <- gsub("..", " ", mat_hemato_rownames, fixed = TRUE)
new_hemato_rownames <- gsub(".", " ", new_hemato_rownames, fixed = TRUE)
new_hemato_rownames <- gsub("Pre HE", "Pre-HE", new_hemato_rownames, fixed = TRUE)
new_hemato_rownames <- gsub("Monocyte Macrophage", "Monocyte/Mac", new_hemato_rownames, fixed = TRUE)
new_hemato_rownames <- gsub(" Score", "", new_hemato_rownames, fixed = TRUE)
new_hemato_rownames <- gsub("Calvanese", "(Calvanese)", new_hemato_rownames, fixed = TRUE)
new_hemato_rownames <- gsub("Zeng", "(Zeng)", new_hemato_rownames, fixed = TRUE)

#Create heatmap matrix for contaminating populations
mat_other <- avg_score_matrix(cds_other, "new_cluster", score_names_other)
mat_other <- mat_other[sort(rownames(mat_other)),]

#Select non-hematopoietic-related score rownames and reformat them
mat_other_rownames <- rownames(mat_other)
new_other_rownames <- gsub("..", " ", mat_other_rownames, fixed = TRUE)
new_other_rownames <- gsub(".", " ", new_other_rownames, fixed = TRUE)
new_other_rownames <- gsub(" Score", "", new_other_rownames, fixed = TRUE)
new_other_rownames <- gsub("Calvanese", "(Calvanese)", new_other_rownames, fixed = TRUE)
new_other_rownames <- gsub("Zeng", "(Zeng)", new_other_rownames, fixed = TRUE)
new_other_rownames <- gsub("DLK", "DLK+", new_other_rownames, fixed = TRUE)
new_other_rownames <- gsub("FBLN5", "FBLN5+", new_other_rownames, fixed = TRUE)
new_other_rownames <- gsub("WT1", "WT1+", new_other_rownames, fixed = TRUE)

#Create heatmap for identification of hematopoietic-related populations
svg(filename = paste0(cds_save_directory, 
                      "Hematopoeitic-Related_Score_Heatmap.svg"),
    width = 10,
    height = 10)
my_pheatmap <- pheatmap(mat_hemato[mat_hemato_rownames,
                                   c(1:5, 7:12, 14, 16:22, 24:25, 
                                     "26a", "26b", 27:28)],
                        legend=T, 
                        show_rownames = T, 
                        show_colnames = T, 
                        cluster_cols=F, 
                        cluster_rows = F,
                        cellheight = 10,
                        cellwidth = 10,
                        labels_row = new_hemato_rownames) 
draw(my_pheatmap)
dev.off()

#Create heatmap for identification of other populations
svg(filename = paste0(cds_save_directory, 
                      "Other_Non-Hemato_Cells_Score_Heatmap.svg"),
    width = 10,
    height = 10)
my_pheatmap <- pheatmap(mat_other[mat_other_rownames,],
                        legend=T, 
                        show_rownames = T, 
                        show_colnames = T, 
                        cluster_cols=F, 
                        cluster_rows = F,
                        cellheight = 10,
                        cellwidth = 10,
                        labels_row = new_other_rownames) 
draw(my_pheatmap)
dev.off()

#Add cell type column to cds; manually called 
colData(cds_scores)$cell_type <- cds_scores$new_cluster
colData(cds_scores)$cell_type <- dplyr::recode(colData(cds_cell_types)$cell_type,
                                                   "1" = "Arterial Endothelium",
                                                   "2" = "Hematopoietic (Erythro-Myelo Precursor)",
                                                   "3" = "Arterial Endothelium",
                                                   "4" = "Venous Endothelium",
                                                   "5" = "Hematopoietic (Erythro-Myelo Precursor)",
                                                   "6" = "Unknown Mesoderm/Stroma",
                                                   "7" = "Venous Endothelium",
                                                   "8" = "Endothelium",
                                                   "9" = "Endothelium",
                                                   "10" = "Hematopoietic (Erythro-Myelo Precursor)",
                                                   "11" = "Endothelium",
                                                   "12" = "Endothelium",
                                                   "13" = "Kidney Fibroblast",
                                                   "14" = "Endothelium",
                                                   "15" = "Unknown Mesoderm/Stroma",
                                                   "16" = "Endothelium",
                                                   "17" = "Endothelium",
                                                   "18" = "Venous Endothelium",
                                                   "19" = "Endothelium",
                                                   "20" = "Venous Endothelium",
                                                   "21" = "Unknown",
                                                   "22" = "Endothelium",
                                                   "23" = "Liver Epithelium",
                                                   "24" = "Endothelium",
                                                   "25" = "Hematopoietic (Myeloid Precursor)",
                                                   "26a" = "Late HE",
                                                   "26b" = "Hematopoietic (Erythro-Myelo Precursor)",
                                                   "27" = "HE",
                                                   "28" = "Endothelium")