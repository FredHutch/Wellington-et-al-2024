#NicheNetR analysis for iPS v. HSC-Competent Human Embryo Comparisons

library(Seurat)
library(dplyr)
library(ggplot2)
library(monocle3)
library(nichenetr)
library(tidyr)
library(tidyverse)
library(reshape2)

#Load Seurat integrated data 
seu_obj <- readRDS("~/RDS_objects/integrated_data_Seurat.RDS")

#Add new column containing HSC-competent info
Idents(seu_obj) <- "dataset"
ident_df <- data.frame(Idents(seu_obj))
ident_df$Idents.seu_obj. <- ident_df$Idents.seu_obj. %>% 
  recode("Human Embryo (Zeng et al)" = "pre-HSC",
         "Human Embryo (Crosse et al)" = "HSC-competent",
         "Human Embryo, Week 4.5 (CS14)" = "HSC-competent",
         "Human Embryo, Week 5 (CS15)" = "HSC-competent",
         "Human Embryo, Week 6 (CS17)" = "HSC-competent",
         "45-iPS Day 8" = "iPS") 
ident_df <- ident_df %>% 
  rename("Idents.seu_obj." = "dataset_type")
seu_obj$dataset_type <- ident_df

#Create seurat subset containing only HSC-Competent and iPS
Idents(seu_obj) <- "dataset_type"
iPS_HSC_comp <- subset(seu_obj, idents = c("iPS","HSC-competent"))

#Create seurat object subsets for HE and EHT (based on Seurat clusters) for only
#iPS and HSC-competent embryo data
Idents(iPS_HSC_comp) <- "seurat_clusters"
HE <- subset(iPS_HSC_comp, idents = c("61"))
EHT <- subset(iPS_HSC_comp, idents = c("61", "59"))

#Load nichenetR ligand target matrices
ligand_target_matrix <- readRDS("~/RDS_objects/NicheNetR_inputs/ligand_target_matrix.rds")
lr_network <- readRDS("~/RDS_objects/NicheNetR_inputs/lr_network.rds")
weighted_networks <- readRDS("~/RDS_objects/NicheNetR_inputs/weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))

######Running NicheNetR to Identify Ligand-Receptors######
#Note: This was done with one of the earliest releases of NicheNetR

#For HE
nichenet_output_HE <- nichenet_seuratobj_aggregate(seurat_obj = HE, 
                                                receiver = "61", 
                                                condition_colname = "dataset_type", 
                                                condition_oi = "iPS", 
                                                condition_reference = "HSC-competent", 
                                                sender = "undefined", 
                                                ligand_target_matrix = ligand_target_matrix, 
                                                lr_network = lr_network, 
                                                weighted_networks = weighted_networks, 
                                                organism = "human", 
                                                assay_oi = "RNA") #saved as RDS


#For EHT
nichenet_output_EHT <- nichenet_seuratobj_aggregate(seurat_obj = EHT, 
                                                   receiver = c("61", "59"), 
                                                   condition_colname = "dataset_type", 
                                                   condition_oi = "iPS", 
                                                   condition_reference = "HSC-competent", 
                                                   sender = "undefined", 
                                                   ligand_target_matrix = ligand_target_matrix, 
                                                   lr_network = lr_network, 
                                                   weighted_networks = weighted_networks, 
                                                   organism = "human", 
                                                   assay_oi = "RNA") #saved as RDS