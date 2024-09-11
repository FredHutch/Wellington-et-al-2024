#Preprocessing and Integration of iPS Day 8 10X data with published data from:
#CS14-17, Calvanese et al.: https://www.nature.com/articles/s41586-022-04571-x
#CS16, Crosse et al.: https://www.cell.com/cell-stem-cell/pdfExtended/S1934-5909(20)30400-8
#CS10-13, Zeng et al.: https://www.nature.com/articles/s41422-019-0228-6

library(stringr)
library(dplyr)
library(biomaRt)
library(Seurat)
library(tidyverse)
library(patchwork)

###Preprocessing of Calvanese et al. 2022 data from GEO csv###

#Convert CSV to usable 10X data
Aorta_4wk <- read.csv("~/csv_files/Downloaded_Calvanese_et_al_2022_GEO/GSM4968831_Aorta-4wk-658.csv.gz", 
                      row.names = 1, 
                      header= TRUE) #CS14
Aorta_5wk_1 <- read.csv("~/csv_files/Downloaded_Calvanese_et_al_2022_GEO/GSM4968832_Aorta-5wk-555.csv.gz", 
                        row.names = 1, 
                        header= TRUE) #CS15
Aorta_5wk_2 <- read.csv("~/csv_files/Downloaded_Calvanese_et_al_2022_GEO/GSM4968833_Aorta-5wk-575.csv.gz", 
                        row.names = 1, 
                        header= TRUE) #CS15
Aorta_6wk <- read.csv("~/csv_files/Downloaded_Calvanese_et_al_2022_GEO/GSM4968834_Aorta-6wk-563.csv.gz", 
                      row.names = 1, 
                      header= TRUE) #CS17

#Replace "." in cell names with "_"
cells_4wk <- str_replace(colnames(Aorta_4wk), ".1", "-3")
colnames(Aorta_4wk) <- cells_4wk
cells_5wk_1 <- str_replace(colnames(Aorta_5wk_1), ".1", "-1")
colnames(Aorta_5wk_1) <- cells_5wk_1
cells_5wk_2 <- str_replace(colnames(Aorta_5wk_2), ".1", "-2")
colnames(Aorta_5wk_2) <- cells_5wk_2

#Change rownames from ENSEMBL IDs to short gene names using biomaRt
#Note: must use the following if you want o use the later mitochondrial count code

#The following function takes a human RNA count table with cells as columns and 
#EMSEMBL gene names as rows. It returns a new RNA count table with the rownames 
#replaced with short gene names instead of ENSEMBL gene names. Note that the row 
#names will stay ENSEMBL IDs if there is no short gene names found by biomaRt or 
#if there are multiple ENSEMBL IDs that map to the same short gene name (this is 
#for clarity purposes).
ENSEMBL_rownames_to_symbol <- function(ENSEMBL_count_table) {
  library(dplyr)
  library(biomaRt)
  ENSEMBL_count_table_copy <- ENSEMBL_count_table
  sgn <- c()
  mart <- useMart(biomart = "ensembl",host = "https://uswest.ensembl.org") #include host argument to avoid access issues
  mart <- useDataset("hsapiens_gene_ensembl", mart = mart)
  genes.table <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"), values= rownames(ENSEMBL_count_table_copy), mart= mart)
  for(gene in rownames(ENSEMBL_count_table_copy)){
    if (gene %in% genes.table$ensembl_gene_id){
      symbol <- genes.table$hgnc_symbol[match(gene, genes.table$ensembl_gene_id)]
      if (symbol != ""){
        sgn <- c(sgn, symbol)}
      else(sgn <- c(sgn, gene))} #must have this line or it will skip entries; can use length() to verify list
    else(sgn <- c(sgn, gene))}
  sgn_dup <- unique(sgn[duplicated(sgn)])
  sgn_update <- c()
  for(gene_name_index in 1:length(sgn)){
    if(sgn[gene_name_index] %in% sgn_dup){
      sgn_update <- c(sgn_update, rownames(ENSEMBL_count_table_copy)[gene_name_index])
    }
    else(sgn_update <- c(sgn_update, sgn[gene_name_index]))}
  row.names(ENSEMBL_count_table_copy) <- sgn_update
  return(ENSEMBL_count_table_copy)}

#Can apply the above function for the imported Calvanese et al. 2022 datasets
Aorta_4wk_sgn <- ENSEMBL_rownames_to_symbol(Aorta_4wk)
Aorta_5wk_1_sgn <- ENSEMBL_rownames_to_symbol(Aorta_5wk_1)
Aorta_5wk_2_sgn <- ENSEMBL_rownames_to_symbol(Aorta_5wk_2)
Aorta_6wk_sgn <- ENSEMBL_rownames_to_symbol(Aorta_6wk)

#Create Seurat objects
wk4.obj <- CreateSeuratObject(counts = Aorta_4wk_sgn, 
                              project = "Calvanese_Aorta_Wk4.5", 
                              min.cells = 3, 
                              min.features = 200)
wk5_1.obj <- CreateSeuratObject(counts = Aorta_5wk_1_sgn, 
                                project = "Calvanese_Aorta_Wk5", 
                                min.cells = 3, 
                                min.features = 200)
wk5_2.obj <- CreateSeuratObject(counts = Aorta_5wk_2_sgn, 
                                project = "Calvanese_Aorta_Wk5", 
                                min.cells = 3, 
                                min.features = 200)
wk6.obj <- CreateSeuratObject(counts = Aorta_6wk_sgn, 
                              project = "Calvanese_Aorta_Wk6", 
                              min.cells = 3, 
                              min.features = 200)

#Add Metadata column containing week only information
wk4.obj <- AddMetaData(object=wk4.obj, 
                       metadata="Human Embryo, Week 4.5 (CS14)", 
                       col.name= "dataset")
wk5_1.obj <- AddMetaData(object=wk5_1.obj, 
                         metadata="Human Embryo, Week 5 (CS15)", 
                         col.name= "dataset")
wk5_2.obj <- AddMetaData(object=wk5_2.obj, 
                         metadata="Human Embryo, Week 5 (CS15)", 
                         col.name= "dataset")
wk6.obj <- AddMetaData(object=wk6.obj, 
                       metadata="Human Embryo, Week 6 (CS16)", 
                       col.name= "dataset")

wk4.obj <- AddMetaData(object=wk4.obj, 
                       metadata="Human Embryo, Week 4.5 (CS14)", 
                       col.name= "crosse_dataset")
wk5_1.obj <- AddMetaData(object=wk5_1.obj, 
                         metadata="Human Embryo, Week 5 (CS15)", 
                         col.name= "crosse_dataset")
wk5_2.obj <- AddMetaData(object=wk5_2.obj, 
                         metadata="Human Embryo, Week 5.5 (CS15)", 
                         col.name= "crosse_dataset")
wk6.obj <- AddMetaData(object=wk6.obj, 
                       metadata="Human Embryo, Week 6 (CS16)", 
                       col.name= "crosse_dataset")

#Add Calvanese et al. 2022 AGM numbering to metadata
wk4.obj <- AddMetaData(object=wk4.obj, 
                       metadata="AGM658", 
                       col.name= "AGM_Sample")
wk5_1.obj <- AddMetaData(object=wk5_1.obj,
                         metadata="AGM555", 
                         col.name= "AGM_Sample")
wk5_2.obj <- AddMetaData(object=wk5_2.obj, 
                         metadata="AGM575", 
                         col.name= "AGM_Sample")
wk6.obj <- AddMetaData(object=wk6.obj, 
                       metadata="AGM563", 
                       col.name= "AGM_Sample")

#Add Calvanese et al. 2022 cluster labels
#To load Mikkola data as downloaded from google drive (find link in Calvanese paper)
load("~/seurat_object.Rdata") #creates Seurat object called "sample"
mikkola.obj <- sample #stores loaded info in a renamed variable
rm(sample) #removes the originally named variable

#Labels only includes AGM555, AGM575, and AGM658 (wk 4.5, 5, 5.5)
mikkola.wk4_clust <- subset(mikkola.obj[[]], 
                            orig.ident == "AGM658", 
                            select=c("seurat_clusters"))
mikkola.wk5_1_clust <- subset(mikkola.obj[[]], 
                              orig.ident == "AGM555", 
                              select=c("seurat_clusters"))
mikkola.wk5_2_clust <- subset(mikkola.obj[[]], 
                              orig.ident == "AGM575", 
                              select=c("seurat_clusters"))

#Add names based on cluster classification by Calvanese et al. 2022 as well
Endo <- c(0, 1, 2, 3)
HSC <- c(4, 5)
Non_HSC <- c(6, 7, 8)

#Week 4 cell type identification
mikkola.wk4_cells <- mikkola.wk4_clust %>% 
  mutate(Calvanese_Cell_Type = case_when((mikkola.wk4_clust$seurat_clusters %in% Endo) ~ "Endothelium", 
                                                                                  (mikkola.wk4_clust$seurat_clusters %in% HSC) ~ "HSC", 
                                                                                  (mikkola.wk4_clust$seurat_clusters %in% Non_HSC) ~ "Non-HSC"))
mikkola.wk4_cells <- mikkola.wk4_cells[2]

#Week 5 cell type identification
mikkola.wk5_1_cells <- mikkola.wk5_1_clust %>% 
  mutate(Calvanese_Cell_Type = case_when((mikkola.wk5_1_clust$seurat_clusters %in% Endo) ~ "Endothelium", 
                                         (mikkola.wk5_1_clust$seurat_clusters %in% HSC) ~ "HSC", 
                                         (mikkola.wk5_1_clust$seurat_clusters %in% Non_HSC) ~ "Non-HSC"))
mikkola.wk5_1_cells <- mikkola.wk5_1_cells[2]

#Week 5.5 cell type identification
mikkola.wk5_2_cells <- mikkola.wk5_2_clust %>% 
  mutate(Calvanese_Cell_Type = case_when((mikkola.wk5_2_clust$seurat_clusters %in% Endo) ~ "Endothelium", 
                                         (mikkola.wk5_2_clust$seurat_clusters %in% HSC) ~ "HSC", 
                                         (mikkola.wk5_2_clust$seurat_clusters %in% Non_HSC) ~ "Non-HSC"))
mikkola.wk5_2_cells <- mikkola.wk5_2_cells[2]

#Addition of the metadata to the Seurat objects
wk4.obj <- AddMetaData(object=wk4.obj, 
                       metadata = mikkola.wk4_clust, 
                       col.name = "Calvanese_Cluster")
wk4.obj <- AddMetaData(object=wk4.obj, 
                       metadata = mikkola.wk4_cells, 
                       col.name = "Calvanese_Cell_Type")
wk5_1.obj <- AddMetaData(object=wk5_1.obj, 
                         metadata = mikkola.wk5_1_clust, 
                         col.name = "Calvanese_Cluster")
wk5_1.obj <- AddMetaData(object=wk5_1.obj, 
                         metadata = mikkola.wk5_1_cells, 
                         col.name = "Calvanese_Cell_Type")
wk5_2.obj <- AddMetaData(object=wk5_2.obj, 
                         metadata = mikkola.wk5_2_clust, 
                         col.name = "Calvanese_Cluster")
wk5_2.obj <- AddMetaData(object=wk5_2.obj, 
                         metadata = mikkola.wk5_2_cells, 
                         col.name = "Calvanese_Cell_Type")


##########Integration of iPS Day 8 10X, Zeng et al 2019, Crosse et al 2020 and 
#Calvanese et al 2022#########

#Load reference human embryo data from Zeng et al. 2019
CS10 <- Read10X(data.dir="~/cellranger/Zeng_et_al/CS10")
CS11 <- Read10X(data.dir="~/cellranger/Zeng_et_al/CS11")
CS13 <- Read10X(data.dir="~/cellranger/Zeng_et_al/CS13")

#Load reference human embryo data from Crosse et al. 2020
CS16_all <- Read10X("~/cellranger/Crosse_et_al/CS16_AoMid/", 
                    strip.suffix = TRUE)
CS16_34 <- Read10X ("~/cellranger/Crosse_et_al/CS16_AoV/", 
                    strip.suffix = TRUE)

#Load iPS Day 8 10X data
iPS_D8 <- Read10X("~/cellranger/iPS_10X/outs/filtered_feeature_bc_matrix")

#Create Seurat Objects
CS10.obj <- CreateSeuratObject(counts = CS10, 
                               project = "CS10", 
                               min.cells = 3, 
                               min.features = 200)
CS11.obj <- CreateSeuratObject(counts = CS11, 
                               project = "CS11", 
                               min.cells = 3, 
                               min.features = 200)
CS13.obj <- CreateSeuratObject(counts = CS13, 
                               project = "CS13", 
                               min.cells = 3, 
                               min.features = 200)
iPS_D8.obj <- CreateSeuratObject(counts = iPS_D8, 
                                 project = "45-iPS", 
                                 min.cells = 3, 
                                 min.features = 200)
CS16_all.obj <- CreateSeuratObject(counts = CS16_all, 
                                   project = "CS16, AoMid, no sort", 
                                   min.cells = 3, 
                                   min.features = 200)
CS16_34.obj <- CreateSeuratObject(counts = CS16_34, 
                                  project = "CS16, AoV, 34+ sorted", 
                                  min.cells = 3, 
                                  min.features = 200)
calvanese_all.obj <- merge(wk4.obj, 
                           y = c(wk5_1.obj, wk5_2.obj, wk6.obj), 
                           add.cell.ids = c("wk4", "wk5_1", "wk5_2", "wk6"), 
                           project = "Calvanese_Human_Embryo")

#Add Zeng et al. labels to each of the objects so will be included in final integrated object
#Note that N/A will be added for all cells in the iPS and Crosse datasets
labels <- read.csv("~/marker_files/Zeng-et-al_CS10-13_cellidents.csv")
labels.df <- data.frame(labels)
labels_10 <- labels.df %>% filter(stage == "CS10")
labels_11 <- labels.df %>% filter(stage == "CS11")
labels_13 <- labels.df %>% filter(stage == "CS13")
new_cells_10 <- sub("[^_]*_[^_]*_(*)", "\\1", labels_10$cell_index) #returns list of cells without CS##_[a-z]_
new_cells_11 <- sub("[^_]*_[^_]*_(*)", "\\1", labels_11$cell_index) 
new_cells_13 <- sub("[^_]*_[^_]*_(*)", "\\1", labels_13$cell_index) 
labels_mod10 <- labels_10["cluster"]
labels_mod11 <- labels_11["cluster"]
labels_mod13 <- labels_13["cluster"]
row.names(labels_mod10) <- new_cells_10
row.names(labels_mod11) <- new_cells_11
row.names(labels_mod13) <- new_cells_13

#Add Crosse et al. labels (acquired from Edie Crosse) to the objects so they will 
#be included in the final integrated object
crosse_ids <- read.csv("marker_files/Crosse-et-al_Cluster_IDs.csv")
labels_ventral1 <- crosse_ids[str_detect(crosse_ids$index, 
                                         "Human10XVentral1"),]
labels_ventral2 <- crosse_ids[str_detect(crosse_ids$index, 
                                         "Human10XVentral2"),]
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
                                   "-1_5")
row.names(labels_all_nodup) <- c(labels_all_nodup$barcode)
labels_all_nodup <- subset(labels_all_nodup, 
                           select=-(barcode))

#Note: do not need to add to the iPS data because column not present; will auto-add 
#NA for this column when integrated

#Add dataset metadata column to Zeng et al. objects
CS10.obj <- AddMetaData(object=CS10.obj, 
                        metadata=labels_mod10, 
                        col.name= "Zeng_Label")
CS11.obj <- AddMetaData(object=CS11.obj, 
                        metadata=labels_mod11, 
                        col.name= "Zeng_Label")
CS13.obj <- AddMetaData(object=CS13.obj, 
                        metadata=labels_mod13, 
                        col.name= "Zeng_Label")

CS10.obj <- AddMetaData(object=CS10.obj, 
                        metadata="Human Embryo (Zeng et al)", 
                        col.name= "dataset")
CS11.obj <- AddMetaData(object=CS11.obj, 
                        metadata="Human Embryo (Zeng et al)", 
                        col.name= "dataset")
CS13.obj <- AddMetaData(object=CS13.obj, 
                        metadata="Human Embryo (Zeng et al)", 
                        col.name= "dataset")

CS10.obj <- AddMetaData(object=CS10.obj, 
                        metadata="CS10", 
                        col.name= "crosse_dataset")
CS11.obj <- AddMetaData(object=CS11.obj, 
                        metadata="CS11", 
                        col.name= "crosse_dataset")
CS13.obj <- AddMetaData(object=CS13.obj, 
                        metadata="CS13", 
                        col.name= "crosse_dataset")

#Add dataset metadata column to Crosse et al objects
CS16_all.obj <- AddMetaData(object=CS16_all.obj, 
                            metadata="Human Embryo (Crosse et al)", 
                            col.name= "dataset")
CS16_34.obj <- AddMetaData(object=CS16_34.obj, 
                           metadata="Human Embryo (Crosse et al)", 
                           col.name= "dataset")

#Edit the iPS dataset so "orig.ident" reflects the assay
iPS_D8.obj[["orig.ident"]] <- "45-iPS Day 8"
iPS_D8.obj <- AddMetaData(object=iPS_D8.obj, 
                            metadata="45-iPS Day 8", 
                            col.name= "dataset")
iPS_D8.obj <- AddMetaData(object=iPS_D8.obj, 
                            metadata="45-iPS Day 8", 
                            col.name= "crosse_dataset")

#Add percent mitochondrial reads to each object
CS10.obj[["percent.mt"]] <- PercentageFeatureSet(CS10.obj, 
                                                 pattern = "^MT-")
CS11.obj[["percent.mt"]] <- PercentageFeatureSet(CS11.obj, 
                                                 pattern = "^MT-")
CS13.obj[["percent.mt"]] <- PercentageFeatureSet(CS13.obj, 
                                                 pattern = "^MT-")
iPS_D8.obj[["percent.mt"]] <- PercentageFeatureSet(iPS_D8.obj, 
                                                   pattern = "^MT-")
CS16_34.obj[["percent.mt"]] <- PercentageFeatureSet(CS16_34.obj, 
                                                    pattern = "^MT-")
CS16_all.obj[["percent.mt"]] <- PercentageFeatureSet(CS16_all.obj, 
                                                     pattern = "^MT-")
calvanese_all.obj[["percent.mt"]] <- PercentageFeatureSet(calvanese_all.obj, 
                                                          pattern = "^MT-")

#Filter each dataset based on the individual cutoffs needed
CS10.obj <- subset(CS10.obj, 
                   subset = nFeature_RNA >100 & nFeature_RNA < 7000 & percent.mt < 20)
CS11.obj <- subset(CS11.obj, 
                   subset = nFeature_RNA >100 & nFeature_RNA < 9000 & percent.mt < 10)
CS13.obj<- subset(CS13.obj, 
                  subset = nFeature_RNA >100 & nFeature_RNA < 7000 & percent.mt < 40)
iPS45_D8.obj <- subset(iPS45_D8.obj, 
                       subset = nFeature_RNA >100 & nFeature_RNA < 100000 & percent.mt < 20)
CS16_34.obj <- subset(CS16_34.obj, 
                      subset = nFeature_RNA >100 & nFeature_RNA < 100000 & percent.mt < 10)
CS16_all.obj <- subset(CS16_all.obj, 
                       subset = nFeature_RNA >100 & nFeature_RNA < 100000 & percent.mt < 10)
calvanese_all.obj <- subset(calvanese_all.obj, 
                            subset = nFeature_RNA >500 & percent.mt < 5)

#Normalize each individual datasets and set FindVariableFeatures to the maximum
#number of features for each dataset (limiting to less features leads to poor
#integration)
CS10.obj <- NormalizeData(CS10.obj)
CS10.mod <- FindVariableFeatures(CS10.obj, 
                                 selection.method = "vst", 
                                 nfeatures = 22880)

CS11.obj <- NormalizeData(CS11.obj)
CS11.mod <- FindVariableFeatures(CS11.obj, 
                                 selection.method = "vst", 
                                 nfeatures = 22362)

CS13.obj <- NormalizeData(CS13.obj)
CS13.mod <- FindVariableFeatures(CS13.obj, 
                                 selection.method = "vst", 
                                 nfeatures = 17301)

iPS45_D8.obj <- NormalizeData(iPS45_D8.obj)
iPS45_D8.mod <- FindVariableFeatures(iPS45_D8.obj, 
                                     selection.method = "vst", 
                                     nfeatures = 24305)

CS16_34.obj <- NormalizeData(CS16_34.obj)
CS16_34.mod <- FindVariableFeatures(CS16_34.obj, 
                                    selection.method = "vst", 
                                    nfeatures = 17375)

CS16_all.obj <- NormalizeData(CS16_all.obj)
CS16_all.mod <- FindVariableFeatures(CS16_all.obj, 
                                     selection.method = "vst", 
                                     nfeatures = 15116)

calvanese_all.obj <- NormalizeData(calvanese_all.obj)
calvanese_all.mod <- FindVariableFeatures(calvanese_all.obj, 
                                          selection.method = "vst", 
                                          nfeatures = 22923)

#Create a list of all of the Seurat objects
obj.list <- c(CS10.mod, CS11.mod, CS13.mod, iPS45_D8.mod, CS16_34.mod, 
              CS16_all.mod, calvanese_all.mod)

#Now merge using CCA
features <- SelectIntegrationFeatures(object.list = obj.list)
#Note, if you set the nfeatures input to something high, it will return a list 
#of all of the common features between the lists, which in this case is supposedly 
#13,567 of them 
#16,966 if the Crosse data is removed
#without setting nfeatures (highly recommend) will use only 2000 top variable 
#genes and leads to poor integration
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, 
                                      anchor.features = features)
obj.combined <- IntegrateData(anchorset = obj.anchors)

#Now can do rest of analysis post-integration
DefaultAssay(obj.combined) <- "integrated"

#Preprocessing of integrated data
obj.combined.scale <- ScaleData(obj.combined, 
                                verbose = FALSE)
obj.combined.pca <- RunPCA(obj.combined.scale, 
                           npcs = 50, 
                           verbose = FALSE)
obj.combined.umap <- RunUMAP(obj.combined.pca, 
                             reduction = "pca", 
                             dims = 1:50) #dims refers to PCA
obj.combined.knn <- FindNeighbors(obj.combined.umap, 
                                  reduction = "umap", 
                                  dims = 1:2) #this matches the k.param to the n.neighbors in UMAP
obj.combined.clust <- FindClusters(obj.combined.knn, 
                                   resolution = 1.2) #can modify resolution up to 1.2

#Can view integrated dataset and clusters in UMAP
DimPlot(obj.combined.clust, 
        reduction = "umap", 
        label = TRUE, 
        repel = TRUE) + 
  NoLegend()