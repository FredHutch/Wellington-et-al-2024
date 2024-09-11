#Aggregate Gene Score Heatmaps for all genes associated with individual GO Terms

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
library(org.Hs.eg.db)

#Colors to be used in heatmaps
mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")

#Fuction for calculating the aggregate gene score for a set of markers
gene_group_scoring <- function(cds, gene_group, column_name) {
  cds_gene_group <- cds[fData(cds)$gene_short_name %in% gene_group,]
  aggregate_expression <- exprs(cds_gene_group)
  aggregate_expression <- Matrix::t(Matrix::t(aggregate_expression) / pData(cds_gene_group)$Size_Factor)
  aggregate_expression <- Matrix::colSums(aggregate_expression)
  pData(cds)[,column_name] <- log(aggregate_expression+1)
  return(cds)
}

#Function for generating a Heatmap matrix
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

#Load integrated data
cds <- readRDS("~/RDS_objects/integrated_data_Monocle_label-transfer.RDS")
colData(cds)$cluster <- clusters(cds) #add monocle cluster labels as a new column

#Generate cds subset containing only HSC-Competent Embryo HE (cluster 75) cells
he_vivo_cds <- cds[,colData(cds)$cluster %in% c("75") & 
                     colData(cds)$dataset %in% c("Human Embryo (Crosse et al)",
                                                 "Human Embryo, Week 4.5 (CS14)",
                                                 "Human Embryo, Week 5 (CS15)",
                                                 "Human Embryo, Week 6 (CS17)")]

#Generate cds subset containing only HSC-Competent Embryo Early HSPC (cluster 74) cells
hspc_vivo_cds <- cds[,colData(cds)$cluster %in% c("74") & 
                       colData(cds)$dataset %in% c("Human Embryo (Crosse et al)",
                                                   "Human Embryo, Week 4.5 (CS14)",
                                                   "Human Embryo, Week 5 (CS15)",
                                                   "Human Embryo, Week 6 (CS16)")]

#Generate cds subset containing only iPS HE (cluster 75) cells
he_vitro_cds <- cds[,colData(cds)$cluster %in% c("75") & 
                      colData(cds)$dataset %in% c("45-iPS Day 8")]

#Generate cds subset containing only iPS Early HSPC (cluster 74) cells
hspc_vitro_cds <- cds[,colData(cds)$cluster %in% c("74") & 
                      colData(cds)$dataset %in% c("45-iPS Day 8")]




######HE-HSPC Downregulated Genes######
hspc_down_go <- read.csv("~/csv_files/Metascape_GO_Analysis/iPS_v_HSC-Competent-Embryo_HSPC/hspc_down-in-vitro_GO.csv")
key_desc <- c("blood vessel development", "positive regulation of cell motility",
              "endothelium development", "regulation of angiogenesis", 
              "cell-cell adhesion", "vasculogenesis", "actin filament-based process",
              "response to wounding", "enzyme-linked receptor protein signaling pathway",
              "eye development", "heart valve morphogenesis", "embryonic morphogenesis",
              "regulation of cell-substrate adhesion", "response to growth factor")
hspc_down_go <- hspc_down_go %>% filter(Description %in% key_desc)
hspc_down_go <- hspc_down_go %>% select(GO, Description)

#create vector of unique GO terms
go_terms <- unique(hspc_down_go$GO)

#generate GO term variables containing list of genes associated with the GO term
#generate vector of column names to be used for storing aggregate gene scores for
#each GO term
go_cols <- c()
for (i in 1:length(go_terms)) {
  go_genes <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys=go_terms[i], columns=c("SYMBOL"))
  go_col <- paste0("GO_", i)
  assign(go_col, unique(go_genes$SYMBOL))
  go_cols <- c(go_cols, go_col)
}

#Add aggregate gene scores for each HE-HSPC Downregulated GO term to the integrated cds
cds_he_hspc <- cds
for (go_col in go_cols) {
  cds_he_hspc <- gene_group_scoring(cds_he_hspc, 
                                    get(go_col), 
                                    go_col)
}

#Subset cds for only iPS and HSC-Competent Embryo HE and HSPC cells
cds_he_hspc <- cds_he_hspc[,colData(cds_he_hspc)$cluster %in% c("75", "74")]
cds_he_hspc <- cds_he_hspc[,!colData(cds_he_hspc)$dataset %in% c("Human Embryo (Zeng et al)")]

#add designation of in vitro or in vivo cell type as a metadata column
colData(cds_he_hspc)$vivo_or_vitro <- colData(cds_he_hspc)$dataset
colData(cds_he_hspc)$vivo_or_vitro <- dplyr::recode(colData(cds_he_hspc)$vivo_or_vitro,
                                                    "Human Embryo (Crosse et al)" = ifelse(colData(cds_he_hspc)$cluster == 75, 
                                                                                           "vivo_HE", 
                                                                                           "vivo_HSPC"),
                                                    "Human Embryo, Week 4.5 (CS14)" = ifelse(colData(cds_he_hspc)$cluster == 75, 
                                                                                             "vivo_HE", 
                                                                                             "vivo_HSPC"),
                                                    "Human Embryo, Week 5 (CS15)" = ifelse(colData(cds_he_hspc)$cluster == 75, 
                                                                                           "vivo_HE", 
                                                                                           "vivo_HSPC"),
                                                    "Human Embryo, Week 6 (CS16)" = ifelse(colData(cds_he_hspc)$cluster == 75, 
                                                                                           "vivo_HE", 
                                                                                           "vivo_HSPC"),
                                                    "45-iPS Day 8" = ifelse(colData(cds_he_hspc)$cluster == 75, 
                                                                            "vitro_HE", 
                                                                            "vitro_HSPC"))



score_cols <- colnames(colData(cds_he_hspc))[c(20:33)]
mat_he_hspc <- avg_score_matrix(cds_he_hspc, "vivo_or_vitro", score_cols)

heatmap_rows <- unique(meta$Description)
rownames(mat_he_hspc) <- heatmap_rows

heatmap_cols <- c("in vitro HE", "in vitro HSPC", "in vivo HE", "in vivo HSPC")
colnames(mat_he_hspc) <- heatmap_cols

save_dir <- r"(C:\Users\rachw\Documents\Paper1_Doulatov\GO_Gene_Heatmaps\)"
save_dir <- gsub("\\\\", "/", save_dir)
svglite::svglite(filename = paste0(save_dir, "GO_Gene_Rel_Exp_Heatmap_HE-HSPC-down.svg"))
my_pheatmap <- pheatmap(mat_he_hspc[row.names(mat_he_hspc),
                                    colnames(mat_he_hspc)], 
                        legend=T,
                        show_rownames = T, 
                        show_colnames = T, 
                        cluster_cols=F, 
                        cluster_rows = F, 
                        cellheight = 10, 
                        cellwidth = 10) 
draw(my_pheatmap)
dev.off()


#GO Term Aggregate Gene Scores for AE-HE Up/Down and HE-HSPC up

meta_ae_up <- r"(C:\Users\rachw\Documents\Paper1_Doulatov\Metascape_Data\AE-HE\batch_metascape\up_genes\all.t510z5jo8\Enrichment_GO\GO_AllLists.csv)"
meta_ae_up <- gsub("\\\\", "/", meta_ae_up)
meta_ae_up <- read.csv(meta_ae_up)

meta_ae_down <- r"(C:\Users\rachw\Documents\Paper1_Doulatov\Metascape_Data\AE-HE\batch_metascape\down_genes\all.tgh4ghi73\Enrichment_GO\GO_AllLists.csv)"
meta_ae_down<- gsub("\\\\", "/", meta_ae_down)
meta_ae_down <- read.csv(meta_ae_down)

meta_he_up <- r"(C:\Users\rachw\Documents\Paper1_Doulatov\Metascape_Data\HE-HSPC\batch_metascape\up_genes\all.to5dctdru\Enrichment_GO\GO_AllLists.csv)"
meta_he_up <- gsub("\\\\", "/", meta_he_up)
meta_he_up <- read.csv(meta_he_up)

key_ae_up <- c("ribosome assembly", "ribosomal large subunit biogenesis",
               "positive regulation of signal transduction by p53 class mediator",
               "carbon dioxide transport", "regulation of proteolysis","protein folding", 
               "regulation of cell activation", "maturation of SSU-rRNA",
               "positive regulation of immune effector process", "regulation of translation",
               "negative regulation of RNA splicing", "maturation of LSU-rRNA") #12
key_ae_down <- c("blood vessel development", "positive regulation of cell motility",
                 "response to wounding", "cellular response to growth factor stimulus",
                 "endothelium development", "tissue morphogenesis", "heart morphogenesis",
                 "cell junction organization", "actin filament-based process",
                 "embryonic morphogenesis", "positive regulation of cell adhesion",
                 "skeletal system development", "multicellular organismal-level homeostasis",
                 "cell population proliferation") #14
key_he_up <- c("regulation of leukocyte degranulation", "positive regulation of leukocyte migration",
               "regulation of immune effector process", 
               "regulation of intrinsic apoptotic signaling pathway by p53 class mediator",
               "regulation of endopeptidase activity",
               "positive regulation of response to external stimulus","liver development",
               "nucleoside triphosphate metabolic process") #8

meta_ae_up <- meta_ae_up %>% filter(Description %in% key_ae_up)
meta_ae_up <- meta_ae_up %>% dplyr::select(GO, Description)
go_terms_ae_up <- unique(meta_ae_up$GO)

meta_ae_down <- meta_ae_down %>% filter(Description %in% key_ae_down)
meta_ae_down <- meta_ae_down %>% dplyr::select(GO, Description)
go_terms_ae_down <- unique(meta_ae_down$GO)

meta_he_up <- meta_he_up %>% filter(Description %in% key_he_up)
meta_he_up <- meta_he_up %>% dplyr::select(GO, Description)
go_terms_he_up <- unique(meta_he_up$GO)


get_go_gene_symbols <- function (cds, organism, go_term_list) {
  cds_copy <- cds
  if (organism == "human") {
    library(org.Hs.eg.db)
    lib <- org.Hs.eg.db
  } else if (organism == "mouse") {
    library(org.Mm.eg.db)
    lib <- org.Mm.eg.db
  }
  go_cols <- c()
  for (i in 1:length(go_term_list)) {
    go_genes <- AnnotationDbi::select(lib, 
                                      keytype="GOALL", 
                                      keys=go_term_list[i], 
                                      columns=c("SYMBOL"))
    go_col <- paste0("GO_", i)
    assign(go_col, unique(go_genes$SYMBOL))
    go_cols <- c(go_cols, go_col)
  }
  for (go_col in go_cols) {
    cds_copy <- gene_group_scoring(cds_copy, get(go_col), go_col)
  }
  return(cds_copy)
}

cds_ae_up <- get_go_gene_symbols(cds, 
                                 organism = "human", 
                                 go_term_list = go_terms_ae_up)
cds_ae_down <- get_go_gene_symbols(cds, 
                                   organism = "human", 
                                   go_term_list = go_terms_ae_down)
cds_he_up <- get_go_gene_symbols(cds, 
                                 organism = "human", 
                                 go_term_list = go_terms_he_up)


colData(cds_ae_up)$cluster <- clusters(cds_ae_up)
cds_ae_up <- cds_ae_up[,colData(cds_ae_up)$cluster %in% c("36", "75") & !colData(cds_ae_up)$dataset %in% c("Human Embryo (Zeng et al)")]
colData(cds_ae_up)$vivo_or_vitro <- colData(cds_ae_up)$dataset
colData(cds_ae_up)$vivo_or_vitro <- dplyr::recode(colData(cds_ae_up)$vivo_or_vitro,
                                                  "Human Embryo (Crosse et al)" = ifelse(colData(cds_ae_up)$cluster == 75, "vivo_HE", "vivo_AE"),
                                                  "Human Embryo, Week 4.5 (CS14)" = ifelse(colData(cds_ae_up)$cluster == 75, "vivo_HE", "vivo_AE"),
                                                  "Human Embryo, Week 5 (CS15)" = ifelse(colData(cds_ae_up)$cluster == 75, "vivo_HE", "vivo_AE"),
                                                  "Human Embryo, Week 6 (CS16)" = ifelse(colData(cds_ae_up)$cluster == 75, "vivo_HE", "vivo_AE"),
                                                  "45-iPS Day 8" = ifelse(colData(cds_ae_up)$cluster == 75, "vitro_HE", "vitro_AE"))

colData(cds_ae_down)$cluster <- clusters(cds_ae_down)
cds_ae_down <- cds_ae_down[,colData(cds_ae_down)$cluster %in% c("36", "75") & !colData(cds_ae_down)$dataset %in% c("Human Embryo (Zeng et al)")]
colData(cds_ae_down)$vivo_or_vitro <- colData(cds_ae_down)$dataset
colData(cds_ae_down)$vivo_or_vitro <- dplyr::recode(colData(cds_ae_down)$vivo_or_vitro,
                                                    "Human Embryo (Crosse et al)" = ifelse(colData(cds_ae_down)$cluster == 75, "vivo_HE", "vivo_AE"),
                                                    "Human Embryo, Week 4.5 (CS14)" = ifelse(colData(cds_ae_down)$cluster == 75, "vivo_HE", "vivo_AE"),
                                                    "Human Embryo, Week 5 (CS15)" = ifelse(colData(cds_ae_down)$cluster == 75, "vivo_HE", "vivo_AE"),
                                                    "Human Embryo, Week 6 (CS16)" = ifelse(colData(cds_ae_down)$cluster == 75, "vivo_HE", "vivo_AE"),
                                                    "45-iPS Day 8" = ifelse(colData(cds_ae_down)$cluster == 75, "vitro_HE", "vitro_AE"))

colData(cds_he_up)$cluster <- clusters(cds_he_up)
cds_he_up <- cds_he_up[,colData(cds_he_up)$cluster %in% c("75", "74") & !colData(cds_he_up)$dataset %in% c("Human Embryo (Zeng et al)")]
colData(cds_he_up)$vivo_or_vitro <- colData(cds_he_up)$dataset
colData(cds_he_up)$vivo_or_vitro <- dplyr::recode(colData(cds_he_up)$vivo_or_vitro,
                                                  "Human Embryo (Crosse et al)" = ifelse(colData(cds_he_up)$cluster == 75, "vivo_HE", "vivo_HSPC"),
                                                  "Human Embryo, Week 4.5 (CS14)" = ifelse(colData(cds_he_up)$cluster == 75, "vivo_HE", "vivo_HSPC"),
                                                  "Human Embryo, Week 5 (CS15)" = ifelse(colData(cds_he_up)$cluster == 75, "vivo_HE", "vivo_HSPC"),
                                                  "Human Embryo, Week 6 (CS16)" = ifelse(colData(cds_he_up)$cluster == 75, "vivo_HE", "vivo_HSPC"),
                                                  "45-iPS Day 8" = ifelse(colData(cds_he_up)$cluster == 75, "vitro_HE", "vitro_HSPC"))

save_dir <- r"(C:\Users\rachw\Documents\Paper1_Doulatov\GO_Gene_Heatmaps\)"
save_dir <- gsub("\\\\", "/", save_dir)

score_ae_up <- colnames(colData(cds_ae_up))[c(20:31)]
mat_ae_up <- avg_score_matrix(cds_ae_up, "vivo_or_vitro", score_ae_up)
heatmap_rows_ae_up <- unique(meta_ae_up$Description)
rownames(mat_ae_up) <- heatmap_rows_ae_up
heatmap_cols_ae_up <- gsub("_", " ", colnames(mat_ae_up))
heatmap_cols_ae_up <- paste0("in ", heatmap_cols_ae_up)
colnames(mat_ae_up) <- heatmap_cols_ae_up
svglite::svglite(filename = paste0(save_dir, "GO_Gene_Rel_Exp_Heatmap_AE-HE-up.svg"))
my_pheatmap <- pheatmap(mat_ae_up[row.names(mat_ae_up),
                                  colnames(mat_ae_up)], 
                        legend=T,
                        show_rownames = T, 
                        show_colnames = T, 
                        cluster_cols=F, 
                        cluster_rows = F, 
                        cellheight = 10, 
                        cellwidth = 10) 
draw(my_pheatmap)
dev.off()


score_ae_down <- colnames(colData(cds_ae_down))[c(20:33)]
mat_ae_down <- avg_score_matrix(cds_ae_down, "vivo_or_vitro", score_ae_down)
heatmap_rows_ae_down <- unique(meta_ae_down$Description)
rownames(mat_ae_down) <- heatmap_rows_ae_down
heatmap_cols_ae_down <- gsub("_", " ", colnames(mat_ae_down))
heatmap_cols_ae_down <- paste0("in ", heatmap_cols_ae_down)
colnames(mat_ae_down) <- heatmap_cols_ae_down
svglite::svglite(filename = paste0(save_dir, "GO_Gene_Rel_Exp_Heatmap_AE-HE-down.svg"))
my_pheatmap <- pheatmap(mat_ae_down[row.names(mat_ae_down),
                                    colnames(mat_ae_down)], 
                        legend=T,
                        show_rownames = T, 
                        show_colnames = T, 
                        cluster_cols=F, 
                        cluster_rows = F, 
                        cellheight = 10, 
                        cellwidth = 10) 
draw(my_pheatmap)
dev.off()


score_he_up <- colnames(colData(cds_he_up))[c(20:27)]
mat_he_up <- avg_score_matrix(cds_he_up, "vivo_or_vitro", score_he_up)
heatmap_rows_he_up <- unique(meta_he_up$Description)
rownames(mat_he_up) <- heatmap_rows_he_up
heatmap_cols_he_up <- gsub("_", " ", colnames(mat_he_up))
heatmap_cols_he_up <- paste0("in ", heatmap_cols_he_up)
colnames(mat_he_up) <- heatmap_cols_he_up
svglite::svglite(filename = paste0(save_dir, "GO_Gene_Rel_Exp_Heatmap_HE-HSPC-up.svg"))
my_pheatmap <- pheatmap(mat_he_up[row.names(mat_he_up),
                                  colnames(mat_he_up)], 
                        legend=T,
                        show_rownames = T, 
                        show_colnames = T, 
                        cluster_cols=F, 
                        cluster_rows = F, 
                        cellheight = 10, 
                        cellwidth = 10) 
draw(my_pheatmap)
dev.off()