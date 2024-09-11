#Creation of the Metascape GO Term FDR plots in Figure 4 and Supplemental Figure 4

library(data.table)
library(RColorBrewer)
library(gplots)
library(dplyr)
library(ComplexHeatmap)
library(svglite)

my_col <- c("white", "orange", "brown")
save_dir <- "~/save_dir/"

#Function for creating iPS v. Human Embryo FDR bar plots
neglogq_bar_graph_meta <- function(data) {
  colnames(data)[colnames(data) == "Log.q.value."] = "log10(FDR)"
  data$"-log(q-value)" <- data$"log10(FDR)" * -1
  data$"q-value" <- 10^(data$"log10(FDR)")
  data_cat <- data %>% 
    filter(!!as.symbol("q-value") < 0.05)
  data_filter <- data_cat %>% 
    select(Description, !!as.symbol("log10(FDR)"), !!as.symbol("-log(q-value)"))
  
  data_top10 <- data_filter %>% 
    slice_min(!!as.symbol("log10(FDR)"), n = 10) %>%
    dplyr::arrange(dplyr::desc(!!as.symbol("log10(FDR)")))
  
  p <- ggplot(data_top10, aes(x = !!as.symbol("-log(q-value)"), 
                              y = reorder(Description, 
                                          !!as.symbol("-log(q-value)")))) +
    geom_bar(stat="identity", 
             fill = "darkred") +
    ylab(NULL)
  return(p)
}




##### HE, UP in vitro v. in vivo ######
he_up_go <- read.csv("~/csv_files/Metascape_GO_Analysis/iPS_v_HSC-Competent-Embryo_HE/HE_up-in-vitro_GO.csv")

#Create and save plot
he_up_plt <- neglogq_bar_graph_meta(he_up_go)
up_terms <- ###
ggsave(plot = he_up_plt, 
       filename = paste0(save_directory, 
                         "HE_upDEG_GO_Terms_-log10qvalue.svg"), 
       width = 5, 
       height = up_terms/4)




##### HE, DOWN in vitro v. in vivo ######
he_down_go <- read.csv("~/csv_files/Metascape_GO_Analysis/iPS_v_HSC-Competent-Embryo_HE/HE_down-in-vitro_GO.csv")

#Create and save plot
he_down_plt <- neglogq_bar_graph_meta(he_down_go)
down_terms <- ###
  ggsave(plot = he_down_plt, 
         filename = paste0(save_directory, 
                           "HE_downDEG_GO_Terms_-log10qvalue.svg"), 
         width = 5, 
         height = down_terms/4)



##### EHT, UP in vitro v. in vivo ######
eht_up_go <- read.csv("~/csv_files/Metascape_GO_Analysis/iPS_v_HSC-Competent-Embryo_EHT/EHT_up-in-vitro_GO.csv")

#Create and save plot
eht_up_plt <- neglogq_bar_graph_meta(eht_up_go)
up_terms <- 7
ggsave(plot = eht_up_plt, 
       filename = paste0(save_directory, 
                         "EHT_upDEG_GO_Terms_-log10qvalue.svg"), 
       width = 5, 
       height = up_terms/4)





##### EHT, DOWN in vitro v. in vivo ######
eht_down_go <- read.csv("~/csv_files/Metascape_GO_Analysis/iPS_v_HSC-Competent-Embryo_EHT/EHT_down-in-vitro_GO.csv")

#Create and save plot
eht_down_plt <- neglogq_bar_graph_meta(eht_down_go)
down_terms <- 10
ggsave(plot = eht_down_plt, 
       filename = paste0(save_directory, 
                         "EHT_downDEG_GO_Terms_-log10qvalue.svg"), 
       width = 5, 
       height = down_terms/4)





##### HSPC, UP in vitro v. in vivo ######
hspc_up_go <- read.csv("~/csv_files/Metascape_GO_Analysis/iPS_v_HSC-Competent-Embryo_HSPC/hspc_up-in-vitro_GO.csv")

#Create and save plot
hspc_up_plt <- neglogq_bar_graph_meta(hspc_up_go) 
up_terms <- 7
ggsave(plot = hspc_up_plt, 
       file = paste0(save_dir, "HSPC_upDEG_GO_Terms_-log10qvalue.png"), 
       width = 5, 
       height = up_terms/4,
       dpi = 600)





##### HSPC, DOWN in vitro v. in vivo ######
hspc_down_go <- read.csv("~/csv_files/Metascape_GO_Analysis/iPS_v_HSC-Competent-Embryo_HSPC/hspc_down-in-vitro_GO.csv")

#Create and save plot
hspc_down_plt <- neglogq_bar_graph_meta(hspc_down_go)
down_terms <- 10
ggsave(plot = hspc_down_plt, 
       file = paste0(save_dir, "HSPC_downDEG_GO_Terms_-log10qvalue.png"), 
       width = 5, 
       height = down_terms/4,
       dpi = 600)




##### AE v. HE, UP in HE v. AE ######
ae_he_up_go <- read.csv("~/csv_files/Metascape_GO_Analysis/AE_v_HE_iPS_HSC-Competent-Embryo_Batch_Analysis/AE-HE_up-in-HE_GO.csv")
ae_he_up_heatmap_terms <- read.csv("~/csv_files/Metascape_GO_Analysis/AE_v_HE_iPS_HSC-Competent-Embryo_Batch_Analysis/AE-HE_up-in-HE_Metascape-Heatmap-Selected-GO.csv")

ae_he_up_go_filt <- ae_he_up_go %>%
  select(Category, GO, Description, "Log.q.value.", GeneList) %>%
  filter(Category == "GO Biological Processes") %>%
  filter(Log.q.value. <= -1.30102999566) #log(0.05)

#select only metascape top terms
ae_he_up_heatmap_terms_list <- c(ae_he_up_heatmap_terms$GO)
ae_he_up_go_filt <- ae_he_up_go_filt[ae_he_up_go_filt$GO %in% ae_he_up_heatmap_terms_list,]

#Calculate -log10(q-value)
ae_he_up_go_filt$"-log(q-value)" <- -1*ae_he_up_go_filt$Log.q.value.

#Limit to Description and -log(q-value)
ae_he_up_go_filt <- ae_he_up_go_filt %>%
  select(Description, GeneList, "-log(q-value)")

#Fix rowname numbering
rownames(ae_he_up_go_filt) <- 1:nrow(ae_he_up_go_filt)

#Unmelt and create matrix from data
ae_he_up_go_unmelt <- dcast(data = ae_he_up_go_filt, 
                            formula = Description~GeneList)
rownames(ae_he_up_go_unmelt) <- ae_he_up_go_unmelt$Description
ae_he_up_go_unmelt <- select(ae_he_up_go_unmelt, 
                             -c("Description"))
ae_he_up_mat <- data.matrix(ae_he_up_go_unmelt)

#Reorder based on in vivo values (for better visualization)
ae_he_up_df <- data.frame(ae_he_up_mat)
ae_he_up_df <- ae_he_up_df[order(-ae_he_up_df$in.vivo),]
ae_he_up_df <- ae_he_up_df %>% select("in.vivo", "in.vitro")
colnames(ae_he_up_df) <- c("in vivo", "in vitro")
ae_he_up_mat <- as.matrix(ae_he_up_df)

#Heatmap
pheatmap_ae_he_up <- pheatmap(ae_he_up_mat,
                              legend = T,
                              show_rownames = T,
                              show_colnames = T,
                              cluster_cols = F,
                              cluster_rows = F,
                              cellheight = 20,
                              cellwidth = 20,
                              color = my_col,
                              main = "GO Biological Process Terms\n(Upregulated DEG, AE-HE)")
svglite(paste0(save_directory,"pheatmap_ae_he_up.svg"),
        width = 7,
        height = 5)
draw(pheatmap_ae_he_up)
dev.off()





##### AE v. HE, DOWN in HE v. AE ######
ae_he_down_go <- read.csv("~/csv_files/Metascape_GO_Analysis/AE_v_HE_iPS_HSC-Competent-Embryo_Batch_Analysis/AE-HE_down-in-HE_GO.csv")
ae_he_down_heatmap_terms <- read.csv("~/csv_files/Metascape_GO_Analysis/AE_v_HE_iPS_HSC-Competent-Embryo_Batch_Analysis/AE-HE_down-in-HE_Metascape-Heatmap-Selected-GO.csv")

ae_he_down_go_filt <- ae_he_down_go %>%
  select(Category, GO, Description, "Log.q.value.", GeneList) %>%
  filter(Category == "GO Biological Processes") %>%
  filter(Log.q.value. <= -1.30102999566) #log(0.05)

#select only metascape top terms
ae_he_down_heatmap_terms_list <- c(ae_he_down_heatmap_terms$GO)
ae_he_down_go_filt <- ae_he_down_go_filt[ae_he_down_go_filt$GO %in% ae_he_down_heatmap_terms_list,]

#Calculate -log10(q-value)
ae_he_down_go_filt$"-log(q-value)" <- -1*ae_he_down_go_filt$Log.q.value.

#Limit to Description and -log(q-value)
ae_he_down_go_filt <- ae_he_down_go_filt %>%
  select(Description, GeneList, "-log(q-value)")

#Fix rowname numbering
rownames(ae_he_down_go_filt) <- 1:nrow(ae_he_down_go_filt)

#Unmelt and create matrix from data
ae_he_down_go_unmelt <- dcast(data = ae_he_down_go_filt, 
                            formula = Description~GeneList)
rownames(ae_he_down_go_unmelt) <- ae_he_down_go_unmelt$Description
ae_he_down_go_unmelt <- select(ae_he_down_go_unmelt, 
                             -c("Description"))
ae_he_down_mat <- data.matrix(ae_he_down_go_unmelt)

#Reorder based on in vivo values (for better visualization)
ae_he_down_df <- data.frame(ae_he_down_mat)
ae_he_down_df <- ae_he_down_df[order(-ae_he_down_df$in.vivo),]
ae_he_down_df <- ae_he_down_df %>% select("in.vivo", "in.vitro")
colnames(ae_he_down_df) <- c("in vivo", "in vitro")
ae_he_down_mat <- as.matrix(ae_he_down_df)

#Heatmap
pheatmap_ae_he_down <- pheatmap(ae_he_down_mat,
                              legend = T,
                              show_rownames = T,
                              show_colnames = T,
                              cluster_cols = F,
                              cluster_rows = F,
                              cellheight = 20,
                              cellwidth = 20,
                              color = my_col,
                              main = "GO Biological Process Terms\n(Downregulated DEG, AE-HE)")
svglite(paste0(save_directory,"pheatmap_ae_he_down.svg"),
        width = 7,
        height = 5)
draw(pheatmap_ae_he_down)
dev.off()




##### HE v. HSPC, UP in HSPC v. HE ######
he_hspc_up_go <- read.csv("~/csv_files/Metascape_GO_Analysis/HE_v_HSPC_iPS_HSC-Competent-Embryo_Batch_Analysis/HE-HSPC_up-in-HSPC_GO.csv")
he_hspc_up_heatmap_terms <- read.csv("~/csv_files/Metascape_GO_Analysis/HE_v_HSPC_iPS_HSC-Competent-Embryo_Batch_Analysis/HE-HSPC_up-in-HSPC_Metascape-Heatmap-Selected-GO.csv")

he_hspc_up_go_filt <- he_hspc_up_go %>%
  select(Category, GO, Description, "Log.q.value.", GeneList) %>%
  filter(Category == "GO Biological Processes") %>%
  filter(Log.q.value. <= -1.30102999566) #log(0.05)

#select only metascape top terms
he_hspc_up_heatmap_terms_list <- c(he_hspc_up_heatmap_terms$GO)
he_hspc_up_go_filt <- he_hspc_up_go_filt[he_hspc_up_go_filt$GO %in% he_hspc_up_heatmap_terms_list,]

#Calculate -log10(q-value)
he_hspc_up_go_filt$"-log(q-value)" <- -1*he_hspc_up_go_filt$Log.q.value.

#Limit to Description and -log(q-value)
he_hspc_up_go_filt <- he_hspc_up_go_filt %>%
  select(Description, GeneList, "-log(q-value)")

#Fix rowname numbering
rownames(he_hspc_up_go_filt) <- 1:nrow(he_hspc_up_go_filt)

#Unmelt and create matrix from data
he_hspc_up_go_unmelt <- dcast(data = he_hspc_up_filt, 
                              formula = Description~GeneList)
rownames(he_hspc_up_go_unmelt) <- he_hspc_up_go_unmelt$Description
he_hspc_up_go_unmelt <- select(he_hspc_up_go_unmelt, 
                               -c("Description"))
he_hspc_up_mat <- data.matrix(he_hspc_up_go_unmelt)

#Reorder based on in vivo values (for better visualization)
he_hspc_up_df <- data.frame(he_hspc_up_mat)
he_hspc_up_df <- he_hspc_up_df[order(he_hspc_up_df$in.vivo),]
he_hspc_up_df <- he_hspc_up_df %>% select("in.vivo", "in.vitro")
colnames(he_hspc_up_df) <- c("in vivo", "in vitro")
he_hspc_up_mat <- as.matrix(he_hspc_up_df)

#Heatmap
pheatmap_he_hspc_up <- pheatmap(he_hspc_up_mat,
                                legend = T,
                                show_rownames = T,
                                show_colnames = T,
                                cluster_cols = F,
                                cluster_rows = F,
                                cellheight = 20,
                                cellwidth = 20,
                                color = my_col,
                                main = "GO Biological Process Terms\n(Upregulated DEG, HE-HSPC)")
svglite(paste0(save_directory,"pheatmap_he_hspc_up.svg"),
        width = 7,
        height = 5)
draw(pheatmap_he_hspc_up)
dev.off()




##### HE v. HSPC, DOWN in HSPC v. HE ######
he_hspc_down_go <- read.csv("~/csv_files/Metascape_GO_Analysis/HE_v_HSPC_iPS_HSC-Competent-Embryo_Batch_Analysis/HE-HSPC_down-in-HSPC_GO.csv")
he_hspc_down_heatmap_terms <- read.csv("~/csv_files/Metascape_GO_Analysis/HE_v_HSPC_iPS_HSC-Competent-Embryo_Batch_Analysis/HE-HSPC_down-in-HSPC_Metascape-Heatmap-Selected-GO.csv")

he_hspc_down_go_filt <- he_hspc_down_go %>%
  select(Category, GO, Description, "Log.q.value.", GeneList) %>%
  filter(Category == "GO Biological Processes") %>%
  filter(Log.q.value. <= -1.30102999566) #log(0.05)

#select only metascape top terms
he_hspc_down_heatmap_terms_list <- c(he_hspc_down_heatmap_terms$GO)
he_hspc_down_go_filt <- he_hspc_down_go_filt[he_hspc_down_go_filt$GO %in% he_hspc_down_heatmap_terms_list,]

#Calculate -log10(q-value)
he_hspc_down_go_filt$"-log(q-value)" <- -1*he_hspc_down_go_filt$Log.q.value.

#Limit to Description and -log(q-value)
he_hspc_down_go_filt <- he_hspc_down_go_filt %>%
  select(Description, GeneList, "-log(q-value)")

#Fix rowname numbering
rownames(he_hspc_down_go_filt) <- 1:nrow(he_hspc_down_go_filt)

#Unmelt and create matrix from data
he_hspc_down_go_unmelt <- dcast(data = he_hspc_down_filt, 
                              formula = Description~GeneList)
rownames(he_hspc_down_go_unmelt) <- he_hspc_down_go_unmelt$Description
he_hspc_down_go_unmelt <- select(he_hspc_down_go_unmelt, 
                               -c("Description"))
he_hspc_down_mat <- data.matrix(he_hspc_down_go_unmelt)

#Reorder based on in vivo values (for better visualization)
he_hspc_down_df <- data.frame(he_hspc_down_mat)
he_hspc_down_df <- he_hspc_down_df[order(he_hspc_up_df$in.vivo),]
he_hspc_down_df <- he_hspc_down_df %>% select("in.vivo", "in.vitro")
colnames(he_hspc_down_df) <- c("in vivo", "in vitro")
he_hspc_down_mat <- as.matrix(he_hspc_down_df)

#Heatmap
pheatmap_he_hspc_down <- pheatmap(he_hspc_down_mat,
                                legend = T,
                                show_rownames = T,
                                show_colnames = T,
                                cluster_cols = F,
                                cluster_rows = F,
                                cellheight = 20,
                                cellwidth = 20,
                                color = my_col,
                                main = "GO Biological Process Terms\n(Downregulated DEG, HE-HSPC)")
svglite(paste0(save_directory,"pheatmap_he_hspc_down.svg"),
        width = 7,
        height = 5)
draw(pheatmap_he_hspc_down)
dev.off()