#NicheNetR Circos Plots and Ligand Activity Tables

library(nichenetr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(svglite)
library(data.table)
library(ggplot2)
library(reshape2)
library(dplyr)

#Load NicheNetR outputs
he <- readRDS("~/RDS_files/NicheNetR_output/nichenetR_output_HE.RDS")
eht <- readRDS("~/RDS_files/NicheNetR_output/nichenetR_output_EHT.RDS")

#Create save directory for plots
save_dir <- "~/save_dir/"

#Numbering function for heatmaps
every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}

############### HE Circos ##################

#Extract ligand-target information
he_ligand_target <- he$ligand_target_matrix
he_ligand_target_df <- data.frame(he_ligand_target)
he_ligand_target_rownames <- rownames(he_ligand_target_df)
he_ligand_activities <- he$ligand_activities
he_ligand_target_df <- he_ligand_target_df[rownames(he_ligand_target_df) %in% he_ligand_activities$test_ligand,]
he_ligand_target_df$ligand <- rownames(he_ligand_target_df)
he_ligand_target_melt <- melt(he_ligand_target_df, id.vars = c("ligand"))
names(he_ligand_target_melt) <- c("ligand", "target", "weight")

#Extract Activity Score information for ligand-receptor pairs
he_aupr <- he_ligand_activities[he_ligand_activities$test_ligand %in% he_ligand_target_rownames,]
he_aupr_df_preord <- data.frame(he_aupr) %>% 
  select(test_ligand, aupr_corrected) %>%
  arrange(-aupr_corrected)
he_aupr_df_preord$aupr_corrected <- round(he_aupr_df_preord$aupr_corrected, digits = 4)
colnames(he_aupr_df_preord) <- c("Ligand","Activity Score (AUPR)")

#HE ligand-receptor values, reordered by AUPR
he_ligand_receptor <- he$ligand_receptor_df
he_ligand_receptor <- as.data.frame(he_ligand_receptor)
he_lr_unmelt <- dcast(data = he_ligand_receptor, formula = ligand~receptor)
he_lr_unmelt <- he_lr_unmelt %>% replace(is.na(.), as.double(0))

#Filter only top ligands
he_top_ligands_filt <- he_aupr_df_preord$Ligand

#reverse list since will be printed top to bottom by geom_tile
he_top_ligands_filt <- rev(he_top_ligands_filt)

#Reorder ligand-receptor dataframe
he_lr_unmelt_reordered <- he_lr_unmelt %>% slice(match(he_top_ligands_filt, ligand))

#remove columns that only have 0 values
he_lr_unmelt_reordered <- he_lr_unmelt_reordered[,colSums(he_lr_unmelt_reordered != 0) > 0]

#Melt the reordered ligand-receptor dataframe
he_lr_remelt <- melt(he_lr_unmelt_reordered)
he_lr_remelt <- he_lr_remelt %>% rename("receptor" = variable, "weight" = value)

#Add factors, add levels and rename columns in the melted ligand-receptor dataframe
he_lr_remelt_structured <- he_lr_remelt
he_lr_remelt_structured$ligand <- factor(he_lr_remelt$ligand, 
                                         levels = unique(he_lr_remelt$ligand))
names(he_lr_remelt_structured) <- c("ligand", "target", "weight")

#Extract circos links, ligand order, receptor order and overall order
circos_links_he <-he_lr_remelt_structured %>% mutate(ligand = paste(ligand," "))
ligand_order <- as.character(unique(circos_links_he$ligand))
receptor_order <- as.character(unique(circos_links_he$target))
order <- c(ligand_order, receptor_order)

#Create circos plot
svglite(paste0(save_dir, "HE_nichenet_ligand_receptor_circos_plt.svg"),
        width = 15, height = 15)
set.seed(42)
chordDiagram(circos_links_he, 
             directional = 1,
             order=order,
             link.sort = TRUE, 
             link.decreasing = FALSE, 
             transparency = 0, 
             diffHeight = 0.005, 
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", 
             link.visible = circos_links_he$weight,
             annotationTrack = c("grid"), 
             preAllocateTracks = list(track.height = 0.075),
             col = grid.col)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)
dev.off()


############### EHT Circos ##################

#Extract ligand-target information
eht_ligand_target <- eht$ligand_target_matrix
eht_ligand_target_df <- data.frame(eht_ligand_target)
eht_ligand_target_rownames <- rownames(eht_ligand_target_df)
eht_ligand_activities <- eht$ligand_activities
eht_ligand_target_df <- eht_ligand_target_df[rownames(eht_ligand_target_df) %in% eht_ligand_activities$test_ligand,]
eht_ligand_target_df$ligand <- rownames(eht_ligand_target_df)
eht_ligand_target_melt <- melt(eht_ligand_target_df, id.vars = c("ligand"))
names(eht_ligand_target_melt) <- c("ligand", "target", "weight")

#Extract Activity Score information for ligand-receptor pairs
eht_aupr <- eht_ligand_activities[eht_ligand_activities$test_ligand %in% eht_ligand_target_rownames,]
eht_aupr_df_preord <- data.frame(eht_aupr) %>% 
  select(test_ligand, aupr_corrected) %>%
  arrange(-aupr_corrected)
eht_aupr_df_preord$aupr_corrected <- round(eht_aupr_df_preord$aupr_corrected, digits = 4)
colnames(eht_aupr_df_preord) <- c("Ligand","Activity Score (AUPR)")

#EHT ligand-receptor values, reordered by AUPR
eht_ligand_receptor <- eht$ligand_receptor_df
eht_ligand_receptor <- as.data.frame(eht_ligand_receptor)
eht_lr_unmelt <- dcast(data = eht_ligand_receptor, formula = ligand~receptor)
eht_lr_unmelt <- eht_lr_unmelt %>% replace(is.na(.), as.double(0))

eht_top_ligands_filt <- eht_aupr_df_preord$Ligand
#reverse list since will be printed top to bottom by geom_tile
eht_top_ligands_filt <- rev(eht_top_ligands_filt)

eht_lr_unmelt_reordered <- eht_lr_unmelt %>% dplyr::slice(match(eht_top_ligands_filt, ligand))

eht_lr_remelt <- melt(eht_lr_unmelt_reordered)
eht_lr_remelt <- eht_lr_remelt %>% rename(variable = "receptor", value = "weight")

eht_lr_remelt_structured <- eht_lr_remelt
eht_lr_remelt_structured$ligand <- factor(eht_lr_remelt$ligand, 
                                          levels = unique(eht_lr_remelt$ligand))

#EHT ligand-receptor values, reordered by AUPR
eht_ligand_receptor <- eht$ligand_receptor_df
eht_ligand_receptor <- as.data.frame(eht_ligand_receptor)
eht_lr_unmelt <- dcast(data = eht_ligand_receptor, formula = ligand~receptor)
eht_lr_unmelt <- eht_lr_unmelt %>% replace(is.na(.), as.double(0))

#Filter only top ligands
eht_top_ligands_filt <- eht_aupr_df_preord$Ligand

#reverse list since will be printed top to bottom by geom_tile
eht_top_ligands_filt <- rev(eht_top_ligands_filt)

#Reorder ligand-receptor dataframe
eht_lr_unmelt_reordered <- eht_lr_unmelt %>% dplyr::slice(match(eht_top_ligands_filt, ligand))

#remove columns that only have 0 values
eht_lr_unmelt_reordered <- eht_lr_unmelt_reordered[,colSums(eht_lr_unmelt_reordered != 0) > 0]

#Melt the reordered ligand-receptor dataframe
eht_lr_remelt <- melt(eht_lr_unmelt_reordered)
eht_lr_remelt <- eht_lr_remelt %>% rename(variable = "receptor", value = "weight")

#Add factors, add levels and rename columns in the melted ligand-receptor dataframe
eht_lr_remelt_structured <- eht_lr_remelt
eht_lr_remelt_structured$ligand <- factor(eht_lr_remelt$ligand, 
                                          levels = unique(eht_lr_remelt$ligand))
names(eht_lr_remelt_structured) <- c("ligand", "target", "weight")

#Extract circos links, ligand order, receptor order and overall order
circos_links_eht <-eht_lr_remelt_structured %>% mutate(ligand = paste(ligand," "))
ligand_order <- as.character(unique(circos_links_eht$ligand))
receptor_order <- as.character(unique(circos_links_eht$target))
order <- c(ligand_order, receptor_order)

#Create EHT circos plot
svglite(paste0(save_dir, "EHT_nichenet_ligand_receptor_circos_plt.svg"),
        width = 11, height = 10)
chordDiagram(circos_links_eht, 
             directional = 1,
             order=order,
             link.sort = TRUE, 
             link.decreasing = FALSE, 
             transparency = 0, 
             diffHeight = 0.005, 
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", 
             link.visible = circos_links_eht$weight,
             annotationTrack = c("grid"), 
             preAllocateTracks = list(track.height = 0.075),
             grid.col = grid_col_named)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)
dev.off()


############### HE Ligand Activity Table ##################

#Extract Matrix and convert to dataframe
he_ligand_target <- he$ligand_target_matrix
he_ligand_target_df <- data.frame(he_ligand_target)

#Extract number of DEGs associated with signaling from each ligand
he_deg_totals <- rowSums(he_ligand_target != 0)
he_deg_totals_df <- data.frame(he_deg_totals)
he_deg_totals_df$test_ligand <- rownames(he_deg_totals_df)

#Merge DEG number into ligand activity table from HE NicheNetR output in prep
#for generating a table of activity scores
he_df <- merge(he_deg_totals_df, he$ligand_activities, by = "test_ligand")
he_df <- he_df %>% 
  select(test_ligand, aupr_corrected, he_deg_totals) %>%
  arrange(-aupr_corrected)

#Round the activity score
he_df$aupr_corrected <- round(he_df$aupr_corrected, digits = 4)

#Rename final ligand activity and DEG dataframe columns
colnames(he_df) <- c("Ligand", "Activity Score (AUPR)", "Regulated DEGs")

#Convert ligand activity and DEG dataframe to tibble
he_tib <- as_tibble(he_df)

#Create ligand activity score and associated DEG table
he_tbl <- gt(he_tib) |>
  tab_header(title = "Ligand Activities and Associated DEGs (HE)",
             subtitle = "Ranked According to Corrected AUPR") |>
  gtsave(filename = paste0(save_dir, "he_ligand_activity_aupr_deg-number.png"))


############### EHT Ligand Activity Table ##################

#Extract Matrix and convert to dataframe
eht_ligand_target <- eht$ligand_target_matrix
eht_ligand_target_df <- data.frame(eht_ligand_target)

#Extract number of DEGs associated with signaling from each ligand
eht_deg_totals <- rowSums(eht_ligand_target != 0)
eht_deg_totals_df <- data.frame(eht_deg_totals)
eht_deg_totals_df$test_ligand <- rownames(eht_deg_totals_df)

#Merge DEG number into ligand activity table from HE NicheNetR output in prep
#for generating a table of activity scores
eht_df <- merge(eht_deg_totals_df, eht$ligand_activities, by = "test_ligand")
eht_df <- eht_df %>% 
  select(test_ligand, aupr_corrected, he_deg_totals) %>%
  arrange(-aupr_corrected)

#Round the activity score
eht_df$aupr_corrected <- round(eht_df$aupr_corrected, digits = 4)

#Rename final ligand activity and DEG dataframe columns
colnames(eht_df) <- c("Ligand", "Activity Score (AUPR)", "Regulated DEGs")

#Convert ligand activity and DEG dataframe to tibble
he_tib <- as_tibble(eht_df)

#Create ligand activity score and associated DEG table
eht_tbl <- gt(eht_tib) |>
  tab_header(title = "Ligand Activities and Associated DEGs (EHT)",
             subtitle = "Ranked According to Corrected AUPR") |>
  gtsave(filename = paste0(save_dir, "eht_ligand_activity_aupr_deg-number.png"))



############### HE Ligand DEG Regulatory Potential Heatmap ##################

#Extract ligand-target matrix
he_ligand_target_matrix <- he$ligand_target_matrix
he_mat <- as.matrix(he_ligand_target_matrix)

#Melt the ligand-target matrix
he_data <- melt(he_mat)

#Transform DEG names into numbers
he_data <- transform(he_data, num = as.numeric(factor(Var2)))

#Add order to numbers to ensure they remain in the correct order
he_data$ord.num <- factor(he_data$num, 
                          ordered = TRUE, 
                          levels = unique(he_data$num))

#Create heatmap (using geom_tile)
he_plt <- ggplot(he_data, 
                 aes(y = ord.num, 
                     x = Var1, 
                     fill = value)) +
  geom_tile(color = "black") +
  labs(x = "Predicted Ligands",
       y = "Prioritized Target Genes",
       fill = "Regulatory Potential") +
  scale_fill_gradient(high = "darkred", 
                      low = "white") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust = 0),
        legend.position = "bottom",
        legend.text = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust = 0.4),
        legend.title = element_text(vjust = 0.9),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  scale_y_discrete(breaks = every_nth(n = 10))



############### EHT Ligand DEG Regulatory Potential Heatmap ##################

#Extract ligand-target matrix
eht_ligand_target_matrix <- eht$ligand_target_matrix
eht_mat <- as.matrix(eht_ligand_target_matrix)

#Melt the ligand-target matrix
eht_data <- melt(eht_mat)

#Transform DEG names into numbers
eht_data <- transform(eht_data, num = as.numeric(factor(Var2)))

#Add order to numbers to ensure they remain in the correct order
eht_data$ord.num <- factor(eht_data$num, 
                          ordered = TRUE, 
                          levels = unique(he_data$num))

#Create heatmap (using geom_tile)
eht_plt <- ggplot(eht_data, 
                 aes(y = ord.num, 
                     x = Var1, 
                     fill = value)) +
  geom_tile(color = "black") +
  labs(x = "Predicted Ligands",
       y = "Prioritized Target Genes",
       fill = "Regulatory Potential") +
  scale_fill_gradient(high = "darkred", 
                      low = "white") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust = 0),
        legend.position = "bottom",
        legend.text = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust = 0.4),
        legend.title = element_text(vjust = 0.9),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  scale_y_discrete(breaks = every_nth(n = 10))