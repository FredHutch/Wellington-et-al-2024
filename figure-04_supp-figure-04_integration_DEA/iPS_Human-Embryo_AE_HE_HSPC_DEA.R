#Comparison of AE and HE, or HE and HSPC, in either iPS or HSC-Competent Human Embyros

library(dplyr)
library(monocle3)
library(Matrix)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(viridis)
library(cowplot)
library(ggplotify)
library(harmony)
library(Seurat)
library(lme4)
library(sf)
library(ggprism)

cds <- readRDS("~/RDS_objects/integrated_data_Monocle_label-transfer.RDS")

#Create column containing monocle clusters
colData(cds)$monocle_clust <- clusters(cds)

#Create is_in_vitro column in metadata
colData(cds)$is_in_vitro <- colData(cds)$orig.ident
colData(cds)$is_in_vitro <- dplyr::recode(colData(cds)$is_in_vitro,
                                          "CS10" = "FALSE",
                                          "CS11" = "FALSE",
                                          "C13" = "FALSE",
                                          "45-iPS Day 8" = "TRUE",
                                          "CS16, AoV, 34+ sorted" = "FALSE",
                                          "CS16, AoMid, no sort" = "FALSE",
                                          "Calvanese_Aorta_Wk4.5" = "FALSE",
                                          "Calvanese_Aorta_Wk5" = "FALSE",
                                          "Calvanese_Aorta_Wk6" = "FALSE")
as.factor(pData(cds)$is_in_vitro)

#Create cds subset containing data only for iPS or HSC-Competent Human Embryo
#datasets
cds_ips_def <- cds[,!colData(cds)$orig.ident %in% c("CS10", "CS11", "CS13")]

###############Differential Expression Analysis###################

####AE-HE, iPS####
#Note: doesn't make sense to use mixed-binomial model since there is only one sample
#positive estimates means up in HE (cluster 75) relative to AE (cluster 36)
#negative estimates means down in HE (cluster 75) relative to AE (cluster 36)
cds_ips_AE_HE <- cds_ips_def[,colData(cds_ips_def)$monocle_clust %in% c("36", "75") &
                                  colData(cds_ips_def)$is_in_vitro %in% c("FALSE")]

gene_fits_ips_AE_HE <- fit_models(cds_ips_AE_HE,
                                  model_formula_str="~cluster") #uses quasipoisson

fit_coefs_ips_AE_HE <- coefficient_table(gene_fits_ips_AE_HE)

iPS_terms_AE_HE <- fit_coefs_ips_AE_HE %>% 
  filter(term == "cluster75") %>% 
  select(gene_short_name, term, q_value, estimate)

iPS_terms_AE_HE_sig <- fit_coefs_ips_AE_HE %>% 
  filter(term == "cluster75") %>% 
  filter(q_value < 0.05) %>%
  filter(estimate >= 0.25 | estimate <= -0.25) %>%
  select(gene_short_name, term, q_value, estimate) #saved to csv

####AE-HE, HSC-Competent Human Embryo####
#positive estimates means up in HE (cluster 75) relative to AE (cluster 36)
#negative estimates means down in HE (cluster 75) relative to AE (cluster 36)
cds_embryo_AE_HE <- cds_ips_def[,colData(cds_ips_def)$monocle_clust %in% c("36", "75") &
                                  colData(cds_ips_def)$is_in_vitro %in% c("FALSE")]

gene_fits_embryo_AE_HE <- fit_models(cds_embryo_AE_HE,
                                     model_formula_str="~cluster + (1|publication)",
                                     expression_family = "mixed-negbinomial")

fit_coefs_embryo_AE_HE <- coefficient_table(gene_fits_embryo_AE_HE)

embryo_terms_AE_HE <- fit_coefs_embryo_AE_HE %>% 
  filter(term == "cluster75") %>% 
  select(gene_short_name, term, q_value, estimate)

embryo_terms_AE_HE_sig <- fit_coefs_embryo_AE_HE %>% 
  filter(term == "cluster75") %>% 
  filter(q_value < 0.05) %>%
  filter(estimate >= 0.25 | estimate <= -0.25) %>%
  select(gene_short_name, term, q_value, estimate) #saved to csv

####HE-HSPC, iPS####
#Note: doesn't make sense to use mixed-binomial model since there is only one sample
#positive estimates means up in Early HSPC (cluster 74) relative to HE (cluster 75)
#negative estimates means down in Early HSPC (cluster 74) relative to HE (cluster 75)
#Note that it makes more sense to consider the transition from HE to HSPC and the
#estimate signs can be reversed to get this info (since the comparison is the same,
#but the analysis just automatically generates results for the highest cluster number
#relative to the lowest cluster number so the signs need to be reversed)

cds_ips_HE_HSPC <- cds_ips_def[,colData(cds_ips_def)$monocle_clust %in% c("75", "74") &
                               colData(cds_ips_def)$is_in_vitro %in% c("FALSE")]

gene_fits_ips_HE_HSPC <- fit_models(cds_ips_HE_HSPC,
                                  model_formula_str="~cluster") #uses quasipoisson

fit_coefs_ips_HE_HSPC <- coefficient_table(gene_fits_ips_HE_HSPC)

iPS_terms_HE_HSPC <- fit_coefs_ips_HE_HSPC %>% 
  filter(term == "cluster74") %>% 
  select(gene_short_name, term, q_value, estimate)

iPS_terms_HE_HSPC_sig <- fit_coefs_ips_HE_HSPC %>% 
  filter(term == "cluster74") %>% 
  filter(q_value < 0.05) %>%
  filter(estimate >= 0.25 | estimate <= -0.25) %>%
  select(gene_short_name, term, q_value, estimate) #saved to csv

####HE-HSPC, HSC-Competent Human Embryo####
#positive estimates means up in Early HSPC (cluster 74) relative to HE (cluster 75)
#negative estimates means down in Early HSPC (cluster 74) relative to HE (cluster 75)
#Note that it makes more sense to consider the transition from HE to HSPC and the
#estimate signs can be reversed to get this info (since the comparison is the same,
#but the analysis just automatically generates results for the lowest cluster number)
#relative to the lowest cluster number so the signs need to be reversed)

cds_embryo_HE_HSPC <- cds_ips_def[,colData(cds_ips_def)$monocle_clust %in% c("75", "74") &
                                  colData(cds_ips_def)$is_in_vitro %in% c("FALSE")]

gene_fits_embryo_HE_HSPC <- fit_models(cds_embryo_HE_HSPC,
                                     model_formula_str="~cluster + (1|publication)",
                                     expression_family = "mixed-negbinomial")

fit_coefs_embryo_HE_HSPC <- coefficient_table(gene_fits_embryo_HE_HSPC)

embryo_terms_HE_HSPC <- fit_coefs_embryo_HE_HSPC %>% 
  filter(term == "cluster75") %>% 
  select(gene_short_name, term, q_value, estimate)

embryo_terms_HE_HSPC_sig <- fit_coefs_embryo_HE_HSPC %>% 
  filter(term == "cluster75") %>% 
  filter(q_value < 0.05) %>%
  filter(estimate >= 0.25 | estimate <= -0.25) %>%
  select(gene_short_name, term, q_value, estimate) #saved to csv
