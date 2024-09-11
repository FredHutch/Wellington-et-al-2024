#Differential Expression Analysis for Integrated Data, iPS v. HSC-Competent
#Human Embryo

library(monocle3)
library(ggplot2)
library(dplyr)
library(ggeasy)
library(garnett)
library(Matrix)
library(tidyr)
library(ggrepel)
library(viridis)
library(cowplot)
library(ggplotify)
library(harmony)
library(Seurat)
library(lme4)
library(sf)
library(xpectr)

cds <- readRDS("~/RDS_objects/integrated_data_Monocle_label-transfer.RDS")

#Create "is_in_vitro" column in cds metadata
colData(cds)$is_in_vitro <- "NA"

#Create cds subset containing only iPS and HSC-competent embryo data and add
#metadata column containing TRUE or FALSE regarding if it is in vitro
cds_iPS <- cds[,colData(cds)$dataset %in% "45-iPS Day 8"]
cds_invivo <- cds[,!colData(cds)$dataset %in% "45-iPS Day 8"]

#If cells are from iPS, put TRUE in the "is_in_vitro" metadata column,
#If cells are from Human Embryo, put FALSE in the "is_in_vitro" metadata column
colData(cds)[colnames(cds_iPS),]$is_in_vitro <- TRUE
colData(cds)[colnames(cds_invivo),]$is_in_vitro <- FALSE

#Create "stage" column in cds metadata
colData(cds)$stage <- "NA"

#Create cds subsets containing individual human embryo stages; note that 
#the cds_ips was already created above
cds_crosse <- cds[,colData(cds)$dataset %in% "Human Embryo (Crosse et al)"]
cds_cs14 <- cds[,colData(cds)$dataset %in% "Human Embryo, Week 4.5 (CS14)"]
cds_cs15 <- cds[,colData(cds)$dataset %in% "Human Embryo, Week 5 (CS15)"]
cds_cs17 <- cds[,colData(cds)$dataset %in% "Human Embryo, Week 6 (CS16)"]

#Add embryo stage information to stage metadata column in cds
#Note that Zeng et al cells will have "NA" since it isn't relevant to this DE analysis
colData(cds)[colnames(cds_iPS),]$stage <- "in_vitro"
colData(cds)[colnames(cds_crosse),]$stage <- "CS16"
colData(cds)[colnames(cds_cs14),]$stage <- "CS14"
colData(cds)[colnames(cds_cs15),]$stage <- "CS15"
colData(cds)[colnames(cds_cs17),]$stage <- "CS17"

#Create cds objects that contain only specific cell population and only 
#HSC-competent Human Embryo or iPS
HE <- cds[,!colData(cds)$dataset %in% "Human Embryo (Zeng et al)" &
            clusters(cds) %in% c("75")]
EHT <- cds[,!colData(cds)$dataset %in% "Human Embryo (Zeng et al)" &
             clusters(cds) %in% c("75", "74")]
Early_HSPC <- cds[,!colData(cds)$dataset %in% "Human Embryo (Zeng et al)" &
                    clusters(cds) %in% c("74")]

####HE####
#Monocle cluster 75
gene_fits_HE <- fit_models(HE, 
                           model_formula_str="~is_in_vitro + stage + (1|dataset)", 
                           expression_family = "mixed-negbinomial")

fit_coefs_HE <- coefficient_table(gene_fits_HE)

in_vitro_terms_HE <- fit_coefs_HE %>% 
  filter(term == "is_in_vitroTRUE")

in_vitro_terms_HE_sig <- in_vitro_terms_HE %>% 
  filter(q_value < 0.05) %>% 
  filter(estimate >= 0.25 | estimate <= -0.25) %>%
  select(gene_short_name, term, q_value, estimate) #saved to csv

####EHT (HE + Early HSPC)####
#Monocle clusters 75 and 74
gene_fits_EHT <- fit_models(EHT, 
                           model_formula_str="~is_in_vitro + stage + (1|dataset)", 
                           expression_family = "mixed-negbinomial")

fit_coefs_EHT <- coefficient_table(gene_fits_EHT)

in_vitro_terms_EHT <- fit_coefs_EHT %>% 
  filter(term == "is_in_vitroTRUE")

in_vitro_terms_EHT_sig <- in_vitro_terms_EHT %>% 
  filter(q_value < 0.05) %>% 
  filter(estimate >= 0.25 | estimate <= -0.25) %>%
  select(gene_short_name, term, q_value, estimate) #saved to csv

####Early HSPC####
#Monocle cluster 74
gene_fits_Early_HSPC <- fit_models(Early_HSPC, 
                                   model_formula_str="~is_in_vitro + (1|dataset)", 
                                   expression_family = "mixed-negbinomial")

fit_coefs_Early_HSPC <- coefficient_table(gene_fits_Early_HSPC)

in_vitro_terms_Early_HSPC <- fit_coefs_Early_HSPC %>% 
  filter(term == "is_in_vitroTRUE")

in_vitro_terms_Early_HSPC_sig <- in_vitro_terms_Early_HSPC %>% 
  filter(q_value < 0.05) %>% 
  filter(estimate >= 0.25 | estimate <= -0.25) %>%
  select(gene_short_name, term, q_value, estimate) #saved to csv