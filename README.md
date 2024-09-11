# Wellington et al. 2024

## Description
This is the code that was used for the data analysis carried out in the Wellington et al. 2024 paper, "Developmental regulation of endothelial-to-hematopoietic transition from iPSCs".

Associated input files for the code can be found in this [google drive](https://drive.google.com/drive/folders/154wYpNe8G2jP0YeaFASthEk0h4uEw_ck?usp=sharing) at the listed path.

### cellranger:
- contains scripts for [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome) processing of raw fastq data for Zeng et al. CS10, CS11 and CS13 and 10X iPS Day 8 datasets
- Crosse et al. fastq processing was done on the 10X cloud and did not require an input script
- Calvanese et al. provided csv files containing raw counts, which were used for analysis
- sci-RNA-seq data was returned by the Brotman Baty Institute as a cell_data_set object post-cellranger

### figure-01_supp-figure-01_sciRNAseq:
- contains the finalized code for the analysis of sci-RNA-seq data as associated with the generation of figure panels in Figure 1 and Figure S1
- sci_endo_hemato_scoring: the code used to identify endothelial and hematopoeitic populations in Figure 1D and Figure S1F
- sci_hemato_related_table: the code used to generate Figure S1E
- sci_preprocessing: the code used for preprocessing the sci-RNA-seq data recieved from the Brotman Baty Institute
- sci_Xu-et-al_marker_scoring: the code used for Figure 1B and Figure S1A-B,C

### figure-02_supp-figure-02_iPS-10X:
- contains the finalized code for the analysis of 10X iPS Day 8 data as associated with the generation of figure panels in Figure 2 and Figure S2
- iPS_10X_26a-v-26b_HSC-genes: the code used to generate Figure 2G
- iPS_10X_marker_scoring: the code used to generate Figure 2D-F
- iPS_10X_preprocessing: the code used for preprocessing the 10X iPS Day 8 data after cellranger

### figure-03_supp-figure-03_integration
- contains the finalized code for the integration of 10X iPS Day 8 data with Zeng et al. 2019, Crosse et al. 2020, and Calvanese et al. 2022 10X Human Embyro datasets as associated with the generation of figure panels in Figure 3 and Figure S3
- iPS_Human-Embryo_Label-Transfer: the code used to generate Figure 3D-M and Figure S3A-R
- iPS_Human-Embryo_Preprocessing_Integration: the code used to integrate iPS Day 8 10X data with Zeng et al. 2019, Crosse et al. 2020, and Calvanese et al. 2022 10X Human Embyro datasets
- iPS_Human-Embryo_Seurat_to_Monocle_Conversion: the code used to convert the integrated Seurat object to a Monocle 3 object for use in label transfer, differential expression analysis, etc.

### figure-04_supp-figure-04_integration_DEA
- contains the finalized code for the differential expression analysis of 10X iPS Day 8 and HSC-Competent Human Embryo data as associated with the generation of figure panels in Figure 4 and Figure S4
- GO-Term-Aggregate-Gene-Scores_AE_HE_HSPC: the code used to generate Figure S4G-J
- iPS_Human-Embryo_AE_HE_HSPC_DEA: the code used for differential expression analysis comparing AE-HE or HE-HSPC for both iPS and HSC-Competent Human Embryo 
- iPS_v_HSC-competent-Human-Embryo: the code used for differential expression analysis comparing iPS and HSC-Competent Human Embryo HE, Early HSPCs or EHT (HE and Early HSPCs)
- Metascape_FDR_Plots: the code used to generate Figure 4B,C,E,F,H,I and Figure S4B,C,E,F

### figure-05_supp-figure-05_nichenetr
- contains the finalized code for the [NicheNetR](https://github.com/saeyslab/nichenetr) ligand-receptor analysis associated with the HE and EHT (HE and Early HSPCs) differential expression profiles from Figure 4 associated with the generation of figure panels in Figure 5 and Figure S5
- NicheNetR_circos-plots_ligand-activity-tables_ligand-DEG-heatmaps: the code used to generate Figure 5A,B and Figure S5A-D
- NicheNetR_iPS_v_HSC-competent-Human-Embryo: the code used to run NicheNetR on the integrated data for HE or EHT (HE and Early HSPCs), comparing iPS Day 8 to HSC-Competent Human Embryos