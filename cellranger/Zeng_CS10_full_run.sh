#!/bin/bash -eu
#SBATCH -N1 -n8
echo "Running on $(host) starting $(date)"
module purge
module load CellRanger/4.0.0
time cellranger count \
     --id=CS10_RCW\
     --fastqs=/fh/fast/hadland_b/Rachel/Human_Embryo_Data_Zeng_etal_2019/CS10 \
     --sample=SRR9875915,SRR9875916,SRR9875917,SRR9875918 \
     --transcriptome=/shared/biodata/ngs/Reference/10X/refdata-gex-GRCh38-2020-A \
     --localcores=8 \
     --localmem=64