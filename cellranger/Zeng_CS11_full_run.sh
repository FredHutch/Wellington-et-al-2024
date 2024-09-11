#!/bin/bash -eu
#SBATCH -N1 -n8
echo "Running on $(host) starting $(date)"
module purge
module load CellRanger/4.0.0
time cellranger count \
     --id=CS11_RCW\
     --fastqs=/fh/fast/hadland_b/Rachel/Human_Embryo_Data_Zeng_etal_2019/CS11 \
     --sample=SRR9875919,SRR9875920,SRR9875921,SRR9875922 \
     --transcriptome=/shared/biodata/ngs/Reference/10X/refdata-gex-GRCh38-2020-A \
     --localcores=8 \
     --localmem=64