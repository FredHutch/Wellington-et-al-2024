#!/bin/bash -eu
#SBATCH -N1 -n8
echo "Running on $(host) starting $(date)"
module purge
module load CellRanger/4.0.0
time cellranger count \
     --id=D_Day8_EB \
     --fastqs=/fh/fast/hadland_b/Rachel/2020-10-23_Sergei_10X/201021_A00613_0190_AHTH3HDRXX/cellranger/mkfastq/HTH3HDRXX/outs/fastq_path/HTH3HDRXX/D_day8_EB_lib \
     --sample=D_day8_EB_lib  \
     --transcriptome=/shared/biodata/ngs/Reference/10X/refdata-gex-GRCh38-2020-A \
     --localcores=8 \
     --localmem=64