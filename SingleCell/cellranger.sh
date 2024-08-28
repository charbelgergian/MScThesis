#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=16:mem=64gb

export PATH=$PATH:/rds/general/user/cg2723/home/cellrangerATAC/cellranger-atac-2.1.0

cd $PBS_O_WORKDIR

cellranger-atac count --id=IGF133743\
 --reference=/rds/general/user/cg2723/home/cellrangerATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--fastqs=/rds/general/user/cg2723/projects/epinott/live/raw_reads/AACHJF2HV/1/8_10X/IGF133743,/rds/general/user/cg2723/projects/epinott/live/raw_reads/AACHJF2HV/2/8_10X/IGF133743 \
--localcores=16 --localmem=64
