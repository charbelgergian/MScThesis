#!/bin/bash
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=200:mem=72gb

#Example of script to run modisco 

module load anaconda3/personal
source activate modisco

cd $PBS_O_WORKDIR

#Run modisco algorithm
modisco motifs -s /rds/general/user/cg2723/home/tfFinal/modelfiles/erg/seq.npy -a /rds/general/user/cg2723/home/tfFinal/modelfiles/erg/shap.npy -w 250 -n 100000 -o erg_modisco_results.h5

#Report motifs in html file
modisco report -i erg_modisco_results.h5 -o /rds/general/user/cg2723/home/tfFinal/modisco/results/erg

#Convert results into meme format
modisco meme -i erg_modisco_results.h5 -o /rds/general/user/cg2723/home/tfFinal/modisco/results/erg/ergPFM.meme -t PFM 

