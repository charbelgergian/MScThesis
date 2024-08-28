#PBS -l walltime=18:00:00
#PBS -l select=1:ncpus=20:mem=100gb

export PATH=$PATH:/rds/general/user/cg2723/home/meme/bin

cd $PBS_O_WORKDIR

fimo --oc '/rds/general/user/cg2723/home/tfFinal/fimo/HomerERG' --verbosity 1 \
'/rds/general/user/cg2723/home/tfFinal/fimo/hocomocoHomerERG.meme'\
'/rds/general/user/cg2723/home/tfFinal/erg_ADvas_AAA_20240422_tf.fa'
