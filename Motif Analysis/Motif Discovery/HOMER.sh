#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=20:mem=96gb

$PBS_O_WORKDIR
PATH=$PATH:/rds/general/user/cg2723/home/Homer/bin

findMotifsGenome.pl\
 /rds/general/user/cg2723/home/tfFinal/erg_ADvas_AAA_20240422_tf.bed\
  hg38\
   /rds/general/user/cg2723/home/tfFinal/Homer/allerg_auto_bg -size 250\
