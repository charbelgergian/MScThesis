#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=10:mem=10gb

$PBS_O_WORKDIR

PATH=$PATH:/rds/general/user/cg2723/home/meme/bin

tomtom -oc /rds/general/user/cg2723/home/tfFinal/tomtom/homer_erg_auto_bg_250 /rds/general/user/cg2723/home/tfFinal/Homer/allerg_auto_bg_250.meme /rds/general/user/cg2723/home/tfFinal/H12CORE_meme_format.meme
