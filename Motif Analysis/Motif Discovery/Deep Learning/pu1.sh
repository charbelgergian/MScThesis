#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=4:mem=24gb:ngpus=1

module load anaconda3/personal
source activate tf2_env

cd $PBS_O_WORKDIR

python /rds/general/user/cg2723/home/tfFinal/Deep/pu1.py
