#!/bin/bash
#SBATCH --job-name neal_netprop
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --nodes 1
#SBATCH -a 1-69
#SBATCH -c 4
#SBATCH -t 8:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/neal_netprop-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/neal_netprop-%j.e
#SBATCH --mail-type END,neal_netprop
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account csd795 




cd /tscc/projects/ps-palmer/brittany/SUD_cross_species/scripts

# Modify the paths accordingly

source activate env-std-py38 

export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH


python neale_netprop.py $SLURM_ARRAY_TASK_ID
