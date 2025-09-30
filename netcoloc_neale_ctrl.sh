#!/bin/bash
#SBATCH --job-name neal_netcoloc
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --nodes 1
#SBATCH -a 1-1
#SBATCH -c 4
#SBATCH -t 1:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/neal_netcoloc-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/neal_netcoloc-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account csd795 


traits=( 'ext' 'loco_final_cf' )
cutoffs=( 'top500' 'FDR' )

t=${traits[$((SLURM_ARRAY_TASK_ID-1))]}
c=${cutoffs[$((SLURM_ARRAY_TASK_ID-1))]}


cd /tscc/projects/ps-palmer/brittany/SUD_cross_species/scripts

# Modify the paths accordingly

source activate env-std-py38 

export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH


python new_neale_ctrl_netcoloc_loco.py $t $c
