#!/bin/bash
#SBATCH -J ext_network_control_permutation
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --nodes 1
#SBATCH -a 1-10
#SBATCH -c 4
#SBATCH -t 8:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH -o /tscc/projects/ps-palmer/brittany/SUD_cross_species/job_run_out/ext_network_subsampling-%j.o
#SBATCH -e /tscc/projects/ps-palmer/brittany/SUD_cross_species/job_run_out/ext_network_subsampling-%j.e  
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account csd795




cd /tscc/projects/ps-palmer/brittany/SUD_cross_species/scripts/

source activate env-std-py38

python /tscc/projects/ps-palmer/brittany/SUD_cross_species/scripts/loco_network_prop.py $SLURM_ARRAY_TASK_ID
