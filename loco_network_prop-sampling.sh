#!/bin/bash 
#SBATCH -J ext_network_subsampling
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu 4Gb
#SBATCH --export ALL
#SBATCH -o /tscc/projects/ps-palmer/brittany/SUD_cross_species/job_run_out/ext_network_subsampling-%j.o
#SBATCH -e /tscc/projects/ps-palmer/brittany/SUD_cross_species/job_run_out/ext_network_subsampling-%j.e  
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH -A csd795    

cd /tscc/projects/ps-palmer/brittany/SUD_cross_species/scripts/

source activate env-std-py38

python /tscc/projects/ps-palmer/brittany/SUD_cross_species/scripts/loco_network_prop-sampling.py
