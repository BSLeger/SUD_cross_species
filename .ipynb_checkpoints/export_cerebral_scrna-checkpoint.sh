#!/bin/bash 
#SBATCH -J scRNA_process_cerebral-b
#SBATCH -p platinum
#SBATCH -q hcp-csd795
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem-per-cpu 16Gb
#SBATCH --export ALL
#SBATCH -o /tscc/projects/ps-palmer/brittany/SUD_cross_species/job_run_out/scRNA_process_cerebral-%j.o
#SBATCH -e /tscc/projects/ps-palmer/brittany/SUD_cross_species/job_run_out/scRNA_process_cerebral-%j.e  
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH -A csd795    

cd /tscc/projects/ps-palmer/brittany/SUD_cross_species/scripts/

source activate env-std-py38

python /tscc/projects/ps-palmer/brittany/SUD_cross_species/scripts/export_cerebral_scrna.py
