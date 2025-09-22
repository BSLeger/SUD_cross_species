#!/bin/bash
#SBATCH --job-name gwas_cat
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --time 6:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 2
#SBATCH --mem-per-cpu 32G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/gwas_cat-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/gwas_cat-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account csd795

source activate env-std-py38

cd /tscc/projects/ps-palmer/brittany/SUD_cross_species/scripts

python EFO_analysis_GWAS_cat.py
