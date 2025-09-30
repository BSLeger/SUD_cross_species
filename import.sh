#!/bin/bash
#SBATCH --job-name import
#SBATCH -p condo
#SBATCH -q condo
#SBATCH --nodes 1
#SBATCH -c 2
#SBATCH -t 2:00:00
#SBATCH --mem-per-cpu 8G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/SUD_cross_species/import-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/SUD_cross_species/import-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account csd795 

cd /tscc/projects/ps-palmer/brittany/SUD_cross_species/neale_ctrl/sumstats
bash wget_commands_rerun.txt