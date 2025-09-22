#!/bin/bash
#SBATCH --job-name nicsa_common_factor
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --time 8:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 8G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/nicsa_common_factor-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/nicsa_common_factor-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account csd795

# Modify the paths accordingly

#activate python2 environment with necessary packages
source activate lzenv
conda list
# Loop through each chromosome
#
#test for amount of memory and time
#for CHR in {1..20}; do
Rscript common_factor.r

