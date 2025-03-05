#!/bin/bash
#SBATCH --job-name tissue_network_format
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --nodes 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 4:00:00
#SBATCH --mem-per-cpu 8G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/tissue_network_format-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/tissue_network_format-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account csd795

source activate env-std-py38

cd /tscc/projects/ps-palmer/brittany/SUD_cross_species/tissue_networks

t='global'

#prefix=https://s3-us-west-2.amazonaws.com/humanbase/networks/

#wget "$prefix$t.gz"

python ../scripts/format_tissue_network.py $t



