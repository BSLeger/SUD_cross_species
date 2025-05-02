#!/bin/bash
#SBATCH --job-name tissue_network_prop_single
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --nodes 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 8:00:00
#SBATCH --mem-per-cpu 8G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/tissue_network_prop_single-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/tissue_network_prop_single-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account csd795

source activate env-std-py38

cd /tscc/projects/ps-palmer/brittany/SUD_cross_species/tissue_networks

t='amygdala'

if [ ! -f "$prefix$t.gz" ]; then
    echo "File not found- downloading from server"
	prefix=https://s3-us-west-2.amazonaws.com/humanbase/networks/
	wget "$prefix$t.gz"

fi


python  ../scripts/loco_network_prop.py $t



