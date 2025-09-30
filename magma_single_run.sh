#!/bin/bash
#SBATCH --job-name magma
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --time 6:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 8G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/magma-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/magma-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account csd795 


# original
#single run- alzheimer's- uses old human genome

cd /tscc/projects/ps-palmer/brittany/SUD_cross_species/magma

file='pau'
n=903147


prefix="magma/"$file
gene_loc_file='/tscc/projects/ps-palmer/brittany/magma_v1/NCBI38/NCBI38.gene.loc'
#gene_loc_file='/tscc/projects/ps-palmer/brittany/magma_v1/NCBI37/NCBI37.3.gene.loc'

suffix=""
dir='/tscc/projects/ps-palmer/brittany/SUD_cross_species/gwas_ctrl_hm/'
bfile_loc='/tscc/projects/ps-palmer/brittany/magma_v1/g1000_eur/g1000_eur'



/tscc/projects/ps-palmer/brittany/magma_v1/magma --annotate window=10 --snp-loc "${dir}${prefix}_pos.tsv" --gene-loc $gene_loc_file --out $dir$prefix$suffix


/tscc/projects/ps-palmer/brittany/magma_v1/magma --bfile $bfile_loc --pval "${dir}${prefix}_pval.tsv" N=$n --gene-annot "${dir}${prefix}${suffix}.genes.annot" --out $dir$prefix$suffix
