#!/bin/bash
#SBATCH --job-name MAGMA_ext_orig
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --time 4:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 8G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/MAGMA_ext_orig-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/MAGMA_ext_orig-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account csd795


# original
prefix="ext_munged"
gene_loc_file='/tscc/projects/ps-palmer/brittany/magma_v1/NCBI38/NCBI38.gene.loc'
suffix=""
dir='/tscc/projects/ps-palmer/brittany/SUD_cross_species/magma/'
bfile_loc='/tscc/projects/ps-palmer/brittany/magma_v1/g1000_eur/g1000_eur'

/tscc/projects/ps-palmer/brittany/magma_v1/magma --annotate window=10 --snp-loc "${dir}${prefix}_pos.tsv" --gene-loc $gene_loc_file --out $dir$prefix$suffix

/tscc/projects/ps-palmer/brittany/magma_v1/magma --bfile $bfile_loc --pval "${dir}${prefix}_pval.tsv" ncol=N --gene-annot "${dir}${prefix}${suffix}.genes.annot" --out $dir$prefix$suffix

