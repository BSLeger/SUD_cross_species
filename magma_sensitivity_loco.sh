#!/bin/bash
#SBATCH --job-name magma_sensitivity_loco_chr10
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --nodes 1
#SBATCH -a 1-3
#SBATCH -c 4
#SBATCH -t 150:00:00
#SBATCH --mem-per-cpu 8G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/SUD_cross_species/magma_ext_chr10-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/SUD_cross_species/magma_ext_chr10-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account csd795


cd /tscc/projects/ps-palmer/brittany/SUD_cross_species/magma

# Modify the paths accordingly
win_ls=( 500 1000)
prefix="loco_meta"
gene_loc_file="rn7.2_annotatedgenes_ncbi/rn7.2_gene_attribute_table_protein_coding_forMAGMA.tsv"
N=8679
rats="loco_meta_gwas_geno/genotypes"
suffix=""

win=${win_ls[$SLURM_ARRAY_TASK_ID-1]}

/tscc/projects/ps-palmer/brittany/magma_v1/magma --annotate nonhuman window=$win --snp-loc "${prefix}_pos.tsv" --gene-loc $gene_loc_file --out $prefix$suffix"_win"$win

/tscc/projects/ps-palmer/brittany/magma_v1/magma --bfile $rats --pval "${prefix}_pval.tsv" N=$N --batch 10 chr --gene-annot  "${prefix}${suffix}"_win"${win}.genes.annot" --out $prefix$suffix"_win"$win
