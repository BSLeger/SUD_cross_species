#!/bin/bash
#SBATCH --job-name magma_hm_psych
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --nodes 1
#SBATCH -a 1-4
#SBATCH -c 4
#SBATCH -t 4:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/magma_hm_psych-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/magma_hm_psych-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account 



# Modify the paths accordingly

#file_ls=('facial_hair' 'age_smkinit' 'antisoc' 'friend_sat' 'hr' 'infant_bw' 'LDL' 'maternal_smok' 'townsend' 'age_menarche' 'neurot')
#file_ls=('age_menarche' '' 'facial_hair' 'friend_sat' 'hr' 'infant_bw' 'maternal_smok')
file_ls=('anxiety' 'panic' 'asd' 'adhd' 'scz' 'bipolar' 'dep' 'ptsd' 'ocd' 'alz' 'park' 'als' 'epilepsy' 'anorexia')
file_ls=('bipolar_euro' 'dep_euro' 'ptsd_euro' 'epilepsy_euro')
file=${file_ls[$SLURM_ARRAY_TASK_ID-1]}



cd /tscc/projects/ps-palmer/brittany/SUD_cross_species/magma

prefix="magma/"$file
gene_loc_file='/tscc/projects/ps-palmer/brittany/magma_v1/NCBI38/NCBI38.gene.loc'
#gene_loc_file='/tscc/projects/ps-palmer/brittany/magma_v1/NCBI37/NCBI37.3.gene.loc'

suffix=""
dir='/tscc/projects/ps-palmer/brittany/SUD_cross_species/gwas_hm_psych/'
bfile_loc='/tscc/projects/ps-palmer/brittany/magma_v1/g1000_eur/g1000_eur'


/tscc/projects/ps-palmer/brittany/magma_v1/magma --annotate window=10 --snp-loc "${dir}${prefix}_pos.tsv" --gene-loc $gene_loc_file --out $dir$prefix$suffix

/tscc/projects/ps-palmer/brittany/magma_v1/magma --bfile $bfile_loc --pval "${dir}${prefix}_pval.tsv" ncol=N --gene-annot "${dir}${prefix}${suffix}.genes.annot" --out $dir$prefix$suffix

