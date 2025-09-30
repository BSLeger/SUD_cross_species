#!/bin/bash
#SBATCH --job-name magma_hm_ctrl
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --nodes 1
#SBATCH -a 1-12
#SBATCH -c 4
#SBATCH -t 8:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/magma_hm_ctrl-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/magma_hm_ctrl-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account csd795 


cd /tscc/projects/ps-palmer/brittany/SUD_cross_species/magma

# Modify the paths accordingly

#SLURM_ARRAY_TASK_ID=1

file_ls=( continuous-30040-both_sexes-irnt continuous-30100-both_sexes-irnt continuous-50-both_sexes-irnt continuous-30120-both_sexes-irnt continuous-30010-both_sexes-irnt continuous-30130-both_sexes-irnt continuous-30080-both_sexes-irnt phecode-594.1-both_sexes phecode-714-both_sexes categorical-4293-both_sexes-3 continuous-30070-both_sexes-irnt categorical-20003-both_sexes-1141194794 )

#file=${file_ls[$SLURM_ARRAY_TASK_ID-1]}

sample_size_ls=( 420473 137837 407995 407993 407993 407992 407987 407265 407265 419596 414965 377454 )

prefix=${file_ls[$SLURM_ARRAY_TASK_ID-1]}
gene_loc_file='/tscc/projects/ps-palmer/brittany/magma_v1/NCBI38/NCBI38.gene.loc'
suffix=""
dir='/tscc/projects/ps-palmer/brittany/SUD_cross_species/neale_ctrl/sumstats/magma/'
bfile_loc='/tscc/projects/ps-palmer/brittany/magma_v1/g1000_eur/g1000_eur'

source activate env-std-py38 

/tscc/projects/ps-palmer/brittany/magma_v1/magma --annotate window=10 --snp-loc "${dir}${prefix}_pos.tsv"  --gene-loc $gene_loc_file --out ${dir}${prefix}

/tscc/projects/ps-palmer/brittany/magma_v1/magma --bfile $bfile_loc --pval "${dir}${prefix}_pval.tsv" N=${sample_size_ls[$SLURM_ARRAY_TASK_ID-1]} --gene-annot "${dir}${prefix}.genes.annot" --out ${dir}${prefix}
