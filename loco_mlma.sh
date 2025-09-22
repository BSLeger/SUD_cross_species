#!/bin/bash 
#SBATCH --job-name mlma_loco_7k_model
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --nodes 1
#SBATCH -a 1-20 
#SBATCH -c 30
#SBATCH -t 8:00:00
#SBATCH --mem-per-cpu 3Gb
#SBATCH -o /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/mlma_loco_7k_model-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/mlma_loco_7k_model-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account 


((CHR = SLURM_ARRAY_TASK_ID))

cd /tscc/projects/ps-palmer/

/tscc/projects/ps-palmer/software/local/src/gcta/gcta64 --mlma --bfile apurva/locomotor_round10.1/genotypes/genotypes.chr${CHR} --grm apurva/locomotor_round10.1/grm/${CHR}chrGRM --pheno brittany/rat_genomics_paper_pipeline_2024/data/phenotype/combined_loco_trait_7k_sorted.fam --out brittany/loco_7k/loco_7k_${CHR} --thread-num 10
