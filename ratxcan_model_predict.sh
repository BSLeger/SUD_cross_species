#!/bin/bash 
#SBATCH --job-name ratxcan_loco_7k_model
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --nodes 1
#SBATCH -a 1-4 
#SBATCH -c 4
#SBATCH -t 4:00:00
#SBATCH --mem-per-cpu 4G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/ratxcan_loco_7k_model-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/ratxcan_loco_7k_model-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account 


#dependent of the dataset
GENO_PREFIX="/tscc/projects/ps-palmer/brittany/SUD_cross_species/magma/loco_meta_gwas_geno/genotypes_7k"
GENO="/tscc/projects/ps-palmer/brittany/SUD_cross_species/magma/loco_meta_gwas_geno"
PREDICT_PREF="loco_meta_7k"
PHENOTYPES="$PRE/data/phenotype/combined_loco_trait_7k.fam"

PRE="/tscc/projects/ps-palmer/brittany/rat_genomics_paper_pipeline_2024"
OUTPUT="$PRE/output"
METAXCAN="$PRE/MetaXcan"
MODEL="$PRE/models"
model_PRE=( "VO" "PL" "LH" "IL" )
# full list model_PRE=( "VO" "PL" "LH" "IL" "AC" )
# already ran AC, so skip- incraese SBATCH array to 1-5 if run all


cd $PRE
source activate imlabtools

m=${model_PRE[$SLURM_ARRAY_TASK_ID-1]}


#locomotor- need to use liftover to map back to rn6 which was used for to develop the model
source activate imlabtools

python ${METAXCAN}/software/Predict.py \
--model_db_path ${MODEL}/${m}-filtered.db \
--model_db_snp_key rsid \
--liftover rn7ToRn6.over.chain.gz \
--vcf_genotypes ${GENO_PREFIX}.vcf.gz \
--vcf_mode genotyped \
--on_the_fly_mapping METADATA "{}_{}_{}_{}" \
--prediction_output $OUTPUT/${m}-filtered-${PREDICT_PREF}_predict.txt  \
--prediction_summary_output $OUTPUT/${m}-filtered-${PREDICT_PREF}_summary.txt \
--throw