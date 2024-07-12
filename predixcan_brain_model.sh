#!/bin/bash
#SBATCH --job-name predixcan_ext
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --nodes 1
#SBATCH -a 1-13
#SBATCH -c 2
#SBATCH -t 2:00:00
#SBATCH --mem-per-cpu 8G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/SUD_cross_species/predixcan_ext-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/SUD_cross_species/predixcan_ext-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account csd795


cd GTEx/brain_models/
db_list=(`ls *.db`)

cd /tscc/projects/ps-palmer/brittany/MetaXcan/
source activate imlabtools

m=${db_list[$SLURM_ARRAY_TASK_ID-1]}
echo $m

software/SPrediXcan.py \
--model_db_path  "GTEx/brain_models/"$m \
--covariance GTEx/gtex_v8_expression_elastic_net_snp_smultixcan_covariance.txt.gz \
--gwas_file /tscc/projects/ps-palmer/brittany/SUD_cross_species/ext_sumstat_2019/FINAL.EXT_COMMON_FACTOR.EXTERNALIZING.20191014.PREPARED.wFREQ.A1.txt.gz \
--snp_column SNP \
--effect_allele_column A1 \
--non_effect_allele_column A2 \
--beta_column BETA.A1 \
--pvalue_column P \
--output_file "results/predixcan_externalizing2019_"${m##*/}".csv"
