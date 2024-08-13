#!/bin/bash
#SBATCH --job-name ext_FUSION
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH -a 1-9
#SBATCH --time 2:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu 4G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/ext_FUSION-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/ext_FUSION-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account csd795



FUS_PATH='/tscc/nfs/home/bsleger/bsl/fusion_twas-master/'
SUD_PATH='/tscc/projects/ps-palmer/brittany/SUD_cross_species/'
#unmunged
DATA_FILE='ext_sumstat_2019/FINAL.EXT_COMMON_FACTOR.EXTERNALIZING.20191014.PREPARED.wFREQ.A1.txt'

OUT_PATH="${SUD_PATH}ext_FUSION/"
OUT_PREF="ext2019"

#munged- recommended to use
#DATA_FILE='ext_sumstat_2019/FINAL.EXT_COMMON_FACTOR.EXTERNALIZING.20191014.sumstats.txt'
#OUT_PREF="ext2019_munged"


cd "${FUS_PATH}WEIGHTS"
db_list=(`ls GTEx.Brain*.pos`)
# echo ${#db_list[*]} 
#db list ls 9 long - make job array that's 1-9

cd $FUS_PATH

source activate lzenv

TISSUE=${db_list[$SLURM_ARRAY_TASK_ID-1]}
#m=${db_list[1]}
echo $TISSUE


for ((CHR = 1; CHR < 23; CHR++));
do
    echo $CHR
    OUT=${OUT_PATH}${OUT_PREF}_${TISSUE}_${CHR}.dat
    echo $OUT
    Rscript FUSION.assoc_test.R \
    --sumstats $SUD_PATH$DATA_FILE \
    --weights "./WEIGHTS/"${TISSUE} \
    --weights_dir ./WEIGHTS/ \
    --ref_ld_chr ./LDREF/1000G.EUR. \
    --chr $CHR \
    --out $OUT
done