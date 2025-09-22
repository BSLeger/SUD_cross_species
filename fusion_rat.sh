#!/bin/bash
#SBATCH --job-name loco_FUSION
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH -a 1-14
#SBATCH --time 2:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu 4G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/loco_FUSION-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/loco_FUSION-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account 

source activate lzenv
db_list=( 'Adipose'  'BLA'  'Brain'  'Eye'  'IL'  'LHb'  'Liver'  'NAcc'  'NAcc1'  'NAcc2'  'OFC'  'PL'  'PL1'  'PL2' )


TISSUE=${db_list[$SLURM_ARRAY_TASK_ID-1]}
VAR_TYPE='expression'

FUS_PATH='/tscc/nfs/home/bsleger/bsl/fusion_twas-master/'

cd $FUS_PATH
SUD_PATH='/tscc/projects/ps-palmer/brittany/SUD_cross_species/'
OUT_PATH="${SUD_PATH}rat_fusion/output/"
#OUT_PREF='loco_final_CF'
OUT_PREF='loco_final_mega'
DATA_FILE='loco_final/mlma_concat/regressedlr_locomotor_mega_chrgwas_zscores.mlma'
#DATA_FILE='loco_final/mlma_concat/regressedlr_gsem_results_commonfactor_F1_common_chrgwas_zscores.mlma'
FUS_REF_PATH=${SUD_PATH}'rat_fusion/'
WEIGHTS_PATH=$FUS_REF_PATH"twas-weights-rn7/"${TISSUE}

for ((CHR = 1; CHR < 21; CHR++)); do
    OUT=${OUT_PATH}${OUT_PREF}_${TISSUE}_${CHR}.dat
    echo $OUT
	if [ ! -f $OUT ]; then
	    Rscript FUSION.assoc_test.R \
	    --sumstats $SUD_PATH$DATA_FILE \
	    --weights ${WEIGHTS_PATH}/${VAR_TYPE}'.pos' \
	    --weights_dir $WEIGHTS_PATH \
	    --ref_ld_chr $FUS_REF_PATH/LDREF/Brain_rn7. \
	    --chr $CHR \
	    --out $OUT
	else
		 echo 'file already calculated, skipping to next file.'
	fi
done