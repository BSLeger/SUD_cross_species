#!/bin/bash
#SBATCH --job-name ldsc_nicsa
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --nodes 1
#SBATCH -a 1-20
#SBATCH -c 1
#SBATCH -t 4:00:00
#SBATCH --mem-per-cpu 4G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/nicsa_ldsc-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/nicsa_ldsc-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account 

# Modify the paths accordingly

#location of genotype files separated by chromosome
BFILE_DIR=/tscc/nfs/home/bsleger/bsl/SUD_cross_species/nicsa_genotypes
#output file for LDSC
OUT_DIR=/tscc/nfs/home/bsleger/bsl/SUD_cross_species/nicsa_ldsc

# Path to ldsc.py
LDSC_EXEC=/tscc/nfs/home/bsleger/bsl/ldsc/ldsc.py

#activate python2 environment with necessary packages
source activate ldsc

#add 2  to array number because running 3-20
#((CHR = SLURM_ARRAY_TASK_ID+2))

((CHR = SLURM_ARRAY_TASK_ID))
echo $CHR
# Define the input and output file paths
BFILE=${BFILE_DIR}/genotypes_${CHR}
OUT_FILE=${OUT_DIR}/ldsc_5Mb_chr${CHR}
echo $BFILE
echo $OUT_FILE
# Create the output directory if it doesn't exist
mkdir -p ${OUT_DIR}

# Execute ldsc.py for each chromosome
python ${LDSC_EXEC} --bfile ${BFILE} \
    --l2 \
    --ld-wind-kb 5000 \
    --out ${OUT_FILE}

conda deactivate
