#!/bin/bash
#SBATCH --job-name ldsc_nicsa
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --time 54:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 8G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/nicsa_ldsc-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/nicsa_ldsc-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account csd795

# Modify the paths accordingly

#location of genotype files separated by chromosome
BFILE_DIR=/tscc/nfs/home/bsleger/bsl/SUD_cross_species/nicsa_genotypes
#output file for LDSC
OUT_DIR=/tscc/nfs/home/bsleger/bsl/SUD_cross_species/nicsa_ldsc

# Path to ldsc.py
LDSC_EXEC=/tscc/nfs/home/bsleger/bsl/ldsc/ldsc.py

#activate python2 environment with necessary packages
source activate ldsc
conda list
# Loop through each chromosome
#
#test for amount of memory and time
#for CHR in {1..20}; do
for CHR in {1..20}; do
    # Define the input and output file paths
    BFILE=${BFILE_DIR}/genotypes_${CHR}
    OUT_FILE=${OUT_DIR}/ldsc_5Mb_chr${CHR}

    # Create the output directory if it doesn't exist
    mkdir -p ${OUT_DIR}

    # Execute ldsc.py for each chromosome
    python ${LDSC_EXEC} --bfile ${BFILE} \
        --l2 \
        --ld-wind-kb 5000 \
        --out ${OUT_FILE}
done
