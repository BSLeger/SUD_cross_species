#!/bin/bash 
#SBATCH -J nicsa_gwasreport_PCA_run2
#SBATCH -p platinum
#SBATCH -q hcp-csd795
#SBATCH -t 30:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 30
#SBATCH --mem-per-cpu 3Gb
#SBATCH --export ALL
#SBATCH -o /tscc/projects/ps-palmer/brittany/SUD_cross_species/job_run_out/nicsa_gwasreport_PCA_run2-%j.o
#SBATCH -e /tscc/projects/ps-palmer/brittany/SUD_cross_species/job_run_out/nicsa_gwasreport_PCA_run2-%j.e  
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH -A csd795    

source activate gwas
cd /tscc/projects/ps-palmer/gwas/GWAS-pipeline
#gwasVersion=$(git describe --tags)
#echo $gwasVersion

python /tscc/projects/ps-palmer/gwas/GWAS-pipeline/gwas_cli.py path=/tscc/projects/ps-palmer/brittany/SUD_cross_species/  genotypes=/tscc/projects/ps-palmer/gwas/databases/rounds/r10.2.1 n_autosome=20 phewas_path=/tscc/projects/ps-palmer/gwas/projects/phewasdb_rn7_g102.parquet.gz round=10.2.1 genome=rn7 founder_genotypes=/tscc/projects/ps-palmer/gwas/databases/founder_genotypes/founders7.2 project=nicsa_gwas threads=30 threshold=5.38  species="rattus norvegicus" subset h2 latent_space clear_directories gwas regressout qtl phewas grm goea eqtl gcorr locuszoom sqtl report gwas_version=v0.2.0-16-gc12000f


#    removed 