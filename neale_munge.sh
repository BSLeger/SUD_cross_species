#!/bin/bash
#SBATCH --job-name neale_munge
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --nodes 1
#SBATCH -a 2-65
#SBATCH -c 4
#SBATCH -t 1:30:00
#SBATCH --mem-per-cpu 16G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/neale_munge-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/neale_munge-%j.e
#SBATCH --mail-type END,neal_netprop
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account csd795 

#export PATH="/tscc/nfs/home/bsleger/miniconda3/bin:$PATH"
conda activate 
INDIR=/tscc/projects/ps-palmer/brittany/SUD_cross_species/neale_ctrl/sumstats
OUTDIR=/tscc/projects/ps-palmer/brittany/SUD_cross_species/munged_sumstats

traits=( 'biomarkers-30610-both_sexes-irnt' 'biomarkers-30670-both_sexes-irnt' 'biomarkers-30700-both_sexes-irnt' 'biomarkers-30710-both_sexes-irnt' 'biomarkers-30720-both_sexes-irnt' 'biomarkers-30750-both_sexes-irnt' 'biomarkers-30760-both_sexes-irnt' 'biomarkers-30770-both_sexes-irnt' 'biomarkers-30810-both_sexes-irnt' 'categorical-1747-both_sexes-5' 'categorical-20002-both_sexes-1065' 'categorical-20002-both_sexes-1111' 'categorical-20002-both_sexes-1226' 'categorical-20002-both_sexes-1474' 'categorical-20003-both_sexes-1141194794' 'categorical-20116-both_sexes-0' 'categorical-20116-both_sexes-2' 'categorical-20160-both_sexes-20160' 'categorical-4293-both_sexes-3' 'categorical-4294-both_sexes-0' 'continuous-102-both_sexes-irnt' 'continuous-20016-both_sexes-irnt' 'continuous-20153-both_sexes-irnt' 'continuous-21002-both_sexes-irnt' 'continuous-2178-both_sexes' 'continuous-23100-both_sexes-irnt' 'continuous-23101-both_sexes-irnt' 'continuous-23106-both_sexes-irnt' 'continuous-23116-both_sexes-irnt' 'continuous-30000-both_sexes-irnt' 'continuous-30010-both_sexes-irnt' 'continuous-30040-both_sexes-irnt' 'continuous-30070-both_sexes-irnt' 'continuous-30080-both_sexes-irnt' 'continuous-30100-both_sexes-irnt' 'continuous-30120-both_sexes-irnt' 'continuous-30130-both_sexes-irnt' 'continuous-30150-both_sexes-irnt' 'continuous-30180-both_sexes-irnt' 'continuous-30200-both_sexes-irnt' 'continuous-30300-both_sexes-irnt' 'continuous-3143-both_sexes-irnt' 'continuous-3148-both_sexes-irnt' 'continuous-4104-both_sexes-irnt' 'continuous-50-both_sexes-irnt' 'continuous-5134-both_sexes-irnt' 'continuous-5257-both_sexes-irnt' 'continuous-AG-both_sexes-irnt' 'continuous-LDLC-both_sexes-medadj_irnt' 'continuous-LDLC-both_sexes-medadj_raw' 'continuous-MAP-both_sexes-manual_medadj_raw' 'continuous-NAP-both_sexes-irnt' 'continuous-PP-both_sexes-combined_medadj_irnt' 'continuous-PP-both_sexes-combined_medadj_raw' 'continuous-SBP-both_sexes-combined_medadj_irnt' 'continuous-SBP-both_sexes-combined_medadj_raw' 'continuous-eGFRcreacys-both_sexes-irnt' 'icd10-E66-both_sexes' 'icd10-M17-both_sexes' 'phecode-250.2-both_sexes' 'phecode-327.3-both_sexes' 'phecode-411.4-both_sexes' 'phecode-496.21-both_sexes' 'phecode-594.1-both_sexes' 'phecode-714-both_sexes' )

Ns=( 400988 400687 400761 400094 400940 400825 367021 398797 366484 419469 420473 420473 420473 420473 420473 418817 418817 418860 137837 137999 396667 135088 137532 419316 418781 412524 413158 413153 413135 407990 407995 407993 407993 407992 407987 407265 407265 407265 407282 407282 401345 242926 242809 132165 419596 74603 87199 368642 398402 398402 35367 368642 417001 417001 417001 417001 401570 420531 420531 418949 419529 405940 391193 414965 377454 )

i=$SLURM_ARRAY_TASK_ID-1
trait=${traits[$i]}

N=${Ns[$i]}

echo "preprocessing $trait with N=$N"

source activate env-std-py38
cd /tscc/projects/ps-palmer/brittany/SUD_cross_species/scripts
python preprocess_munge_neale.py $trait



echo "preprocessing complete"
conda activate ldsc

cd /tscc/projects/ps-palmer/brittany/ldsc

python munge_sumstats.py \
--signed-sumstats beta_EUR,0 \
--out /tscc/projects/ps-palmer/brittany/SUD_cross_species/munged_sumstats/${trait} \
--merge-alleles w_hm3.snplist \
--N $N \
--a1 alt \
--a2 ref \
--snp SNP \
--frq af_EUR \
--sumstats /tscc/projects/ps-palmer/brittany/SUD_cross_species/munged_sumstats/${trait}_preprocessed.tsv.gz  \
--p P_EURO
