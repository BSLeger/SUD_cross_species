# Genomic SEM for NICSA

This exploratory project uses the NICSA GWAS data from the Chen lab. More detailed instructions and information for using the GenomicSem package can be found here: https://github.com/GenomicSEM/GenomicSEM
Outlined here are the steps taken for this particular project and notes specific to this dataset. Montana data in /tscc/projects/ps-palmer/mklara/genomic_sem/

## Contents

- [Genotypes](#genotypes)
- [LD scores](#ld-scores)
- [Munge GWAS](#munge-gwas)
- [LDSC](#ldsc)
- [Sumstats](#sumstats)
- [commonfactorGWAS](#commonfactorGWAS)

## Genotypes

Genotype files (.bed, .bim, .fam) for HS rats (n=629) from the NICSA project can be found here: /tscc/projects/ps-palmer/gwas/projects/p50_hao_chen_nicsa/
These genotype files contain autosomes, sex chromosomes, and mitochondria. For this analysis, only autosomes are used. Additionally, these files have to be separated into individual chromosome files. This was done using PLINK. For more detailed information on installing and using PLINK, go here: https://www.cog-genomics.org/plink/

### base code
'''bash
for chr in {1..20}; do \
plink --bfile genotypes --chr $chr --make-bed --out genotypes_${chr}; \
done
'''
### code ran
'''bash
for chr in {1..20}; do \
/tscc/projects/ps-palmer/software/local/src/plink2
 --bfile genotypes --chr $chr --make-bed --out genotypes_${chr}; \
done
'''

The output of this was updated genotype files in the plink format (.bed, .bim, .fam) separated by chromosome. Find in /tscc/nfs/home/bsleger/bsl/SUD_cross_species/nicsa_genotypes

## LD scores

Using the updated genotype files separated by chromosome in the plink format (.bed, .bim, .fam), LD scores were calculated for each SNP using LDSC. For more detailed information on installing and using ldsc, go here: https://github.com/bulik/ldsc
For this project, use script 'ldsc_bychr_array.sh' to run a job array over each chromosome. LDSC used to estimate the LD scores from the genotypes.

'''bash
bash ~/bsl/SUD_cross_species/scripts/ldsc_bychr_array.sh
'''
Of note, this only worked with python2 for me. 

ldsc_bychr_array.sh script:
#------------------------------------------------
'''bash
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

#add 2  to array number because running 3-20
#((CHR = SLURM_ARRAY_TASK_ID+2))

((CHR = SLURM_ARRAY_TASK_ID+2))
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
'''
#------------------------------------------------


## Munge GWAS

GWAS summary statistics for NICSA were munged using munge() in r. For this project, use r markdown 'munge_sumstats.ipynb'.

Files required for this step include gwas summary statistics, LD scores, and a reference file. The reference file was created from one of the .mlma files using the 'makeref.ipynb' notebook. (this actually happened in the gwas folder, and it's a little rough, but it did the trick - there is a conversation to be had about frequency, A1, and A2, but for now use the "notmaf" ref). 

The r markdown also formats the gwas and reference files to work with munge() by appending "chr" to the SNP column. Importantly, downstream functions will NOT work without this "chr". Additionally, take note of the order of the traits, which have to be consistent throughout all directories and inputs. Here, a simple alphabetical format was used in front of the trait names for all steps. 

## LDSC

To run multivariable LD-score regression, use r markdown 'ldsc.ipynb'. 

Files required for this step include munged sumstats and ld scores from steps above. Also note that 'ldsc_addchrtoSNP.ipynb' was used to append "chr" to the SNP column for all the LD score files. The r markdown also formats the ld score file names to be able to work with ldsc(). They must be named by the chromosome number only. 

The output is saved as LDSCoutput_dd.RData
Outputs:
S   estimated genetic covariance matrix

V   variance covariance matrix of the parameter estimates in S

I   matrix containing the "cross trait intercepts", or the error covariance between traits induced by overlap, in terms of subjects, between the GWASes on which the analyses are based

N   a vector contsaining the sample size (for the genetic variances) and the geometric mean of sample sizes (i.e. sqrt(N1,N2)) between two samples for the covariances

m   number of SNPs used to compute the LD scores with.

Warning message in ldsc(traits = traits, sample.prev = sample.prev, population.prev = population.prev, :
“Your genetic covariance matrix includes traits estimated to have a negative heritability.”
Genetic correlation results could not be computed due to negative heritability estimates.
LDSC finished running at 2024-03-18 14:25:59.206074
Running LDSC for all files took 23 minutes and 28 seconds
output of run saved in nicsa_LD_score_regression/gsem_ldsc_runoutput.txt
 /tscc/projects/ps-palmer/brittany/SUD_cross_species/nicsa_munge/nicsa_munge_wchr/regressedlr_nicsa_day1_infusion.sumstats.gz has negative h2
 rerun without day1_infusion- no error code and genetic correlation calculated.

## Sumstats

According to GenomicSEM wiki, in order to run a multivariate GWAS, the SNPs, and corresponding SNP effects, need to be coded across phenotypes so that the same allele is the reference allele in all cases. 

For this step, use scripts/gen_sumstats.ipynb. Files required include the gwas .mlma summary statistics with "chr" in the SNP column (created for the munge() step) and the reference file. 

The output is saved as SUMSTATSoutput.RData

## commonfactorGWAS

Next, combine the summary statistics and LDSC output and run the common factor GWAS. For this project, use markdown 'commonfactor.ipynb'. 

Files required for this step are the LDSCoutput.RData and SUMSTATSoutput.RData.

Note this is an expensive step and had to be broken up using Sarah Colbert's pipeline: https://github.com/sarahcolbert/gsemGWAS 

in directory gsemGWAS

subset sumstats:
gsemGWAS/scripts/split_sumstats.R
    #ran in R module

concat all subsets together:
awk 'FNR>1 || NR==1' *.csv > CF_sumstats.csv

for single factor (will try to combine all factors- NOT a good idea for this data idk why I did this.): 
gsemGWAS/scripts/commonfactor_GWAS.R
gsemGWAS/scripts/commonfactor_GWAS.sh
    #had to run twice, first for array 1-600, then for 601-1047

## comparison of LDSC and GCTA genetic correlation and heritability
analyze_LDSC_vs_GCTA.ipynb - our LDSC findings are horrendous

FROM AUTHOR:
    I have not used LDSC (or Genomic SEM) with nonhuman data, so I can't provide much guidance on that front. I assume you are following best practices for estimating rG in HS rats and have an appropriate set of lds scores. the Ns certainly seem small by human gwas standards so that is a concern.
Below is how I would go about inspecting the ldsc output to see whether it is sensible to move forward with a Genomic SEM model with these data.
we can get a sense of the rGs by converting the S matrix to an R matrix using cov2cor (you could have requrested standardized ldsc output when you ran it, but it appears that you didn't, so we'll create the R matrix after the fact):
cov2cor(LDSCoutput$S)
     p50_david_dietz_loco_distance p50_hao_chen_open_field_first_15 p50_shelly_flagel_2014_d1_total_dist
[1,]                     1.0000000                         47.10806                            0.4141109
[2,]                    47.1080599                          1.00000                           60.4659724
[3,]                     0.4141109                         60.46597                            1.0000000
[4,]                     1.4610404                         97.90025                            1.1392513
[5,]                     1.2685243                        160.74589                            1.1426221
[6,]                     0.7605677                         37.00778                            0.8759522
     u01_peter_kalivas_oft_distance_1_first_15 u01_suzanne_mitchell_locomotor_t1_total_distance u01_tom_jhou_locomotor1
[1,]                                 1.4610404                                        1.2685243               0.7605677
[2,]                                97.9002515                                      160.7458944              37.0077836
[3,]                                 1.1392513                                        1.1426221               0.8759522
[4,]                                 1.0000000                                        2.5493044               0.3537181
[5,]                                 2.5493044                                        1.0000000               0.2220304
[6,]                                 0.3537181                                        0.2220304               1.0000000
There appear to be many out-of-bounds estimates, which can stem from very low h2s and/or low power.
The h2s are the diag of S:
diag(LDSCoutput$S)
[1] 0.21967473984 0.00002519852 0.15855032707 0.02835159616 0.07015541286 0.11715178105
The 2nd and 4th variable are particularly concerning to me here. They also appear to be the biggest offenders with respect to rGs>>1.
Let's look at the SEs of those h2s...
k<-nrow(LDSCoutput$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSCoutput$V))
diag(SE)
[1] 0.05769571 0.05163793 0.06338283 0.08752559 0.14561862 0.09478580
Looks like some are rather large, indicating possibly low power. Let's see whether the h2s are significant by inspecting their Z stats.
diag(LDSCoutput$S)/diag(SE)
[1] 3.8074709510 0.0004879847 2.5014714039 0.3239234991 0.4817750062 1.2359634498
Looks like only phenotypes 1 and 3 have significant h2. So I think power is certainly an issue here. It's hard to model genetic covariance structure if you can't detect significant genetic signal to start with.
Finally, let's look at the LDSC intercepts.
diag(LDSCoutput$I)
[1] 1.180277 1.366183 1.077194 1.104384 1.167949 1.060990
Looks like there is some inflation for variable 2 in particular. But many of the others are still a good bit above one too...
To conclude, there are certainly some red flags with respect to the data. It's possible that the implementation of LDSC in HS rats was suboptimal, but also possible that the N just isn't large enough (or that some of the phonetypes aren't appreciably heritable (which can happen if they are measured with lots of error). One thing that could be tried to increase power is HDL. Michel provides a tutorial on that on our github wiki. However, this would require using appropriate LD matrices for the HS rats.
I hope that helps.
Elliot



# Genetic Correlation Methods------------------------

we use REML, which uses individual-level data rather than sumstats that LDSC and HDL use

HDL and LDSC both recommended methods for implementation with genomic SEM.

bash(f'''{self.gcta} --reml-bivar {d_[trait1]} {d_[trait2]} {self.thrflag} \
                --grm {self.autoGRM} --pheno {self.path}data/allpheno.txt --reml-maxit 200 \
                --reml-bivar-lrt-rg 0 --out {self.path}temp/rG/gencorr.temp{randomid}''', print_call=False)




 '''
        Generates a genetic correlation matrix using GCTA.
    
        Parameters
        ----------
        traitlist : list, optional
            List of traits for which genetic correlations will be computed.
            If not provided, the function uses the traits from the object's attribute.
        
        print_call : bool, optional
            Whether the command line call for each trait is printed.
        
        Returns
        -------
        pd.DataFrame
            A DataFrame representing the genetic correlation matrix.
    
        Design
        ------
        This function performs the following steps:
        1. Checks for or creates a Dask client for parallel computation.
        2. Prepares data and directories.
        3. Computes genetic correlations and phenotypic correlations for all trait pairs.
        4. Processes the results and creates a melted table.
        5. Saves the genetic correlation matrix with hierarchical clustering.
        6. Adds heritability information to the matrix.
        7. Generates and saves a clustered heatmap using Seaborn.
        8. Returns the final genetic correlation matrix.
        '''

