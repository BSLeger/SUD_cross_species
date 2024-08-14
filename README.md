# SUD_cross_species

# purpose: 
This respository contains the analysis notebooks used studying the cross-species genetic overlap for locomotor activity and externalizing and locomotor activity, and between nicotine self-administration (NICSA) and smoking traits. This includes testing genomic SEM and PCA analysis as methods for combining signal from NICSA traits.

# reference tables
- nicsa_traits.csv
    * list of NICSA traits with their parallel human traits recommended by Hao Chen
- data_dict_nicsa_gwas_PCA.csv
    * full data dictionary of Hao Chen P50 NICSA traits. 
# Notebooks testing genomic SEM for NICSA (Hao Chen recommended traits):
- README_NICSA_gSEM.md
    * describes the pipeline for running genomic SEM- based on Montana Lara's pipeline. It describes the use of the below notebooks and scripts
- common_factor.ipynb
    * purpose: run common factor analysis- this is the notebook version with comments- HOWEVER, this step is VERY computationally intensive- to the point that it had to be broken down into severally steps. common_factor.r is the R code that runs, and common_factor.sh is the shell script for running it.
    - common_factor.r
    - common_factor.sh
- gen_sumstats.ipynb
    * purpose: in order to run a multivariate GWAS using genomicSEM, the SNPs, and corresponding SNP effects, need to be coded across phenotypes so that the same allele is the reference allele in all cases. This notebook is for this purpose.
- ldsc.ipynb
    * run LDSC on NICSA traits: output is saved in 
- ldsc_addchrtoSNP.ipynb
    * purpose: to append "chr" to the SNP column for all the LD score files. The r markdown also formats the ld score file names to be able to work with ldsc(). They must be named by the chromosome number only.   
- ldsc_bychr_array.sh
    *  purpose: calculate LD scores for each SNP using LDSC as a batch array by chromosome.
- analyze_LCSD_vs_GCTA.ipynb
    * purpose: to compare the genetic correlations and heritability calculated by LDSC bs GCTA. LDSC is used for genomicSEM, however, the results were... questionable so this notebook was for checking how they compared to GCTA.
- munge_sumstats.ipynb
    * purpose: standardize format of NISCA sumstats for running LDSC on.
- makeref.ipynb
    * purpose: make the reference file was from one of the .mlma files using the  for munging. used in the NICSA genomic SEM pipeline.
# Notebooks for PCA analysis of NICSA traits:
- PCA_phenotypes_nicsa.ipynb
    * purpose: run PCA analysis NICSA traits that Hao Chen recommended for parallels for smoking initiation and cigarettes per day
- PCA_phenotypes_run2.ipynb
    * purpose: follow up to PCA_phenotypes_nicsa. As those phenotypes did not work well in making a unified phenotype, this time I only used traits that have qtls, as well as looking at genetic correlation matrix to define phenotypes, using heirarchical clustering.
- format_PCA_forgwas.ipynb
    * purpose: format the datasets from the first run of PCA analysis on Hao Chen's recomended parellels for smokinit and cigarettes per day into datasets of phenotypes and genotypes for running GWAS using Thiago's pipeline.
- genomic_PCA.ipynb
    * purpose: to try running PCA analysis on the summary statistic level GWAS output for NICSA. This DID NOT WORK.
- nicsa_gwasreport_PCA.sh
    * purpose: shellscript for running Thiago's GWAS pipeline for the NICSA traits- probably has to be corrected to run now because the pipeline changes fairly regularly.
    - nicsa_gwasreport_secondhalf_PCA.sh: for running just the second half of the GWAS pipeline- used when the primary one failed.
    - nicsa_gwasreport_PCA_prt2.sh: for running the GWAS pipeline on the second half of the PCA analysis dataset.

# SNP to Gene mapping (TWAS, proximity)
## TWAS
- predixcan_brain_model-munged.sh
    * purpose: runscript for running S-PrediXcan on externalizing 1.0 for all GTEX brain tissues
- predixcan_brain_model.sh
    * purpose: runscript for running S-PrediXcan on externalizing 1.0 for all GTEX brain tissues
- ratxcan_locomotor.ipynb
    * purpose: to run RatXcan on lomotor meta analysis output. this notebook follows the ratXcan tutorial (jupyter notebook in the ratXcan folder). follow up analysis done in python in the TWAS_FUSION_predixcan_comparison_rat_human notebook
- ratxcan_model_predict.sh
    * batch shellscript for running the prediction portion of ratXcan on locomotor activity for all tissues. part of ratxcan_locomotor.ipynb notebook.
- ratxcan_tutorial.ipynb
    * purpose: this is the tutorial for using ratXcan, translated over from an R markdown file, written by Sabrina Mi
- QC_loco_phenotypes.ipynb
    * purpose: to check the data distributions for the locomotor datasets that were used for meta-analysis. Done primarily at Abe's recommendation to see if it's impacting the ratXcan (if not normalized)- they are in fact normalized.
- FUSION_TWAS_s-predixcan_ext.ipynb
    * purpose: run predixcan and FUSION TWAS on externalizing1.0 data
- fusion_brain_tissue_ext.sh
    * purpose: runscript for running FUSION for all brain region from GTEX. Run using unmunged FUSION, though T-SEM recommends using munged results.
- TWAS_FUSION_predixcan_comparison_human_rat.ipynb
    * purpose: analyze FUSION vs PREDIXCAN/RatXcan for externalizing/locomotor activity. The output from Daniel's TWAS (FUSION) pipeline from TWAS (as of 20 June 2024) for all tissues, RatXcan for all 5 models, and S-PrediXcan and FUSION TWAS for externalizing- for determining whether they can be used for tsem and compare to the output from S-PrediXcan
- TSEM_tutorial.ipynb
    * purpose: run through the tutorial for transcriptome-wide structural equation modeling. check it using the example human data before trying with human-rat data to see if compatible.
## Proximity
- analyze_SNPvsGene_magma.ipynb
    * purpose: to look at MAGMA S2G mapping for rats and humans and compare the output and signal 
    - manhattan plots
    - qq plots
    - histograms of pvalues
    datasets looked at:
     - rat:
         - run 1 traits:
             - locomotor activity
                 - window size optimization for locomotor (using N=8.9k, which I think actually should have been 7.7k, but I won't rerun until I have the final dataset for loco)
                 - rerun 10kb window using 7.7k rats (actualy number of phenotyped rats, not genotyped, as Hao Chen's dataset was removed)
                 - as of 1 August 2024, all results have been calculated with 8.9k rats. However, checking the seed genes (below) shows that FDR and bonferroni cutoffs did not lead to different genes based on cutoff so all network analysis does not need to be rerun.
             - PCA analysis of NICSA traits using traits recommended by hao chen
         - PCA analysis run 2 for NICSA traits using QTLs and with and all traits with QTLs
     - human:
         - externalizing1.0 (2019)
- analyze_magma_output_allrat_vs_nicsarat.ipynb
    * purpose: check how the if using all of the round 10.2 rats as the reference population for MAGMA vs using just the rats from the NICSA study substantially changes the results. All rats is ~10.2k, so it took 6x as long.
- magma_ext_munged.sh
    * purpose: runscript for running MAGMA on the externalizing 1.0 munged sumstats.
- magma_ext_orig.sh
    * purpose: runscript for running MAGMA on the externalizing 1.0 unmunged sumstats.
- magma_nicsa_qtl.sh
    * purpose: runscript for running MAGMA on the run 2 PCA analysis using QTL's traits.
- magma_sensitivity_loco.sh
    * purpose: runscript for running MAGMA on locomotor activity using different window sizes (batch)
- sensitivity_loco_singlerun.sh
    * purpose: runscript for running MAGMA on locomotor activity using different window sizes
- run_magma_human.ipynb
    * purpose: run MAGMA on the externalizing sumstats from doi: 10.1038/s41593-021-00908-3. Epub 2021 Aug 26. PMID: 34446935; PMCID: PMC8484054. They used european reference panel so i will also use european (all samples from european background). 2 versions provided- MUNGED and not MUNGED.
- run_magma_rat.ipynb
    * purpose: to run magma on the gwas output for SMKINIT and CIGDAY PC1s (from input traits recommended by Hao Chen- MCs made using PCA_phenotypes_nicsa.ipynb, format_PCA_forgwas.ipynb, and gwas pipeline) and Locomotor phenotypes
- run_snpeff.ipynb
    * purpose: run SNPEff on rat datasets, PCA run1 cigday, smkinit, and locomotor meta-analysis.
# Network Analysis
- define_seed_genes_orthology_mapping.ipynb
    * purpose: translate Rat GWAS into human orthologs (and maybe human to rat once a co-expression network made). After ortholog mapping, define seed gene sets. For ortho mapping- using bestmatch.
- loco_network_colocalization.ipynb
    * purpose: generate network colocalization values for Locomotor-Externalizing network and compare NPS scores
- loco_network_prop-sampling.ipynb
    * purpose: run network propagation for externalizing using the 500 gene subsampling that we used for BMI. It is computationally intensive, so the python script is for running it as a job, and the bash script is the job shellscript.
    - loco_network_prop-sampling.py
    - loco_network_prop-sampling.sh
- loco_network_prop.ipynb
    * purpose: run network propagation for a given dataset, in this case used for locomotor activity and externalizing.
- loco_network_single_species_specific.ipynb
    * purpose: identify rat and human only subnetworks networks for the loco-ext network and compare to the BMI paper's combined species-enriched approach.
# Post-network analysis
- validation_loco_ext.ipynb
    * purpose: identify validation datasets for Locomotor-Exploratory integration. for model organisms, use datasets that are in the alliance of genome resources. from GWAS catalog, use EFO terms. 