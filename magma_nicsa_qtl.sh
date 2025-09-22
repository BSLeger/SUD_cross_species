#!/bin/bash
#SBATCH --job-name nicsa_qtl
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --nodes 1
#SBATCH -a 1-37
#SBATCH -c 1
#SBATCH -t 8:00:00
#SBATCH --mem-per-cpu 8G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/SUD_cross_species/magma_nicsa-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/SUD_cross_species/magma_nicsa-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account 


cd /tscc/projects/ps-palmer/brittany/SUD_cross_species/magma

#just the rats in the dataset used
prefixes=('nicsa_day8_infusion' 'nicsa_day6_activelick' 'nicsa_day2_activelick' 'nicsa_total_activelick_10days' 'nicsa_day11_activelick' 'nicsa_total_inactivelick_10days' 'nicsa_day4_inactivelick' 'nicsa_day8_activelick' 'nicsa_day5_infusion' 'nicsa_active_inactive_ratio_all_days' 'nicsa_day10_inactivelick' 'nicsa_day3_active_inactive_ratio' 'nicsa_day6_infusion' 'nicsa_day3_inactivelick' 'nicsa_day1_inactivelick' 'nicsa_slope_nicotine_infusion' 'nicsa_day4_infusion' 'nicsa_first_three_days_infusion_median' 'nicsa_day1_infusion' 'nicsa_day8_active_inactive_ratio' 'nicsa_active_inactive_ratio_last_three_median' 'nicsa_day6_active_inactive_ratio' 'nicsa_day4_active_inactive_ratio' 'nicsa_day3_activelick' 'nicsa_day11_infusion' 'nicsa_last_three_days_infusion_median' 'nicsa_total_infusion_10days' 'nicsa_day11_active_inactive_ratio' 'nicsa_day2_infusion' 'nicsa_day9_inactivelick' 'nicsa_day1_activelick' 'nicsa_day4_activelick' 'pc1_smkinit_qtl' 'pc1_cigday1a_qtl' 'pc1_smkinit1a_qtl' 'pc1_cigday2a_qtl' 'pc1_cigday_qtl')
Ns=(1785 1783 1804 1811 1769 1811 1802 1785 1802 1811 1811 1774 1784 1804 1811 1811 1794 1811 1811 1738 1811 1746 1766 1799 1771 1811 1811 1624 1805 1799 1811 1794 1811 1750 1811 1750 1750)
gene_loc_file="rn7.2_annotatedgenes_ncbi/rn7.2_gene_attribute_table_protein_coding_forMAGMA.tsv"
nicsa_rats="nicsa_gwas_geno/genotypes"
suffix="_nicsa_geno"    
prefix=${prefixes[$SLURM_ARRAY_TASK_ID-1]}
N=${Ns[$SLURM_ARRAY_TASK_ID-1]}


/tscc/projects/ps-palmer/brittany/magma_v1/magma --annotate nonhuman window=10 --snp-loc "${prefix}_pos.tsv" --gene-loc $gene_loc_file --out $prefix$suffix"_win"$win

/tscc/projects/ps-palmer/brittany/magma_v1/magma --bfile $nicsa_rats --pval "${prefix}_pval.tsv" N=$N --gene-annot "${prefix}${suffix}"_win"${win}.genes.annot" --out $prefix$suffix"_win"$win
