import os
import pandas as pd
import numpy as np
import sys

os.chdir('/tscc/projects/ps-palmer/brittany/SUD_cross_species/')

i=sys.argv[1]

bim = pd.read_csv("/tscc/projects/ps-palmer/brittany/magma_v1/g1000_eur/g1000_eur.bim", sep="\t", header=None,names=["CHR", "SNP", "CM", "POS", "A1", "A2"])

bim["CHRPOS"] = "chr" + bim["CHR"].astype(str) + ":" + bim["POS"].astype(str)

tbl=pd.read_csv(f'neale_ctrl/sumstats/{i}.tsv.bgz',compression='gzip',sep='\t',low_memory=False)

tbl['CHRPOS']='chr'+tbl['chr'].astype(str)+':'+tbl['pos'].astype(str)
tbl=tbl[['chr','pos','ref','alt','beta_EUR','CHRPOS','neglog10_pval_EUR','af_EUR']]
tbl=tbl.merge(bim[["SNP","CHRPOS"]], left_on="CHRPOS", right_on="CHRPOS")
tbl['neglog10_pval_EUR']=tbl['neglog10_pval_EUR'].apply(lambda x: 10 ** (-x)).astype(float)


tbl.rename(columns={'neglog10_pval_EUR': 'P_EURO'},inplace=True)

tbl=tbl.drop('CHRPOS',axis=1)

tbl=tbl[['SNP','chr','pos','alt','ref','af_EUR','beta_EUR','P_EURO']]
print(tbl.head(5))
tbl.to_csv(f'/tscc/projects/ps-palmer/brittany/SUD_cross_species/munged_sumstats/{i}_preprocessed.tsv.gz',index=False,sep='\t',compression='gzip')
print(f'file written to /munged_sumstats/{i}_preprocessed.tsv.gz ')
