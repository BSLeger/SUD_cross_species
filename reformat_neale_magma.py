import os
import pandas as pd
import numpy as np
import sys

os.chdir('/tscc/projects/ps-palmer/brittany/SUD_cross_species/')

bim = pd.read_csv("/tscc/projects/ps-palmer/brittany/magma_v1/g1000_eur/g1000_eur.bim", sep="\t", header=None, names=["CHR", "SNP", "CM", "POS", "A1", "A2"])

bim["CHRPOS"] = "chr" + bim["CHR"].astype(str) + ":" + bim["POS"].astype(str)

print('bim imported')
i=sys.argv[1]

print(i)
tbl=pd.read_csv(f'neale_ctrl/sumstats/{i}.tsv.bgz',compression='gzip',sep='\t',low_memory=False)

tbl['CHRPOS']='chr'+tbl['chr'].astype(str)+':'+tbl['pos'].astype(str)

print('tbl imported')
tbl=tbl[['CHRPOS','chr','pos','neglog10_pval_EUR']]
tbl=tbl.dropna()

tbl=tbl.merge(bim[["SNP","CHRPOS"]], left_on="CHRPOS", right_on="CHRPOS")

tbl['neglog10_pval_EUR']=tbl['neglog10_pval_EUR'].apply(lambda x: 10 ** (-x)).astype(float)

tbl=tbl.drop('CHRPOS',axis=1)
tbl.columns=['CHR','POS','P','SNP']

print('t formatted')
tbl[['SNP','CHR','POS']].to_csv(f'neale_ctrl/sumstats/magma/{i}_pos.tsv',index=False,sep='\t')
print('position file formatted')
tbl[['SNP','P']].to_csv(f'neale_ctrl/sumstats/magma/{i}_pval.tsv',index=False,sep='\t')
print('pval file formatted')


