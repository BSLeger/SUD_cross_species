import os
import h5py
import scanpy as sc
import pandas as pd
import numpy as np
import json
import math
import nexusformat.nexus as nx
import anndata as ad
os.chdir('/tscc/projects/ps-palmer/brittany/SUD_cross_species/')

print('analyzing Cerebral cortex')
roi='Cerebral cortex'
adata=sc.read_h5ad(f"scRNA_seq/{roi.replace(' ','-')}_combined.h5ad",backed='r')


groups = adata.obs['group'].unique()
groups=set(groups)

for group in groups:
	print(f'analyzing {group}')		
	group_data = adata[adata.obs['group'] == group, :].to_memory()
	sc.pp.log1p(group_data)

	group_data=group_data.to_df()
	if (group!='Cerebral-cortex_Neuron'):
		print('exporting mean')
		group_data.mean().T.to_csv(f'scRNA_seq/processed/{group}_pseudobulk_mean.csv')
		print('exporting sum')
		group_data.sum().T.to_csv(f'scRNA_seq/processed/{group}_pseudobulk_sum.csv')
		print('exporting stdev')
		group_data.std().T.to_csv(f'scRNA_seq/processed/{group}_pseudobulk_stdev.csv')
		print('exporting count')
		group_data.count().T.to_csv(f'scRNA_seq/processed/{group}_pseudobulk_count.csv')
		print('exporting non-zero count')
		group_data.apply(lambda df: (df>0).sum()).T.to_csv(f'scRNA_seq/processed/{group}_pseudobulk_count_nonzero.csv')
	else:
		print('exporting non-zero count')
		group_data.apply(lambda df: (df>0).sum()).T.to_csv(f'scRNA_seq/processed/{group}_pseudobulk_count_nonzero.csv')

	del(group_data)
del(adata)
print('analyzing Cerebral nuclei')
roi='Cerebral nuclei'

adata=sc.read_h5ad(f"scRNA_seq/{roi.replace(' ','-')}_combined.h5ad")
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_counts=10)
sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=False)
sc.pp.log1p(adata)
print('writing roi data')
adata.to_df().groupby(adata.obs['group']).std().T.to_csv(f'scRNA_seq/processed/{roi}_pseudobulk_stdev.csv')
adata.to_df().groupby(adata.obs['group']).apply(lambda df: (df > 0).sum()).T.to_csv(f'scRNA_seq/processed/{roi}_pseudobulk_count_nonzero.csv')
adata.to_df().groupby(adata.obs['group']).count().T.to_csv(f'scRNA_seq/processed/{roi}_pseudobulk_count.csv')



