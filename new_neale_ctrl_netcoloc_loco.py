import os
import pandas as pd
import ndex2
import networkx as nx
from netcoloc import netprop_zscore
from netcoloc import netprop
from netcoloc import network_colocalization
import sys
import random

os.chdir('/tscc/projects/ps-palmer/brittany/SUD_cross_species/scripts')
from network_functions import *
from plotting_functions import *
os.chdir('/tscc/projects/ps-palmer/brittany/SUD_cross_species/')
random_seed=random.seed(211)

save_fig=True

trait_r= sys.argv[1]
cut_r=sys.argv[2]
print(f'{trait_r}_{cut_r}')

#create a file called environ_ndex_meta.py where you save variables 'ndex_user' and 'ndex_password'
#otherwise will prompt you to define those within the notebooks
if os.path.isfile('/tscc/projects/ps-palmer/brittany/environ_ndex_meta.py'):
    print ('NDEx credentials imported from meta file')
    sys.path.insert(1, '../')
    from environ_ndex_meta import *
    sys.path.pop(1)
else:
    # Prompt the user for a username
    ndex_user = input("Enter your NDEx username: ")
    # Prompt the user for a password
    ndex_password = input("Enter your NDEx password: ")


plt.rcParams.update({'font.size': 16})

interactome_name='PCNet2.0'
interactome=import_interactome(UUIDs=UUIDs,interactome_name=interactome_name)
all_nodes=list(interactome.nodes())
# pre calculate the matricies used for network propagation
print('\ncalculating w_prime')
w_prime = netprop.get_normalized_adjacency_matrix(interactome, conserve_heat=True)

print('\ncalculating w_double_prime')
w_double_prime = netprop.get_individual_heats_matrix(w_prime, .5)

magma=True

seed_dict=import_seed_dict(mag_dir,file_dict,ctrl_traits,ctrl_traits_rat,psych_traits,bonf_dict,gene_col_dict,all_nodes)

NPS_dict,NPS_dict_series=import_NPS_scores(seed_dict,interactome_name)

ls=[w for w in list(seed_dict.keys()) if any(label in w for label in psych_traits)]
print(len(ls))

zlist=[cut_comb]
z12list=[cut_single]

trait_h=None
cut_h=None
overwrite=False
#for cut_r in ['FDR','bonf','top500']:
_,label_r,_,seed_r,_,NPS_r,_=return_analysis_datasets(trait_r,cut_r,trait_h,cut_h,seed_dict,NPS_dict,interactome_name)
print(trait_r)
for label_h in ls:
		print(label_h)
		coloc_filename=f'colocalization_scores/colocScore_{label_r}_{label_h}_{interactome_name}.tsv'
		if not (os.path.exists(coloc_filename)and overwrite==False):
			print('running analysis')
			seed_h=seed_dict[label_h]
			print(f'{len(seed_h)} seed genes for this dataset.')
			npsh_label=label_h+'_'+interactome_name
			if not (npsh_label in NPS_dict.keys()):
				print(f'{npsh_label} does not have scores, possibly due to number of seed genes.')
			else:
				NPS_h=NPS_dict[npsh_label]
				netcoloc_enrichment_df = network_colocalization.calculate_network_enrichment(NPS_r,NPS_h,
																							 zthresh_list = zlist,
																							 z12thresh_list=z12list,
																							 verbose=False)
				#netcoloc_enrichment_df=netcoloc_enrichment_df[netcoloc_enrichment_df['z_comb']>=netcoloc_enrichment_df['NPS_single']]
				#print(netcoloc_enrichment_df)
				netcoloc_enrichment_df['rat_dataset']=label_r
				netcoloc_enrichment_df['human_dataset']=label_h
				if save_fig:
					netcoloc_enrichment_df.to_csv('colocalization_scores/colocScore_'+label_r+'_'+label_h+'_'+interactome_name+'.tsv',sep='\t',index=False)
		else:
			print('file already exists')