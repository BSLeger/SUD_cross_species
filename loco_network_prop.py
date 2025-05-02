#!/usr/bin/env python
# coding: utf-8

# purpose: run network propagation for a given dataset, in this case used for locomotor activity and externalizing.
import os
import pandas as pd
import ndex2
import networkx as nx
from netcoloc import netprop_zscore
from netcoloc import netprop
from netcoloc import network_colocalization
import sys
import random
#os.chdir('/tscc/projects/ps-palmer/brittany/rare_common_alcohol/rare_common_alcohol_comparison/notebooks/')
#from rca_functions import *
os.chdir('/tscc/projects/ps-palmer/brittany/SUD_cross_species/scripts')
from network_functions import *
#from network_validation_functions import *
from plotting_functions import *
os.chdir('/tscc/projects/ps-palmer/brittany/SUD_cross_species/')

random_seed=random.seed(211)

save_file=True

#create a file called environ_ndex_meta.py where you save variables 'ndex_user' and 'ndex_password'
#otherwise will prompt you to define those within the notebooks
if os.path.isfile('../environ_ndex_meta.py'):
    print ('NDEx credentials imported from meta file')
    sys.path.insert(1, '../')
    from environ_ndex_meta import *
    sys.path.pop(1)
else:
    # Prompt the user for a username
    ndex_user = input("Enter your NDEx username: ")
    # Prompt the user for a password
    ndex_password = input("Enter your NDEx password: ")


# # Interactome Set-up

# pcnet2- versions 
# from wright et al. 2024 preprint:
# PCNet 2.0= best-performing ranked composite (top 15 interactomes, 3.85M interactions)
# PCNet 2.1= top 8 interactomes, 1.75M interactions
# PCNet 2.2= top 10 co-citation-free interactomes, 3.32M interactions 


tissue_network=True

if tissue_network:
    #tissue='global'
	tissue=sys.argv[1]
else:
    interactome_name='PCNet2.0'


if not tissue_network:
    interactome=import_interactome(UUIDs=UUIDs,interactome_name=interactome_name)
    all_nodes=list(interactome.nodes())
    # pre calculate the matricies used for network propagation
    print('\ncalculating w_prime')
    w_prime = netprop.get_normalized_adjacency_matrix(interactome, conserve_heat=True)
    print('\ncalculating w_double_prime')
    w_double_prime = netprop.get_individual_heats_matrix(w_prime, .5)
    edges=list(interactome.edges())
    all_nodes=list(interactome.nodes())
    degree=interactome.degree()
else:
    netdir='tissue_networks/intermediate/'
    
    print(f'available files for this tissue: {[x for x in os.listdir(netdir) if tissue in x]}')
    interactome_name=f'hb_tissue_{tissue}_top'
    #import node list
    with open(f'{netdir}node_list_{tissue}_top.txt', 'r') as file:
        lines = file.readlines()
    # Remove newline characters from each line
    all_nodes=[line.strip() for line in lines]
    #import degrees
    degree=pd.read_csv(f'{netdir}degree_{tissue}_top.csv', header=None)
    degree.index=degree[0].astype(str)
    degree=degree[1].to_dict()

    print(f'importing files for {tissue}_top')
    #import w_double_prime
    if os.path.exists(f'{netdir}w_double_prime_{tissue}_top.npy'):
        print('importing w_double_prime from file')
        w_double_prime=np.load(f'{netdir}w_double_prime_{tissue}_top.npy')
    elif os.path.exists(f'{netdir}normalized_adjacency_{tissue}_top.npz'):
        print('w_double_prime not calculated previously- importing w_prime, from which w_double_prime will be calculated')
        w_prime=sp.load_npz(f'{netdir}normalized_adjacency_{tissue}.npz')
        w_double_prime=netprop.get_individual_heats_matrix(w_prime, .5)
        np.save(f'{outdir}w_double_prime_{tissue}',w_double_prime)
    else:
        print('files not found- interactome files must be calculated before continuing')


# # calculate gwas NPS

if not tissue_network:
    seed_dict=import_seed_dict(mag_dir,file_dict,bonf_dict,gene_col_dict,all_nodes)
else:
    hgnc=pd.read_csv('hgnc_complete_set.txt',sep='\t',low_memory=False)
    hgnc=hgnc[['symbol','entrez_id']].dropna()
    hgnc['entrez_id']=hgnc['entrez_id'].astype(int).astype(str)
    seed_dict=import_seed_dict(mag_dir,file_dict,bonf_dict,gene_col_dict,hgnc[hgnc.entrez_id.isin(all_nodes)]['symbol']) 
seed_dict.keys()


#dictionary of human control traits
ctrl_dict={}
ctrl_traits=['facial_hair', 'age_smkinit', 'antisoc', 'friend_sat', 'hr', 'infant_bw', 'LDL', 'maternal_smok', 'townsend', 'age_menarche', 'neurot','addict-rf']
for t in ctrl_traits:
    ctrl_dict[t]=pd.read_csv('gwas_ctrl_hm/magma/seed_genes/'+t+'_annot.tsv',sep='\t')
for t in ctrl_traits:
    seed_dict[t+'_FDR']=(set(ctrl_dict[t][ctrl_dict[t]['Q']<0.05]['GENE']))
    seed_dict[t+'_bonf']=(set(ctrl_dict[t][ctrl_dict[t]['P']<0.05/len(ctrl_dict[t])]['GENE']))
    if not tissue_network:
        seed_dict[t+'_top500']=set(ctrl_dict[t][(ctrl_dict[t]['GENE'].isin(all_nodes))].nsmallest(500,'P')['GENE'])
    else:
        seed_dict[t+'_top500']=set(ctrl_dict[t][(ctrl_dict[t]['GENE'].isin(hgnc[hgnc.entrez_id.isin(all_nodes)]['symbol']))].nsmallest(500,'P')['GENE'])



if tissue_network:
    seed_dict={k: set(hgnc[hgnc.symbol.isin(v)]['entrez_id']) for k, v in seed_dict.items()}


NPS_dict,NPS_dict_series=import_NPS_scores(seed_dict,interactome_name)

overwrite=False

#loop for only subset of traits genes- define ls. otherwise use ls=seed_dict.keys()
#ls=[i for i in seed_dict.keys() if 'final' in i]
ls=seed_dict.keys()
for k in ls:  
    seed_genes = list(seed_dict[k].intersection(all_nodes))
    print(f'analyzing {k}')
    if (len(seed_genes)>0):
        file_path='network_scores/'+k+'_'+interactome_name+'_zscore.tsv'
        if ((os.path.exists(file_path))&(overwrite==False)):
            print('File already exists. If you would like to overwrite this file, set overwrite=True,and rerun')
        else:
            NPSc, Fnew_score, Fnew_rand_score = netprop_zscore.calculate_heat_zscores(
                w_double_prime,  
                list(all_nodes),
                dict(degree), 
                seed_genes, num_reps=1000,
                minimum_bin_size=100,
                random_seed=random_seed)
            print(NPSc.head())
            if save_file:
                file_path='network_scores/'+k+'_'+interactome_name+'_zscore.tsv'
                print(f'saving to path: {file_path}')
                NPSc.to_csv(file_path,sep='\t',header=False)
    else:
        print('not enough seed genes for propagation (n=0)')
