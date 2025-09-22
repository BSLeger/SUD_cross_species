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
from network_validation_functions import *
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


print(f'value passed={(sys.argv[1])}')
min_loop=100*int(sys.argv[1])
max_loop=100*(int(sys.argv[1])+1)
print(f'{min_loop}: {max_loop}')
#choose which interactome you want to import
tissue_network=False

if tissue_network:
    tissue='basal_ganglion'
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

if not tissue_network:
    seed_dict=import_seed_dict(mag_dir,file_dict,ctrl_traits,ctrl_traits_rat,psych_traits,bonf_dict,gene_col_dict,all_nodes)
else:
    hgnc=pd.read_csv('hgnc_complete_set.txt',sep='\t',low_memory=False)
    hgnc=hgnc[['symbol','entrez_id']].dropna()
    hgnc['entrez_id']=hgnc['entrez_id'].astype(int).astype(str)
    seed_dict=import_seed_dict(mag_dir,file_dict,ctrl_traits,ctrl_traits_rat,bonf_dict,gene_col_dict,hgnc[hgnc.entrez_id.isin(all_nodes)]['symbol'])

rat_TWAS_tissue_label={'BLA':'Basolateral amygdala',
'Brain':'Brain hemisphere',
'IL':'Infralimbic cortex',
'LHb':'Lateral habenula',
'NAcc':'Nucleus accumbens core',
'NAcc1':'Nucleus accumbens core 1',
'NAcc2':'Nucleus accumbens core 2',
'OFC':'Orbitofrontal cortex',
'PL':'Prelimbic cortex',
'PL1':'Prelimbic cortex 1',
'PL2':'Prelimbic cortex 2',
'Adipose':'Adipose',
'Eye':'Eye',
'Liver':'Liver'}

overwrite=True
prefix='loco_final_CF'
if not tissue_network: 
    for t in rat_TWAS_tissue_label.keys():
        path=f'rat_fusion/output/FUSION_concat/{prefix}_{t}_seed_genes.dat'
        if os.path.exists(path):
            tbl=pd.read_csv(path,low_memory=False, sep='\t')
            seed_dict[f'{prefix.lower()}_{t}_bonf']=set(tbl[tbl['TWAS.P']<0.05/len(tbl)]['HM_ORTHO'].dropna())
            seed_dict[f'{prefix.lower()}_{t}_FDR']=set(tbl[tbl['Q']<0.05]['HM_ORTHO'].dropna())
            tbl=tbl[~tbl.HM_ORTHO.isna()]
            seed_dict[f'{prefix.lower()}_{t}_top500']=set(tbl[tbl['HM_ORTHO'].isin(all_nodes)].nsmallest(500,'TWAS.P')['HM_ORTHO'])
            print(f'{t} (NORTHO={len(set(tbl.HM_ORTHO))}) BONF: {len(seed_dict[f"{prefix.lower()}_{t}_bonf"])} FDR: {len(seed_dict[f"{prefix.lower()}_{t}_FDR"])} top500: {len(seed_dict[f"{prefix.lower()}_{t}_top500"])}')
        else:
           print(f'{path} does not exist')  

if tissue_network:
    seed_dict={k: set(hgnc[hgnc.symbol.isin(v)]['entrez_id']) for k, v in seed_dict.items()}

k='loco_final_cf_FDR'
n=len(seed_dict['loco_final_cf_FDR'].intersection(all_nodes))


print('running propagations')
for i in range(min_loop,max_loop):
    seed_genes = random.sample(all_nodes, n)
    NPSc, Fnew_score, Fnew_rand_score = netprop_zscore.calculate_heat_zscores(
        w_double_prime,  
        list(all_nodes),
        dict(degree), 
        seed_genes, num_reps=1000,
        minimum_bin_size=100,
        random_seed=random_seed)
    print(NPSc.head())
    if save_file:
        file_path=f'network_scores/permuted_control/{k}_permuted_control_{i}_{interactome_name}_zscore.tsv'
        pd.DataFrame(seed_genes).to_csv(f'network_scores/permuted_control/{k}_permuted_control_{i}_{interactome_name}_seedgenes.tsv',index=False)
        if ((overwrite==False)&(os.path.exists(file_path))):
            print('File already exists. If you would like to overwrite this file, set overwrite=True, and rerun')
        else:
            NPSc.to_csv(file_path,sep='\t',header=False)
