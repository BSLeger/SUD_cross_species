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


interactome_name='PCNet2.0'
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

seed_dict=import_seed_dict(mag_dir,file_dict,ctrl_traits,ctrl_traits_rat,psych_traits,bonf_dict,gene_col_dict,all_nodes)
NPS_dict,NPS_dict_series=import_NPS_scores(seed_dict,interactome_name)

'''traits=['oud_deak_2022','cud_dbGAP','cigday','pau','biomarkers-30610-both_sexes-irnt',
 'biomarkers-30670-both_sexes-irnt',
 'biomarkers-30700-both_sexes-irnt',
 'biomarkers-30710-both_sexes-irnt',
 'biomarkers-30720-both_sexes-irnt',
 'biomarkers-30750-both_sexes-irnt',
 'biomarkers-30760-both_sexes-irnt',
 'biomarkers-30770-both_sexes-irnt',
 'biomarkers-30810-both_sexes-irnt',
 'categorical-1747-both_sexes-5',
 'categorical-20002-both_sexes-1065',
 'categorical-20002-both_sexes-1111',
 'categorical-20002-both_sexes-1226',
 'categorical-20002-both_sexes-1474',
 'categorical-20003-both_sexes-1141194794',
 'categorical-20116-both_sexes-0',
 'categorical-20116-both_sexes-2',
 'categorical-20160-both_sexes-20160',
 'categorical-4293-both_sexes-3',
 'categorical-4294-both_sexes-0',
 'continuous-102-both_sexes-irnt',
 'continuous-20016-both_sexes-irnt',
 'continuous-20153-both_sexes-irnt',
 'continuous-21002-both_sexes-irnt',
 'continuous-2178-both_sexes',
 'continuous-23100-both_sexes-irnt',
 'continuous-23101-both_sexes-irnt',
 'continuous-23106-both_sexes-irnt',
 'continuous-23116-both_sexes-irnt',
 'continuous-30000-both_sexes-irnt',
 'continuous-30010-both_sexes-irnt',
 'continuous-30040-both_sexes-irnt',
 'continuous-30070-both_sexes-irnt',
 'continuous-30080-both_sexes-irnt',
 'continuous-30100-both_sexes-irnt',
 'continuous-30120-both_sexes-irnt',
 'continuous-30130-both_sexes-irnt',
 'continuous-30150-both_sexes-irnt',
 'continuous-30180-both_sexes-irnt',
 'continuous-30200-both_sexes-irnt',
 'continuous-30300-both_sexes-irnt',
 'continuous-3143-both_sexes-irnt',
 'continuous-3148-both_sexes-irnt',
 'continuous-4104-both_sexes-irnt',
 'continuous-5134-both_sexes-irnt',
 'continuous-5257-both_sexes-irnt',
 'continuous-AG-both_sexes-irnt',
 'continuous-LDLC-both_sexes-medadj_irnt',
 'continuous-LDLC-both_sexes-medadj_raw',
 'continuous-MAP-both_sexes-manual_medadj_raw',
 'continuous-NAP-both_sexes-irnt',
 'continuous-PP-both_sexes-combined_medadj_irnt',
 'continuous-PP-both_sexes-combined_medadj_raw',
 'continuous-SBP-both_sexes-combined_medadj_irnt',
 'continuous-SBP-both_sexes-combined_medadj_raw',
 'continuous-eGFRcreacys-both_sexes-irnt',
 'icd10-E66-both_sexes',
 'icd10-M17-both_sexes',
 'phecode-250.2-both_sexes',
 'phecode-327.3-both_sexes',
 'phecode-411.4-both_sexes',
 'phecode-496.21-both_sexes',
 'phecode-594.1-both_sexes',
 'phecode-714-both_sexes','continuous-50-both_sexes-irnt']'''

traits=['categorical-20002-both_sexes-1226',
'continuous-30010-both_sexes-irnt',
'continuous-30040-both_sexes-irnt',
'continuous-30070-both_sexes-irnt',
'continuous-30080-both_sexes-irnt',
'continuous-4104-both_sexes-irnt']

traitn=int(sys.argv[1])-1
print(traitn)
trait=traits[traitn]

# modified for rerun
#ls=[w for w in list(seed_dict.keys()) if trait in w]
ls=[f'{trait}_top500',f'{trait}_bonf',f'{trait}_FDR']
print(ls)

overwrite=False

#ls=seed_dict.keys()
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
