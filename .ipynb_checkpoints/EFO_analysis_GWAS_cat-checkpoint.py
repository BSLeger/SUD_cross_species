#!/usr/bin/env python
# coding: utf-8

# purpose: query the EFO for parent and children terms to use for defining parallel traits in gwas catalog
# based on code from obonet's tutorial:
#     https://github.com/dhimmel/obonet/blob/main/examples/go-obonet.ipynb

# # set-up

# In[1]:


import os
import pandas as pd
import obonet as obo
import networkx as nx
import ndex2
import numpy as np
import matplotlib.pyplot as plt
from upsetplot import plot as upplot
from upsetplot import from_contents
from upsetplot import UpSet
import scipy.stats as stats
from adjustText import adjust_text
import matplotlib.patches as mpatches
import sys


# In[2]:

# Add the directory containing rca_functions to the system path
sys.path.append('/tscc/projects/ps-palmer/brittany/rare_common_alcohol/rare_common_alcohol_comparison/notebooks/')
from rca_functions import *
sys.path.remove('/tscc/projects/ps-palmer/brittany/rare_common_alcohol/rare_common_alcohol_comparison/notebooks/')

# Add the scripts directory to the system path
sys.path.append('/tscc/projects/ps-palmer/brittany/SUD_cross_species/scripts')
from network_functions import *
from network_validation_functions import *
from plotting_functions import *


# In[3]:


os.chdir('/tscc/projects/ps-palmer/brittany/SUD_cross_species/')


# In[4]:


save_fig=True


# In[5]:


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


# # functions

# # import network genelists

# In[6]:


remove_ext=True # remove GWAS results for externalizing and component GWAS from GWAS catalog validation set


# ## Interactome Set-up

# pcnet2- versions 
# from wright et al. 2024 preprint:
# PCNet 2.0= best-performing ranked composite (top 15 interactomes, 3.85M interactions)
# PCNet 2.1= top 8 interactomes, 1.75M interactions
# PCNet 2.2= top 10 co-citation-free interactomes, 3.32M interactions 

# In[7]:


interactome_name='PCNet2.0'


# In[8]:


if (interactome_name=='PCNet2.0'):
    all_nodes=list(pd.read_csv('PCNET2.0_allNodes.tsv',header=None)[0])
else:
    interactome=import_interactome(UUIDs=UUIDs,interactome_name=interactome_name)
    all_nodes=list(interactome.nodes())


# ## import NPS scores and seed genes

# In[9]:


seed_dict=import_seed_dict(mag_dir,file_dict,ctrl_traits,ctrl_traits_rat,bonf_dict,gene_col_dict,all_nodes)
seed_dict.keys()


# In[10]:


NPS_dict,NPS_dict_series=import_NPS_scores(seed_dict,interactome_name)


# # import the network of EFO terms

# In[11]:


#graphical representation of the EFO ontology database
graph = obo.read_obo('https://www.ebi.ac.uk/efo/efo.obo')


# In[12]:


#get id from name or name from ID
id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}


# In[13]:


nodes=list(graph.nodes())


# In[14]:


nx.is_directed_acyclic_graph(graph) #undirected network??? HOW???


# # import EFO table from GWAS catalog

# In[15]:


#ontology mapping table for GWAS catalog
onto=pd.read_csv('validation_datasets/trait_mappings.txt',sep='\t')


# In[16]:


onto['EFO_term']=onto['EFO URI'].apply(lambda x: x.split('/')[len(x.split('/'))-1])
onto['EFO_term']=onto['EFO_term'].str.replace('_',':')


# In[17]:


onto['EFO_parent']=onto['Parent URI'].apply(lambda x: x.split('/')[len(x.split('/'))-1])
onto['EFO_parent']=onto['EFO_parent'].str.replace('_',':')


# In[18]:


onto['Disease trait']=onto['Disease trait'].str.lower()
o=set(onto['Disease trait'])
e=(set(onto['EFO_term']))


# In[19]:


onto['EFO_parent_idtoname']=onto['EFO_parent'].map(id_to_name)


# In[20]:


onto[['Parent term','EFO_parent_idtoname','EFO_parent']].drop_duplicates()


# ## import gwas catalog associations

# In[21]:


gwas_catalog=pd.read_csv('validation_datasets/gwas-catalog_v1.0.2_20240807.tsv',sep='\t',low_memory=False)
ext_studies=['GCST90061435','GCST007543','GCST009886','GCST006716','GCST006421','GCST007327','GCST007474']
if remove_ext:
    gwas_catalog=gwas_catalog[~gwas_catalog['STUDY ACCESSION'].isin( ext_studies)]
cat=format_catalog(gwas_catalog)


# In[22]:


cat=cat.merge(onto, left_on='DISEASE/TRAIT',right_on='Disease trait')


# # redefine parent terms in GWAS catalog ontology table

# parent terms in the GWAS catalog are at random depths within the catalog, which is confusing for the enrichment analysis. To fix this, I'm going to redefine them in the following terms:
# 
# disease
# measurement
# biological function
# experimental factor (if not present in the above 3)
# NR- not in EFO
# 
# this way, most of the parents will be in one of the above terms, THEN, we can look at children in a more informative way

# In[23]:


sorted(set(onto['EFO_parent_idtoname'].dropna()))


# ## identify shared parent terms

# In[24]:


parents=set(cat['EFO_parent'])
parents.discard('NR')
for k in sorted(parents):
    #print(id_to_name[k])
    if k in nodes:
        #return_efo_ancestors_name(graph,id_to_name,k)
        #print(': \n')
        print((list((id_to_name[supterm] for supterm in nx.descendants(graph, k))))) #descendents get superterms- not sure why)
        print('\n')
    else:
        print(' not in EFO. \n')


# ## define new parent catagories

# In[25]:


par=['disease','measurement','biological_process','experimental factor']


# In[26]:


par_children={}
for x in par:
    par_children[x]=nx.ancestors(graph,name_to_id[x])


# In[27]:


par_children['experimental factor']=par_children['experimental factor']-par_children['biological_process']-par_children['measurement']-par_children['disease']


# In[28]:


par_children_terms={}
for j in par_children.keys():
    print(j)
    par_children_terms[j]=[id_to_name[k] for k in list(par_children[j])]


# In[30]:


#children=return_descendents_name(graph,id_to_name, name_to_id['measurement'])


# In[31]:


for p in par_children_terms.keys():
    my_list=list(return_descendents_name(graph,id_to_name, name_to_id[p]))
    with open(f'validation_datasets/efo_parent_children_terms_names_{p}.txt', "w") as file:
        for item in my_list:
            file.write(item + "\n")


# ## make upset plot showing overlap of new parent traits

# In[32]:


par_children_cat={}
for x in par:
    par_children_cat[x]=par_children[x].intersection(set(onto.EFO_term))


# In[33]:


par_children_cat


# In[34]:


#check overlap for only traits in the GWAS catalog
table=from_contents(par_children_cat)  
UpSet(table, subset_size='count',show_counts=True).plot()
if save_fig:
    plt.savefig('figures/' + 'gwas_catalog_parent_upset.pdf')


# In[35]:


#check for all traits
table=from_contents(par_children)  
UpSet(table, subset_size='count',show_counts=True).plot()
if save_fig:
    plt.savefig('figures/' + 'gwas_catalog_parent_upset.pdf')


# ## calculate number of genes per trait

# In[37]:


if not remove_ext:
    gpt_path='validation_datasets/GWAS-CAT-EFO_genes_per_trait.csv'
else:
    gpt_path='validation_datasets/GWAS-CAT-EFO_genes_per_trait_ext_removed.csv'


if (os.path.isfile(gpt_path)):
    genes_per_trait=pd.read_csv(gpt_path)
    print('reading in file')
else:
    print('calculating genes per trait')
    genes_per_trait=count_genes_per_trait([name_to_id[x] for x in par],cat,'GENE','EFO_term',id_to_name,name_to_id,graph)
    genes_per_trait.to_csv(gpt_path,index=False)


# # choose datasets for analysis- put gene lists into dictionary

# In[ ]:


#modify for correct genesets
trait_h='ext'
cut_h='top500'

trait_r='loco_final_cf'
cut_r= 'FDR'



#choose which community to check enrichment for
#must be keys from coloc_dict
sub_community='seed_r'

label_h,label_r,seed_h,seed_r,NPS_h,NPS_r,NPS=return_analysis_datasets(trait_r,cut_r,trait_h,cut_h,seed_dict,NPS_dict,interactome_name)


# In[39]:


filter_traits=False
gpt=genes_per_trait[(genes_per_trait.ngenes>=5)]


# In[40]:


#get rest of datasets based on what was specified above
coloc_dict=def_coloc_dict(seed_r,seed_h,NPS,all_nodes,cut_single,cut_comb,cut_rat_specific,cut_hm_specific)
val_lab_dict=def_val_label_dict(label_h,label_r,interactome_name,cut_single,cut_comb)


# In[44]:


#sub_community is the community being tested for enrichment (i.e. the ext-loco network, 'net')
#whole_community is the greater pool of genes to test against (i.e. PCNET ('graph'))

if remove_ext:
    rmv='_remove_ext'
else:
    rmv=''
if sub_community in val_lab_dict.keys():
    whole_community=coloc_dict_ref[sub_community]
    outpath=f'validation_output/GWAS-CAT-EFO_enr_{val_lab_dict[sub_community]}_enr{rmv}.csv'
    #outpath=f'validation_output/GWAS-CAT-EFO_enr_{val_lab_dict[sub_community]}_enr_original.csv'
else:
    print('sub_community not in val_lab_dict- using graph as the whole_community')
    whole_community='graph'
    outpath=f'validation_output/GWAS-CAT-EFO_enr_temp{rmv}.csv'
print(f'path for this output file: {outpath}')
if (outpath==f'validation_output/GWAS-CAT-EFO_enr_temp{rmv}.csv'):
    print('sub_community not in dictionary- will be saved as a temporary file (may overwrite previous temporary file).')
elif (os.path.isfile(outpath)):
    run_analysis=False
    print('this analysis has been run previously- importing from file. If want to rerun, set run_analysis to True')
    tbl=pd.read_csv(outpath)
else:
    run_analysis=True
    print('this analysis has not been run- run_analysis set to True')


# # validate whole geneset at all depths of the GWAS catalog

# In[45]:
run_analysis=True

if run_analysis:
    coloc_dict_cat={}
    for k in coloc_dict.keys():
        coloc_dict_cat[k]=set(coloc_dict[k]).intersection(list(cat['GENE'].dropna()))


# In[46]:


print(f'WARNING: YOU ARE CURRENTLY ANALYZING THE ENRICHMENT OF {sub_community} RELATIVE TO {whole_community}.\n\tTHIS IS CALCULATED FROM {label_r}, {label_h}, AND {interactome_name}.\n\tIF THIS IS CORRECT, CONTINUE :)')


# In[ ]:

# In[ ]:


if run_analysis:
    tbl=recurse_enrichment(par,graph,id_to_name, 
        name_to_id,'exploratory factor',cat,'EFO_term','GENE',
        coloc_dict_cat,sub_community,
        whole_community,outpath,depth=0,depth_term=5,verbose=False,enr_concat=None)
    #tbl.to_csv(outpath,index=False)
    #print(f'table saved as {outpath}')


# In[ ]:


if ('net' in sub_community):
    if not 'n_seed' in tbl.columns:
        if 'community_genes' in tbl.columns:
            col='community_genes'
        else:
            col='network_genes'
        if type(tbl[col].iloc[0])!=set:
            tbl['seed_h']=tbl[col].apply(lambda x: set(x.replace("'","").replace("{","").replace("}","").split(',')).intersection(coloc_dict['seed_h']))
            tbl['seed_r']=tbl[col].apply(lambda x: set(x.replace("'","").replace("{","").replace("}","").split(',')).intersection(coloc_dict['seed_h']))
        else:
            tbl['seed_h']=tbl[col].apply(lambda x: x.intersection(coloc_dict['seed_h']))
            tbl['seed_r']=tbl[col].apply(lambda x: x.intersection(coloc_dict['seed_r']))
        tbl['n_seed_h']=tbl.seed_h.apply(lambda x: len(x))
        tbl['n_seed_r']=tbl.seed_r.apply(lambda x: len(x))
        tbl['n_seed']=tbl['n_seed_h']+tbl['n_seed_r']
        if not filter_traits:
            tbl.to_csv(outpath,index=False)
        tbl.sort_values('n_seed',ascending=False).head()
if not('ngenes' in tbl.columns):
    tbl=tbl.merge(gpt,left_on='trait',right_on='trait_name').drop('trait_name',axis=1)
tbl=tbl.drop_duplicates()
tbl.to_csv(outpath,index=False)


# In[ ]:


len_t=len(set(tbl.trait))


# In[46]:


filter_traits=True


# In[47]:


print(f'filter traits={filter_traits}')


# In[48]:


if filter_traits:
    print('filtering traits')
    tbl=tbl[tbl.trait.isin(gpt.trait_name)]


# In[66]:


tbl=tbl.drop_duplicates()


# In[71]:


tbl[tbl.p_intersect<(0.05/len_t)].sort_values('ngenes')


# # plot enrichment

# In[72]:


save_fig=True


# In[73]:


if save_fig==True:
    outpath_dir=outpath[:len(outpath)-4]
    if not os.path.exists(outpath_dir):
        os.makedirs(outpath_dir)


# In[74]:


colormap=plt.colormaps.get_cmap('tab20b')


# In[75]:


for p in par:
    children=return_descendents_name(graph,id_to_name, name_to_id[p])
    #make the terms somewhat mutually exclusive- otherwise this is unintelligible
    '''    if p=='disease':
            bp_child=return_descendents_name(graph,id_to_name, name_to_id['biological_process'])
            meas_child=return_descendents_name(graph,id_to_name, name_to_id['measurement'])
            children=set(children).difference(meas_child+bp_child)
        if p=='biological_process':
            meas_child=return_descendents_name(graph,id_to_name, name_to_id['measurement'])
            disease_child=return_descendents_name(graph,id_to_name, name_to_id['disease'])
            children=set(children).difference(meas_child)
        if p=='experimental factor':
            meas_child=return_descendents_name(graph,id_to_name, name_to_id['measurement'])
            disease_child=return_descendents_name(graph,id_to_name, name_to_id['disease'])
            bp_child=return_descendents_name(graph,id_to_name, name_to_id['biological_process'])
            children=set(children).difference(meas_child+disease_child+bp_child)'''
    tb=tbl[tbl.trait.isin(children)]    
    depth_set=set(tb['depth'])
    for d in set(tb['depth']):
        t=tb[tb.depth==d]
        t=t[t['p_intersect']<(0.05/len_t)]
        if len(t)>0:
            t=t.sort_values('p_intersect',ascending=False)
            # Assign a unique color for each parent_trait based on its index
            unique_traits = t.parent_trait.unique()
            color_mapping = {trait: colormap(i / len(unique_traits)) for i, trait in enumerate(unique_traits)}
            
            # Plotting
            fig = plt.figure(figsize=(10, 15))
            colors = [color_mapping[trait] for trait in t.parent_trait]  # Color assignment based on parent_trait
            
            plt.barh(y=t.trait, width=-np.log10(t.p_intersect*len_t),label=t.parent_trait, color=colors)
            
            plt.ylabel('Experimental Phenotype')
            plt.xlabel('corrected -log10(p-value)')
            plt.axvline(-np.log10(0.05), color='black', ls=':')
            plt.axvline(0, color='black', ls='-')
        
            #plt.xticks(rotation=90)
            #plt.title(d)
            plt.title(f'GWAS Cat ({p}) val {val_lab_dict[sub_community]}, depth={d}')
            # Create legend patches (only one for each unique parent_trait)
            legend_patches = [mpatches.Patch(color=color_mapping[trait], label=trait) for trait in unique_traits]
            
            # Add legend to the plot
            plt.legend(handles=legend_patches, title='Parent Trait', loc='upper right',bbox_to_anchor=(2, 1))
            #plt.tight_layout()
            plt.margins(y=0)
            # Show the plot
            if save_fig:
                if not filter_traits:
                    plt.savefig((f'{outpath_dir}/GWAS_Cat_{p}_val_bar_depth-{d}.svg'), bbox_inches = "tight")
                else:
                    plt.savefig((f'{outpath_dir}/GWAS_Cat_{p}_val_filtered_bar_depth-{d}.svg'), bbox_inches = "tight")

            plt.show()
            
            # Clear the figure after showing
            plt.clf()


# In[56]:


for p in par:
    tb=tbl[(tbl.parent_trait==p)&(tbl.p_intersect<(0.05/len_t))].sort_values('p_intersect')
    display(tb.head(10))

