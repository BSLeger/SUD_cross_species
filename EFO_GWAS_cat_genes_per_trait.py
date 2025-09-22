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

os.chdir('/tscc/projects/ps-palmer/brittany/SUD_cross_species/scripts')
from network_functions import *
from network_validation_functions import *
from plotting_functions import *

os.chdir('/tscc/projects/ps-palmer/brittany/SUD_cross_species/')


remove_ext=True # remove GWAS results for externalizing and component GWAS from GWAS catalog validation set

interactome_name='PCNet2.0'

if (interactome_name=='PCNet2.0'):
    all_nodes=list(pd.read_csv('PCNET2.0_allNodes.tsv',header=None)[0])
else:
    interactome=import_interactome(UUIDs=UUIDs,interactome_name=interactome_name)
    all_nodes=list(interactome.nodes())


#graphical representation of the EFO ontology database
graph = obo.read_obo('https://www.ebi.ac.uk/efo/efo.obo')


#get id from name or name from ID
id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}

nodes=list(graph.nodes())
print('import of EFO network complete')

#ontology mapping table for GWAS catalog
onto=pd.read_csv('validation_datasets/trait_mappings.txt',sep='\t')

onto['EFO_term']=onto['EFO URI'].apply(lambda x: x.split('/')[len(x.split('/'))-1])
onto['EFO_term']=onto['EFO_term'].str.replace('_',':')

onto['EFO_parent']=onto['Parent URI'].apply(lambda x: x.split('/')[len(x.split('/'))-1])
onto['EFO_parent']=onto['EFO_parent'].str.replace('_',':')


onto['Disease trait']=onto['Disease trait'].str.lower()
o=set(onto['Disease trait'])
e=(set(onto['EFO_term']))

onto['EFO_parent_idtoname']=onto['EFO_parent'].map(id_to_name)
print('import of trait mapping complete')

gwas_catalog=pd.read_csv('validation_datasets/gwas-catalog_v1.0.2_20240807.tsv',sep='\t',low_memory=False)
ext_studies=['GCST90061435','GCST007543','GCST009886','GCST006716','GCST006421','GCST007327','GCST007474']
if remove_ext:
    gwas_catalog=gwas_catalog[~gwas_catalog['STUDY ACCESSION'].isin( ext_studies)]
cat=format_catalog(gwas_catalog)

cat=cat.merge(onto, left_on='DISEASE/TRAIT',right_on='Disease trait')
print('import of GWAS catalog complete')

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

par=['disease','measurement','biological_process','experimental factor']
par_children={}
for x in par:
    par_children[x]=nx.ancestors(graph,name_to_id[x])

par_children['experimental factor']=par_children['experimental factor']-par_children['biological_process']-par_children['measurement']-par_children['disease']


par_children_terms={}
for j in par_children.keys():
    print(j)
    par_children_terms[j]=[id_to_name[k] for k in list(par_children[j])]

for p in par_children_terms.keys():
    my_list=list(return_descendents_name(graph,id_to_name, name_to_id[p]))
    with open(f'validation_datasets/efo_parent_children_terms_names_{p}.txt', "w") as file:
        for item in my_list:
            file.write(item + "\n")
par_children_cat={}
for x in par:
    par_children_cat[x]=par_children[x].intersection(set(onto.EFO_term))

if not remove_ext:
    gpt_path='validation_datasets/GWAS-CAT-EFO_genes_per_trait.csv'
else:
    gpt_path='validation_datasets/GWAS-CAT-EFO_genes_per_trait_ext_removed.csv'

print('restructuring of parent terms complete')
if (os.path.isfile(gpt_path)):
    genes_per_trait=pd.read_csv(gpt_path)
    print('reading in file')
else:
    print('calculating genes per trait')
    genes_per_trait=count_genes_per_trait([name_to_id[x] for x in par],cat,'GENE','EFO_term',id_to_name,name_to_id,graph)
    genes_per_trait.to_csv(gpt_path,index=False)