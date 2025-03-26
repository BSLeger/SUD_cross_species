import os
import pandas as pd
import requests
os.chdir('/tscc/projects/ps-palmer/brittany/ddot')
import ddot
import obonet as obo
import networkx as nx
from scipy.stats import fisher_exact, hypergeom
import scipy.stats as stats
import numpy as np
import itertools

coloc_dict_ref={
    'seed_r':'magma_rat_ref',
    'seed_h':'magma_hm_ref',
    'net':'graph',
    'seed_hr':'magma_hm_rat_overlap_ref',
    'hm_net':'graph',
    'rat_net':'graph'
}

#dataset import functions-----------
def def_coloc_dict(seed_r,seed_h,NPS,all_nodes,cut_single,cut_comb,cut_rat_specific,cut_hm_specific):
    #human magma reference dataframe
    ref=pd.read_csv('/tscc/projects/ps-palmer/brittany/magma_v1/NCBI38/NCBI38.gene.loc',sep='\t',header=None)
    #rat magma reference- only those in magma reference file that have human orthologs- hm ortho must be used
    ortho=pd.read_csv('/tscc/projects/ps-palmer/brittany/orthology_ref_tbls/ORTHOLOGY-ALLIANCE_COMBINED_2024.tsv',sep='\t',skiprows=15)
    #downloaded from https://www.alliancegenome.org/downloads#orthology on 11 June 2024
    #filter for rat-human
    ortho=ortho[(ortho['Gene1SpeciesName']=='Rattus norvegicus')&(ortho['Gene2SpeciesName']=='Homo sapiens')]
    #filter for best match
    ortho=ortho[ortho['IsBestScore']=='Yes']
    gene_loc_file=pd.read_csv("magma/rn7.2_annotatedgenes_ncbi/rn7.2_gene_attribute_table_protein_coding_forMAGMA.tsv",sep='\t',header=None)
    
    coloc_dict={
        'seed_r':seed_r,
        'seed_h':seed_h,
        'seed_hr':list(set(seed_r).intersection(seed_h)),
        'net':list(NPS[(NPS.zh>cut_single)&(NPS.zr>cut_single)&(NPS.zhr>cut_comb)].index),
        'graph':all_nodes,
        'magma_hm_ref':set(ref[5]),
        'magma_rat_ref':set(ortho[ortho['Gene1Symbol'].isin(gene_loc_file[0])]['Gene2Symbol']),
        'hm_net':set(NPS[(NPS['zh'] > cut_hm_specific['zh']) & (NPS['zr'] < cut_hm_specific['zr']) &(NPS['zhr']<cut_hm_specific['zhr'])].index),
        'rat_net':set(NPS[(NPS['zr'] > cut_rat_specific['zr']) & (NPS['zh'] < cut_rat_specific['zh']) &(NPS['zhr']<cut_rat_specific['zhr'])].index)
    }
    coloc_dict['magma_hm_rat_overlap_ref']=coloc_dict['magma_hm_ref'].intersection(coloc_dict['magma_rat_ref'])
    return coloc_dict


def def_val_label_dict(label_h,label_r,interactome_name,cut_single,cut_comb):
	val_lab_dict={
	    'seed_h':f'seed_genes-{label_h}_{interactome_name}',
	    'seed_r':f'seed_genes-{label_r}_{interactome_name}',
	    'seed_hr':f'seed_genes-{label_r}-intersection-{label_h}_{interactome_name}',
	    'net':f'network-{label_h}-{label_r}_NPS-{cut_single}-{cut_comb}_{interactome_name}',
	    'hm_net':f'network-{label_h}-{label_r}_human_specific_{interactome_name}',
	    'rat_net':f'network-{label_h}-{label_r}_rat_specific_{interactome_name}'
	}
	return(val_lab_dict)



def count_genes_per_trait(head,graph_df,gene_col,term_col,id_to_name,name_to_id,graph):
    trait_list=list(itertools.chain.from_iterable([return_descendents_name(graph,id_to_name,x) for x in head]))
    genes_per_trait={}
    for trait in trait_list:
        children=return_descendents_name(graph,id_to_name, name_to_id[trait])
        t=set(graph_df[graph_df[term_col].isin([name_to_id[x] for x in children])][gene_col].dropna())
        genes_per_trait[trait]=len(t)
    genes_per_trait=pd.DataFrame(genes_per_trait.items(), columns=['trait_name', 'ngenes'])
    genes_per_trait['trait_ID']=[name_to_id[x] for x in genes_per_trait['trait_name']]
    return genes_per_trait
# ontology import functions--------------------------------------
def import_MPO_description(url='http://www.informatics.jax.org/downloads/reports/MPheno_OBO.ontology'):
    """
    Function to parse and load mouse phenotype ontology, using DDOT's ontology module
    Modified from :py:func:`netcoloc.validation.load_MPO`

    :param url: URL containing MPO ontology file
    :type url: str
    :return: MPO parsed using DDOT
    :rtype: :py:class:`pd.DataFrame`
    """

    # download the mammalian phenotype ontology, parse with ddot
    r = requests.get(url,allow_redirects=True)
    open('MPheno_OBO.ontology','wb').write(r.content)

    ddot.parse_obo('MPheno_OBO.ontology',
                   'parsed_mp.txt',
                  'id2name_mp.txt',
                  'id2namespace_mp.txt',
                  'altID_mp.txt')


    MP2desc = pd.read_csv('id2name_mp.txt',sep='\t',
                          names=['MP','description'],index_col='MP')

    MP2desc=MP2desc.loc[MP2desc.index.dropna()] # drop NAN from index
    print(len(MP2desc))
    return(MP2desc)
	
def return_descendents_name(graph,id_to_name, term):
	#return list set of traits in  ontology including parent- reported by name
	#modified from DDOT example
    l=list(sorted(id_to_name[subterm] for subterm in nx.ancestors(graph, term))) #descendents get subterms- not sure why but I tested it and OBONET says so as well
    l.append(id_to_name[term])
    return list(set(l))

def return_ancestors_name(graph,id_to_name, term):
	#return list set of traits in  ontology including parent- reported by name
	#modified from DDOT example
    print('getting ancestors for '+term+' : '+id_to_name[term])
    l=list((id_to_name[supterm] for supterm in nx.descendants(graph, term))) #descendents get superterms- not sure why
    #l.append(id_to_name[term])
    return list(set(l))



# validation functions------------------------------------------

def calculate_enrichment(t, coloc_dict_cat, k, sub='net', total='graph', verbose=True):
    # modified from rare_common_alcohol
    # Calculate enrichment for a group of genes (sub) versus all in larger group (total)
    # Calculate values for the contingency table
    M = len(coloc_dict_cat[total])  # Population size: genes in PCNet annotated in the GWAS catalog
    n = len(coloc_dict_cat[total].intersection(t))  # Genes in PCNet annotated for the trait of interest
    N = len(coloc_dict_cat[sub])  # Genes in the network annotated in the GWAS catalog
    x = len(coloc_dict_cat[sub].intersection(t))  # Genes in network annotated for the trait of interest

    # Build contingency table
    contingency_table = [
        [x, N - x],
        [n - x, M - N - (n - x)]
    ]
    
    # Perform Fisher's exact test
    odds_ratio, p_intersect = stats.fisher_exact(contingency_table, alternative='greater')
    gene_list = t.intersection(coloc_dict_cat[sub])
    
    try:
        se = np.sqrt(1/x + 1/(N-x) + 1/(n-x) + 1/(M - N - (n - x)))
    except ZeroDivisionError:
        se = None
    
    # Create the results dictionary
    if verbose:
        print(f"Enrichment of network nodes in genes in the GWAS catalog annotated for {k}: p={p_intersect}")
        p_value_hypergeom = stats.hypergeom.sf(x-1, M, n, N)
        print(f'Enrichment calculated using hypergeom.sf for {k}: p={p_value_hypergeom}')
        print(f"Odds ratio: OD={odds_ratio}")
        print(f"Number of annotated genes in {total}: {len(t.intersection(coloc_dict_cat[total]))}")
        print(f"Number of annotated genes in {sub}: {len(t.intersection(coloc_dict_cat[sub]))}\n")   
    return odds_ratio, se, p_intersect, gene_list


def recurse_enrichment(par,graph,id_to_name, name_to_id,heirarchy_name,graph_df,term_col,gene_col,coloc_dict_cat,sub_community,whole_community,outpath=None,depth=0,depth_term=None,verbose=False,enr_concat=None):
    '''
    par- list of starting parent traits
    graph- graph structure of the network
    id_to_name- dictionary that transforms names to ID#s
    name_to_id- dictionary that transforms IDs to names
    heirarchy_name- string- name of the heirarchical structure
    graph_df- dataframe of the graph that shows relationship between trait and gene
    term_col- string- column of dataframe that has the term names of the communities in the heriarchy. For MGI it's 'MP' for GWAS catalog it's 'EFO_term'
    gene_col- string- column of dataframe that has the human gene names. For MGI it's 'human_ortholog', for GWAS catalog it's 'GENE'
    coloc_dict_cat- dictionary of genesets (sourced from coloc_dict) that are annoated in graph
    sub_community- string- community you're testing enrichment for. should be a key in coloc_dict_cat
    whole_community- string- community you're testing enrichment against. should be a key in coloc_dict_cat
    outpath- string- path to where the file is saved
    depth- int- variable that tracks what depth of the heirarchical structure you are in
    depth_term- int- overwrite depth to terminate at. Use None unless you want to terminate early. You will probably want to terminate early for EFO recommend depth_term=5.
    verbose- string- determines whether to print the enrichments as they are calculated
    enr_concat- dataframe- initially pass None 

    returns: dataframe-enr_concat
    '''
    if depth==None:
        depth=0
    enr_tbl=pd.DataFrame(columns=['trait','parent_trait','community_genes','n_community_genes','odds_ratio','log_se_or','p_intersect','depth'])

    print(f"analyzing structure depth={depth}")
    '''
    this code was written like this for the following reason- If you just started with the top of the heirarchy's single trait
    this would be unnecessary- HOWEVER, I wanted to be able to input a starting list of parent terms, as opposed to a singular term, mostly because of the EFO.
    so this takes the list of parent terms- then after that it goes through the structure following the standard way. 
    if you did it the other way, you would just want to use the else statement part of this code.
    '''
    if (depth==0):        
        for p in par:
            #get all children terms
            children=return_descendents_name(graph,id_to_name, name_to_id[p])
            #get all genes annotated for children terms 
            t=set(graph_df[graph_df[term_col].isin([name_to_id[x] for x in children])][gene_col].dropna())
            odds_ratio, log_se_or, p_intersect, gene_list= calculate_enrichment(t,coloc_dict_cat,p,sub_community,whole_community,verbose)
            if (len(gene_list)>0):
                enr_tbl = pd.concat([pd.DataFrame([[p, heirarchy_name, gene_list, len(list(gene_list)),odds_ratio, log_se_or, p_intersect,depth]], columns=enr_tbl.columns), enr_tbl], ignore_index=True)
                enr_tbl['depth']=depth
                if len(enr_tbl)>0:
                    enr_tbl.to_csv(outpath,index=False)

    else:
       for p in par:
            #select parent term
            #get all children terms (not all descendents, just first level children)- loop over
            for c in list(map(lambda x: id_to_name[x],list(graph.predecessors(name_to_id[p])))):
                # get all grand children- use all terms as subterm
                children=return_descendents_name(graph,id_to_name, name_to_id[c])
                t=set(graph_df[graph_df[term_col].isin([name_to_id[x] for x in children])][gene_col].dropna())
                odds_ratio, log_se_or, p_intersect, gene_list= calculate_enrichment(t,coloc_dict_cat,c,sub_community,whole_community,verbose)
                if (len(gene_list)>0):
                    enr_tbl = pd.concat([pd.DataFrame([[c, p, gene_list,len(list(gene_list)), odds_ratio, log_se_or, p_intersect,depth]], columns=enr_tbl.columns), enr_tbl], ignore_index=True)
                    enr_tbl['depth']=depth
                if len(enr_tbl)>0:
                    enr_tbl.to_csv(outpath,index=False,header=False,mode='a')
             
    children=set(enr_tbl['trait'])
    if enr_concat is None:
        enr_concat = pd.DataFrame(columns=['trait', 'parent_trait', 'community_genes', 'n_community_genes','odds_ratio', 'log_se_or', 'p_intersect', 'depth'])
    if not (enr_tbl.empty):
        enr_concat = pd.concat([enr_concat, enr_tbl])
    
    print(f"\tlength of enrichment table={len(enr_tbl)}")
    print(f"\tlength of concatenated enrichment table={len(enr_concat)}")
    print(f"\tlength of children={len(children)}")
    if (depth==depth_term):
        print('reached overwritten depth-returning concatenated enrichment table')
        print(enr_concat.head())
        return enr_concat
    else:
        if (len(children)!=0):
            return recurse_enrichment(children,graph,id_to_name, name_to_id,heirarchy_name,graph_df,term_col,gene_col,coloc_dict_cat,sub_community,whole_community,outpath,depth+1,depth_term,verbose,enr_concat)
        else:
            print('returning concatenated enrichment table')
            print(enr_concat.head())
            return enr_concat
