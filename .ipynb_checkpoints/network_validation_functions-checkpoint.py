import os
import pandas as pd
import requests
import sys
sys.path.append('/tscc/projects/ps-palmer/brittany/ddot')
import ddot
sys.path.remove('/tscc/projects/ps-palmer/brittany/ddot')
import obonet as obo
import networkx as nx
from scipy.stats import fisher_exact, hypergeom
import scipy.stats as stats
import numpy as np
import itertools
from scipy.stats import norm
from statsmodels.stats import contingency_tables


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
    net=list(NPS[(NPS.zh>cut_single)&(NPS.zr>cut_single)&(NPS.zhr>cut_comb)].index)
    outNet=NPS[~NPS.index.isin(net)]
    coloc_dict={
        'seed_r':seed_r,
        'seed_h':seed_h,
        'seed_hr':list(set(seed_r).intersection(seed_h)),
        'net':net,
        'graph':all_nodes,
        'magma_hm_ref':set(ref[5]),
        'magma_rat_ref':set(ortho[ortho['Gene1Symbol'].isin(gene_loc_file[0])]['Gene2Symbol']),
		'hm_net':set(NPS[NPS['zh'] > cut_hm_specific['zh']].index).difference(net),
        'rat_net':set(NPS[NPS['zr'] > cut_rat_specific['zr']].index).difference(net),
        'hm_net_alt':set(outNet.sort_values('zh',ascending=False).head(len(net)).index),
        'rat_net_alt':set(outNet.sort_values('zr',ascending=False).head(len(net)).index)}

    '''	'hm_net':set(NPS[(NPS['zh'] > cut_hm_specific['zh']) & (NPS['zr'] < cut_hm_specific['zr']) &(NPS['zhr']  <cut_hm_specific['zhr'])].index),'rat_net':set(NPS[(NPS['zr'] > cut_rat_specific['zr']) & (NPS['zh'] < cut_rat_specific['zh']) &(NPS['zhr']<cut_rat_specific['zhr'])].index)
    '''
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
    CT = contingency_tables.Table2x2(contingency_table)

    OR_p_temp = CT.oddsratio_pvalue()
    OR_CI_temp = CT.oddsratio_confint()
    OR = CT.oddsratio

    try:
        se = np.sqrt(1/x + 1/(N-x) + 1/(n-x) + 1/(M - N - (n - x)))
    except ZeroDivisionError:
        se = None
    
    # Create the results dictionary
    if verbose:
        print(f"Enrichment of network nodes in genes in the GWAS catalog annotated for {k}: p={p_intersect}")
        p_value_hypergeom = stats.hypergeom.sf(x-1, M, n, N)
        print(f'Enrichment calculated using hypergeom.sf for {k}: p={p_value_hypergeom}, p_contingency_table={OR_p_temp}')
        print(f"Odds ratio: OD={odds_ratio}, via_contingency OR={OR}")
        print(f"Number of annotated genes in {total}: {len(t.intersection(coloc_dict_cat[total]))}")
        print(f"Number of annotated genes in {sub}: {len(t.intersection(coloc_dict_cat[sub]))}\n")
        print(f"Number of genes in interactome annotated for trait: {N}")
    return odds_ratio, se, p_intersect, gene_list


'''def calculate_enrichment(t, coloc_dict_cat, k, sub='net', total='graph', verbose=True):
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
    return odds_ratio, se, p_intersect, gene_list'''


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

## Extensions to NetColoc from CrossSpeciesBMI github: https://github.com/sarah-n-wright/CrossSpeciesBMI --------------------------------------------------------------------
def calculate_mean_z_score_distribution(z1, z2, num_reps=1000, zero_double_negatives=True, 
                                        overlap_control="remove", seed1=[], seed2=[]):
    """Determines size of expected mean combined `z=z1*z2` by randomly shuffling gene names

    Args:
        z1 (pd.Series, pd.DataFrame): Vector of z-scores from network propagation of trait 1
        z2 (pd.Series, pd.DataFrame): Vector of z-scores from network propagation of trait 2
        num_reps (int): Number of perumation analyses to perform. Defaults to 1000
        zero_double_negatives (bool, optional): Should genes that have a negative score in both `z1` and `z2` be ignored? Defaults to True.
        overlap_control (str, optional): 'bin' to permute overlapping seed genes separately, 'remove' to not consider overlapping seed genes. Any other value will do nothing. Defaults to "remove".
        seed1 (list, optional): List of seed genes used to generate `z1`. Required if `overlap_control!=None`. Defaults to [].
        seed2 (list, optional): List of seed genes used to generate `z2`. Required if `overlap_control!=None`. Defaults to [].

    Returns:
        float: The observed mean combined z-score from network colocalization
        list: List of permuted mean combined z-scores
    """
    #convert to correct format
    if isinstance(z1, pd.Series):
        z1 = pd.DataFrame(z1, columns=["z"])
    if isinstance(z2, pd.Series):
        z2 = pd.DataFrame(z2, columns=["z"])
    #combine table
    z1z2 = z1.join(z2, lsuffix="1", rsuffix="2")
    z1z2 = z1z2.assign(zz=z1z2.z1 * z1z2.z2)
    #print(z1z2.head())
    if overlap_control == "remove":
        seed_overlap = list(set(seed1).intersection(set(seed2)))
        print("Overlap seed genes:", len(seed_overlap))
        z1z2.drop(seed_overlap, axis=0, inplace=True)
        
    elif overlap_control == "bin":
        seed_overlap = list(set(seed1).intersection(set(seed2)))
        print("Overlap seed genes:", len(seed_overlap))
        overlap_z1z2 = z1z2.loc[seed_overlap]
        overlap_z1 = np.array(overlap_z1z2.z1)
        z1z2.drop(seed_overlap, axis=0, inplace=True)
    z1 = np.array(z1z2.z1)
    z2 = np.array(z1z2.z2)

    if z1z2.empty:
        raise ValueError("All genes removed after overlap control. Cannot proceed.")

    if zero_double_negatives:
        for node in z1z2.index:
            if (z1z2.loc[node].z1 < 0 and z1z2.loc[node].z2 < 0):
                z1z2.loc[node, 'zz'] = 0


    permutation_means = np.zeros(num_reps)
    i = 0
    while i < num_reps:
        include=True
        perm_z1z2 = np.zeros(len(z1))

        #shuffle rat seed genes
        np.random.shuffle(z1)
        #calculate NPScombined after rat seed shuffle
        for node in range(len(z1)):
            if not zero_double_negatives or not (z1[node] < 0 and z2[node] < 0):
                perm_z1z2[node] = z1[node] * z2[node]
            else:
                perm_z1z2[node] = 0
        if (np.isnan(perm_z1z2).any()):
            include=False
        if overlap_control == "bin":
            overlap_perm_z1z2 = np.zeros(len(overlap_z1))
            np.random.shuffle(overlap_z1) 
            for node in range(len(overlap_z1)):
                if zero_double_negatives and (overlap_z1[node] < 0 and z2[node] < 0):
                    overlap_perm_z1z2[node] = 0
                else:
                    overlap_perm_z1z2[node] = overlap_z1[node] * z2[node]
            perm_z1z2 = np.concatenate([perm_z1z2, overlap_perm_z1z2])
            if (np.isnan(perm_z1z2).any() or np.isnan(overlap_perm_z1z2).any()):
                include=False
                print('iteration skipped, value is NA')
        if include:
            permutation_means[i] = np.mean(perm_z1z2)
            i+=1   
                            
    return np.mean(z1z2.zz), permutation_means
def filter_go_annotations(go_df, term_min=10, term_max=5000, p_th=1e-5, min_intersection=3):
    """Filters available annotations for a community based on specificity and significance.
    Args:
        go_df (pandas.DataFrame): All available significant GO terms for each community
        term_min (int, optional):   The minimum size of a term to keep. Defaults to 50.
        term_max (int, optional): The maximum size of a term to keep. Defaults to 1000.
        p_th (float, optional): The significance threshold. Defaults to 1e-4.
        min_intersection (int, optional): Minimum number of community terms annotated to the GO term. Defaults to 3.

    Returns:
        pandas.DataFrame: A filter dataframe of GO annotations per community, sorted by sum of precision and recall.
    """
    go_df = go_df[(go_df['term_size'] <= term_max) & (go_df['term_size'] >= term_min)]
    go_df = go_df[go_df['intersection_size'] >= min_intersection]
    go_df = go_df[go_df['p_value'] < p_th] # set a stringent pvalue threshold
    go_df['sum_PR'] = go_df['recall'] + go_df['precision']
    go_df = go_df.sort_values('sum_PR',ascending=False)
    return go_df
def format_catalog(catalog=None):
	"""
	From RCA FUNCTIONS
	Processes and formats the GWAS catalog associations DataFrame by standardizing gene and trait names, 
	filtering relevant entries, and organizing data for easier querying and analysis.

	This function performs several operations to prepare GWAS catalog data for analysis:
	1. Converts 'MAPPED_TRAIT' and 'DISEASE/TRAIT' to lowercase for consistent querying.
	2. Filters out entries without mapped genes or traits.
	3. Splits gene entries that contain multiple genes listed together.
	4. Removes intergenic regions and entries labeled as 'mapped'.
	5. Combines trait information into a single column with PubMed ID references.

	Parameters:
	- catalog (DataFrame): A pandas DataFrame containing GWAS catalog data. If not provided, the function attempts to process, but will likely fail silently within the try-except block.

	Returns:
	DataFrame: A formatted DataFrame with each gene associated with its traits and citations.

	Raises:
	- Prints an error message if the input catalog is None or processing fails due to other issues.

	Notes:
	- This function assumes the input DataFrame contains specific columns: 'MAPPED_GENE', 'REPORTED GENE(S)', 'MAPPED_TRAIT', 'DISEASE/TRAIT', and 'PUBMEDID'.
	- The output DataFrame consolidates trait information into a single 'TRAIT' column and normalizes gene names.
	"""
	try:
		#make all annotations lowercase for consistency for querying
		catalog['MAPPED_TRAIT']=catalog['MAPPED_TRAIT'].str.lower()
		catalog['DISEASE/TRAIT']=catalog['DISEASE/TRAIT'].str.lower()
		#filter for genes that were mapped
		mapped=catalog[~catalog['MAPPED_GENE'].isna()]
		mapped=mapped[~mapped['MAPPED_TRAIT'].isna()]
		mapped=mapped[['MAPPED_GENE','MAPPED_TRAIT','DISEASE/TRAIT','PUBMEDID']]
		mapped.columns=['GENE','MAPPED_TRAIT','DISEASE/TRAIT','PUBMEDID']
		#filter for genes that were reported
		rep=catalog[~catalog['REPORTED GENE(S)'].isna()]
		rep=rep[~rep['MAPPED_TRAIT'].isna()]
		rep=rep[~rep['REPORTED GENE(S)'].str.contains('Intergenic')]
		rep=rep[['REPORTED GENE(S)','MAPPED_TRAIT','DISEASE/TRAIT','PUBMEDID']]
		rep.columns=['GENE','MAPPED_TRAIT','DISEASE/TRAIT','PUBMEDID']
		cat=pd.concat([rep, mapped])
		cat['GENE']=cat['GENE'].str.split('; ')
		cat=cat.explode('GENE')
		cat=cat[~(cat['GENE'].str.contains('mapped'))]
		cat['GENE']=cat['GENE'].str.split(', ')
		cat=cat.explode('GENE')
		cat['GENE']=cat['GENE'].str.split(' - ')
		cat=cat.explode('GENE')
		cat['GENE']=cat['GENE'].astype('str')
		cat=cat[~(cat['GENE'].str.contains('intergenic'))]
		cat['TRAIT']=cat['MAPPED_TRAIT'] + ": " +cat['DISEASE/TRAIT']+ " (PMID: "+(cat['PUBMEDID'].astype(str))+")"
		cat=cat.dropna()
		return(cat)
	except:
		print('please add gwas catalog file.')

########################################