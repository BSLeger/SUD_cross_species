
import os
import pandas as pd
import ndex2
import networkx as nx
from netcoloc import netprop_zscore
from netcoloc import netprop
from netcoloc import network_colocalization
import sys
import random


sys.path.insert(0, '/tscc/projects/ps-palmer/brittany/rare_common_alcohol/rare_common_alcohol_comparison/notebooks/')
#os.chdir('/tscc/projects/ps-palmer/brittany/rare_common_alcohol/rare_common_alcohol_comparison/notebooks/')
print(os.getcwd())
from rca_functions import *


os.chdir('/tscc/projects/ps-palmer/brittany/SUD_cross_species/')

random_seed=random.seed(211)

save_file=True


UUIDs={
    'PCNet2.0':'d73d6357-e87b-11ee-9621-005056ae23aa',
    'PCNet2.1':'e9c574f4-e87a-11ee-9621-005056ae23aa',
    'PCNet2.2':'8b4b54fa-e87d-11ee-9621-005056ae23aa'
}

# # functions

def import_seedgenes(path,pcol='P',gene_col='GENE NAME',delim='comma', cutoff=None):
    if delim=='comma':
        df=pd.read_csv(path,sep=',')
    else:
        df=pd.read_csv(path,sep='\t')
    if pcol==None:
        print('pvalue column not specified- all genes will be used')
        cutoff=None
    if cutoff=='bonferroni':
        df=df[df[pcol]<0.05/len(df)]
    elif cutoff=='FDR_05':
        df=df[df[pcol]<0.05]
    else:
        print('cutoff not defined/custom- using all genes ')
        df=df
    print(df.head())
    return(df)


def import_interactome(interactome_name=None, UUIDs=UUIDs,ndex_user=None, ndex_password=None,UUID=None):
    """
    Imports a gene interactome from the NDEx database and returns it as a NetworkX graph object. Optionally,
    the function allows for importing using a unique identifier (UUID) or by an interactome name.

    The function checks if the interactome name provided corresponds to a predefined dictionary of UUIDs. If it does, it
    retrieves the network using the specified credentials. If an interactome name is not provided but a UUID is,
    it retrieves the network using the provided UUID. The nodes of pcnet_v14 are relabelled by their gene name rather than ID number.

    Parameters:
    - interactome_name (str, optional): The name of the interactome as defined in the UUIDs dictionary. If not provided
      but a UUID is, the interactome associated with the UUID is imported instead.
    - ndex_user (str, optional): The NDEx account username for accessing private networks.
    - ndex_password (str, optional): The NDEx account password for accessing private networks.
    - UUID (str, optional): A specific UUID to directly download an interactome from NDEx if the interactome name is not used.

    Returns:
    networkx.Graph: A graph object representing the interactome. Nodes and edges represent genes and their interactions, respectively.

    Notes:
    - The function uses the NDEx2 Python client and requires Internet access to NDEx's servers.
    - Depending on the access rights of the NDEx account, private or public interactomes can be retrieved.
    - The function prints the number of nodes and edges of the imported graph for diagnostic purposes.

    Raises:
    - ValueError: If neither an interactome name nor a UUID is provided.
    """    
    interactome_uuid=UUIDs[interactome_name]
    print(interactome_name)
    ndex_server='public.ndexbio.org'
    #import network based on provided interactome key
    if (interactome_name in UUIDs.keys()):
        graph = ndex2.create_nice_cx_from_server(
                    ndex_server, 
                    username=ndex_user, 
                    password=ndex_password, 
                    uuid=interactome_uuid
                ).to_networkx()
        if (interactome_name=='pcnet_v14'):
            graph=nx.relabel_nodes(graph, nx.get_node_attributes(graph, 'HGNC Symbol'))
        # print out interactome num nodes and edges for diagnostic purposes
        print('number of nodes:')
        print(len(graph.nodes))
        print('\nnumber of edges:')
        print(len(graph.edges))
        return(graph)
    elif(interactome_name==None & UUID!=None):
        print('using novel UUID. For UUIDs used in this study, see UUID_dict')
        graph = ndex2.create_nice_cx_from_server(
            ndex_server, 
            username=ndex_user, 
            password=ndex_password, 
            uuid=UUID
        ).to_networkx()
        # print out interactome num nodes and edges for diagnostic purposes
        print('number of nodes:')
        print(len(graph.nodes))
        print('\nnumber of edges:')
        print(len(graph.edges))
        return(graph)
    else:
        print('UUID/interactome name not provided- please provide either to import interactome.')


#from rat bmi notebooks not netcoloc
def calculate_heat_zscores_with_sampling(data, nodes, individual_heats, G_PC, trait="BMI", max_genes=500, num_samples=100,
                                        nominal_sig=0.05, num_reps=1000, out_path="", minimum_bin_size=10,outfile='sample.tsv'):
    """Takes a set of summary statistics and a molecular interaction and performs sampling of the significant genes.
    For each sample a random selection of seed genes is chosen, weighted by the p-value of each gene in the summary
    statistics. Network propagation with zscore calculation is performed for each sample to generate a distribution
    of z-scores for each gene in the seed_gene set.

    Args:
        data (pd.DataFrame): Gene level summary statistics
        nodes (list): list of nodes in the interaction network
        individual_heats (np.array): Heat matrix calculated by `netprop_zscore.get_individual_heats_matrix()`
        G_PC (nx.Graph): molecular interaction network
        trait (str, optional): name of trait being investigated. Defaults to "BMI".
        max_genes (int, optional): Maximum number of seed genes to include in each sample (maximum=500). Defaults to 500.
        num_samples (int, optional): Number of times to perform sampling. Defaults to 100.
        nominal_sig (float, optional): Significance cutoff for keeping genes in data (Note: this value will be Bonferroni corrected). Defaults to 0.05.
        num_reps (int, optional): Number of repetitions of randomization for generating null distribution for z_scores. Defaults to 1000.
        out_path (str, optional): File path prefix for saving results of sampling. Defaults to "".
        minimum_bin_size (int, optional): minimum number of genes that should be in each degree matching bin. Defaults to 10.

    Returns:
        pd.DataFrame: Gene x sampling run dataframe of sampled z-scores
    """
    assert max_genes <= 500, "NetColoc is only valid for sets of 500 or less genes so maximum number of genes for sampling must be <= 500"
    #outfile = out_path + trait + "sampling_" + str(max_genes) + "_" + str(num_samples) + ".tsv"
    data = data.loc[data.gene_symbol.isin(nodes)]  # subset to genes in interaction network
    all_seeds = data.loc[data.pvalue <= nominal_sig / len(data)]  # Bonferroni correction
    all_seeds = all_seeds.assign(log10p=-1 * np.log10(all_seeds.pvalue))  # get -log10p for weighted sampling
    sampling_results = []
    for i in tqdm(range(num_samples)):
        # perform propagation for sample
        sample_seeds = random.choices(population=all_seeds.gene_symbol.values, weights=all_seeds.log10p.values, k=max_genes)
        sample_results = netprop_zscore.calculate_heat_zscores(individual_heats, nodes=list(G_PC.nodes), degrees=dict(G_PC.degree),
                                                seed_genes=sample_seeds, num_reps=num_reps,
                                                minimum_bin_size=minimum_bin_size, random_seed=i)[0]
        sample_z = pd.DataFrame(sample_results, columns=["z" + str(i)])
        # save running results of sampling
        if i == 0:
            sample_z.to_csv(outfile, sep="\t")
        else:
            existing = pd.read_csv(outfile, sep="\t", index_col=0)
            existing = existing.join(sample_z)
            existing.to_csv(outfile, sep="\t")
        sampling_results.append(sample_z)
    return pd.concat(sampling_results, axis=1)

## import datasets
os.getcwd()


mag_dir='magma/seed_genes/'
file_dict={
    'loco':mag_dir+'loco_win10_annot.tsv',
    'ext_munged':mag_dir+'ext_munged_annot.tsv',
    'ext':mag_dir+'ext_orig_annot.tsv',
    'ext_st22':mag_dir+'all_tests_ext1_st22_genes.csv'
}


ext=pd.read_csv(file_dict['ext'],sep='\t')
ext_sub=ext[['GENE','P']]
ext_sub.columns=['gene_symbol','pvalue']
k='ext_subsample'

print(ext.head())


# # Interactome Set-up

# pcnet2- versions 
# from wright et al. 2024 preprint:
# PCNet 2.0= best-performing ranked composite (top 15 interactomes, 3.85M interactions)
# PCNet 2.1= top 8 interactomes, 1.75M interactions
# PCNet 2.2= top 10 co-citation-free interactomes, 3.32M interactions 

interactome_name='PCNet2.0'


graph=import_interactome(interactome_name)


all_nodes=list(graph.nodes())


# pre calculate the matricies used for network propagation
print('\ncalculating w_prime')
w_prime = netprop.get_normalized_adjacency_matrix(graph, conserve_heat=True)

print('\ncalculating w_double_prime')
w_double_prime = netprop.get_individual_heats_matrix(w_prime, .5)


# # calculate gwas NPS



print('starting propagation')
NPSc, Fnew_score, Fnew_rand_score = calculate_heat_zscores_with_sampling(
    data=ext_sub,
    nodes=all_nodes,
    individual_heats=w_double_prime,
    G_PC=graph,
    trait=k,
    outfile='network_scores/'+k+'_'+interactome_name+'_zscore.tsv'
)

