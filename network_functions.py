import pandas as pd
import numpy as np
import os
import ndex2
import networkx as nx


UUIDs={
    'PCNet2.0':'d73d6357-e87b-11ee-9621-005056ae23aa',
    'PCNet2.1':'e9c574f4-e87a-11ee-9621-005056ae23aa',
    'PCNet2.2':'8b4b54fa-e87d-11ee-9621-005056ae23aa',
	'signor_rat':'76be57cd-afe8-11e9-8bb4-0ac135e8bacf',
	'signor_human':'523fff27-afe8-11e9-8bb4-0ac135e8bacf',
	'signor_mouse':'656370fa-afe8-11e9-8bb4-0ac135e8bacf'
}
mag_dir='magma/seed_genes/'
file_dict={
    'loco':mag_dir+'loco_win10_annot.tsv',
    'loco_gsem':mag_dir+'loco_gsem_annot.tsv',
    'ext':mag_dir+'ext_orig_annot.tsv',
    'ext_st22':mag_dir+'all_tests_ext1_st22_genes.csv',
    'loco_mega_fus_naac':'loco_twas_dan/loco_fusion_NACC_seed.tsv',
    'ext_fus_naac':'ext_FUSION/ext_fusion_NACC_seed.tsv'
}
bonf_dict={
    'loco_gsem':2.650129856362962e-06,
    'loco':2.6389402016150313e-06,
	'loco_mega_fus_naac':9.338812103100487e-06,
}
gene_col_dict={
    'loco':'HM_ORTHO',
    'loco_gsem':'HM_ORTHO',
	'loco_mega_fus_naac':'human_ortholog',
    'ext':'GENE',
	'ext_fus_naac':'ID',
	'ext_st22':'GENE NAME'
}
# define network cutoffs
cut_single=1.5
cut_comb=3
cut_rat_specific={
    'zr':3,
    'zh':0,
    'zhr':0
}
cut_hm_specific={
    'zr':0,
    'zh':3,
    'zhr':0
}
def import_interactome(interactome_name=None, UUIDs=UUIDs,ndex_user=None, ndex_password=None, UUID=None):
    """
    Imports a gene interactome from the NDEx database and returns it as a NetworkX graph object.
    Optionally, the function allows for importing using a unique identifier (UUID) or by an interactome name.
    
    Parameters:
    - interactome_name (str, optional)
	- UUIDs (str,optional): dictionary of UUIDs
    - ndex_user (str, optional)
    - ndex_password (str, optional)
    - UUID (str, optional)

    Returns:
    networkx.Graph: A graph object representing the interactome.
    
    Raises:
    - ValueError: If neither an interactome name nor a UUID is provided.
    """
    
    ndex_server = 'public.ndexbio.org'
    # import network based on provided interactome key
    if interactome_name in UUIDs.keys():
        interactome_uuid = UUIDs[interactome_name]
        print(interactome_name)
        g = ndex2.create_nice_cx_from_server(
            ndex_server, 
            username=ndex_user, 
            password=ndex_password, 
            uuid=interactome_uuid
        )
        g.print_summary()
        graph=g.to_networkx()
        if interactome_name == 'pcnet_v14':
            graph = nx.relabel_nodes(graph, nx.get_node_attributes(graph, 'HGNC Symbol'))
        
        # print out interactome num nodes and edges for diagnostic purposes
        print('number of nodes:')
        print(len(graph.nodes))
        print('\nnumber of edges:')
        print(len(graph.edges))
        return graph

    elif interactome_name is None and UUID is not None:
        print('using novel UUID. For UUIDs used in this study, see UUIDs')
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
        return graph

    else:
        print('UUID/interactome name not provided- please provide either to import interactome.')


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



def import_seed_dict(mag_dir,file_dict,bonf_dict,gene_col_dict,all_nodes):
    #written for MAGMA output- need to rewrite for fusion or ratXcan
    seed_dict={}
    for f in file_dict.keys():
        t=pd.read_csv(file_dict[f],sep='\t')
        gene_col=gene_col_dict[f]
        if f in bonf_dict.keys():
            bonf_cutoff=bonf_dict[f]
        else:
            bonf_cutoff=0.05/len(t)
        if ('fus' in f):
            Pcol='TWAS.P'
        else:
            Pcol='P'
        if (f=='ext_st22'):
            seed_dict[f]=(set(t[gene_col]))
        else:
            try:
                seed_dict[f'{f}_bonf']=(set(t[t[Pcol]<bonf_cutoff][gene_col]))
                seed_dict[f'{f}_top500']=set(t[(t[gene_col].isin(all_nodes))].nsmallest(500,Pcol)[gene_col])
                seed_dict[f'{f}_FDR']=(set(t[t['Q']<0.05][gene_col]))
            except:
                print(f'error occurred importing {f}')
    return seed_dict


def import_NPS_scores(seed_dict,UUIDs):
    NPS_dict_series={}
    for k in seed_dict.keys():
        for u in UUIDs.keys():
            p=('network_scores/'+k+'_'+u+'_zscore.tsv')
            if os.path.isfile(p):
                t=pd.read_csv('network_scores/'+k+'_'+u+'_zscore.tsv',header=None, sep='\t')
                t.index=t[0]
                t=t.drop(columns=[0])
                #t=t[1].squeeze()
                #t = pd.DataFrame({'z':t})
                NPS_dict_series[k+'_'+u]=t
    NPS_dict={}
    for k in NPS_dict_series.keys():
        t=NPS_dict_series[k]
        t=t[1].squeeze()
        t = pd.DataFrame({'z':t})
        NPS_dict[k]=t
    return NPS_dict, NPS_dict_series


def return_analysis_datasets(trait_r,cut_r,trait_h,cut_h,seed_dict,NPS_dict,interactome_name):
    #labels
    if cut_h==None:
        label_h=trait_h
    else:
        label_h=trait_h+'_'+cut_h
    if cut_r==None:
        label_r=trait_r
    else:
        label_r=trait_r+'_'+cut_r
    #seed genes
    seed_r=seed_dict[label_r]
    seed_h=seed_dict[label_h]
    #NPS scores
    NPS_r=NPS_dict[label_r+'_'+interactome_name]
    NPS_h=NPS_dict[label_h+'_'+interactome_name]
    NPS = NPS_h.join(NPS_r, lsuffix="h", rsuffix="r")
    NPS = NPS.assign(zhr=NPS.zh * NPS.zr)
    return label_h,label_r,seed_h,seed_r,NPS_h,NPS_r,NPS

## Utilities for systems map- from human rat bmi -------------------------------------------------------------
def get_seed_gene_fractions(hier_df, seeds1, seeds2, seed1_name='h_seed', seed2_name='r_seed'):
    """Assess the number of genes in each community that were seed genes from the orginal inputs

    Args:
        hier_df (pd.DataFrame): The hierarchy information of gene communities
        seeds1 (list): List of genes from input 1
        seeds2 (list): List of genes from input 2
        seed1_name (str, optional): Name to give the seed genes from input 1. Defaults to 'h_seed'.
        seed2_name (str, optional): Name to give the seed genes from input 2. Defaults to 'r_seed'.

    Returns:
        pd.DataFrame: The fraction of genes in each community that are seeds in both seeds1 and seeds2, the fraction just in seeds1, the fraction just in seeds2
    """
    hier_df["CD_MemberList"] = hier_df.CD_MemberList.apply(lambda x: x if type(x)==list else x.split(" "))
    comm_genes = hier_df.explode("CD_MemberList")
    comm_genes[seed1_name] = [1 if x in seeds1 else 0 for x in comm_genes.CD_MemberList]
    comm_genes[seed2_name] = [1 if x in seeds2 else 0 for x in comm_genes.CD_MemberList]
    comm_genes["overlap"] = comm_genes.apply(lambda x: x[seed1_name] * x[seed2_name], axis=1)
    a = comm_genes.groupby(level=0).overlap.sum()
    b = comm_genes[comm_genes.overlap != 1].groupby(level=0)[seed1_name].sum()
    c = comm_genes[comm_genes.overlap != 1].groupby(level=0)[seed2_name].sum()
    d = comm_genes.groupby(level=0).CD_MemberList.count()
    counts = pd.concat([a,b,c,d], axis=1)
    counts["network"] = counts.apply(lambda x: x.CD_MemberList - x.overlap - x[seed1_name] - x[seed2_name], axis=1)
    fracs = counts.div(counts.CD_MemberList, axis=0)
    return fracs