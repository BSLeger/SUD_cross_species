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
	'signor_mouse':'656370fa-afe8-11e9-8bb4-0ac135e8bacf',
	'signor_rat_protein_only':'aa57a66c-0842-11f0-9806-005056ae3c32'
}
ctrl_traits=['facial_hair', 'age_smkinit', 'antisoc', 'hr', 'infant_bw', 'LDL', 'maternal_smok', 'age_menarche','addict-rf','adhd', 'dpw', 'risk', 'auto_speed', 'nsex', 'bmi', 'height']

ctrl_traits_rat=['bmi_rn6','body_length_rn6']

mag_dir='magma/seed_genes/'

ctrl_mag_dir='gwas_ctrl_hm/magma/seed_genes/'
file_dict={
    'loco':mag_dir+'loco_win10_annot.tsv',
    'loco_gsem':mag_dir+'loco_gsem_annot.tsv',
    'ext':mag_dir+'ext_orig_annot.tsv',
    'ext_st22':mag_dir+'all_tests_ext1_st22_genes.csv',
    'loco_mega_fus_naac':'loco_twas_dan/loco_fusion_NACC_seed.tsv',
    'ext_fus_naac':'ext_FUSION/ext_fusion_NACC_seed.tsv',
	'loco_final_cf':mag_dir+'loco_final_cf_win10_annot.tsv',
	'loco_final_mega':mag_dir+'loco_final_mega_win10_annot.tsv',
	'ext_rat':mag_dir+'ext_orig_annot_rat_ortho.tsv',
	'loco_final_cf_rat':mag_dir+'loco_final_cf_win10_annot.tsv',
	'ext_db':'ext_2factor/MAGMA_v108_DB_bonf.tsv',
	'ext_rtb':'ext_2factor/MAGMA_v108_RTB_bonf.tsv',
	'loco_final_cf_25':mag_dir+'loco_final_cf_win25_annot.tsv'
}
bonf_dict={
    'loco_gsem':2.650129856362962e-06,
    'loco':2.6389402016150313e-06,
	'loco_mega_fus_naac':9.338812103100487e-06,
	'loco_final_cf':2.635601707869907e-06,
	'loco_final_mega':2.6467630088401888e-06,
	'ext_rat':2.7003672499459928e-06,
	'ext_db':2.75e-06,
	'ext_rtb':2.75e-06,
	'loco_final_cf_25':2.6467630088401888e-06,
	'body_length_rn6':2.8370403994552885e-06,
	'bmi_rn6':2.8370403994552885e-06
}
gene_col_dict={
    'loco':'HM_ORTHO',
    'loco_gsem':'HM_ORTHO',
	'loco_mega_fus_naac':'human_ortholog',
    'ext':'GENE',
	'ext_fus_naac':'ID',
	'ext_st22':'GENE NAME',
	'loco_final_cf':'HM_ORTHO',
	'loco_final_mega':'HM_ORTHO',
	'loco_final_cf_rat':'GENE',
	'ext_rat':'Gene1Symbol',
	'ext_db':'Gene symbol',
	'ext_rtb':'Gene symbol',
	'loco_final_cf_25':'HM_ORTHO',
	'hm_ctrl':'GENE',
	'rat_ctrl':'HM_ORTHO'
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
#from the original tissue dictionary
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
'Liver':'Liver',
'AC':'Nucleus accumbens core',
'LH':'Lateral habenula',
'VO':'Orbitofrontal cortex'}

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



def import_seed_dict(mag_dir,file_dict,ctrl_traits,ctrl_traits_rat,bonf_dict,gene_col_dict,all_nodes):
    #written for MAGMA output- need to rewrite for fusion or ratXcan
    seed_dict={}
    for f in file_dict.keys():
        t=pd.read_csv(file_dict[f],sep='\t')
        gene_col=gene_col_dict[f]
        #print('successfully read in file')
        #print(t.head())
        if f in bonf_dict.keys():
            bonf_cutoff=bonf_dict[f]
        else:
            bonf_cutoff=0.05/len(t)
        #print(f'bonferroni cuttoff ={bonf_cutoff}')
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
    for f in ctrl_traits:
        gene_col=gene_col_dict['hm_ctrl']
        t=pd.read_csv(f'{ctrl_mag_dir}{f}_annot.tsv',sep='\t')
        seed_dict[f'{f}_FDR']=(set(t[t['Q']<0.05][gene_col]))
        seed_dict[f'{f}_bonf']=(set(t[t['P']<(0.05/len(t))][gene_col]))
        seed_dict[f'{f}_top500']=set(t[(t['GENE'].isin(all_nodes))].nsmallest(500,'P')[gene_col])
    for f in ctrl_traits_rat:
        gene_col=gene_col_dict['rat_ctrl']
        t=pd.read_csv(f'{mag_dir}{f}_annot.tsv',sep='\t')
        seed_dict[f'{f}_FDR']=(set(t[t['Q']<0.05][gene_col]))
        seed_dict[f'{f}_bonf']=(set(t[t['P']<bonf_dict[f]][gene_col]))
        seed_dict[f'{f}_top500']=set(t[(t[gene_col].isin(all_nodes))].nsmallest(500,'P')[gene_col])

    return seed_dict


'''def import_NPS_scores(seed_dict,UUIDs):
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
    return NPS_dict, NPS_dict_series'''
def import_NPS_scores(seed_dict,interactome_name):
    NPS_dict_series={}
    for k in seed_dict.keys():
        p=('network_scores/'+k+'_'+interactome_name+'_zscore.tsv')
        if os.path.isfile(p):
            t=pd.read_csv('network_scores/'+k+'_'+interactome_name+'_zscore.tsv',header=None, sep='\t')
            t.index=t[0]
            t=t.drop(columns=[0])
			#t=t[1].squeeze()
			#t = pd.DataFrame({'z':t})
            NPS_dict_series[k+'_'+interactome_name]=t
    NPS_dict={}
    for k in NPS_dict_series.keys():
        t=NPS_dict_series[k]
        t=t[1].squeeze()
        t = pd.DataFrame({'z':t})
        NPS_dict[k]=t
    return NPS_dict, NPS_dict_series

def return_analysis_datasets(trait_r,cut_r,trait_h,cut_h,seed_dict,NPS_dict,interactome_name):
    #labels
	if not (trait_h==None):
	    if cut_h==None:
	        label_h=trait_h
	    else:
	        label_h=trait_h+'_'+cut_h
	    seed_h=seed_dict[label_h]
	    NPS_h=NPS_dict[label_h+'_'+interactome_name]
	else:
		trait_h=None
		label_h=None
		cut_h=None
		seed_h=None
		NPS_h=None

	if not (trait_r==None):
	    if cut_r==None:
	        label_r=trait_r
	    else:
	        label_r=trait_r+'_'+cut_r
	    seed_r=seed_dict[label_r]
	    NPS_r=NPS_dict[label_r+'_'+interactome_name]

	else:
		trait_r=None
		label_r=None
		cut_r=None
		seed_r=None
		NPS_r=None
    #seed genes
    #NPS scores
	if ((trait_h!=None) and (trait_r!=None)):
	    NPS = NPS_h.join(NPS_r, lsuffix="h", rsuffix="r")
	    NPS = NPS.assign(zhr=NPS.zh * NPS.zr)
	else:
	    NPS=None
	return label_h,label_r,seed_h,seed_r,NPS_h,NPS_r,NPS
	
def format_network(network, traitr, traitc, seedr, seedc,zr, zc):
    """
    Formats the colocalized network for easy secondary analysis in Cytoscape
    This function takes a network (graph) and updates its nodes with several attributes:
    - Seed gene indicators for two traits with secondary "color scheme" indicator for use in cytoscape color scheme assignment
    - Z-scores for each trait and their combination

    Parameters:
    - network (NetworkX graph): The graph to be formatted, where nodes represent genes.
    - traitr (str): Identifier for the rare trait.
    - traitc (str): Identifier for the common trait.
    - seedr (list): List of seed genes associated with the rare trait.
    - seedc (list): List of seed genes associated with the common trait.
    - zr (dict): Dictionary mapping genes to their z-scores for the rare trait.
    - zc (dict): Dictionary mapping genes to their z-scores for the common trait.

    Returns:
    - NetworkX graph: The original network updated with node attributes for seed gene status and z-scores for traits and their combination.

    Note:
    The function utilizes pandas DataFrames for intermediate data manipulation and requires 
    NetworkX for working with the network. It expects `network.nodes()` to return a list-like 
    object of genes.
    """
    nodes_df=pd.DataFrame(network.nodes())
    nodes_df.columns=['Gene']
    #make node seed gene dataframe from which to make dictionaries
    nodes_df[('seed_'+traitr)]=nodes_df['Gene'].isin(seedr)
    nodes_df[('seed_'+traitc)]=nodes_df['Gene'].isin(seedc)
    nodes_df['seed_both']=(nodes_df['Gene'].isin(seedr) & nodes_df['Gene'].isin(seedc))
    
    nodes_df['color_scheme']=0
    nodes_df['color_scheme'] = np.where(nodes_df[('seed_'+traitr)] == True, 1, nodes_df['color_scheme'])
    nodes_df['color_scheme'] = np.where(nodes_df[('seed_'+traitc)] == True, 2, nodes_df['color_scheme'])
    nodes_df['color_scheme'] = np.where((nodes_df["seed_both"] == True), 3, nodes_df['color_scheme'])    
    nodes_df.index=nodes_df['Gene']
    #set zscores as node attributes
    nx.set_node_attributes(network, dict(zr), ('z_'+traitr))
    nx.set_node_attributes(network, dict(zc),('z_'+traitc))
    nx.set_node_attributes(network, dict(zr*zc), 'z_comb')
    #add seed genes as node attributes
    nx.set_node_attributes(network,dict(zip(nodes_df['Gene'], nodes_df[('seed_'+traitr)])), ('seed_'+traitr))
    nx.set_node_attributes(network,dict(zip(nodes_df['Gene'], nodes_df[('seed_'+traitc)])), ('seed_'+traitc))
    nx.set_node_attributes(network,dict(zip(nodes_df['Gene'], nodes_df['seed_both'])), ('seed_both'))
    nx.set_node_attributes(network,dict(zip(nodes_df['Gene'], nodes_df['color_scheme'])), ('seed_color_scheme'))
    return(network)
    
def export_network(network, name, user, password, ndex_server='public.ndexbio.org'):
    '''
    shell for net_cx upload network function, that creates nicecx network, then exports to NDEx in the CX format.

    This function converts a NetworkX graph to a NiceCXNetwork object using the ndex2 Python package,
    sets the network's name, and uploads it to the specified NDEx server. Upon successful upload, 
    the function returns the UUID of the network in the NDEx platform.

    Parameters:
    - network (NetworkX graph): The graph to be exported.
    - name (str): The name to assign to the network in NDEx.
    - user (str): NDEx account username.
    - password (str): NDEx account password.
    - ndex_server (str, optional): The URL of the NDEx server to which the network is to be uploaded. Defaults to 'public.ndexbio.org'.

    Returns:
    - exports network to NDEx
    - str: The UUID of the uploaded network on NDEx.
    Notes:
    Tequires the ndex2 package. Ensure that you have a valid NDEx account and that the specified server URL is correct.
    """
    '''
    print(user)
    print(password)
    if ((user==None) | (password==None)):
        print('please provide a NDEx username and password.')
    else:
        net_cx = ndex2.create_nice_cx_from_networkx(network)
        net_cx.set_name(name)
        network_uuid = net_cx.upload_to(ndex_server, user, password)

## Utilities for systems map- modified from human rat bmi -------------------------------------------------------------
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
    comm_genes[f'{seed1_name}_ratio'] = [1 if x in seeds1 else 0 for x in comm_genes.CD_MemberList]
    comm_genes[f'{seed2_name}_ratio'] = [1 if x in seeds2 else 0 for x in comm_genes.CD_MemberList]
    comm_genes["seed_hr_ratio"] = comm_genes.apply(lambda x: x[f'{seed1_name}_ratio'] * x[f'{seed2_name}_ratio'], axis=1)
    a = comm_genes.groupby(level=0).seed_hr_ratio.sum()
    b = comm_genes[comm_genes.seed_hr_ratio != 1].groupby(level=0)[f'{seed1_name}_ratio'].sum()
    c = comm_genes[comm_genes.seed_hr_ratio != 1].groupby(level=0)[f'{seed2_name}_ratio'].sum()
    d = comm_genes.groupby(level=0).CD_MemberList.count()
    counts = pd.concat([a,b,c,d], axis=1)
    counts["network_ratio"] = counts.apply(lambda x: x.CD_MemberList - x.seed_hr_ratio - x[f'{seed1_name}_ratio'] - x[f'{seed2_name}_ratio'], axis=1)
    fracs = counts.div(counts.CD_MemberList, axis=0)    
    return fracs
def get_gene_NPS_scores(hier_df,NPS_r,NPS_h,cut_r=1.5,cut_h=1.5,cut_hr=3):
    hier_df["CD_MemberList"] = hier_df.CD_MemberList.apply(lambda x: x if type(x)==list else x.split(" "))
    comm_genes = hier_df.explode("CD_MemberList")
    comm_genes['NPSr']=comm_genes.CD_MemberList.apply(lambda x: NPS_r.z[x])
    comm_genes['NPSh']=comm_genes.CD_MemberList.apply(lambda x: NPS_h.z[x])
    comm_genes['NPShr']=comm_genes['NPSh']*comm_genes['NPSr']
    
    NPS=NPS_r.merge(NPS_h,left_index=True,right_index=True)
    NPS.columns=['NPSr','NPSh']
    NPS['NPShr']=NPS.NPSh*NPS.NPSr
    netgenes=list(NPS[(NPS.NPSh>cut_h)&(NPS.NPSr>cut_r)&(NPS.NPShr>cut_hr)].index)
    
    #median
    medr=comm_genes.groupby(level=0)['NPSr'].median()
    medh=comm_genes.groupby(level=0)['NPSh'].median()
    medhr=comm_genes.groupby(level=0)['NPShr'].median()
    
    #median_normalized
    medr_n=comm_genes.groupby(level=0)['NPSr'].median()/NPS.NPSr[netgenes].median()
    medh_n=comm_genes.groupby(level=0)['NPSh'].median()/NPS.NPSh[netgenes].median()
    medhr_n=comm_genes.groupby(level=0)['NPShr'].median()/NPS.NPShr[netgenes].median()
    
    
    meanr=comm_genes.groupby(level=0)['NPSr'].mean()
    meanh=comm_genes.groupby(level=0)['NPSh'].mean()
    meanhr=comm_genes.groupby(level=0)['NPShr'].mean()
    
    #mean_normalized
    meanr_n=comm_genes.groupby(level=0)['NPSr'].mean()/NPS.NPSr[netgenes].mean()
    meanh_n=comm_genes.groupby(level=0)['NPSh'].mean()/NPS.NPSh[netgenes].mean()
    meanhr_n=comm_genes.groupby(level=0)['NPShr'].mean()/NPS.NPShr[netgenes].mean()
    
    c = pd.concat([medr,medh,medhr,medr_n,medh_n,medhr_n,meanr,meanh,meanhr,meanr_n,meanh_n,meanhr_n], axis=1)
    c.columns=['medr','medh','medhr','medr_n','medh_n','medhr_n','meanr','meanh','meanhr','meanr_n','meanh_n','meanhr_n']
    return c
def get_expanded_network_gene_fractions(hier_df, hm_genes, rat_genes, coloc_genes, hm_genes_name='hm_net', rat_genes_name='rat_net', coloc_genes_name='coloc_net'):
    """Assess the number of genes in each community that were from species-specific networks
    Args:
        hier_df (pd.DataFrame): The hierarchy information of gene communities
        hm_genes (list): List of genes from input 1 (typically human network)
        rat_genes (list): List of genes from input 2 (typically rat network)
		coloc_genes (list): List of genes from input 2 (typically colocalized network)
        hm_genes_name (str, optional): Name to give the seed genes from input 1. Defaults to 'hm_net'.
        rat_genes_name (str, optional): Name to give the seed genes from input 2. Defaults to 'rat_net'.
		coloc_genes_name (str,optional): Name to give the seed genes from input 3. Defaults to 'coloc_net'.

    Returns:
        pd.DataFrame: The fraction of genes in each community that are seeds in both hm_genes and rat_genes, the fraction just in hm_genes, the fraction just in rat_genes
    """
    hier_df["CD_MemberList"] = hier_df.CD_MemberList.apply(lambda x: x if type(x)==list else x.split(" "))
    comm_genes = hier_df.explode("CD_MemberList")
    
    comm_genes[f'{hm_genes_name}_ratio'] = [1 if x in hm_genes else 0 for x in comm_genes.CD_MemberList]
    comm_genes[f'{rat_genes_name}_ratio'] = [1 if x in rat_genes else 0 for x in comm_genes.CD_MemberList]
    comm_genes[f'{coloc_genes_name}_ratio'] = [1 if x in coloc_genes else 0 for x in comm_genes.CD_MemberList]    
    
    a=comm_genes.groupby(level=0)[f'{hm_genes_name}_ratio'].sum()
    b=comm_genes.groupby(level=0)[f'{rat_genes_name}_ratio'].sum()
    c=comm_genes.groupby(level=0)[f'{coloc_genes_name}_ratio'].sum()
    d=comm_genes.groupby(level=0).CD_MemberList.count()
    counts=pd.concat([a,b,c,d], axis=1)
    fracs=counts.div(counts.CD_MemberList, axis=0)    
    return fracs