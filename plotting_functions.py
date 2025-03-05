import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, hypergeom
import matplotlib.pyplot as plt
import scipy.stats as stats
from matplotlib_venn import venn2
import scipy.stats as stats
import networkx as nx


IBM=['#648FFF','#785EF0','#DC267F','#FE6100','#FFB000']
colour_dict={
    'ext':IBM[1],
    'ext_alt':IBM[0],
    'ext_munged':IBM[1],
    'ext_munged_alt':IBM[0],
    'loco':IBM[4],
    'loco_alt':IBM[3],
    'loco_gsem':IBM[4],
    'loco_gsem_alt':IBM[3],
    'shared':IBM[2],
    'other':'#CCCCCC',
    'shared_alt':'#780534',
	'addict-rf':'#A8E5A0',
	'addict-rf_alt':'#63AC00',
	'shared_addict-rf':'#31739C',
	'shared_alt_addict-rf':'#56B4E9'
}

def plt_scatter_NPS(tblr, tblc, tblr_label, tblc_label, tblr_seed, tblc_seed,
                    colour_r=colour_dict['loco'], colour_c=colour_dict['ext'], colour_shared=colour_dict['shared'],colour_nonseed=colour_dict['shared'],
                    tblr_lim=1.5, tblc_lim=1.5, comb_lim=3, savefig=False, filename='scatter_NPS'):
    """
    Visualizes the NPScommon and NPSrare scores as a scatter plot with and without seed genes

    The function generates four plots in a single figure: two histograms of NPScommon and NPS rare with and without seed genes,
    a combined histogram of NPScommon-rare, and a scatter plot of combined scores with network threshold lines,
    where NPSrare is plotted on the x-axis and NPScommon is plotted on the y-axis.

    Parameters:
    - tblr (DataFrame): DataFrame containing NPSrare.
    - tblc (DataFrame): DataFrame containing NPScommon.
    - tblr_label (str): Label for rare trait, used in the Venn diagram.
    - tblc_label (str): Label for common trait, used in the Venn diagram.
    - tblr_seed (list of str): A list of rare seed genes.
    - tblc_seed (list of str): A list of common seed genes.
    - tblr_lim (float, optional): The NPSrare cutoff. Defaults to 1.5.
    - tblc_lim (float, optional): The NPScommon cutoff. Defaults to 1.5.
    - comb_lim (float, optional): The NPScommon-rare. Defaults to 3.
    - savefig (bool, optional): If True, saves the plot as an SVG file in a predefined directory. Defaults to False.

    Returns:
    None
    """
    fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(12.5, 5))

    # Combine zscore tables
    tbl_z = pd.concat([tblr, tblc], axis=1)
    tbl_z.columns = ('z1', 'z2')
    tbl_z['z_comb'] = tbl_z['z1'] * tbl_z['z2']

    # Separate in-network and out-network data
    inNetwork = tbl_z[(tbl_z['z1'] > tblr_lim) & (tbl_z['z2'] > tblc_lim) & (tbl_z['z_comb'] > comb_lim)]
    outNetwork = tbl_z[(tbl_z['z1'] <= tblr_lim) | (tbl_z['z2'] <= tblc_lim) | (tbl_z['z_comb'] <= comb_lim)]
    seedInNetwork = inNetwork[(~inNetwork.index.isin(tblr_seed)) & (~inNetwork.index.isin(tblc_seed))]

    # Plot with seed genes
    ax1.scatter(x=outNetwork['z1'], y=outNetwork['z2'], s=1, color=colour_dict['other'],label='non-network genes')
    ax1.scatter(x=inNetwork['z1'], y=inNetwork['z2'], s=1, color=colour_nonseed,label='network seed genes')
    ax1.scatter(x=seedInNetwork['z1'], y=seedInNetwork['z2'], s=1, color=colour_shared,label='network genes')
    # Set labels and lines for ax1
    ax1.set_xlabel(tblr_label)
    ax1.set_ylabel(tblc_label)
    ax1.axvline(x=tblr_lim, color=colour_r, linestyle='dashed', linewidth=1,label='NPSr')
    ax1.axhline(y=tblc_lim, color=colour_c, linestyle='dashed', linewidth=1,label='NPSh')
    x_points = [(i + 0.0001) / 10 for i in range(-50, 250)]
    combo_line = [comb_lim / x for x in x_points if x > comb_lim / 50]
    ax1.plot([x for x in x_points if x > comb_lim / 40], combo_line, color=colour_shared, linestyle='dashed', linewidth=1,label='NPShr')
    ax1.axvline(x=0, color='black', linestyle='solid', linewidth=1)
    ax1.axhline(y=0, color='black', linestyle='solid', linewidth=1)
    ax1.legend(loc='center right', bbox_to_anchor=(2.75, 0.5))

    
    # Filter out seed genes
    tbl_z = tbl_z[(~tbl_z.index.isin(tblr_seed)) & (~tbl_z.index.isin(tblc_seed))]
    inNetwork = tbl_z[(tbl_z['z1'] > tblr_lim) & (tbl_z['z2'] > tblc_lim) & (tbl_z['z_comb'] > comb_lim)]
    outNetwork = tbl_z[(tbl_z['z1'] <= tblr_lim) | (tbl_z['z2'] <= tblc_lim) | (tbl_z['z_comb'] <= comb_lim)]

    # Plot without seed genes
    ax2.scatter(x=outNetwork['z1'], y=outNetwork['z2'], s=1, color=colour_dict['other'])
    ax2.scatter(x=inNetwork['z1'], y=inNetwork['z2'], s=1, color=colour_shared)

    # Set labels and lines for ax2
    ax2.set_xlabel(tblr_label)
    ax2.set_ylabel(tblc_label)
    ax2.axvline(x=tblr_lim, color=colour_r, linestyle='dashed', linewidth=1)
    ax2.axhline(y=tblc_lim, color=colour_c, linestyle='dashed', linewidth=1)
    x_points = [(i + 0.0001) / 10 for i in range(-50, 250)]
    combo_line = [comb_lim / x for x in x_points if x > comb_lim / 50]
    ax2.plot([x for x in x_points if x > comb_lim / 40], combo_line, color=colour_shared, linestyle='dashed', linewidth=1)
    ax2.axvline(x=0, color='black', linestyle='solid', linewidth=1)
    ax2.axhline(y=0, color='black', linestyle='solid', linewidth=1)
    ax2.set_xlim([min(tbl_z['z1'])-.5, max(tbl_z['z1'])+.5])
    ax2.set_ylim([min(tbl_z['z2'])-.5, max(tbl_z['z2'])+.5])

    # Save the figure if requested
    if savefig:
        plt.savefig('figures/' + filename + '.svg', bbox_inches='tight')

    # Display the plot
    plt.show()
def venn_seeds(tblr_seed, tblc_seed, tblr_label, tblc_label, colour_r, colour_c,all_nodes,interactome_name,savefig=False):
    """
    Generates and displays a Venn diagram visualizing the overlap between seed genes from two lists within a given set of nodes.

    This function filters seed genes from two lists to ensure they are within a specified set of nodes, calculates the overlap between these filtered lists, and visualizes this overlap in a Venn diagram. The significance of the overlap is calculated using a hypergeometric test, similarly to methodologies used in specific scientific literature. Optionally, the diagram can be saved as an SVG file.

    Parameters:
    - tblr_seed (list of str): A list of rare seed genes.
    - tblc_seed (list of str): A list of common seed genes.
    - tblr_label (str): Label for rare trait, used in the Venn diagram.
    - tblc_label (str): Label for common trait, used in the Venn diagram.
    - colour_r (str): colour used for rare trait
    - colour_c (str): colour used for common trait
    - all_nodes (set): A set of all possible nodes within which the seeds should be filtered. Typically this is all nodes in PCNet. 
    - savefig (bool, optional): If True, saves the plot as an SVG file in a predefined directory. Defaults to False.

    Returns:
    None
    """
    tblr_seed=list(set(tblr_seed).intersection(all_nodes))
    tblc_seed=list(set(tblc_seed).intersection(all_nodes))  
    #define overlap for seed genes plot
    seed_overlap=set(tblr_seed).intersection(set(tblc_seed))
    print(seed_overlap)
    #compute significance of seed genes overlap- same test as used in BMI paper
    hyper = hypergeom(M=len(all_nodes), n=len(tblr_seed), N=len(tblc_seed))
    p_intersect_seed = hyper.sf(len(seed_overlap))
    
    venn2((len(tblr_seed)-len(seed_overlap), len(tblc_seed)-len(seed_overlap), len(seed_overlap)), 
          set_labels=(tblr_label, tblc_label), 
          set_colors=(colour_r, colour_c), alpha = 0.7)
    plt.title(' Seed Gene Overlap, p='+str(p_intersect_seed))
    if (savefig):
        plt.savefig('figures/seed_venn_'+tblr_label+'_'+tblc_label+'_'+interactome_name+'.svg',bbox_inches='tight')
    plt.show()
def plt_histogram (tblr, tblc, tblr_label, tblc_label, tblr_seed, tblc_seed,colour_r=colour_dict['loco'],colour_c=colour_dict['ext'],colour_shared=colour_dict['shared'], tblr_lim=1.5, tblc_lim=1.5, comb_lim=3, savefig=False,filename=('histogram_NPS')):
    """
    Visualizes the NPScommon and NPSrare scores as histograms and as a scatter plot

    The function generates four plots in a single figure: two histograms of NPScommon and NPS rare with and without seed genes, a combined histogram of NPScommon-rare, and a scatter plot of combined scores with network threshold lines, where NPSrare is plotted on the x-axis and NPScommon is plotted on the y-axis.

    Parameters:
    - tblr (DataFrame): DataFrame containing NPSrare.
    - tblc (DataFrame): DataFrame containing NPScommon.
    - tblr_label (str): Label for rare trait, used in the Venn diagram.
    - tblc_label (str): Label for common trait, used in the Venn diagram.
    - tblr_seed (list of str): A list of rare seed genes.
    - tblc_seed (list of str): A list of common seed genes.
    - tblr_lim (float, optional): The NPSrare cutoff. Defaults to 1.5.
    - tblc_lim (float, optional): The NPScommon cutoff. Defaults to 1.5.
    - comb_lim (float, optional): The NPScommon-rare. Defaults to 3.
    - savefig (bool, optional): If True, saves the plot as an SVG file in a predefined directory. Defaults to False.

    Returns:
    None
    """
    fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(18.75, 5))
    _, bins, _ = ax1.hist(tblr, bins=100, alpha=0.7, density=True, label=tblr_label, color=colour_r)
    _ = ax1.hist(tblc, bins=bins, alpha=0.7, density=True, label=tblc_label, color=colour_c)
    ax1.set_ylabel("density")
    ax1.set_xlabel("proximity zscore")
    ax1.legend()

    _, bins, _ = ax2.hist(tblr[~tblr.index.isin(tblr_seed)], bins=100, alpha=0.7, density=True, label=tblr_label,color=colour_r )
    _ = ax2.hist(tblc[~tblc.index.isin(tblc_seed)], bins=bins, alpha=0.7, density=True, label=tblc_label, color=colour_c)
    ax2.set_ylabel("density")
    ax2.set_xlabel("proximity zscore (no seed genes)")
    ax2.legend()
    
    _, bins, _ = ax3.hist(tblc['z']*tblr['z'], bins=bins, alpha=0.7, density=True, label='combined score', color=colour_shared)
    ax3.set_xlabel("Combined proximity zscore")
    ax3.set_ylabel("density")
    
    if (savefig):
        plt.savefig('figures/'+filename+'.svg',bbox_inches='tight')
    plt.show()
def combine_nps_table(tblr, tblc):
    tbl_z=pd.concat([tblr, tblc], axis=1)
    tbl_z.columns=('z1','z2')
    tbl_z['z_comb']=tbl_z['z1']*tbl_z['z2']
    tbl_z.columns=['NPSr','NPSc','NPScr']
    return(tbl_z)
def venn_net(tblr, tblc, tblr_label, tblc_label, p_net_overlap,colour_r=colour_dict['loco'],colour_c=colour_dict['ext'],colour_shared=colour_dict['shared'],tblr_lim=1.5, tblc_lim=1.5, comb_lim=3, savefig=False):
    """
    Generates and displays a Venn diagram visualizing the overlap between network implicated genes based calculated NPSrare and NPScommon.

    This function combines two tables of z-scores, filters entries based on network limits, and then visualizes the overlap between them in a Venn diagram. Additional metadata and network parameters are displayed in the diagram's title. Optionally, the diagram can be saved as an SVG file.

    Parameters:
    - tblr (DataFrame): A pandas DataFrame containing NPSrare.
    - tblc (DataFrame): A pandas DataFrame containing NPScommon.
    - tblr_label (str): Label for the rare trait, used in the Venn diagram.
    - tblc_label (str): Label for the common trait, used in the Venn diagram.
    - p_net_overlap (float): A probability or statistical value associated with the network overlap, shown in the plot title.
    - tblr_lim (float, optional): The NPSrare cutoff. Defaults to 1.5.
    - tblc_lim (float, optional): The NPScommon cutoff. Defaults to 1.5.
    - comb_lim (float, optional): The NPScommon-rare. Defaults to 3.
    - savefig (bool, optional): If True, saves the plot as an SVG file in a predefined directory. Defaults to False.

    Returns:
    None
    """
    print(tblr_lim)
    #combine zscore tables
    tbl_z=combine_nps_table(tblr, tblc)
    #subset table to those within network limit parameters
    inNetwork=tbl_z[(tbl_z['NPSr']>tblr_lim) & (tbl_z['NPSc']>tblc_lim) & (tbl_z['NPScr']>comb_lim)]
    print(len(inNetwork))
    #plot venn diagram
    Nr=(len(tbl_z[tbl_z['NPSr']>tblr_lim])-len(inNetwork))
    Nc=(len(tbl_z[tbl_z['NPSc']>tblc_lim])-len(inNetwork))
    Nboth=len(inNetwork)
    venn2((Nr,Nc,Nboth), 
		  set_labels=(tblr_label, tblc_label),
      set_colors=(colour_r, colour_c), alpha = 0.7)
    plt.title('p='+str(p_net_overlap)+ ', single cut='+str(tblr_lim)+', comb cut='+str(comb_lim))
    if savefig:
        plt.savefig('figures/network_venn_'+tblr_label+'_'+tblc_label+'.svg',bbox_inches='tight')
    plt.show()