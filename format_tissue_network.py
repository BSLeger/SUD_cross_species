import pandas as pd
import os
import networkx as nx
import gzip
import ndex2
import sys
import pickle
import gzip
import numpy as np
import scipy.sparse as sp
from collections import defaultdict

from netcoloc import netprop_zscore
from netcoloc import netprop
from netcoloc import network_colocalization



print('setting up environment')
os.chdir('/tscc/projects/ps-palmer/brittany/SUD_cross_species/')

tissue=sys.argv[1]
print(f'tissue is {tissue}')
file_path=f'tissue_networks/{tissue}.gz'
print(os.path.exists(file_path)) #check if input data exists

#read in network from edge file
print('reading in network')

def process_network_in_chunks(file_path, chunk_size=100000):
    node_set = set()
    edges = []
    node_index = {}
    
    print('First pass: determine unique nodes and create an index')
    with gzip.open(file_path, 'rt') as file:
        while True:
            lines = [file.readline() for _ in range(chunk_size)]
            if not lines or lines[0] == '':
                break  # Stop at EOF
            
            for line in lines:
                parts = line.strip().split()
                if len(parts) < 2:
                    continue  # Skip malformed lines
                
                node1, node2 = (parts[0]), (parts[1])
                
                node_set.add(node1)
                node_set.add(node2)
                edges.append((node1, node2))
                edges.append((node2, node1))  # Ensure symmetry
    
    print('Create node index mapping')
    node_list = sorted(node_set)  # Ensure consistent order- MUST REPLACE WHEN FIXED
    node_index = {node: idx for idx, node in enumerate(node_list)}
    
    print('Second pass: construct sparse adjacency matrix')
    row_idx, col_idx = [], []
    degree = defaultdict(int)
    
    with gzip.open(file_path, 'rt') as file:
        while True:
            lines = [file.readline() for _ in range(chunk_size)]
            if not lines or lines[0] == '':
                break  # Stop at EOF
            
            for line in lines:
                parts = line.strip().split()
                if len(parts) < 2:
                    continue  # Skip malformed lines
                
                node1, node2 = (parts[0]), (parts[1])
    
                i, j = node_index[node1], node_index[node2]
                
                # Store indices (all weights are 1)
                row_idx.append(i)
                col_idx.append(j)
                row_idx.append(j)
                col_idx.append(i)  # Ensure symmetric adjacency matrix
                
                # Update node degree
                degree[node1] += 1
                degree[node2] += 1

    n = len(node_list)
    adj_matrix = sp.csr_matrix((np.ones(len(row_idx)), (row_idx, col_idx)), shape=(n, n))
    return adj_matrix, node_index, degree, node_list

def normalize_adjacency_matrix(adj_matrix, node_index, degree, node_list, conserve_heat=True):
    print('normalize adjacency matrix')
    row_idx, col_idx = adj_matrix.nonzero()
    weight_vals = adj_matrix.data.copy()  # Copy original weights (all 1s)
    
    if conserve_heat:
        for k in range(len(weight_vals)):
            i, j = row_idx[k], col_idx[k]
            weight_vals[k] = 1/degree[node_list[i]]  # Normalize by the degree of the destination node
    
    else:
        for k in range(len(weight_vals)):
            i, j = row_idx[k], col_idx[k]
            weight_vals[k] = 1/np.sqrt(degree[node_list[i]] * degree[node_list[j]])
    
    x=sp.csr_matrix((weight_vals, (row_idx, col_idx)), shape=adj_matrix.shape)
    return x


print('importing file')
file_path=f'tissue_networks/{tissue}.gz'
outdir='tissue_networks/intermediate/'
print(f'file exists={os.path.exists(file_path)}')

adj_matrix, node_index, degree, node_list = process_network_in_chunks(file_path)
wp = normalize_adjacency_matrix(adj_matrix, node_index, degree, node_list, conserve_heat=True)

print('writing the node list, adjacency matrix, and degree dictionary to file')
#outdir='tissue_networks/intermediate/'

# Save adjacency matrix as a sparse file
sp.save_npz(f'{outdir}normalized_adjacency_{tissue}.npz', wp)

# Save the node list in the same order as adjacency matrix
with open(f'{outdir}node_list_{tissue}.txt', 'w') as f:
    for node in node_list:
        f.write(f"{node}\n")

pd.DataFrame.from_dict(degree,orient='index').to_csv(f'{outdir}degree_{tissue}.csv',header=None)

print('calculate w_double_prime and write it to file')
wdp = netprop.get_individual_heats_matrix(wp, .5)
np.save(f'{outdir}w_double_prime_{tissue}',wdp)
