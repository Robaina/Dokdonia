#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import os
import json
import copy
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import Dokdonia_code as Dc


# In[2]:


import pickle

def saveToPickleFile(python_object, path_to_file='object.pkl'):
    """
    Save python object to pickle file
    """
    out_file = open(path_to_file,'wb')
    pickle.dump(python_object, out_file)
    out_file.close()
    
def readFromPickleFile(path_to_file='object.pkl'):
    """
    Load python object from pickle file.
    Returns python object.
    """
    in_file = open(path_to_file,'rb')
    python_object = pickle.load(in_file)
    return python_object


# In[3]:


# Parsing GBK file and KEGG pathways
with open('Data/Function_Annotations/KEGG/kegg_pathways.json') as json_file:
    kegg_pathways = json.load(json_file)['children']
kegg_dict = Dc.assignSystemsToEnzymes(kegg_pathways)
    
gbk = Dc.GenomeGBK('Data/DokdoniaMED134.gbk')


# In[4]:


# Loading counts and removing genes with low read counts across samples
min_count = 10

counts = pd.read_csv('Data/DokdoniaCounts.csv', index_col=0)
counts = counts[counts.filter(regex='^[^T]+$').columns]
conditions = [name.split('.sam')[0] for name in counts.columns]
counts.columns = conditions
counts = counts[(counts > min_count).all(1)]
counts.reset_index(level=0, inplace=True)


p_value_cutoff = 1e-2
fold_cutoff = k = 0.5

# DE light vs dark across temperatures
L_D_res, L_D_stats = {}, {}
L_D_res['all'], L_D_stats['all'] = Dc.runDEtest(counts, test='Wald', alpha=p_value_cutoff,
                                                formula='~ lighting', log2fold_cutoff=k)
# L vs D for each temperature
for T in ['10', '18', '25', '34']:
    counts_T = counts[counts.filter(regex=f'{T}|index').columns]    
    L_D_res[T], L_D_stats[T] = Dc.runDEtest(counts_T, test='Wald', alpha=p_value_cutoff,
                                            formula='~ lighting', log2fold_cutoff=k)


# Do not discriminate between Light and Dark conditions
T_res, T_stats = {}, {}
T_res['all'], T_stats['all'] = Dc.runDEtest(counts, test='LRT', alpha=p_value_cutoff,
                                            formula='~ temperature', reduced_formula='~ 1')

# Discriminate between Light and Dark conditions
for L in ['L', 'D']:
    counts_L = counts[counts.filter(regex=f'{L}|index').columns]    
    T_res[L], T_stats[L] = Dc.runDEtest(counts_L, test='LRT', alpha=p_value_cutoff,
                                        formula='~ temperature', reduced_formula='~ 1')
    


# Get TPM values
clust_data = pd.read_csv('Data/tpm_counts.csv', index_col=0)
clust_data = clust_data[clust_data.filter(regex='^[^T]+$').columns]
clust_data.index.name = 'ID'
clust_data.columns = conditions


# In[14]:


# Cluster only DE genes across temperatures (n_counts)
res_id = 'CLUSTER_DE_GENES'
gene_ids = T_res['all'].index
workdir = os.path.join(os.getcwd(),'Results')
DE_Data = clust_data[clust_data.index.isin(gene_ids)]
clusters = Dc.getGeneClusters(DE_Data, path_to_wd=workdir, 
                              out_dir=os.path.join(workdir, res_id),
                              cluster_tightness=5)

pdata = pd.read_csv(os.path.join(
    os.getcwd(),f'Results/{res_id}/Processed_Data/clust_input.tsv_processed.tsv'),
                    sep='\t', index_col='Genes')
# Add fake column to separate datasets
pdata.insert(4, '', [np.nan for n in range(pdata.shape[0])])

# Obtain DE genes without cluster
genes_in_cluster = [gene for gene_ids in clusters.values() for gene in gene_ids]
genes_without_cluster = np.setdiff1d(gene_ids, genes_in_cluster)




# In[15]:


print(f'There a total of {len(genes_in_cluster)} DE genes assigned to a cluster')
print(f'There a total of {len(genes_without_cluster)} DE genes not assigned to any cluster')


# ## DE genes not assigned to any cluster
# Let's run cluster again using the subset of genes not assigned to any cluster.

# In[16]:


# Cluster only DE genes across temperatures (n_counts)
res_id = 'CLUSTER_DE_GENES_NOT_ASSIGNED'
gene_ids = genes_without_cluster
workdir = os.path.join(os.getcwd(),'Results')
DE_Data = clust_data[clust_data.index.isin(gene_ids)]
na_clusters = Dc.getGeneClusters(DE_Data, path_to_wd=workdir, 
                              out_dir=os.path.join(workdir, res_id),
                              cluster_tightness=5)

pdata = pd.read_csv(os.path.join(
    os.getcwd(),f'Results/{res_id}/Processed_Data/clust_input.tsv_processed.tsv'),
                    sep='\t', index_col='Genes')
# Add fake column to separate datasets
pdata.insert(4, '', [np.nan for n in range(pdata.shape[0])])

# Obtain DE genes without cluster
genes_in_cluster = [gene for gene_ids in na_clusters.values() for gene in gene_ids]
genes_without_cluster = np.setdiff1d(gene_ids, genes_in_cluster)

print(f'There a total of {len(genes_in_cluster)} DE genes assigned to a cluster')
print(f'There a total of {len(genes_without_cluster)} DE genes not assigned to any cluster')



# These clusters show patterns that are similar to the some of the first. Specifically, we could join C0 here with C2 in the first cluster set, C1 with C1, C2 with C6 and C3 with C5.

# In[17]:


# Joined clusters
joined_clusters = copy.deepcopy(clusters)
joined_clusters['C2'].extend(na_clusters['C0'])
joined_clusters['C1'].extend(na_clusters['C1'])
joined_clusters['C6'].extend(na_clusters['C2'])
joined_clusters['C5'].extend(na_clusters['C3'])



# Read eggNOG - Mapper (http://eggnog-mapper.embl.de/) results 
eggNOG = pd.read_excel('Data/Function_Annotations/KEGG/result_eggNOGMapper.xlsx', header=2)
ko_pathway_dict = Dc.getKEGGPathwayDict(kegg_pathways)
gene_ko_dict= Dc.getGeneKOs(eggNOG)
gene_list = T_res['all'].index
gene_list = list(gene_ko_dict.keys())
KEGG_pathway_counts = Dc.computeKEGGpathwaySize(gene_list, 
                                                gene_ko_dict, ko_pathway_dict)



# Trying this out         
gene_list = np.array(T_res['all'].index) 
n_permutations = [100, 200, 1000, 10000, 20000, 40000, 60000, 80000, 100000, 200000]

for N in n_permutations:
    p_KEGG_paths = Dc.runClusterPathwayEnrichmentAnalysis(
        gene_list, joined_clusters, KEGG_pathway_counts,
        ko_pathway_dict, gene_ko_dict, n_permutations=N)
    
    saveToPickleFile(p_KEGG_paths, path_to_file=f'p_KEGG_paths_{N}.pkl')
    

os.system("shutdown /s /t 1")
