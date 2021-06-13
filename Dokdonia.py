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
# eggNOG = pd.read_excel('Data/Function_Annotations/KEGG/result_eggNOGMapper.xlsx', header=2)
# ko_pathway_dict = Dc.getKEGGPathwayDict(kegg_pathways)
# gene_ko_dict= Dc.getGeneKOs(eggNOG)
# gene_list = T_res['all'].index
# gene_list = list(gene_ko_dict.keys())
# KEGG_pathway_counts = Dc.computeKEGGpathwaySize(gene_list, 
#                                                 gene_ko_dict, ko_pathway_dict)



# # Trying this out         
# gene_list = np.array(T_res['all'].index) 
# n_permutations = [100, 200, 1000, 10000, 20000, 40000, 60000, 80000, 100000, 200000, 500000]

# for N in n_permutations:
#     p_KEGG_paths = Dc.runClusterPathwayEnrichmentAnalysis(
#         gene_list, joined_clusters, KEGG_pathway_counts,
#         ko_pathway_dict, gene_ko_dict, n_permutations=N)
    
#     saveToPickleFile(p_KEGG_paths, path_to_file=f'p_KEGG_paths_{N}.pkl')


# PATRIC

patric_features = pd.read_csv('Data/Function_Annotations/PATRIC/Dokdonia_MED134_Craig_PATRIC_genome_feature.csv')
patric_pathways = pd.read_csv('Data/Function_Annotations/PATRIC/Dokdonia_MED134_PATRIC_pathways.csv')
patric_pathways_genes = pd.read_csv('Data/Function_Annotations/PATRIC/Dokdonia_MED134_Craig_PATRIC_pathways_genes.csv')

# Perform pathway analysis using PATRIC pathways
def locusTag2PatricID(locus_tag, patric_features):
    return patric_features['PATRIC ID'][patric_features['RefSeq Locus Tag'] == locus_tag].item()


def getPatricPathway(patric_id, patric_pathways_genes, patric_pathways):
    path_name = patric_pathways_genes['Pathway Name'][patric_pathways_genes['PATRIC ID'] == patric_id].item()
    path_class = patric_pathways['Pathway Class'][patric_pathways['Pathway Name'] == path_name].item()
    return {'subsystem': path_name, 'system': path_class}


def getPatricPathwaysForLocusTag(locus_tag, patric_features,
                                patric_pathways_genes, patric_pathways):
    try:
        patric_id = locusTag2PatricID(locus_tag, patric_features)
        pathway = getPatricPathway(patric_id, patric_pathways_genes, patric_pathways)
        return pathway
    except Exception:
        return {'system':'', 'subsystem': ''}
    
    
def getPathwayCountsInGeneList(gene_list, gene_pathways):
    
    pathway_counts = {'system': [], 'subsystem': []}
    systems, subsystems = [], []
    for gene_id in gene_list:
        pathway = gene_pathways[gene_id]
        if pathway['system'] != '':
            systems.append(pathway['system'])
            
        if pathway['subsystem'] != '':
            subsystems.append(pathway['subsystem'])
    
    pathway_counts['system'] = Dc.getCounts(systems)
    pathway_counts['subsystem'] = Dc.getCounts(subsystems)
    
    return pathway_counts
    
    
def getPathwayRepresentationInGeneList(gene_list, total_pathway_counts, gene_pathways):
    
    pathway_rep = {'system': {}, 'subsystem': {}}
    pathway_counts = getPathwayCountsInGeneList(gene_list, gene_pathways)
    
    pathway_rep['system'] = {k: v/total_pathway_counts['system'][k] for k,v in pathway_counts['system'].items()}
    pathway_rep['subsystem'] = {k: v/total_pathway_counts['subsystem'][k] for k,v in pathway_counts['subsystem'].items()}
    
    return pathway_rep


# Modify to accommodate frequencies within pathways in each cluster
def permuteGenesInClusters(total_pathway_counts, gene_list, clusters,
                           gene_pathways, n_permutations=10):
    """
    Obtain frequencies of KEGG pathways in permuted clusters
        
    Note: since we are randomly shuffling genes in bins, we can't ensure
    that each permutation will produce the same set of pathways. Thus,
    some pathways may get empty frequency value in a permutation. Filling
    with 0s to ammend this issue.
    """
    
    bin_sizes = [len(v) for v in clusters.values()]
    cluster_ids = [k for k in clusters.keys()]
    bin_sizes.append(len(gene_list) - sum(bin_sizes))
    
    # Initialize result dict
    res = {
        k: {
            'system': {p: [] for p in total_pathway_counts['system']},
            'subsystem': {p: [] for p in total_pathway_counts['subsystem']}
        } for k in cluster_ids
    }
    # Run permutation
    for i in range(n_permutations):
        partition = Dc.randomPartition(gene_list, bin_sizes)

        for cluster_id, rand_bin in zip(cluster_ids, partition):

            pathway_representation = getPathwayRepresentationInGeneList(
                rand_bin, total_pathway_counts, gene_pathways)

            for k, v in pathway_representation['system'].items():
                res[cluster_id]['system'][k].append(v)
            for k, v in pathway_representation['subsystem'].items():
                res[cluster_id]['subsystem'][k].append(v)
                
    # Add 0s if sample size is smaller than n_permutations 
    # (because pathway not represented in random samples)
    for cluster_id in cluster_ids:
        for sys_type in ('system', 'subsystem'):
            for k in res[cluster_id][sys_type].keys():
                sample_size = len(res[cluster_id][sys_type][k])
                size_diff = n_permutations - sample_size
                if size_diff > 0:
                    res[cluster_id][sys_type][k].extend([0.0 for _ in range(size_diff)])

    return res


def runClusterPathwayEnrichmentAnalysis(gene_list, clusters, total_pathway_counts, 
                                        gene_pathways, n_permutations=10, sort_by_pvalue=True):
    """
    Run permutation analysis
    """
    def computeSamplePvalue(sample, value):
        return len(np.where(np.array(sample) >= value)[0]) / len(sample)
    
    def computePathwayPvalue(pathways_freq, pathways_permuted_freq):
        p_pathways = {}
        for pathway, freq in pathways_freq.items():
            pvalue = computeSamplePvalue(pathways_permuted_freq[pathway], freq)
            p_pathways[pathway] = (freq, pvalue)
        return p_pathways
        
    cluster_path_rep = {}
    for cluster_id, cluster in clusters.items():
        cluster_path_rep[cluster_id] = getPathwayRepresentationInGeneList(cluster,
                                                                          total_pathway_counts,
                                                                          gene_pathways)
        
    permuted_path_rep = permuteGenesInClusters(total_pathway_counts, gene_list, clusters,
                                                gene_pathways,
                                                n_permutations=n_permutations)
    p_paths = {}
    for cluster_id in clusters.keys():
    
        systems = computePathwayPvalue(cluster_path_rep[cluster_id]['system'],
                                                     permuted_path_rep[cluster_id]['system'])
        subsystems = computePathwayPvalue(cluster_path_rep[cluster_id]['subsystem'],
                                                     permuted_path_rep[cluster_id]['subsystem'])
    
        if sort_by_pvalue:
            sorted_keys = np.array(list(systems.keys()))[np.argsort([pvalue for f, pvalue in systems.values()])]
            systems = {k: systems[k] for k in sorted_keys}

            sorted_keys = np.array(list(subsystems.keys()))[np.argsort([pvalue for f, pvalue in subsystems.values()])]
            subsystems = {k: subsystems[k] for k in sorted_keys}
        
        p_paths[cluster_id] = {
                'system': systems,
                'subsystem': subsystems
            }
        
    return p_paths


# Get dict of patric pathways for each locus tag (gene id)
gene_list = counts['index'].values
gene_pathways = {}
for gene_id in gene_list:
    gene_pathways[gene_id] = getPatricPathwaysForLocusTag(gene_id, patric_features,
                                                          patric_pathways_genes, patric_pathways)
    
    
n_permutations = [100, 200, 1000, 10000, 20000, 40000, 60000, 80000, 100000, 200000, 500000]
total_pathway_counts = getPathwayCountsInGeneList(gene_list, gene_pathways)

for N in n_permutations:
    print(f'Running {N} permutations')
    p_PATRIC_paths = runClusterPathwayEnrichmentAnalysis(gene_list, clusters, total_pathway_counts,
                                                         gene_pathways,
                                                         n_permutations=N,
                                                         sort_by_pvalue=True)
    Dc.saveToPickleFile(p_PATRIC_paths, path_to_file=f'p_PATRIC_paths_{N}.pkl')
    

os.system("shutdown /s /t 1")
