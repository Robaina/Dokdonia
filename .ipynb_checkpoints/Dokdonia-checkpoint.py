#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import os
import json
import copy
import pickle
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import Dokdonia_code as Dc


# Parsing GBK file and KEGG pathways
with open('Data/Function_Annotations/KEGG/kegg_pathways.json') as json_file:
    kegg_pathways = json.load(json_file)['children']
kegg_dict = Dc.assignSystemsToEnzymes(kegg_pathways)
    
gbk = Dc.GenomeGBK('Data/DokdoniaMED134.gbk')


# Loading counts and removing genes with low read counts across samples
min_count = 10

counts = pd.read_csv('Data/DokdoniaCounts.csv', index_col=0)
counts = counts[counts.filter(regex='^[^T]+$').columns]
conditions = [name.split('.sam')[0] for name in counts.columns]
counts.columns = conditions
counts = counts[(counts > min_count).all(1)]
counts.reset_index(level=0, inplace=True)

clust_tightness = 5 # Maximizes number of clusters in both cases

########################################################################
# Cluster ALL genes across temperatures using TPM values
########################################################################

# # Get TPM values
# clust_data_TPM = pd.read_csv('Data/tpm_counts.csv', index_col=0)
# clust_data_TPM = clust_data_TPM[clust_data_TPM.filter(regex='^[^T]+$').columns]
# clust_data_TPM.index.name = 'ID'
# clust_data_TPM.columns = conditions

# # Cluster genes
# res_id = 'CLUSTER_ALL_GENES_TPM'
# print(f'Running condition: {res_id}')
# workdir = os.path.join(os.getcwd(),'Results')
# clusters = Dc.getGeneClusters(clust_data_TPM, path_to_wd=workdir, 
#                               out_dir=os.path.join(workdir, res_id),
#                               cluster_tightness=clust_tightness,
#                               normalization_file='clust_TPM_normalization.txt')


# # Obtain genes without assigned cluster
# gene_ids = list(clust_data_TPM.index)
# genes_in_cluster = [gene for gene_ids in clusters.values() for gene in gene_ids]
# genes_without_cluster = np.setdiff1d(gene_ids, genes_in_cluster).tolist()
# clusters['No_cluster_assigned'] = genes_without_cluster

# # Save clusters
# Dc.saveToPickleFile(clusters, path_to_file=f'Results/Clusters_{res_id}.pkl')

# print(f'There a total of {len(genes_in_cluster)} genes assigned to a cluster')
# print(f'There a total of {len(genes_without_cluster)} genes not assigned to any cluster')

########################################################################
# Cluster ALL genes across temperatures using Transcript/cell values
########################################################################

# # Loading Transcript/Cell data and scale up data
# n_counts = pd.read_csv('Data/Dokdonia_transcripts_cell.csv', index_col=0)
# n_counts.iloc[:,1:] = 1e4 * n_counts.iloc[:,1:]

# clust_data_TC = n_counts
# clust_data_TC = clust_data_TC[clust_data_TC.filter(regex='^[^T]+$').columns]
# clust_data_TC = clust_data_TC.set_index('index')
# clust_data_TC.index.name = 'ID'

# # Cluster genes
# res_id = 'CLUSTER_ALL_GENES_TRANSCRIPT_CELL'
# print(f'Running condition: {res_id}')
# workdir = os.path.join(os.getcwd(),'Results')

# clusters = Dc.getGeneClusters(clust_data_TC, path_to_wd=workdir, 
#                               out_dir=os.path.join(workdir, res_id),
#                               cluster_tightness=clust_tightness,
#                               normalization_file='clust_TPM_normalization.txt')


# # Obtain genes without assigned cluster
# gene_ids = list(clust_data_TC.index)
# genes_in_cluster = [gene for gene_ids in clusters.values() for gene in gene_ids]
# genes_without_cluster = np.setdiff1d(gene_ids, genes_in_cluster).tolist()
# clusters['No_cluster_assigned'] = genes_without_cluster

# # Save clusters
# Dc.saveToPickleFile(clusters, path_to_file=f'Results/Clusters_{res_id}.pkl')

# print(f'There a total of {len(genes_in_cluster)} genes assigned to a cluster')
# print(f'There a total of {len(genes_without_cluster)} genes not assigned to any cluster')


##############################################################################
# Permutation tests
##############################################################################

N = 500000
res_ids = ['CLUSTER_ALL_GENES_TPM', 'CLUSTER_ALL_GENES_TRANSCRIPT_CELL']

##############################################################################
# Permutation test with KEGG annotations
##############################################################################

# # Read eggNOG - Mapper (http://eggnog-mapper.embl.de/) results 
# eggNOG = pd.read_excel('Data/Function_Annotations/KEGG/result_eggNOGMapper.xlsx', header=2)
# ko_pathway_dict = Dc.getKEGGPathwayDict(kegg_pathways)
# gene_ko_dict= Dc.getGeneKOs(eggNOG)
# gene_list = list(gene_ko_dict.keys())
# KEGG_pathway_counts = Dc.computeKEGGpathwaySize(gene_list, 
#                                                 gene_ko_dict, ko_pathway_dict)
        
# for res_id in res_ids:
#     clusters = Dc.readFromPickleFile(f'Results/Clusters_{res_id}.pkl')
#     clusters = {k: v for k,v in clusters.items() if k != 'No_cluster_assigned'}
#     p_KEGG_paths = Dc.runClusterPathwayEnrichmentAnalysis(gene_list, clusters, KEGG_pathway_counts,
#                                                           ko_pathway_dict, gene_ko_dict, n_permutations=N)

#     Dc.saveToPickleFile(p_KEGG_paths, path_to_file=f'Results/p_KEGG_paths_{N}_{res_id}.pkl')


##############################################################################
# Permutation test with PATRIC annotations
##############################################################################

patric_features = pd.read_csv('Data/Function_Annotations/PATRIC/Dokdonia_MED134_Craig_PATRIC_genome_feature.csv')
patric_pathways = pd.read_csv('Data/Function_Annotations/PATRIC/Dokdonia_MED134_PATRIC_pathways.csv')
patric_pathways_genes = pd.read_csv('Data/Function_Annotations/PATRIC/Dokdonia_MED134_Craig_PATRIC_pathways_genes.csv')

# Get dict of patric pathways for each locus tag (gene id)
gene_list = counts['index'].values
gene_pathways = {}
for gene_id in gene_list:
    gene_pathways[gene_id] = Dc.getPatricPathwaysForLocusTag(gene_id, patric_features,
                                                             patric_pathways_genes, patric_pathways)
total_pathway_counts = Dc.getPathwayCountsInGeneList(gene_list, gene_pathways)

for res_id in res_ids:
    clusters = Dc.readFromPickleFile(f'Results/Clusters_{res_id}.pkl')
    clusters = {k: v for k,v in clusters.items() if k != 'No_cluster_assigned'}
    p_PATRIC_paths = Dc.runClusterPathwayEnrichmentAnalysisPatric(gene_list, clusters, total_pathway_counts,
                                                                  gene_pathways,
                                                                  n_permutations=N,
                                                                  sort_by_pvalue=True)

    Dc.saveToPickleFile(p_PATRIC_paths, path_to_file=f'Results/p_PATRIC_paths_{N}_{res_id}.pkl')


# os.system("shutdown /s /t 1")
