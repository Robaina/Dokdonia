from Bio import SeqIO
import pandas as pd
import numpy as np
import random
import os
import re
from subprocess import call
from sklearn.metrics import silhouette_samples
from diffexpr.py_deseq import py_DESeq2
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
from rpy2.robjects import Formula
from matplotlib import pyplot as plt
import logging

from staticinteract import StaticInteract, DropDownWidget, RangeWidget
rpy2_logger.setLevel(logging.ERROR)


class GenomeGBK:

    def __init__(self, path_to_gbk):
        self._gbk = list(SeqIO.parse(path_to_gbk, 'genbank'))[0]

    @property
    def meta(self):
        return dict(self._gbk.features[0].qualifiers)

    @property
    def features(self):
        return [f for f in self._gbk.features[1:]]

    def getGeneInfo(self, gene_id: str):
        try:
            gene, cds = [f for f in self._gbk.features[1:]
                         if gene_id.lower() in f.qualifiers['locus_tag'][0].lower()]
            res = dict(cds.qualifiers)
            res.update({'location': gene.location})
            return res
        except Exception:
            raise ValueError(f'Gene {gene_id} not found in GBK')

    def has_EC_number(self, gene_id: str):
        return 'EC_number' in self.getGeneInfo(gene_id).keys()

    def getEnzymeGene(self, ec_number: str):
        try:
            return [f.qualifiers['locus_tag'][0] for f in self._gbk.features[1:]
                    if ('EC_number' in f.qualifiers.keys()
                        and ec_number in f.qualifiers['EC_number'][0])]
        except Exception:
            raise ValueError(f'Enzyme {ec_number} not found in GBK')


def getMetaMatrix(counts):
    return pd.DataFrame({
        'lighting': [s[0] for s in counts.columns[1:]],
        'temperature': [int(s[2:4]) for s in counts.columns[1:]],
        'replicate': [s[5:] for s in counts.columns[1:]]
    }, index=counts.columns[1:])


def runDEtest(counts, test='Wald', alpha=1e-2,
              formula='~ lighting', reduced_formula=None,
              log2fold_cutoff=0):
    '''
    Runs DeSeq2
    reduced_formula only for LRT test
    log2fold_cutoff: threshold to consider genes as DE when
    pair-wise comparisons with Wald test
    (Schurch et al., 2016 recommends 0.5 for 3 replicates)
    '''
    meta = getMetaMatrix(counts)
    dds = py_DESeq2(count_matrix=counts,
                    design_matrix=meta,
                    design_formula=formula,
                    gene_column='index')

    if test == 'LRT':
        dds.run_deseq(test=test, reduced=Formula(reduced_formula))
        dds.get_deseq_result(alpha=alpha)
    else:
        dds.run_deseq(test=test)
        dds.get_deseq_result(alpha=alpha)  #, lfcThreshold=log2fold_cutoff)
    res = dds.deseq_result
    res = res[res.padj < alpha]
    stats = {
        'DE+': res.log2FoldChange.where(
            res.log2FoldChange >= log2fold_cutoff).dropna().shape[0],
        'DE-': res.log2FoldChange.where(
            res.log2FoldChange <= -log2fold_cutoff).dropna().shape[0]
    }
    return (res, stats)


def addProteinNamesToDeseqResult(gbk, deseq_res):
    proteins = [gbk.getGeneInfo(gene_id)['product'][0] for gene_id in deseq_res.index]
    deseq_res['proteins'] = proteins
    return deseq_res


def writeClustInputFiles(clust_data, path_to_wd='Data'):
    clust_data.to_csv(os.path.join(path_to_wd, 'clust_input.tsv'), sep='\t')
    conds = np.unique([s[:4] for s in clust_data.columns])
    open(os.path.join(path_to_wd, 'clust_replicates.txt'), 'w').close()
    with open(os.path.join(path_to_wd, 'clust_replicates.txt'), 'a+') as file:
        for cond in conds:
            reps = ",".join(list(clust_data.filter(regex=f'{cond}').columns))
            txt_s = f'clust_input.tsv, {cond}, {reps}\n'
            file.write(txt_s)
          
    open(os.path.join(path_to_wd, 'clust_no_normalization.txt'), 'w').close()
    with open(os.path.join(path_to_wd, 'clust_no_normalization.txt'), 'a+') as file:
        file.write('clust_input.tsv 0') # no normalization
        
    open(os.path.join(path_to_wd, 'clust_transcript_cell_normalization.txt'), 'w').close()
    with open(os.path.join(path_to_wd, 'clust_transcript_cell_normalization.txt'), 'a+') as file:
        file.write('clust_input.tsv 101 4') # Quantile + z-score
        
    open(os.path.join(path_to_wd, 'clust_TPM_normalization.txt'), 'w').close()
    with open(os.path.join(path_to_wd, 'clust_TPM_normalization.txt'), 'a+') as file:
        file.write('clust_input.tsv 101 3 4') # Quantile + log2 + z-score


def runClust(path_to_wd, out_dir, cluster_tightness=1,
             replicates_file=None, normalization_file=None):
    """
    Compute clusters with clust
    clust_data: pandas DataFrame.
    """
    call_list = ['clust', os.path.join(path_to_wd, 'clust_input.tsv'),
        '-t', f'{cluster_tightness}',
        '-o', f'{out_dir}']
    if replicates_file is not None:
        call_list.append('-r')
        call_list.append(os.path.join(path_to_wd, replicates_file))
    if normalization_file is None:
        call_list.append('-n')
        call_list.append(os.path.join(path_to_wd, 'clust_no_normalization.txt'))
    else:
        call_list.append('-n')
        call_list.append(os.path.join(path_to_wd, normalization_file))
        
    call(call_list, cwd=path_to_wd)


def getGeneClusters(clust_data, path_to_wd, out_dir, cluster_tightness=1,
                    replicates_file=None, normalization_file=None):
    "Returns dict with Clust gene clusters"
    writeClustInputFiles(clust_data, path_to_wd)
    runClust(path_to_wd=path_to_wd,
             out_dir=out_dir, cluster_tightness=cluster_tightness, 
             replicates_file=replicates_file,
             normalization_file=normalization_file)

    clusters = pd.read_csv(
        os.path.join(out_dir, 'Clusters_Objects.tsv'), sep='\t', header=1)
    return {f'C{i}': clusters.iloc[:, i].dropna().values.tolist()
            for i in range(clusters.shape[1])}


def plotClusters(pdata, clusters):
    n_rows = int(np.ceil(len(clusters) / 2))
    fig, axes = plt.subplots(nrows=n_rows, ncols=2)
    plt.subplots_adjust(hspace=0.3)
    coords = list(np.ndindex((n_rows, 2)))
    for n, cluster_id in enumerate(clusters):
        i, j = coords[n]
        try:
            axis = axes[i, j]
        except Exception:
            axis = axes[j]
        cluster = clusters[cluster_id]
        pdata[pdata.index.isin(cluster)].transpose().plot(
            legend=False, figsize=(15, 18), title=f'{cluster_id}, size={len(cluster)}',
            ax=axis, color='#9a9a9a', linewidth=0.8,
            marker='.', markerfacecolor='#ee9929', markersize=12)
    plt.show()


def getAverageStandardRatio(IS_counts, standards_data):
    """
    Compute average standard ratios from standard counts and metadata
    """
    cond_out = 'D_25_R1'
    IS = IS_counts['index'].values
    conditions = IS_counts.columns.values.tolist()
    conditions.remove('index')
    conditions.remove(cond_out)

    avg_st_ratios = {}
    for cond_id in conditions:
        st_ratios = []
        for st_id in IS:

            st_copies = standards_data[
                (standards_data['Sample ID'] == cond_id) & (standards_data['Standard'] == st_id)
            ]['Standard added (copias)'].values[0]

            st_counts = IS_counts[IS_counts['index'] == st_id][cond_id].values[0]

            st_ratios.append(st_copies / st_counts)

        avg_st_ratios[cond_id] = {'average': np.mean(st_ratios),
                                  'std': np.std(st_ratios),
                                  'cv': np.std(st_ratios) / np.mean(st_ratios)}
    return avg_st_ratios


def getTranscriptsPerCell(counts, avg_st_ratios, abundance_meta):
    """
    Normalize counts by internal standards and cell abundances
    transcripts/cell = (counts * avg_st_ratio) / total_cell_abundance

    """
    cond_out = 'D_25_R1'
    conditions = abundance_meta['Sample'].values.tolist()
    conditions.remove(cond_out)

    n_counts = counts[counts.columns.intersection(conditions + ['index'])].copy()
    for cond_id in conditions:
        n_cells = abundance_meta[abundance_meta['Sample'] == cond_id]['Total_cell_abundance'].values[0]
        avg_ratio = avg_st_ratios[cond_id]['average']
        n_counts[cond_id] = (n_counts[cond_id] * avg_ratio) / n_cells

    return n_counts


def getECnumber(rxn_str):
    m = re.search('\[EC:(.+?)\]', rxn_str)
    if m:
        return m.group(1)
    else:
        return False


def assignSystemsToEnzymes(kegg_pathways):
    """
    Supersystem 4 to organismal systems
    Superystem 5 to human diseases
    Supersystem 6 corresponds to BRITE pathway database
    Supersystem 7 corresponds to pathwways not in BRITE or KEGG
    """
    kegg_dict = {}
    for supersystem in kegg_pathways[:6] + kegg_pathways[7:8]:
        for system in supersystem['children']:
            for subsystem in system['children']:
                for rxn in subsystem['children']:
                    ec = getECnumber(rxn['name'])
                    if ec and ec not in kegg_dict.keys():
                        kegg_dict[ec] = {
                            'supersystem': supersystem['name'],
                            'system': system['name'],
                            'subsystem': subsystem['name']
                        }
    return kegg_dict


def getCounts(array_like, sort_by_value=True):
    """
    Get dict of array item counts. Dict sorted
    by keys unless sort_by_value set to True.
    """
    unique_cond, counts_cond = np.unique(array_like, return_counts=True)
    counts = dict(zip(unique_cond, counts_cond))
    if sort_by_value:
        return dict(sorted(
            counts.items(), key=lambda item: item[1], reverse=True))
    else:
        return counts


def getRankedSystems(kegg_dict, EC_numbers, system_type='system'):
    systems = []
    for ec in EC_numbers:
        try:
            systems.append(kegg_dict[ec][system_type])
        except:
            pass
    return getCounts(systems)


# KEGG pathways analysis

def extractKoPathwayName(Ko_str):
    return re.sub('\[.*?\]', '', Ko_str).strip()


def extractKoID(Ko_str):
    return 'ko' + re.search('ko(.*?)\]', Ko_str).group(1)


def getKEGGPathwayDict(kegg_pathways):
    """
    Obtain dictionary containing KEGG pathway classification for each
    KEGG pathway id (koXXXXX)
    """
    kegg_dict = {}

    for supersystem in kegg_pathways[:4]:
        for system in supersystem['children']:
            for subsystem in system['children']:
                try:
                    ko_id = extractKoID(subsystem['name'])
                    ko_path_name = subsystem['name']
                except Exception:
                    pass
                kegg_dict[ko_id] = {
                    'subsystem': ko_path_name,
                    'system': system['name'],
                    'supersystem': supersystem['name']
                }
    return kegg_dict


def getGeneKOs(eggNOGresult):
    """
    Obtain dictionary containing KEGG pathway identifiers (koXXXX)
    for each gene as retrieved by EggNOG mapper
    """
    gene_Kos = {}
    for i, row in eggNOGresult.iterrows():
        gene_id = row['query']
        Kos = row['KEGG_Pathway']
        if type(Kos) == str:
            gene_Kos[gene_id] = [Ko for Ko in Kos.split(',') if 'ko' in Ko]
        else:
            gene_Kos[gene_id] = []
    return gene_Kos


def computeKEGGpathwaySize(gene_list, gene_ko_dict, ko_pathway_dict):
    """
    Count Dokdonia genes assigned to each KEGG pathway
    """
    pathway_counts = {'system': {}, 'subsystem': {}}
    for gene_id in np.intersect1d(gene_list, list(gene_ko_dict.keys())):
        kos = gene_ko_dict[gene_id]
        for ko in kos:
            try:
                system = ko_pathway_dict[ko]['system']
                subsystem = ko_pathway_dict[ko]['subsystem']

                if subsystem not in pathway_counts['subsystem'].keys():
                    pathway_counts['subsystem'][subsystem] = 1
                else:
                    pathway_counts['subsystem'][subsystem] += 1

                if system not in pathway_counts['system'].keys():
                    pathway_counts['system'][system] = 1
                else:
                    pathway_counts['system'][system] += 1
            except Exception:
                pass     
            
    return pathway_counts


def getGenesInKEGGsystem(ko_pathway_dict, gene_ko_dict, system:str, system_type:str) -> list:
    """
    Retrieve all genes assigned to selected system 
    system_type: system or subsystem.
    """
    kos = [k for k,v in ko_pathway_dict.items() if system.lower() in v[system_type].lower()]
    genes = [g for g,v in gene_ko_dict.items() if any(x in kos for x in v)]
    return genes


def summarizeKEGGpathwaysForGene(ko_pathway_dict, gene_ko_dict, gene_id, unique=False):
    """
    Summary of KEGG pathway classification for gene (may have multiple repeated ones)
    """
    res = {'subsystem': [], 'system': []}
    try:
        kos = gene_ko_dict[gene_id]
        if len(kos) > 0:
            for ko in kos:
                if ko in ko_pathway_dict.keys():
                    ko_path = ko_pathway_dict[ko]
                    res['subsystem'].append(ko_path['subsystem'])
                    res['system'].append(ko_path['system'])
         
        if unique:
            for k, v in res.items():
                res[k] = np.unique(v).tolist()
        return res

    except Exception:
        return res

    
def getKEGGpathwaysForGeneList(ko_pathway_dict, gene_ko_dict, gene_list, unique=False):
    """
    Get KEGG pathway classification for genes in gene_list
    if unique is set to True, it returns only the list of (unique)
    pathways found in the gene list
    """
    res = {'subsystem': [], 'system': []}
    for gene_id in gene_list:
        gene_res = summarizeKEGGpathwaysForGene(ko_pathway_dict, gene_ko_dict, gene_id)
        for k, v in gene_res.items():
            res[k].extend(v)
    if unique:
        for k, v in res.items():
            res[k] = np.unique(v).tolist()
    return res


def getKEGGPathwaysForLocusTag(gene_id, gene_ko_dict, ko_pathway_dict):
    res = {}
    try:
        kos = gene_ko_dict[gene_id]
        return {
            'system': [ko_pathway_dict[ko]['system'] for ko in kos if ko in ko_pathway_dict.keys()],
            'subsystem': [ko_pathway_dict[ko]['subsystem'] for ko in kos if ko in ko_pathway_dict.keys()]
        }
    except Exception:
        return {'system':'', 'subsystem': ''}


def getElementsFrequency(array_like, ranked=True):
    """
    Return dictionary with keys equal to unique elements in
    array_like and values equal to their frequencies.
    If ranked then dictionary is sorted
    """
    elems, counts = np.unique(array_like, return_counts=True)
    res = dict(zip(elems, counts / sum(counts)))
    if ranked:
        return dict(sorted(res.items(), key=lambda item: item[1], reverse=True))
    else:
        return res
    
    
def computeKEGGPathwayRepresentation(gene_list_pathway_counts, pathway_counts):
    """
    Compute proportion of genes assigned to pathway that
    are present in gene set.
    """
    res = {
        'system': {k: 0 for k in pathway_counts['system'].keys()},
        'subsystem': {k: 0 for k in pathway_counts['subsystem'].keys()}
    }
    
    for system_type, pathways in gene_list_pathway_counts.items():
        for pathway, counts in pathways.items():
            res[system_type][pathway] = counts / pathway_counts[system_type][pathway]
        
    return res
    

def extractSubsystemsFromSystem(subsystems_list, system_name, ko_pathway_dict):
    """
    Filter subsystems from list that belong to specified KEGG system name
    """
    return [s for s in subsystems_list if ko_pathway_dict[extractKoID(s)]['system'] == system_name]


def getEggNOGInputFile(gbk_file):
    "First line cannot be blank"
    with open('eggNOG_Input.fasta', 'a') as file:
        for rec in SeqIO.parse(gbk_file, "genbank"):
            for feature in rec.features:
                if feature.type == 'CDS':
                    gene = feature.qualifiers["locus_tag"][0].replace("'", "")
                    aas = feature.qualifiers["translation"][0].replace("'", "")
                    file.write(f'\n>{gene}\n{aas}')
                    

                                      
                    
                    
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
    gene_list = np.intersect1d(gene_list, list(gene_pathways.keys()))
    for gene_id in gene_list:
        pathway = gene_pathways[gene_id]
        if pathway['system'] != '':
            systems.append(pathway['system'])
            
        if pathway['subsystem'] != '':
            subsystems.append(pathway['subsystem'])
    
    pathway_counts['system'] = getCounts(systems)
    pathway_counts['subsystem'] = getCounts(subsystems)
    
    return pathway_counts


def getPathwayRepresentationInGeneList(gene_list, total_pathway_counts, gene_pathways):
    """
    Need to add pathways with representation of 0
    """
    pathway_rep = {'system': {}, 'subsystem': {}}
    pathway_counts = getPathwayCountsInGeneList(gene_list, gene_pathways)
    
    for system in total_pathway_counts['system'].keys():
        if system not in pathway_counts['system'].keys():
            pathway_rep['system'][system] = 0
        else:
            system_counts = pathway_counts['system'][system]
            system_total_counts = total_pathway_counts['system'][system]
            pathway_rep['system'][system] = system_counts / system_total_counts
            
    for subsystem in total_pathway_counts['subsystem'].keys():
        if subsystem not in pathway_counts['subsystem'].keys():
            pathway_rep['subsystem'][subsystem] = 0
        else:
            system_counts = pathway_counts['subsystem'][subsystem]
            system_total_counts = total_pathway_counts['subsystem'][subsystem]
            pathway_rep['subsystem'][subsystem] = system_counts / system_total_counts
    
    return pathway_rep


def permutePatricGenesInClusters(total_pathway_counts, gene_list, clusters,
                                 gene_pathways, n_permutations=10):
    """
    Obtain frequencies of PATRIC pathways in permuted clusters
        
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
        partition = randomPartition(gene_list, bin_sizes)

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


def runClusterPathwayEnrichmentAnalysisPatric(gene_list, clusters, total_pathway_counts, 
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
        
    permuted_path_rep = permutePatricGenesInClusters(total_pathway_counts, gene_list, clusters,
                                                     gene_pathways,n_permutations=n_permutations)
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


# Custom enrichment analysis based on permutation (randomization)
def randomPartition(elems:list, bin_sizes:list):
    """
    Randomly partition list elemnts into
    bin of given sizes
    """
    def shuffleList(elems:list):
        random.shuffle(elems)
        return elems

    elems = shuffleList(elems)
    partition = []
    start, end = 0, 0
    for bin_size in bin_sizes:
        end += bin_size
        partition.append(elems[start:end])
        start += bin_size
    return partition

# Modify to accommodate frequencies within pathways in each cluster
def permuteGenesInClusters(KEGG_pathway_counts, ko_pathway_dict, gene_ko_dict,
                           gene_list, clusters, n_permutations=10):
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
            'system': {p: [] for p in KEGG_pathway_counts['system']},
            'subsystem': {p: [] for p in KEGG_pathway_counts['subsystem']}
        } for k in cluster_ids
    }
    # Run permutation
    for i in range(n_permutations):
        partition = randomPartition(gene_list, bin_sizes)

        for cluster_id, rand_bin in zip(cluster_ids, partition):

            data = getKEGGpathwaysForGeneList(
                ko_pathway_dict, gene_ko_dict, rand_bin)
            
            data_counts = {k: getCounts(v, sort_by_value=False) for k,v in data.items()}
            
            pathway_representation = computeKEGGPathwayRepresentation(
                data_counts, KEGG_pathway_counts)

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


def runClusterPathwayEnrichmentAnalysis(gene_list, clusters, KEGG_pathway_counts, ko_pathway_dict,
                                        gene_ko_dict, n_permutations=10, sort_by_pvalue=True):
    """
    Run permutation analysis
    """
    def getKEGGrepresentationInClusters(clusters):
        res = {k: {} for k in clusters.keys()}
        for k, v in clusters.items():
            data = getKEGGpathwaysForGeneList(
                ko_pathway_dict, gene_ko_dict, v)
            data_counts = {k: getCounts(v, sort_by_value=False) for k,v in data.items()}
            pathway_representation = computeKEGGPathwayRepresentation(
                data_counts, KEGG_pathway_counts)      
            res[k]['system'] = pathway_representation['system']
            res[k]['subsystem'] = pathway_representation['subsystem']
        return res
    
    def computeSamplePvalue(sample, value):
        return len(np.where(np.array(sample) >= value)[0]) / len(sample)
    
    def computePathwayPvalue(pathways_freq, pathways_permuted_freq):
        p_pathways = {}
        for pathway, freq in pathways_freq.items():
            pvalue = computeSamplePvalue(pathways_permuted_freq[pathway], freq)
            p_pathways[pathway] = (freq, pvalue)
        return p_pathways
        
    
    cluster_path_freq = getKEGGrepresentationInClusters(clusters)
    permuted_path_freq = permuteGenesInClusters(KEGG_pathway_counts, ko_pathway_dict, gene_ko_dict,
                                                gene_list, clusters, n_permutations)
    p_KEGG_paths = {}
    for cluster_id in clusters.keys():
    
        systems = computePathwayPvalue(cluster_path_freq[cluster_id]['system'],
                                                     permuted_path_freq[cluster_id]['system'])
        subsystems = computePathwayPvalue(cluster_path_freq[cluster_id]['subsystem'],
                                                     permuted_path_freq[cluster_id]['subsystem'])
    
        if sort_by_pvalue:
            sorted_keys = np.array(list(systems.keys()))[np.argsort([pvalue for f, pvalue in systems.values()])]
            systems = {k: systems[k] for k in sorted_keys}

            sorted_keys = np.array(list(subsystems.keys()))[np.argsort([pvalue for f, pvalue in subsystems.values()])]
            subsystems = {k: subsystems[k] for k in sorted_keys}
        
        p_KEGG_paths[cluster_id] = {
                'system': systems,
                'subsystem': subsystems
            }
        
    return p_KEGG_paths


# Plotting functions
def plotClusterData(pdata, cluster, ax=None, cluster_id=None):
    pdata[pdata.index.isin(cluster)].transpose().plot(
        legend=False, title=f'{cluster_id}, size={len(cluster)}',
        ax=ax, color='#9a9a9a', linewidth=0.8,
        marker='.', markerfacecolor='#ee9929', markersize=12)


def plotKEGGFrequencies(data: dict, color=None, axis=None):
    """
    Bar plot of sorted KEGG systems or subsystems
    data: dictionary wth keys being system or subsystem names and
    values, their frequencies/representation.
    """
    if color is None:
        color = 'C0'
    pvalues = [v[1] for v in data.values()]
    clean_name_data = 100 * pd.Series(
        {extractKoPathwayName(k): data[k][0] for k in data.keys()}
    )
    ax = clean_name_data.plot.bar(figsize=(12, 8), color=color, ax=axis)
    for i, p in enumerate(ax.patches):
        ax.annotate(f'({pvalues[i]:.4f})', (p.get_x() * 1.005, p.get_height() * 1.006))
    

def plotSystemsAndSubsystemsWebPage(clusters, cluster_data, p_Data_paths,
                                    plot_first_N=10, color=None, 
                                    img_folder_name=None):
    """
    p_Data_paths is a dict with keys equal to database name and values equal
    to results.
    NOTE: Increase tick labels font size!
    """

    plt.rcParams.update({'figure.max_open_warning': 0})
    if color is None:
        color = 'C0'
    if img_folder_name is None:
        img_folder_name = 'iplot'

    system_types = ['system', 'subsystem']
    data_normalizations = ['TPM', 'TC']
    databases = ['KEGG', 'PATRIC']
    cluster_ids = list(np.unique(
    [c for n in data_normalizations for c in clusters[n].keys()])
                      )

    def plot_fun(database, data_normalization, system_type, cluster_id):
        
        fig, ax = plt.subplot_mosaic(
            """
            A
            B
            """,
            gridspec_kw = {"height_ratios": [0.7, 1]}
        )

        ax['B'].set_ylabel('Pathway representation (%)')
        ax['B'].set_title(f'{database} {system_type}s (sample p-value)')
        
        if cluster_id in clusters[data_normalization].keys():
            p_Data = p_Data_paths[data_normalization][database]
            kdata = {k: v 
                     for k,v in p_Data[cluster_id][system_type].items()}
            if len(kdata) > plot_first_N:
                kdata = {k: kdata[k] for k in list(kdata.keys())[:10]}
            plotClusterData(cluster_data[data_normalization], clusters[data_normalization],
                            ax['A'], cluster_id)
            plotKEGGFrequencies(kdata,
                                color=color, axis=ax['B'])
        fig.set_figwidth(20)
        fig.set_figheight(20)
        return fig

    i_fig = StaticInteract(plot_fun,
                           database=DropDownWidget(databases,
                                        description='Database'),
                           data_normalization=DropDownWidget(data_normalizations,
                                        description='Normalization'),
                           system_type=DropDownWidget(system_types,
                                        description='Level'),
                           cluster_id=DropDownWidget(cluster_ids,
                                        description='Cluster ID'),
                           interact_name=img_folder_name)
    return i_fig


def getRepresentationStackedPlotData(p_paths, pvalue_cutoff=1):
    """
    Make pandas dataframe of pathway representation across
    clusters
    """
    unique_systems = np.unique([
        system for c in p_paths.values()
        for system in c['system'].keys() 
    ])
    unique_subsystems = np.unique([
        subsystem
        for c in p_paths.values()
        for subsystem in c['subsystem'].keys()
    ])
    
    rep_systems, rep_subsystems = {}, {}
    for system in unique_systems:
        rep_systems[system] = []
        for cluster_id, values in p_paths.items():
            rep, pvalue = values['system'][system]
            if pvalue < pvalue_cutoff:
                rep_systems[system].append(rep)
            else:
                rep_systems[system].append(0.0)
        # Add rep outside clusters
        rep_systems[system].append(1 - sum(rep_systems[system]))
            
    for subsystem in unique_subsystems:
        rep_subsystems[subsystem] = []
        for cluster_id, values in p_paths.items():
            rep, pvalue = values['subsystem'][subsystem]
            if pvalue < pvalue_cutoff:
                rep_subsystems[subsystem].append(rep)
            else:
                rep_subsystems[subsystem].append(0.0)
        # Add rep outside clusters
        rep_subsystems[subsystem].append(1 - sum(rep_subsystems[subsystem]))
        
        rep_subsystems = {extractKoPathwayName(k): v
                          for k,v in rep_subsystems.items()}
            
    return {'system': rep_systems, 'subsystem': rep_subsystems}


def plotSystemsAndSubsystemsStacked(p_Data_paths, cluster_colors, img_folder_name):
    """
    """
    plt.rcParams.update({'figure.max_open_warning': 0})
    if img_folder_name is None:
        img_folder_name = 'iplot'
        
    system_types = ['system', 'subsystem']
    data_normalizations = ['TPM', 'TC']
    databases = ['KEGG', 'PATRIC']
        
    def plot_fun(database, data_normalization, system_type, tpvalue_cutoff):
        fig = plt.figure(figsize=(10, 8))
        ax = plt.gca()
        p_Data = p_Data_paths[data_normalization][database]
        rep_data = getRepresentationStackedPlotData(p_Data, tpvalue_cutoff)

        df = pd.DataFrame(rep_data[system_type], index=list(p_Data.keys()) + ['No cluster'])
        df = df.loc[:, df.loc['No cluster'] < 1.0].sort_values('No cluster', axis=1)
        ax = (100 * df).transpose().plot.bar(stacked=True, legend=True, figsize=(10,8),
                                             color=list(cluster_colors[data_normalization].values()),
                                             ax=ax)
        return fig
    
    i_fig = StaticInteract(plot_fun,
                           database=DropDownWidget(databases,
                                        description='Database'),
                           data_normalization=DropDownWidget(data_normalizations,
                                        description='Normalization'),
                           system_type=DropDownWidget(system_types,
                                        description='Level'),
                           tpvalue_cutoff=RangeWidget(min=0.05, max=1, 
                                                     step=0.05,
                                                     default=0.1,         
                                                     description='p-value cutoff'),
                           interact_name=img_folder_name)
    return i_fig


#####################################################
# Cluster analysis
#####################################################

def writeExcelOfClusterGenes(clusters, out_path, gbk,
                             patric_features, patric_pathways_genes, patric_pathways,
                             gene_ko_dict, ko_pathway_dict):
    
    silhouette = type(clusters[list(clusters.keys())[0]]) == dict
    writer = pd.ExcelWriter(out_path, engine='xlsxwriter')
    
    for cluster_id, cluster in clusters.items():
        gene_pathways = {}
        for gene_id in cluster:
            PATRIC_gene_pathways = getPatricPathwaysForLocusTag(gene_id, patric_features,
                                                                   patric_pathways_genes, patric_pathways)
            KEGG_gene_pathways = getKEGGPathwaysForLocusTag(gene_id, gene_ko_dict, ko_pathway_dict)

            if silhouette:
                gene_pathways[gene_id] = {
                    'Gene silhouette': cluster[gene_id],
                    'Gene product': gbk.getGeneInfo(gene_id)['product'][0],
                    'KEGG system': ', '.join(np.unique(np.array(KEGG_gene_pathways['system']))),
                    'KEGG subsystem': ', '.join(np.unique(np.array(KEGG_gene_pathways['subsystem']))),
                    'PATRIC system': ', '.join(np.unique(np.array(PATRIC_gene_pathways['system']))),
                    'PATRIC subsystem': ', '.join(np.unique(np.array(PATRIC_gene_pathways['subsystem'])))
                }
            else:
                gene_pathways[gene_id] = {
                    'Gene product': gbk.getGeneInfo(gene_id)['product'][0],
                    'KEGG system': ', '.join(np.unique(np.array(KEGG_gene_pathways['system']))),
                    'KEGG subsystem': ', '.join(np.unique(np.array(KEGG_gene_pathways['subsystem']))),
                    'PATRIC system': ', '.join(np.unique(np.array(PATRIC_gene_pathways['system']))),
                    'PATRIC subsystem': ', '.join(np.unique(np.array(PATRIC_gene_pathways['subsystem'])))
                }

        pd.DataFrame(gene_pathways).transpose().to_excel(writer, sheet_name=f'Cluster {cluster_id}')
    writer.save()
    
    
def computeGeneSilhouettes(clusters, data):
    """
    Compute gene silhouettes for each gene in clusters
    """
    gene_sil = {}
    genes_in_clusters, cluster_labels = [], []
    for cluster_id, cluster in clusters.items():
        genes_in_clusters.extend(cluster)
        cluster_labels.extend([cluster_id for _ in range(len(cluster))])

    X = data.loc[genes_in_clusters,:].values
    sil_values = silhouette_samples(X, cluster_labels)
    
    return {gene_id: sil_values[i] for i, gene_id in enumerate(genes_in_clusters)}


def rankGenesWithinClusters(clusters, data):
    """
    Rank genes within each cluster based on their silhouette
    """
    
    gene_sil = computeGeneSilhouettes(clusters, data)
    
    ranked_clusters = {}
    for cluster_id, cluster in clusters.items():
        sil_dict = {gene_id: gene_sil[gene_id] for gene_id in cluster}
        ranked_dict = dict(
            sorted(sil_dict.items(), key=lambda item: item[1], reverse=True)
        )
        ranked_clusters[cluster_id] = ranked_dict
        
    return ranked_clusters




# These functions were meant to do hypothesis testing, abandoning this idea.
# import statsmodels
# def correctPvalues(pvalues, FWER=0.05, method='fdr_bh'): 
#     """
#     Control FWER in reported p-values
#     """
#     reject, qvalues, alphaSidak, alphaBon = statsmodels.stats.multitest.multipletests(
#         pvalues, alpha=FWER, method=method, is_sorted=False, returnsorted=False)
#     return qvalues


# def extractSignificantPathways(p_KEGG_pathways, pvalue_cutoff):
#     """
#     Extract significantly enriched or depleted pathways
#     """
#     sig_pathways = {k: {'system': {'enriched': {}, 'depleted': {}},
#                        'subsystem': {'enriched': {}, 'depleted': {}}
#                        } for k in p_KEGG_pathways.keys()}
    
#     for cluster_id, pathways in p_KEGG_pathways.items():
#         for pathway, (freq, pvalue) in pathways['system'].items():
#             if pvalue < pvalue_cutoff:
#                 sig_pathways[cluster_id]['system']['enriched'][pathway] = (freq, pvalue)
#             elif pvalue > (1 - pvalue_cutoff): 
#                 sig_pathways[cluster_id]['system']['depleted'][pathway] = (freq, pvalue)
                
#         for pathway, (freq, pvalue) in pathways['subsystem'].items():
#             if pvalue < pvalue_cutoff:
#                 sig_pathways[cluster_id]['subsystem']['enriched'][pathway] = (freq, pvalue)
#             elif pvalue > (1 - pvalue_cutoff): 
#                 sig_pathways[cluster_id]['subsystem']['depleted'][pathway] = (freq, pvalue)
    
#     return sig_pathways