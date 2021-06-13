from Bio import SeqIO
import pandas as pd
import numpy as np
import random
import pickle
import os
import re
from subprocess import call
from diffexpr.py_deseq import py_DESeq2
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
from rpy2.robjects import Formula
from matplotlib import pyplot as plt
import logging
from ipywidgets_New.interact import StaticInteract
from ipywidgets_New.widgets import DropDownWidget
rpy2_logger.setLevel(logging.ERROR)


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
        dds.get_deseq_result(alpha=alpha, lfcThreshold=log2fold_cutoff)
    res = dds.deseq_result
    res = res[res.padj < alpha]
    stats = {'DE+': res.log2FoldChange.where(res.log2FoldChange > 0).dropna().shape[0],
             'DE-': res.log2FoldChange.where(res.log2FoldChange < -0).dropna().shape[0]}
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
        file.write('clust_input.tsv 0')


def runClust(path_to_wd, out_dir, cluster_tightness=1, normalization=True):
    """
    Compute clusters with clust
    clust_data: pandas DataFrame.
    """
    call_list = ['clust', os.path.join(path_to_wd, 'clust_input.tsv'),
        '-r', os.path.join(path_to_wd, 'clust_replicates.txt'),
        '-t', f'{cluster_tightness}',
        '-o', f'{out_dir}']
    if not normalization:
        call_list.append('-n')
        call_list.append(os.path.join(path_to_wd, 'clust_no_normalization.txt'))
    call(call_list, cwd=path_to_wd)


def getGeneClusters(clust_data, path_to_wd, out_dir, cluster_tightness=1,
                    input_data_normalization=True):
    "Returns dict with Clust gene clusters"
    writeClustInputFiles(clust_data, path_to_wd)
    runClust(path_to_wd=path_to_wd,
             out_dir=out_dir, cluster_tightness=cluster_tightness, 
             normalization=input_data_normalization)

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
    def getKEGGfrequenciesInClusters(clusters):
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
        
    
    cluster_path_freq = getKEGGfrequenciesInClusters(clusters)
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
def plotCluster(pdata, clusters, cluster_id, ax):
    cluster = clusters[cluster_id]
    pdata[pdata.index.isin(cluster)].transpose().plot(
        legend=False, figsize=(15, 18), title=f'{cluster_id}, size={len(cluster)}',
        ax=ax, color='#9a9a9a', linewidth=0.8,
        marker='.', markerfacecolor='#ee9929', markersize=12)


def plotKEGGFrequencies(data, color=None, axis=None):
    """
    Bar plot of sorted KEGG systems or subsystems
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
    

def plotSystemsAndSubsystemsWebPage(clusters, pdata, p_Data_paths,
                                    plot_first_N=10, color=None, 
                                    img_folder_name=None):
    """
    p_Data_paths is a dict with keys equal to database name and values equal
    to results.
    """

    plt.rcParams.update({'figure.max_open_warning': 0})
    if color is None:
        color = 'C0'
    if img_folder_name is None:
        img_folder_name = 'iplot'

    cluster_ids = list(clusters.keys())
    system_names = np.unique(
        [k for v in p_KEGG_paths.values() for k in v['system'].keys()]
    ).tolist()
    
    system_types = ['system', 'subsystem']
    databases = list(p_Data_paths.keys())

    def plot_fun(database, system_type, cluster_id):

        fig, ax = plt.subplot_mosaic(
            """
            A
            B
            """,
            gridspec_kw = {"height_ratios": [0.7, 1]}
        )

        ax['B'].set_ylabel('Pathway representation (%)')
        ax['B'].set_title(f'{database} {system_type}s (sample p-value)')

        kdata = {k: v 
                 for k,v in p_Data_paths[database][cluster_id][system_type].items()}
        if len(kdata) > plot_first_N:
            kdata = {k: kdata[k] for k in list(kdata.keys())[:10]}
        plotCluster(pdata, clusters, cluster_id, ax['A'])
        plotKEGGFrequencies(kdata,
                            color=color, axis=ax['B'])
        fig.set_figwidth(20)
        fig.set_figheight(20)
        return fig

    i_fig = StaticInteract(plot_fun,
                           database=DropDownWidget(databases,
                                        description='Database'),
                           system_type=DropDownWidget(system_types,
                                        description='Level'),
                           cluster_id=DropDownWidget(cluster_ids,
                                        description='Cluster ID'),
                           interact_name=img_folder_name)
    return i_fig


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
