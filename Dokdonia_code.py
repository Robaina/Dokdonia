from Bio import SeqIO
import pandas as pd
import numpy as np
import os
import re
from subprocess import call
from diffexpr.py_deseq import py_DESeq2
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
from rpy2.robjects import Formula
from matplotlib import pyplot as plt
import logging
from ipywidgets import widgets, interact
from ipywidgets_New.interact import StaticInteract
from ipywidgets_New.widgets import DropDownWidget
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
        dds.get_deseq_result(alpha=alpha, lfcThreshold=log2fold_cutoff)
    res = dds.deseq_result
    res = res[res.padj < alpha]
    stats = {'DE+': res.log2FoldChange.where(res.log2FoldChange > 0).dropna().shape[0],
             'DE-': res.log2FoldChange.where(res.log2FoldChange < -0).dropna().shape[0]}
    return (res, stats)


def writeClustInputFiles(DE_TPM, path_to_wd='Data'):
    DE_TPM.to_csv(os.path.join(path_to_wd, 'clust_input.tsv'), sep='\t')
    conds = np.unique([s[:4] for s in DE_TPM.columns])
    open(os.path.join(path_to_wd, 'clust_replicates.txt'), 'w').close()
    with open(os.path.join(path_to_wd, 'clust_replicates.txt'), 'a+') as file:
        for cond in conds:
            reps = ",".join(list(DE_TPM.filter(regex=f'{cond}').columns))
            txt_s = f'clust_input.tsv, {cond}, {reps}\n'
            file.write(txt_s)


def runClust(DE_TPM, path_to_wd, out_dir, cluster_tightness=1):
    """
    Compute clusters with clust
    DE_TPM: pandas DataFrame.
    """
    call([
        'clust', os.path.join(path_to_wd, 'clust_input.tsv'),
        '-r', os.path.join(path_to_wd, 'clust_replicates.txt'),
        f'-t {cluster_tightness}',
        '-o', f'{out_dir}'
    ], cwd=path_to_wd)


def getGeneClusters(DE_TPM, path_to_wd, out_dir, cluster_tightness=1):
    "Returns dict with Clust gene clusters"
    writeClustInputFiles(DE_TPM, path_to_wd)
    runClust(DE_TPM, path_to_wd=path_to_wd,
             out_dir=out_dir, cluster_tightness=cluster_tightness)

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
        cluster = clusters[cluster_id]
        pdata[pdata.index.isin(cluster)].transpose().plot(
            legend=False, figsize=(15, 18), title=f'{cluster_id}, size={len(cluster)}',
            ax=axes[i, j], color='#9a9a9a', linewidth=0.8,
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
        gene_id = row['Query']
        Kos = row['KEGG pathway']
        if type(Kos) == str:
            gene_Kos[gene_id] = [Ko for Ko in Kos.split(',') if 'ko' in Ko]
        else:
            gene_Kos[gene_id] = []
    return gene_Kos

def summarizeKEGGpathwaysForGene(ko_pathway_dict, gene_ko_dict, gene_id):
    """
    Summary of KEGG pathway classification for gene (may have multiple repeated ones)
    """
    res = {'subsystem': [], 'system': [], 'supersystem': []}
    try:
        kos = gene_ko_dict[gene_id]
        if len(kos) > 0:
            for ko in kos:
                if ko in ko_pathway_dict.keys():
                    ko_path = ko_pathway_dict[ko]
                    res['subsystem'].append(ko_path['subsystem'])
                    res['system'].append(ko_path['system'])
                    res['supersystem'].append(ko_path['supersystem'])

        for k, v in res.items():
            res[k] = np.unique(v).tolist()
        return res

    except Exception:
        return res

def getKEGGpathwaysForGeneList(ko_pathway_dict, gene_ko_dict, gene_list):
    """
    Get KEGG pathway classification for genes in gene_list
    """
    res = {'subsystem': [], 'system': [], 'supersystem': []}
    for gene_id in gene_list:
        gene_res = summarizeKEGGpathwaysForGene(ko_pathway_dict, gene_ko_dict, gene_id)
        for k, v in gene_res.items():
            res[k].extend(v)
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

def extractSubsystemsFromSystem(subsystems_list, system_name, ko_pathway_dict):
    """
    Filter subsystems from list that belong to specified KEGG system name
    """
    return [s for s in subsystems_list if ko_pathway_dict[extractKoID(s)]['system'] == system_name]


def plotKEGGFrequencies(data, color=None, axis=None):
    """
    Bar plot of sorted KEGG systems or subsystems
    """
    if color is None:
        color = 'C0'
    clean_name_data = {extractKoPathwayName(k): data[k] for k in data.keys()}
    pd.Series(clean_name_data).plot.bar(figsize=(12, 8), color=color, ax=axis)


def plotSystemsAndSubsystems(data, ko_pathway_dict, color=None):
    if color is None:
        color = 'C0'

    system_freqs = getElementsFrequency(data['system'])

    def plot_fun(system_name):
        subsystems = extractSubsystemsFromSystem(data['subsystem'],
                                                      system_name,
                                                      ko_pathway_dict)
        subsystem_freqs = getElementsFrequency(subsystems)
        colors = ['grey' for _ in range(len(system_freqs))]
        colors[list(system_freqs.keys()).index(system_name)] = color
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
        ax1.set_ylabel('freq')
        ax2.set_ylabel('freq')
        ax1.set_title('KEGG systems')
        ax2.set_title(f'KEGG pathways of {system_name}')
        fig.set_figwidth(25, forward=True)
        fig.set_figheight(10, forward=True)
        plotKEGGFrequencies(system_freqs, color=colors, axis=ax1)
        plotKEGGFrequencies(subsystem_freqs, color=color, axis=ax2)
        plt.show()

    interact(plot_fun,
             system_name=widgets.Dropdown(options=list(system_freqs.keys()),
                                             description='KEGG pathway')
            )

def plotCluster(pdata, clusters, cluster_id, ax):
    cluster = clusters[cluster_id]
    pdata[pdata.index.isin(cluster)].transpose().plot(
        legend=False, figsize=(15, 18), title=f'{cluster_id}, size={len(cluster)}',
        ax=ax, color='#9a9a9a', linewidth=0.8,
        marker='.', markerfacecolor='#ee9929', markersize=12)

def plotSystemsAndSubsystemsWebPage(clusters, pdata, ko_pathway_dict, gene_ko_dict,
                             color=None, img_folder_name=None):

    plt.rcParams.update({'figure.max_open_warning': 0})
    if color is None:
        color = 'C0'
    if img_folder_name is None:
        img_folder_name = 'iplot'

    cluster_ids = list(clusters.keys())
    system_names = np.unique([v['system'] for v in ko_pathway_dict.values()]).tolist()

    def plot_fun(system_name, cluster_id):

        data = getKEGGpathwaysForGeneList(
            ko_pathway_dict, gene_ko_dict, clusters[cluster_id])

        system_freqs = getElementsFrequency(data['system'])

        colors = ['grey' for _ in range(len(system_freqs))]

        fig, ax = plt.subplot_mosaic(
            """
            AA
            BC
            """,
            gridspec_kw={"height_ratios": [0.7, 1]}
        )

        ax['B'].set_ylabel('freq')
        ax['C'].set_ylabel('freq')
        ax['B'].set_title('KEGG systems')
        ax['C'].set_title(f'KEGG pathways of {system_name}')

        try:
            subsystems = extractSubsystemsFromSystem(data['subsystem'],
                                                     system_name,
                                                     ko_pathway_dict)

            subsystem_freqs = getElementsFrequency(subsystems)
            plotKEGGFrequencies(subsystem_freqs, color=color, axis=ax['C'])
            colors[list(system_freqs.keys()).index(system_name)] = color
        except Exception:
            pass


        plotKEGGFrequencies(system_freqs, color=colors, axis=ax['B'])
        plotCluster(pdata, clusters, cluster_id, ax['A'])
        fig.set_figwidth(20)
        fig.set_figheight(20)
        return fig

    i_fig = StaticInteract(plot_fun,
                           system_name=DropDownWidget(system_names,
                                        description='KEGG pathway'),
                           cluster_id=DropDownWidget(cluster_ids,
                                        description='Cluster ID'),
                           interact_name=img_folder_name)
    return i_fig

def getEggNOGInputFile(gbk_file):
    "First line cannot be blank"
    with open('eggNOG_Input.fasta', 'a') as file:
        for rec in SeqIO.parse(gbk_file, "genbank"):
            for feature in rec.features:
                if feature.type == 'CDS':
                    gene = feature.qualifiers["locus_tag"][0].replace("'", "")
                    aas = feature.qualifiers["translation"][0].replace("'", "")
                    file.write(f'\n>{gene}\n{aas}')
