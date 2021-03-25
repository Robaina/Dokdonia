from Bio import SeqIO
import pandas as pd
import numpy as np
import os
from subprocess import call
from diffexpr.py_deseq import py_DESeq2
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
from rpy2.robjects import Formula
import logging
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
    return {f'C{i}': clusters.iloc[:, i].dropna().values
            for i in range(clusters.shape[1])}


def plotClusters(pdata, clusters):
    n_rows = int(np.ceil(len(clusters) / 2))
    fig, axes = plt.subplots(nrows=n_rows, ncols=2)
    plt.subplots_adjust(hspace=0.3)
    coords = list(np.ndindex((n_rows, 2)))
    for n, cluster_id in enumerate(clusters):
        i, j = coords[n]
        cluster = clusters[cluster_id]
        ax = pdata[pdata.index.isin(cluster)].transpose().plot(
            legend=False, figsize=(15, 18), title=f'{cluster_id}, size={len(cluster)}',
            ax=axes[i,j])
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