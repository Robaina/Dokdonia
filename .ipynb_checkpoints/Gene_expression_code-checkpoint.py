import pandas as pd
import numpy as np
import os
from subprocess import call
from diffexpr.py_deseq import py_DESeq2
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
from rpy2.robjects import Formula
from matplotlib import pyplot as plt
import logging
rpy2_logger.setLevel(logging.ERROR)


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


def writeClustInputFiles(DE_TPM, path_to_wd='Results'):
    DE_TPM.to_csv(os.path.join(path_to_wd, 'clust_input.tsv'), sep='\t')
    conds = np.unique([s[:4] for s in DE_TPM.filter(regex='^[^index]+$').columns])
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
        pdata[pdata.index.isin(cluster)].transpose().plot(
            legend=False, figsize=(15, 18), title=f'{cluster_id}, size={len(cluster)}',
            ax=axes[i, j])
    plt.show()

    
def tpmNormalize(counts, patric):
    gene_lengths = {}
    for gene_id in counts['index']:
        try:
            gene_lengths[gene_id] = patric.loc[patric['PATRIC ID'] == gene_id]['Length'].item()
        except Exception:
            gene_lengths[gene_id] = np.nan

    rpk = counts.iloc[:, 1:].divide(list(gene_lengths.values()), axis=0)
    tpm = rpk.divide(rpk.sum(axis=0).values / 1e6, axis=1)
    tpm.insert(0, 'index', counts['index'])
    return tpm
