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
              formula='~ lighting', reduced_formula=None):
    '''
    Runs DeSeq2
    reduced_formula only for LRT test
    '''
    meta = getMetaMatrix(counts)
    dds = py_DESeq2(count_matrix=counts,
                    design_matrix=meta,
                    design_formula=formula,
                    gene_column='index')

    if test == 'LRT':
        dds.run_deseq(test=test, reduced=Formula(reduced_formula))
    else:
        dds.run_deseq(test=test)
    dds.get_deseq_result(alpha=alpha)
    res = dds.deseq_result
    res = res[res.padj < alpha]
    stats = {'DE+': res.log2FoldChange.where(res.log2FoldChange > 0).dropna().shape[0],
             'DE-': res.log2FoldChange.where(res.log2FoldChange < 0).dropna().shape[0]}
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
