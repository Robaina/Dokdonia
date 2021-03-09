# Count mapped fragments to genes in SAM file
import re
import os
import numpy as np

sam_dir = ''
sam_files = []


def findAll(string, pattern):
    return [m.start() for m in re.finditer(pattern, string)]


def countMappedFragments(sam_str, gene_list):
    return {gene: findAll(sam_str, gene) for gene in gene_list}


def getGeneList(gtf_path):
    gene_list = []
    gtf = open(gtf_path, "r")
    lines = gtf.readlines()
    gtf.close()

    for line in lines:
        idx = line.find('MED134_')
        gene_id = line[idx:idx+12]
        gene_list.append(gene_id)

    return np.unique(gene_list).tolist()


gene_list = getGeneList("../Data/DokdoniaMED134.gtf")
counts = {}

for file in sam_files:
    with open(os.path.join(sam_dir, file), "r") as text_file:
        sam_str = text_file.read()
    counts[file.split('.sam')[0]] = countMappedFragments(sam_str, gene_list)
