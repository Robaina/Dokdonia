# Python script to convert Internal Standard Fasta (with sequence identifier as "Chromosomes")
# to a GTF file to be input in countReads() which contains these same IDs and the loci which
# are just the total length of each sequence.

import os
from Bio import SeqIO
import Dokdonia_code as Dc

def getGeneStrand(gbk, gene_id):
    strand = 0
    for f in gbk.features:
        try:
            if gene_id.lower() in f.qualifiers['gene'][0].lower():
                strand = f.strand
        except Exception:
            pass
    return strand

path_to_wd = os.getcwd()
input_file = 'Internal_Standard_sequences.fasta'
output_file = 'Internal_Standard_sequences.gtf'

# Parse GBK
gbk = Dc.GenomeGBK('SulfolobussolfataricusP2.gbk')

# Parse FASTA
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
fasta_dic = {fasta.id: (1, len(fasta.seq)) for fasta in fasta_sequences}

# Write GTF
strand = ''
open(os.path.join(path_to_wd, output_file), 'w').close()
with open(os.path.join(path_to_wd, output_file), 'a+') as file:
    for gene, locus in fasta_dic.items():
        
        gene_id = gene.split('_')[1]
        s = getGeneStrand(gbk, gene_id)
        if s == 1:
            strand = '+'
        elif s == -1:
            strand = '-'
        else:
            strand = '+'
        txt_s = f'{gene}\tGenbank\tgene\t{locus[0]}\t{locus[1]}\t.\t{strand}\t0\tgene_id "{gene}" ; transcript_id "{gene}"\n'
        file.write(txt_s)
