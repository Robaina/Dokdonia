# Python script to convert Internal Standard Fasta (with sequence identifier as "Chromosomes") 
# to a GTF file to be input in countReads() which contains these same IDs and the loci which 
# are just the total length of each sequence.

import os
from Bio import SeqIO

path_to_wd = os.getcwd()
input_file = 'Internal_Standard_sequences.fasta'
output_file = 'Internal_Standard_sequences.gtf'

# Parse FASTA
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
fasta_dic = {fasta.id: (1, len(fasta.seq)) for fasta in fasta_sequences}

# Write GTF
open(os.path.join(path_to_wd, output_file), 'w').close()
with open(os.path.join(path_to_wd, output_file), 'a+') as file:
    for gene, locus in fasta_dic.items():
        txt_s = f'{gene} Genbank gene {locus[0]} {locus[1]} . + 0 gene_id "{gene}" ; transcript_id "{gene}"\n'
        file.write(txt_s)