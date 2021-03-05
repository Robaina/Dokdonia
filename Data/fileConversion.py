from Bio import SeqIO

file = 'DokdoniaMED134.gbk'
SeqIO.convert(file, 'genbank', file.split('.')[0] + '_full.fasta', 'fasta')
