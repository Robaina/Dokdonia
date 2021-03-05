#!/usr/bin/env python

# conda activate salmon
# conda deactivate

import os, subprocess
from Bio import SeqIO

ppath = "/home/robaina/Documents/Proyecto_rodopsina/Data/LauraDokdoniaReadsCleaned"
ppatern1="_1."
ppatern2="_2."

gbk_filename = "/home/robaina/Documents/Proyecto_rodopsina/Data/DokdoniaMED134.gbk"
faa_filename = "/home/robaina/Documents/Proyecto_rodopsina/Data/DokdoniaMED134aa.fasta"
fnt_filename = "/home/robaina/Documents/Proyecto_rodopsina/Data/DokdoniaMED134.fasta"

input_handle  = open(gbk_filename, "r")
output_handleaa = open(faa_filename, "w")
output_handlent = open(fnt_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank"):
	print ("Dealing with GenBank record %s" % seq_record.id)
	for seq_feature in seq_record.features:
		if seq_feature.type=="CDS":
			#print (seq_feature)
			try:
				output_handleaa.write(">%s %s\n%s\n" % (
						seq_feature.qualifiers['locus_tag'][0],
						seq_feature.qualifiers['product'][0],
						seq_feature.qualifiers['translation'][0]))

				output_handlent.write(">%s\n%s\n" % (
						seq_feature.qualifiers['locus_tag'][0],
						seq_feature.extract(seq_record.seq)))

			except:
				pass

output_handleaa.close()
output_handlent.close()
input_handle.close()


print("readseq -format=gff -o " + gbk_filename.replace(".gbk", "") + ".gff " + gbk_filename.replace(".gbk", "") + ".gbk")
os.system("readseq -format=gff -o " + gbk_filename.replace(".gbk", "") + ".gff " + gbk_filename.replace(".gbk", "") + ".gbk")

print ("python GffToGtf.py " + gbk_filename.replace(".gbk", ".gff") + " > " + gbk_filename.replace(".gbk", ".gtf"))
os.system("python GffToGtf.py " + gbk_filename.replace(".gbk", ".gff") + " > " + gbk_filename.replace(".gbk", ".gtf"))

cmd = ["rm -r /home/robaina/Documents/Proyecto_rodopsina/Data/ + fnt_filename[:fnt_filename.find('.')]"]
pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["rm -r /home/robaina/Documents/Proyecto_rodopsina/Data/quants"]
pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
p_status = pipe.wait(); out, err = pipe.communicate()

# Salgae.fasta has the orfs
print ("salmon index -t " + fnt_filename+" -i " + fnt_filename[:fnt_filename.find(".")])
os.system("salmon index -t " + fnt_filename+" -i " + fnt_filename[:fnt_filename.find(".")])

dirList = os.listdir(ppath)
ffiles = []
for fname in dirList:
	if fname.find(ppatern1) > -1:
		ffiles.append(fname)

ffiles.sort()
print(ffiles)
print(len(ffiles))

for i in range(len(ffiles)):
	print(ffiles[i])

	a = ffiles[i].replace("1.fastq.gz", "")

	# Paired-end:
	print("salmon quant -i ./" + fnt_filename[:fnt_filename.find(".")] + " -l A -1 " + ppath+ffiles[i] + " -2 " + ppath + ffiles[i].replace(ppatern1, ppatern2) + " -p 3 --gcBias --validateMappings -o ./quants/" + a)
	os.system("salmon quant -i ./" + fnt_filename[:fnt_filename.find(".")] + " -l A -1 " + ppath + ffiles[i] + " -2 " + ppath + ffiles[i].replace(ppatern1, ppatern2) + " -p 3 --gcBias --validateMappings -o ./quants/" + a)

    # Single-end
	# print ("salmon quant -i ./"+fnt_filename[:fnt_filename.find(".")]+" -l A -r "+ppath+ffiles[i]+" -p 3 --validateMappings -o ./quants/"+a)
	# os.system("salmon quant -i ./"+fnt_filename[:fnt_filename.find(".")]+" -l A -r "+ppath+ffiles[i]+" -p 3 --validateMappings -o ./quants/"+a)

	#exit()

print ("rm -r ./" + fnt_filename[:fnt_filename.find(".")])
cmd = ["rm -r ./" + fnt_filename[:fnt_filename.find(".")]]
pipe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
p_status = pipe.wait(); out, err=pipe.communicate()

exit()
