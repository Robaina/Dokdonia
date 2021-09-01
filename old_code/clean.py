#!/usr/bin/env python
# source activate salmon
# afterqc only works with python 2
# Notice that some pair reads don't match at the end. is.py fixes this.
#


ppath="/usr/gonzalez/transcriptomes/LauraDokdoniaReads/"
ppatern1="_1."
ppatern2="_2."
seconddb="dbDokdonia.fasta"

import os, subprocess


cmd = ["rm -r ./temp"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir ./temp"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

dirList=os.listdir(ppath)
ffiles=[]
ffiles2=[]
for fname in dirList:
	if fname.find(ppatern1)>-1:
		ffiles.append(fname)
		ffiles2.append("")

ffiles.sort()
print (ffiles)
print (len(ffiles))

for i in range(len(ffiles)):
	print (ffiles[i])
	
	print ("grep '"+ffiles[i]+"' correspondencia.txt")
	cmd = ["grep '"+ffiles[i]+"' correspondencia.txt"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.strip(), err.strip())

	m=out.strip().split()

	print("cp "+ppath+ffiles[i]+" ./temp/"+m[5]+"_1.fastq.gz")
	cmd = ["cp "+ppath+ffiles[i]+" ./temp/"+m[5]+"_1.fastq.gz"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.strip(), err.strip())
	
	ffiles2[i]=m[5]+"_1.fastq.gz"

	print("cp "+ppath+ffiles[i].replace(ppatern1,ppatern2)+" ./temp/"+m[5]+"_2.fastq.gz")
	cmd = ["cp "+ppath+ffiles[i].replace(ppatern1,ppatern2)+" ./temp/"+m[5]+"_2.fastq.gz"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.strip(), err.strip())

###########################################################################################################################
# Quality trimming:

cmd = ["rm -r ./temp2"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir ./temp2"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

for i in range(len(ffiles2)):
	print (ffiles2[i])
	print ("after.py -1 ./temp/"+ffiles2[i]+" -2 ./temp/"+ffiles2[i].replace(ppatern1,ppatern2))
	os.system("after.py -1 ./temp/"+ffiles2[i]+" -2 ./temp/"+ffiles2[i].replace(ppatern1,ppatern2))

	print("cp ./good/"+ffiles2[i].replace("_1.fastq.gz","_1.good.fq.gz")+" ./temp2/"+ffiles2[i].replace("_1.good.fq.gz","_1.fastq.gz"))
	cmd = ["cp ./good/"+ffiles2[i].replace("_1.fastq.gz","_1.good.fq.gz")+" ./temp2/"+ffiles2[i].replace("_1.good.fq.gz","_1.fastq.gz")]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))
	
	print("cp ./good/"+ffiles2[i].replace(ppatern1,ppatern2).replace("_2.fastq.gz","_2.good.fq.gz")+" ./temp2/"+ffiles2[i].replace(ppatern1,ppatern2).replace("_2.good.fq.gz","_2.fastq.gz"))
	cmd = ["cp ./good/"+ffiles2[i].replace(ppatern1,ppatern2).replace("_2.fastq.gz","_2.good.fq.gz")+" ./temp2/"+ffiles2[i].replace(ppatern1,ppatern2).replace("_2.good.fq.gz","_2.fastq.gz",)]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out.decode('ascii'), err.decode('ascii'))

	print ("zcat ./temp/"+ffiles2[i]+" | wc -l")
	cmd = ["zcat ./temp/"+ffiles2[i]+" | wc -l"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	
	nlines1=int(out.decode('ascii'))
	
	print ("zcat ./temp2/"+ffiles2[i].replace(ppatern1,ppatern2).replace("_2.good.fq.gz","_2.fastq.gz")+" | wc -l")
	cmd = ["zcat ./temp2/"+ffiles2[i].replace(ppatern1,ppatern2).replace("_2.good.fq.gz","_2.fastq.gz")+" | wc -l"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

	nlines2=int(out.decode('ascii'))
	
	print (nlines1)
	print (nlines2)
	print (100*float(nlines2)/float(nlines1))

	cmd = ["rm -r ./bad"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["rm -r ./QC"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["rm -r ./good"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

###########################################################################################################################


###########################################################################################################################
# To remove phiX174:

cmd = ["rm *.fna.*"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm *.sam"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm *.bam"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["bwa index NC_001422.fna"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["rm -r ./temp3"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir ./temp3"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

for i in range(len(ffiles2)):
	print (ffiles2[i])

	print ("bwa mem -M -t 20 NC_001422.fna ./temp2/"+ffiles2[i]+" ./temp2/"+ffiles2[i].replace(ppatern1,ppatern2)+" > file.sam")
	os.system("bwa mem -M -t 20 NC_001422.fna ./temp2/"+ffiles2[i]+" ./temp2/"+ffiles2[i].replace(ppatern1,ppatern2)+" > file.sam")
	
	print ("samtools view -S -b file.sam > file.bam")
	os.system("samtools view -S -b file.sam > file.bam")
	
	print ("samtools view -b -f 4 file.bam > tmp.bam")
	os.system("samtools view -b -f 4 file.bam > tmp.bam") # unmapped reads
	
	print ("samtools sort tmp.bam -o filesorted.bam")
	os.system("samtools sort tmp.bam -o filesorted.bam")

	print ("bedtools bamtofastq -i filesorted.bam -fq ./temp3/"+ffiles2[i].replace(".gz","")+" -fq2 ./temp3/"+ffiles2[i].replace(ppatern1,ppatern2).replace(".gz",""))
	os.system("bedtools bamtofastq -i filesorted.bam -fq ./temp3/"+ffiles2[i].replace(".gz","")+" -fq2 ./temp3/"+ffiles2[i].replace(ppatern1,ppatern2).replace(".gz",""))

	print ("zcat ./temp2/"+ffiles2[i]+" | wc -l")
	cmd = ["zcat ./temp2/"+ffiles2[i]+" | wc -l"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	
	nlines1=int(out.split()[0])
	print (nlines1)
	
	print ("wc -l ./temp3/"+ffiles2[i].replace(".gz",""))
	cmd = ["wc -l ./temp3/"+ffiles2[i].replace(".gz","")]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out, err)

	nlines2=int(out.split()[0])

	print (nlines2)
	print (100*float(nlines2)/float(nlines1))
	
	
	print ("wc -l ./temp3/"+ffiles2[i].replace(".gz",""))
	cmd = ["wc -l ./temp3/"+ffiles2[i].replace(".gz","")]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out, err)
	
	nseq1=int(out.decode('ascii').split()[0])
	print(); print ("No. sequences: "+str(nseq1/4))

	print ("wc -l ./temp3/"+ffiles2[i].replace(".gz","").replace(ppatern1,ppatern2))
	cmd = ["wc -l ./temp3/"+ffiles2[i].replace(".gz","").replace(ppatern1,ppatern2)]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out, err)
	
	nseq2=int(out.decode('ascii').split()[0])
	print(); print ("No. sequences: "+str(nseq2/4))
	
	if nseq1!=nseq2:
		print ("Different number of sequences")
		exit()

	nlines=0
	with open("./temp3/"+ffiles2[i].replace(".gz","")) as file1, open("./temp3/"+ffiles2[i].replace(".gz","").replace(ppatern1,ppatern2)) as file2:
		for line1, line2 in zip(file1, file2):
			if line1.startswith("@") and line2.startswith("@"):
				if not(line1.endswith("1\n") and line2.endswith("2\n")):
					print (line1)
					print (line2)
					print ("Problem!")
					exit()
			nlines+=1	
			   
	print (nlines/4)

	if nseq1!=nlines:
		print (nlines/4)
		print ("Something is wrong!")
		exit()
    
cmd = ["rm *.fna.*"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm *.sam"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["rm *.bam"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

cmd = ["rm NC_001422.fna.*"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
###########################################################################################################################


###########################################################################################################################
# Removes rRNA

cmd = ["rm ./db/db.idx.*"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

os.system("indexdb_rna --ref ./db/db.fasta,./db/db.idx --fast -v")

cmd = ["rm -r ./temp4"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir ./temp4"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

from Bio.SeqIO.QualityIO import FastqGeneralIterator

for i in range(len(ffiles2)):
	print (ffiles2[i])

	print ("bash ./merge-paired-reads.sh ./temp3/"+ffiles2[i].replace(".gz","")+" ./temp3/"+ffiles2[i].replace(".gz","").replace(ppatern1,ppatern2)+" merged.fastq")
	cmd = ["bash ./merge-paired-reads.sh ./temp3/"+ffiles2[i].replace(".gz","")+" ./temp3/"+ffiles2[i].replace(".gz","").replace(ppatern1,ppatern2)+" merged.fastq"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out, err)
    
	print ("cat merged.fastq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > file_sorted.fastq")
	os.system("cat merged.fastq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > file_sorted.fastq")
	os.system("mv file_sorted.fastq merged.fastq")

	print ("sortmerna --ref ./db/db.fasta,./db/db.idx --reads merged.fastq --log --num_alignments 1 --blast 0 -a 20 --fastx --aligned ./rrna --other ./norrna")
	os.system("sortmerna --ref ./db/db.fasta,./db/db.idx --reads merged.fastq --log --num_alignments 1 --blast 0 -a 20 --fastx --aligned ./rrna --other ./norrna")
	#################################################################################################
	# Removes unpaired reads
	mixed_file = "norrna.fastq"
	paired_file = "paired.fastq"
	
	out_handle = open(paired_file, "w")
	prev = None
	for curr in FastqGeneralIterator(open(mixed_file, "rU")):
	    if curr[0].split()[0].endswith("/1"):
	        prev = curr
	    elif not curr[0].split()[0].endswith("/2"):
	        raise ValueError("Expect IDs to end /1 and /2,\n%s" % curr[0])
	    elif prev and prev[0].split()[0] == curr[0].split()[0][:-2] + "/1":
	        out_handle.write("@%s\n%s\n+\n%s\n" % prev)
	        out_handle.write("@%s\n%s\n+\n%s\n" % curr)
	        prev = None
	out_handle.close()
	
	print ("cat norrna.fastq | wc -l")
	cmd = ["cat norrna.fastq | wc -l"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out, err)
	nseq1=float(out)/4
	
	print ("Sequences before removing unpairs: "+str(nseq1))

	print ("cat paired.fastq | wc -l")
	cmd = ["cat paired.fastq | wc -l"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out, err)
	nseq2=float(out)/4
	
	print ("Sequences after removing unpairs: "+str(nseq2))
	print ("Difference: "+str(nseq1-nseq2))
	print ("Percent difference: "+str(100*(nseq1-nseq2)/nseq1))

	#################################################################################################
	print ("bash ./unmerge-paired-reads.sh paired.fastq ./temp4/"+ffiles2[i].replace(".gz","")+" ./temp4/"+ffiles2[i].replace(ppatern1,ppatern2).replace(".gz",""))
	os.system("bash ./unmerge-paired-reads.sh paired.fastq ./temp4/"+ffiles2[i].replace(".gz","")+" ./temp4/"+ffiles2[i].replace(ppatern1,ppatern2).replace(".gz",""))

	print ("wc -l ./temp3/"+ffiles2[i].replace(".gz",""))
	cmd = ["wc -l ./temp3/"+ffiles2[i].replace(".gz","")]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out, err)

	nlines1=int(out.split()[0])
	print (nlines1)

	print ("wc -l ./temp4/"+ffiles2[i].replace(".gz",""))
	cmd = ["wc -l ./temp4/"+ffiles2[i].replace(".gz","")]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out, err)

	nlines2=int(out.split()[0])
	print (nlines2)
	print (100*float(nlines2)/float(nlines1))

	cmd = ["rm *.fastq"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["rm *.blast"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["rm *.log"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

	print ("wc -l ./temp4/"+ffiles2[i].replace(".gz",""))
	cmd = ["wc -l ./temp4/"+ffiles2[i].replace(".gz","")]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out, err)
	
	nseq1=int(out.decode('ascii').split()[0])
	print(); print ("No. sequences: "+str(nseq1/4))

	print ("wc -l ./temp4/"+ffiles2[i].replace(".gz","").replace(ppatern1,ppatern2))
	cmd = ["wc -l ./temp4/"+ffiles2[i].replace(".gz","").replace(ppatern1,ppatern2)]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out, err)
	
	nseq2=int(out.decode('ascii').split()[0])
	print(); print ("No. sequences: "+str(nseq2/4))
	
	if nseq1!=nseq2:
		print ("Different number of sequences")
		exit()

	nlines=0
	with open("./temp4/"+ffiles2[i].replace(".gz","")) as file1, open("./temp4/"+ffiles2[i].replace(".gz","").replace(ppatern1,ppatern2)) as file2:
		for line1, line2 in zip(file1, file2):
			if line1.startswith("@") and line2.startswith("@"):
				if not(line1.endswith("/1\n") and line2.endswith("/2\n")):
					print (line1)
					print (line2)
					print ("Problem!")
					exit()
			nlines+=1	
			   
	print (nlines/4)

	if nseq1!=nlines:
		print (nlines/4)
		print ("Something is wrong!")
		exit()

cmd = ["rm ./db/db.idx.*"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

###########################################################################################################################

# Removes Dokdonia 16s, 23s and 5s rRNA

cmd = ["rm ./db/db.idx.*"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

os.system("indexdb_rna --ref ./db/"+seconddb+",./db/db.idx --fast -v")

cmd = ["rm -r ./temp5"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
cmd = ["mkdir ./temp5"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

from Bio.SeqIO.QualityIO import FastqGeneralIterator

for i in range(len(ffiles2)):
	print (ffiles2[i])

	print ("bash ./merge-paired-reads.sh ./temp4/"+ffiles2[i].replace(".gz","")+" ./temp4/"+ffiles2[i].replace(".gz","").replace(ppatern1,ppatern2)+" merged.fastq")
	cmd = ["bash ./merge-paired-reads.sh ./temp4/"+ffiles2[i].replace(".gz","")+" ./temp4/"+ffiles2[i].replace(".gz","").replace(ppatern1,ppatern2)+" merged.fastq"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out, err)

	print ("cat merged.fastq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > file_sorted.fastq")
	os.system("cat merged.fastq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > file_sorted.fastq")
	os.system("mv file_sorted.fastq merged.fastq")

	print ("sortmerna --ref ./db/"+seconddb+",./db/db.idx --reads merged.fastq --log --num_alignments 1 --blast 0 -a 20 --fastx --aligned ./rrna --other ./norrna")
	os.system("sortmerna --ref ./db/"+seconddb+",./db/db.idx --reads merged.fastq --log --num_alignments 1 --blast 0 -a 20 --fastx --aligned ./rrna --other ./norrna")
	#################################################################################################
	# Removes unpaired reads
	mixed_file = "norrna.fastq"
	paired_file = "paired.fastq"
	
	out_handle = open(paired_file, "w")
	prev = None
	for curr in FastqGeneralIterator(open(mixed_file, "rU")):
	    if curr[0].split()[0].endswith("/1"):
	        prev = curr
	    elif not curr[0].split()[0].endswith("/2"):
	        raise ValueError("Expect IDs to end /1 and /2,\n%s" % curr[0])
	    elif prev and prev[0].split()[0] == curr[0].split()[0][:-2] + "/1":
	        out_handle.write("@%s\n%s\n+\n%s\n" % prev)
	        out_handle.write("@%s\n%s\n+\n%s\n" % curr)
	        prev = None
	out_handle.close()
	
	print ("cat norrna.fastq | wc -l")
	cmd = ["cat norrna.fastq | wc -l"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out, err)
	nseq1=float(out)/4
	
	print ("Sequences before removing unpairs: "+str(nseq1))

	print ("cat paired.fastq | wc -l")
	cmd = ["cat paired.fastq | wc -l"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out, err)
	nseq2=float(out)/4
	
	print ("Sequences after removing unpairs: "+str(nseq2))
	print ("Difference: "+str(nseq1-nseq2))
	print ("Percent difference: "+str(100*(nseq1-nseq2)/nseq1))

	#################################################################################################
	print ("bash ./unmerge-paired-reads.sh paired.fastq ./temp5/"+ffiles2[i].replace(".gz","")+" ./temp5/"+ffiles2[i].replace(ppatern1,ppatern2).replace(".gz",""))
	os.system("bash ./unmerge-paired-reads.sh paired.fastq ./temp5/"+ffiles2[i].replace(".gz","")+" ./temp5/"+ffiles2[i].replace(ppatern1,ppatern2).replace(".gz",""))

	print ("wc -l ./temp4/"+ffiles2[i].replace(".gz",""))
	cmd = ["wc -l ./temp4/"+ffiles2[i].replace(".gz","")]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out, err)

	nlines1=int(out.split()[0])
	print (nlines1)

	print ("wc -l ./temp5/"+ffiles2[i].replace(".gz",""))
	cmd = ["wc -l ./temp5/"+ffiles2[i].replace(".gz","")]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out, err)

	nlines2=int(out.split()[0])
	print (nlines2)
	print (100*float(nlines2)/float(nlines1))

	cmd = ["rm *.fastq"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["rm *.blast"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	cmd = ["rm *.log"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

	print ("wc -l ./temp5/"+ffiles2[i].replace(".gz",""))
	cmd = ["wc -l ./temp5/"+ffiles2[i].replace(".gz","")]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out, err)
	
	nseq1=int(out.split()[0])
	print(); print ("No. sequences: "+str(nseq1/4))

	print ("wc -l ./temp5/"+ffiles2[i].replace(".gz","").replace(ppatern1,ppatern2))
	cmd = ["wc -l ./temp5/"+ffiles2[i].replace(".gz","").replace(ppatern1,ppatern2)]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()
	print (out, err)
	
	nseq2=int(out.split()[0])
	print(); print ("No. sequences: "+str(nseq2/4))
	
	if nseq1!=nseq2:
		print ("Different number of sequences")
		exit()

	nlines=0
	with open("./temp5/"+ffiles2[i].replace(".gz","")) as file1, open("./temp5/"+ffiles2[i].replace(".gz","").replace(ppatern1,ppatern2)) as file2:
		for line1, line2 in zip(file1, file2):
			if line1.startswith("@") and line2.startswith("@"):
				if not(line1.endswith("/1\n") and line2.endswith("/2\n")):
					print (line1)
					print (line2)
					print ("Problem!")
					exit()
			nlines+=1	
			   
	print (nlines/4)

	if nseq1!=nlines:
		print (nlines/4)
		print ("Something is wrong!")
		exit()

cmd = ["rm ./db/db.idx.*"]; pipe = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE); p_status = pipe.wait(); out, err = pipe.communicate()

###########################################################################################################################


exit()



