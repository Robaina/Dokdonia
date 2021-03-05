library(Rsubread)
rm(list=ls())

# define shared directory for RNAseq data
RNAseqDATADIR <- "Data/LauraDokdoniaReadsCleaned"
# list the fastq files in the raw data directory
fastq_files <- dir(RNAseqDATADIR)

# define the reference genome fasta file
REF_GENOME <- "Data/DokdoniaMED134.fasta"
# define the output directory for the Rsubread index
# (admin note: requires data/ref_data/download_hg19.sh to be run first)
RSUBREAD_INDEX_PATH <- "Data/ref_data"
# define the basename for the index
RSUBREAD_INDEX_BASE <- "MED134"

# build the index
buildindex(basename=file.path(RSUBREAD_INDEX_PATH, RSUBREAD_INDEX_BASE), reference=REF_GENOME)

# define the fastq file with forward reads
forward_pattern = "_1.fastq.gz"
reverses_pattern = "_1.fastq.gz"

inputfilefwd <- file.path(RNAseqDATADIR,"D_10_R1_1.fastq.gz")
# define the fastq file with reverse reads
inputfilervs <- file.path(RNAseqDATADIR,"D_10_R1_2.fastq.gz")

# run the align command to map the reads
#align(index=file.path(RSUBREAD_INDEX_PATH, RSUBREAD_INDEX_BASE), readfile1=inputfilefwd, readfile2=inputfilervs, output_file="ERR420388.sam", output_format="SAM")
