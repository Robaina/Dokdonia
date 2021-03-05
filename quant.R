# Workflow: http://monashbioinformaticsplatform.github.io/RNAseq-DE-analysis-with-R/RNAseq_DE_analysis_with_R.html

"
FASTA file provided by Jose is too large to be used in buildindex.
An alternative could be to employ a FASTA with the full genome without annotations
and use the gtf file after mapping to find annotated sequences. This is described
in the second web page.
"
# Another one: http://www.compbio.dundee.ac.uk/user/pschofield/Projects/teaching_pg/workshops/biocNGS.html
# Grant all permissions to working directory: sudo chmod -R a+rwx /path/to/folder
"
Rsubread publication: The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads
"

library(Rsubread)
# library(BiocParallel)
rm(list=ls())

# define shared directory for RNAseq data
RNAseqDATADIR <- "Data/LauraDokdoniaReadsCleaned"
# list the fastq files in the raw data directory
fastq_files <- dir(RNAseqDATADIR)

# define the reference genome fasta file
REF_GENOME <- "Data/DokdoniaMED134_full.fasta"
Annotated_GTF <- "Data/DokdoniaMED134.gtf"
# define the output directory for the Rsubread index
# (admin note: requires data/ref_data/download_hg19.sh to be run first)
RSUBREAD_INDEX_PATH <- "Data/ref_data"
# define the basename for the index
RSUBREAD_INDEX_BASE <- "MED134"

# build the index
buildindex(basename=file.path(RSUBREAD_INDEX_PATH, RSUBREAD_INDEX_BASE), reference=REF_GENOME)

# Alignment
# define the fastq file with forward reads
forward_pattern <- "_1.fastq.gz"
reverses_pattern <- "_2.fastq.gz"
OUTPUT_PATH <- file.path("Data/BAM_files", "")

files <- dir(RNAseqDATADIR)
conditions <- vector()
for (file in files) {
  if (grepl(forward_pattern, file, fixed=TRUE) == TRUE) {
     conditions <- c(conditions, strsplit(file, forward_pattern)[[1]])
  }
}

# This loop can be parallelized, see above publication
for (condition in conditions) {

  print(paste0("Aligning condition: ", condition))
  # Define paired-end fasta files
  inputfilefwd <- file.path(RNAseqDATADIR, paste0(condition, "_1.fastq.gz"))
  inputfilervs <- file.path(RNAseqDATADIR, paste0(condition, "_2.fastq.gz"))

  # run the align command to map the reads
  align(index=file.path(RSUBREAD_INDEX_PATH, RSUBREAD_INDEX_BASE),
   readfile1=inputfilefwd, readfile2=inputfilervs,
   output_file=file.path(OUTPUT_PATH, paste0(condition, ".sam")),
   output_format="SAM")
}
