# This script computes reads from cleaned Fastaq files. It follows the workflow
# described in: http://www.compbio.dundee.ac.uk/user/pschofield/Projects/teaching_pg/workshops/biocNGS.html
#
# Rsubread publication: The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads
#
# A few notes:
# 1. Run script through 'ssh: nohup Rscript quant.R &' to let process continue
#    after ending remote session. Terminal output logged in nohup.out.
# 2. buildIndex cannot handle very large fils, that's why I'm not using the
#    annotated fasta file in the server
# 3. Grant all permissions to workdir: 'sudo chmod -R a+rwx /path/to/folder'
# 4. An additional workflow based in Rsubread can be found at: http://monashbioinformaticsplatform.github.io/RNAseq-DE-analysis-with-R/RNAseq_DE_analysis_with_R.html
# 5. Alignment can be parallelized with 'BiocParallel' following directions in the
#    first workflow.
#
# Semidán Robaina Estévez, March 2021.

library(Rsubread)

rm(list=ls())
RNAseqDATADIR <- "/Data/LauraDokdoniaReadsCleaned"
fastq_files <- dir(RNAseqDATADIR)
REF_GENOME <- "/Data/DokdoniaMED134_full.fasta"
Annotated_GTF <- "/Data/DokdoniaMED134.gff"
RSUBREAD_INDEX_PATH <- "/Data/ref_data"
RSUBREAD_INDEX_BASE <- "MED134"
forward_pattern <- "_1.fastq.gz"
reverses_pattern <- "_2.fastq.gz"
BAM_OUTPUT_PATH <- "/Data/BAM_files"

getDataIDs <- function() {
  files <- dir(RNAseqDATADIR)
  conditions <- vector()
  for (file in files) {
    if (grepl(forward_pattern, file, fixed=TRUE) == TRUE) {
       conditions <- c(conditions, strsplit(file, forward_pattern)[[1]])
    }
  }
  return(conditions)
}

alignSequences <- function(conditions) {
  # This loop can be parallelized, see above publication
  for (condition in conditions) {
    print(paste0("Aligning condition: ", condition))
    # Define paired-end fasta files
    inputfilefwd <- file.path(RNAseqDATADIR, paste0(condition, "_1.fastq.gz"))
    inputfilervs <- file.path(RNAseqDATADIR, paste0(condition, "_2.fastq.gz"))

    # run the align command to map the reads
    align(index=file.path(RSUBREAD_INDEX_PATH, RSUBREAD_INDEX_BASE),
     readfile1=inputfilefwd, readfile2=inputfilervs,
     output_file=file.path(BAM_OUTPUT_PATH, paste0(condition, ".sam")),
     output_format="SAM")
  }
}

countReads <- function() {
  bamFiles <- list.files(BAM_OUTPUT_PATH)
  curdir <- getwd()
  setwd(file.path(curdir, BAM_OUTPUT_PATH))
  print("Counting reads that match overlap exons and grouping exons by gene_id")
  fcLim <- featureCounts(files=bamFiles,
    GTF.featureType="exon", GTF.attrType="gene_id",
    annot.ext=Annotated_GTF, isGTFAnnotationFile=TRUE)
  setwd(curdir)
  print("Saving results")
  knitr::kable(fcLim$stat) # Print stats
  save(fcLim, file="/Data/LauraDokdoniaCounts.RData")
}

# buildindex(basename=file.path(RSUBREAD_INDEX_PATH, RSUBREAD_INDEX_BASE), reference=REF_GENOME)
# conditions <- getDataIDs()
# alignSequences(conditions)
countReads()
