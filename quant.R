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
# 5. Alignment can be parallelized with 'BiocParallel' following directions in the first workflow.
# 6. featureCounts options: https://rdrr.io/bioc/Rsubread/man/featureCounts.html
# 7. SAM format: http://www.metagenomics.wiki/tools/samtools/bam-sam-file-format
# Semidán Robaina Estévez, March 2021.

library(Rsubread)

rm(list = ls())
RNAseqDATADIR <- "Data/LauraDokdoniaReadsCleaned"
fastq_files <- dir(RNAseqDATADIR)
REF_GENOME <- "Data/DokdoniaMED134_full.fasta" # "Data/DokdoniaMED134_80char.fasta"
Annotated_GTF <- "Data/DokdoniaMED134.gtf" # "Data/DokdoniaMED134_geneIDs.gtf"
RSUBREAD_INDEX_PATH <- "Data/ref_data"
RSUBREAD_INDEX_BASE <- "MED134"
forward_pattern <- "_1.fastq.gz"
reverses_pattern <- "_2.fastq.gz"
SAM_OUTPUT_PATH <- "Data/SAM_files"

getDataIDs <- function() {
  files <- dir(RNAseqDATADIR)
  conditions <- vector()
  for (file in files) {
    if (grepl(forward_pattern, file, fixed = TRUE) == TRUE) {
       conditions <- c(conditions, strsplit(file, forward_pattern)[[1]])
    }
  }
  return(conditions)
}

alignSequences <- function(conditions, ncores = 14, paired = TRUE) {
  inputfilesfwd <- vector()
  inputfilesrev <- vector()
  output_files <- vector()
  for (condition in conditions) {
    inputfilesfwd <- c(inputfilesfwd, file.path(
      RNAseqDATADIR, paste0(condition, "_1.fastq.gz")))
    if (paired == TRUE) {
      inputfilesrev <- c(inputfilesrev, file.path(
       RNAseqDATADIR, paste0(condition, "_2.fastq.gz")))
    }
    output_files <- c(output_files, file.path(
      SAM_OUTPUT_PATH, paste0(condition, ".sam")))
  }
  if (paired == TRUE) {
    align(index = file.path(RSUBREAD_INDEX_PATH, RSUBREAD_INDEX_BASE),
       readfile1 = inputfilesfwd,
       readfile2 = inputfilesrev,
       type = "rna",
       nthreads = ncores,
       output_file = output_files,
       output_format = "SAM"
     )
  } else {
    align(index = file.path(RSUBREAD_INDEX_PATH, RSUBREAD_INDEX_BASE),
       readfile1 = inputfilesfwd,
       readfile2 = NULL,
       type = "rna",
       nthreads = ncores,
       output_file = output_files,
       output_format = "SAM"
     )
  }
 }

alignSequencesBioParallel <- function(conditions, ncores = 8) {
  # Rsubread doesn't parallelize very well... only 14% CPU usage!
  # I'll use here an external parallelization tool
  # Carefull, runs out of RAM with 14 threads...
  library(BiocParallel)
  bpparam <- MulticoreParam(ncores)

  myAlign <- function(condition){
             align(index=file.path(RSUBREAD_INDEX_PATH, RSUBREAD_INDEX_BASE),
                   readfile1=file.path(
                     RNAseqDATADIR, paste0(condition, "_1.fastq.gz")
                   ),
                   readfile2=file.path(
                     RNAseqDATADIR, paste0(condition, "_2.fastq.gz")
                   ),,
                   type = "rna",
                   output_file=file.path(
                     SAM_OUTPUT_PATH, paste0(condition, ".sam")
                   ),
                   output_format = "SAM"
                 )
            }

  bplapply(conditions, myAlign, BPPARAM = bpparam)
}

countReads <- function(paired_end = TRUE, ncores = 8) {
  samFiles <- list.files(SAM_OUTPUT_PATH, pattern = ".*sam$")
  curdir <- getwd()
  setwd(file.path(curdir, SAM_OUTPUT_PATH))
  print("Counting reads that match overlap exons and grouping exons by gene_id")

  fcLim <- featureCounts(files = samFiles,
    # annotation
    annot.ext = file.path(curdir, Annotated_GTF),
    isGTFAnnotationFile = TRUE,
    GTF.featureType = "gene",
    GTF.attrType = "gene_id",
    GTF.attrType.extra = NULL,
    chrAliases = NULL,

    # level of summarization
    useMetaFeatures = TRUE,

    # overlap between reads and features
    allowMultiOverlap = TRUE,
    minOverlap = 1,
    fracOverlap = 0,
    fracOverlapFeature = 0,
    largestOverlap = FALSE,
    nonOverlap = NULL,
    nonOverlapFeature = NULL,

    # Read shift, extension and reduction
    readShiftType = "upstream",
    readShiftSize = 0,
    readExtension5 = 0,
    readExtension3 = 0,
    read2pos = NULL,

    # multi-mapping reads
    countMultiMappingReads = TRUE,

    # fractional counting
    fraction = FALSE,

    # long reads
    isLongRead = FALSE,

    # read filtering
    minMQS = 0,
    splitOnly = FALSE,
    nonSplitOnly = FALSE,
    primaryOnly = FALSE,
    ignoreDup = FALSE,

    # strandness
    strandSpecific = 0,

    # exon-exon junctions
    juncCounts = FALSE,
    genome = NULL,

    # parameters specific to paired end reads
    isPairedEnd = paired_end,
    countReadPairs = TRUE,
    requireBothEndsMapped = FALSE,
    checkFragLength = FALSE,
    minFragLength = 50,
    maxFragLength = 600,
    countChimericFragments = TRUE,
    autosort = TRUE,

    # number of CPU threads
    nthreads = ncores,

    # read group
    byReadGroup = FALSE,

    # report assignment result for each read
    reportReads = NULL,
    reportReadsPath = NULL,

    # miscellaneous
    maxMOp = 10,
    tmpDir = ".",
    verbose = FALSE
  )

  setwd(curdir)
  print("Saving results")
  save(fcLim, file="Data/LauraDokdoniaCounts.RData")
}

# buildindex(basename=file.path(RSUBREAD_INDEX_PATH, RSUBREAD_INDEX_BASE), reference=REF_GENOME)
conditions <- getDataIDs()
alignSequences(conditions, ncores = 10, paired = TRUE)
countReads(paired_end = TRUE, ncores = 10)
