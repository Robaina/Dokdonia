# Bioinformatic Methods
All code generated in this study is available in the following repository: https://github.com/Robaina/Dokdonia and has been organized in a Jupyter Notebook, with Python 3 as default programming environment.

## RNAseq raw data preprocessing and read counts
Pair-end reads in the FastQ files were first preprocessed to remove low quality sequences, adapters and contaminating RNA. To this end, we first employed AfterQC [REF], using deftault settings, to remove low quality sequences as well as to trim adapters. We then filtered out Dokdonia 5s, 16s and 23s rRNA sequences and phage ΦX174 sequences with a custom Python script available in the GitHub repository indicated above..

After preprocessing, pair-end reads were mapped to the Dokdonia MED134 genome (accession number: CP009301) using Rsubread’s function align [REFs] with default settings (allowing, at most, one mismatch). Read counts were obtained with Rsubread’s function featureCount using default settings for pair-end reads.

## Differential expression analysis
Differentially expressed genes were obtained with DeSeq2 [REF], a tool developed to analyze differential expression directly from cout data via shrinkage estimation for dispersions and fold changes. To analyze differential expression between light and dark conditions for each temeperature, we employed the Wald test within DeSeq2 using default settings and with a p-value cuttof of 0.01 and fold-change cuttof value of 0.5 to compensate for the small number of replicates (≤4) in our experiment [REF to Schurch et al)]. To obtain the sets of differentially expressed genes across temperatures --- under light, dark and light or dark conditions --- we applied the Likelihood Ratio Test within DeSeq2 with default settings and a p-value cuttof of 0.01.

## Gene clustering
Gene clusters were obtained with Clust, a clustering algorithm specifically designed to cluster gene expression data [REF]. Contrary to other methods, such as those based on hierarchical or K-means approaches, Clust minimizes the variance within each cluster and does not attempt to assign a cluster to every gene in the sample. As a result, Clust minimizes spurios assignments of genes that are not trully correlated. We employed Clust with default parameters with the exception of the cluster tighness parameter, which was set to 5 --- smaller values decreased the resolution of the clusters, producing aggregated clusters of genes with upward or downward trends across temperatures).

To use Clust, we first applied a TPM normalization to the raw count data using a custom function, and then selected only the genes for which the Likelihood Ratio Test perfromed by DeSeq2 was positive, i.e., genes which were differentially expressed across temperatures. Finally, we applied the default normalization of Clust for TPM-normalized data, i.e., quantile normalization, followed by log2 and z-score transformations.

## Gene function annotation
Gene function annotations for Dokdonia MED134 were generated using the online version of the EggNOG mapper [REF], to this end, we prepared a FASTA file containing all protein sequences from the GeneBank file with version CP009301.1. After annotation, we selected the KEGG pathway Ko annotations returned by EggNOG mapper and employed a custom Python function, getKEGGPathwayDict, to assign KEGG classes and superclasses to each pathway.

# Results & Discussion

## Analysis of differential expression
Two main environmental factors potentially affecting gene expression varied in the culture conditons employed in our study: light and temperature. Thus, were were interested in assessing to which extent gene expression was influenced by each factor. To this end, we first performed an analysis of differential expression between light and dark conditions for each culture temperature, followed by an analysis of differential expression across temperatures.

Employing a fold-change cutoff value of 0.5 (Methods), we found few differentially expressed (DE) genes between light and dark conditions. Specifically, we found a set of 17 genes which were DE at 10°C and 18°C, with greater expression in light than in dark. Among these genes, we found a subset directly involved in the proteorhodopsin-mediated light-haversting system (Table X), a system previously described in Dokdonia and which is ubiquitous in marine microbes [REF]

## Gene clusters
