#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
condition <- args[1]
clusters_path <- args[2]
results_dir <- args[3]

library(clusterProfiler)
# options(clusterProfiler.download.method = "wget")


clusters <- read.delim(clusters_path, sep = "\t", header = FALSE)
cluster_genes <- lapply(clusters[-(1:2), ], function(cluster) {
  cluster[cluster != ""]
})
universe <- unlist(cluster_genes)

# Tutorial:
#  https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
#
# Notes:
#  geneList is a named vector containing log2 fold changes
#  red genes in KEGG are DE genes in our dataset

# TODO:
# """
# 1. Create nested list of clusters
# 2. Add one cluster with all genes not included in any other cluster
# 3. Ensure clusters form a partition (no duplicates)
# 4. Try below function to run ORA (hypergeometric + BH)
# """

cp_ora <- compareCluster(
  geneClusters = cluster_genes,
  fun = "enrichKEGG", # ORA function to apply to each cluster
  organism = "dok",
  pvalueCutoff = 0.1,
  qvalueCutoff = 1,
  universe = universe,
  minGSSize = 5,
  pAdjustMethod = "BH", # p-values are adjusted within clusters
)

# Export results for all clusters
write.csv(
  cp_ora@compareClusterResult,
  paste0(results_dir, "/results_", condition, ".csv"),
  row.names = FALSE
)
