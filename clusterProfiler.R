library(clusterProfiler)
options(clusterProfiler.download.method = "wget")
# library(enrichplot)
# library(ggplot2)


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


work_dir <- "/home/robaina/Documents/Aquifex/Dokdonia/"


# counts_path <- paste0(work_dir, "data/DokdoniaCounts.csv")
# total_genes <- read.delim(counts_path, sep=",", header=FALSE)[-1, 1]

# clusters_path <- paste0(work_dir, "results/CLUSTER_ALL_GENES_TRANSCRIPT_CELL/Clusters_Objects.tsv")
clusters_path <- paste0(work_dir, "results/CLUSTER_NONDE_GENES_TRANSCRIPT_CELL/Clusters_Objects.tsv")
# clusters_path <- paste0(work_dir, "results/CLUSTER_ALL_GENES_TPM/Clusters_Objects.tsv")
clusters <- read.delim(clusters_path, sep="\t", header=FALSE)

cluster_genes <- lapply(clusters[-(1:2),], function(cluster) {cluster[cluster != ""];})

# unclustered_genes <- setdiff(total_genes, unlist(cluster_genes))
# cluster_genes <- append(cluster_genes, list(unclustered_genes))

universe <- unlist(cluster_genes)


cp_ora <- compareCluster(
  geneClusters = cluster_genes, 
  fun = "enrichKEGG", # ORA function to apply to each cluster
  organism = 'dok',
  pvalueCutoff = 0.1,
  qvalueCutoff = 1,
  universe = universe,
  minGSSize = 5,
  pAdjustMethod = "BH", # p-values are adjusted within clusters
)

# Export results for all clusters
write.csv(cp_ora@compareClusterResult, "enrichment_results/results_non_DE_genes.csv", row.names=FALSE)

# Visualization
# urls <- vector()
# i <- 1
# for (pathway in kk$ID) {
#     url <- browseKEGG(kk, pathway)
#     urls <- c(urls, paste0(as.character(i), ". ", url))
#     i <- i + 1
#     print(pathway, url)
# }
# outfilename <- "enrichment_results/CLUSTER_ALL_GENES_TRANSCRIPT_CELL.txt"
# outfile <- file.create(outfilename)
# outfile <- file(outfilename, open="w")
# writeLines(as.character(urls), outfile)
# close(outfile)


# # geneList <- DEgenes[,2]  # DE values
# # names(geneList) <- DEgenes[,1]  # Gene IDs
# # genes <- names(geneList)
# kk <- enrichKEGG(gene = cluster_genes,
#                  organism = 'dok',
#                  pvalueCutoff = 1,
#                  qvalueCutoff = 1,
#                  universe = total_genes,
#                  minGSSize = 10)


# # Gene set enrichment
# kk2 <- gseKEGG(geneList = geneList,
#                organism = 'pel',
#                minGSSize = 3,
#                maxGSSize = 800,
#                pAdjustMethod = "BH",
#                pvalueCutoff = 0.05,
#                verbose = FALSE)

# # K-module gene set enrichment
# mkk2 <- gseMKEGG(geneList = geneList,
#                  organism = 'pel',
#                  pvalueCutoff = 0.05)

# emapplot(kk2, showCategory = 10)