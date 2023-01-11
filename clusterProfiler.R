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


work_dir <- "/home/robaina/Documents/Aquifex/Dokdonia/"


counts_path <- paste0(work_dir, "data/DokdoniaCounts.csv")
total_genes <- read.delim(counts_path, sep=",", header=FALSE)[-1, 1]

clusters_path <- paste0(work_dir, "results/CLUSTER_ALL_GENES_TRANSCRIPT_CELL/Clusters_Objects.tsv")
clusters <- read.delim(clusters_path, sep="\t", header=FALSE)

cluster_genes <- clusters[-(1:2),1]


# geneList <- DEgenes[,2]  # DE values
# names(geneList) <- DEgenes[,1]  # Gene IDs
# genes <- names(geneList)
kk <- enrichKEGG(gene = cluster_genes,
                 organism = 'dok',
                 pvalueCutoff = 1,
                 qvalueCutoff = 1,
                 universe = total_genes,
                 minGSSize = 10)

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




# # Some extra functions:
# mkk <- enrichMKEGG(gene = genes,
#                    organism = 'pel',
#                    pvalueCutoff = 0.05)

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