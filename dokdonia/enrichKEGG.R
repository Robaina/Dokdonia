#!/usr/bin/env Rscript

library(clusterProfiler)
options(clusterProfiler.download.method = "wget")

work_dir <- "/home/robaina/Documents/Aquifex/Dokdonia/"

all_genes <- read.csv(
    paste0(work_dir, "results/transcript_cell_clust_input.csv")
)
temp_independent <- read.csv(
    paste0(work_dir, "results/results_temperature_independent/temp_independent_genes_avg_expression.csv")
)

enrich <- enrichKEGG(
    gene = unlist(temp_independent$gene_id),
    organism = "dok",
    pvalueCutoff = 0.05,
    keyType = "kegg",
    pAdjustMethod = "BH",
    universe = unlist(all_genes$ID),
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    use_internal_data = FALSE
)
