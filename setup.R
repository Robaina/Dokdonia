if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2", update = FALSE, ask = FALSE, force=TRUE)
BiocManager::install("Rsubread", update = FALSE, ask = FALSE, force=TRUE)
BiocManager::install("clusterProfiler", update = FALSE, ask = FALSE, force=TRUE)