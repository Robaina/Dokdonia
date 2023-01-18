if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Rsubread", update = FALSE, ask = FALSE)
BiocManager::install("clusterProfiler", update = FALSE, ask = FALSE)