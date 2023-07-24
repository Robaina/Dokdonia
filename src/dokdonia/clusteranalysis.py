from sklearn.metrics import silhouette_samples
import pandas as pd
import numpy as np
import os
from subprocess import call


def writeClustInputFiles(clust_data, path_to_wd="Data"):
    clust_data.to_csv(os.path.join(path_to_wd, "clust_input.tsv"), sep="\t")
    conds = np.unique([s[:4] for s in clust_data.columns])
    open(os.path.join(path_to_wd, "clust_replicates.txt"), "w").close()
    with open(os.path.join(path_to_wd, "clust_replicates.txt"), "a+") as file:
        for cond in conds:
            reps = ",".join(list(clust_data.filter(regex=f"{cond}").columns))
            txt_s = f"clust_input.tsv, {cond}, {reps}\n"
            file.write(txt_s)

    open(os.path.join(path_to_wd, "clust_no_normalization.txt"), "w").close()
    with open(os.path.join(path_to_wd, "clust_no_normalization.txt"), "a+") as file:
        file.write("clust_input.tsv 0")  # no normalization

    open(
        os.path.join(path_to_wd, "clust_transcript_cell_normalization.txt"), "w"
    ).close()
    with open(
        os.path.join(path_to_wd, "clust_transcript_cell_normalization.txt"), "a+"
    ) as file:
        file.write("clust_input.tsv 101 4")  # Quantile + z-score

    open(os.path.join(path_to_wd, "clust_TPM_normalization.txt"), "w").close()
    with open(os.path.join(path_to_wd, "clust_TPM_normalization.txt"), "a+") as file:
        file.write("clust_input.tsv 101 3 4")  # Quantile + log2 + z-score


def runClust(
    path_to_wd,
    out_dir,
    cluster_tightness=1,
    replicates_file=None,
    normalization_file=None,
):
    """
    Compute clusters with clust
    clust_data: pandas DataFrame.
    """
    call_list = [
        "clust",
        os.path.join(path_to_wd, "clust_input.tsv"),
        "-t",
        f"{cluster_tightness}",
        "-o",
        f"{out_dir}",
    ]
    if replicates_file is not None:
        call_list.append("-r")
        call_list.append(os.path.join(path_to_wd, replicates_file))
    if normalization_file is None:
        call_list.append("-n")
        call_list.append(os.path.join(path_to_wd, "clust_no_normalization.txt"))
    elif normalization_file is "auto":
        pass
    else:
        call_list.append("-n")
        call_list.append(os.path.join(path_to_wd, normalization_file))

    call(call_list, cwd=path_to_wd)


def getGeneClusters(
    clust_data,
    path_to_wd,
    out_dir,
    cluster_tightness=1,
    replicates_file=None,
    normalization_file=None,
    scaling_factor=1e4,
):
    "Returns dict with Clust gene clusters"
    clust_data = clust_data.multiply(scaling_factor)  # scale to avoid numerical issues
    writeClustInputFiles(clust_data, path_to_wd)
    runClust(
        path_to_wd=path_to_wd,
        out_dir=out_dir,
        cluster_tightness=cluster_tightness,
        replicates_file=replicates_file,
        normalization_file=normalization_file,
    )

    clusters = pd.read_csv(
        os.path.join(out_dir, "Clusters_Objects.tsv"), sep="\t", header=1
    )
    return {
        f"C{i}": clusters.iloc[:, i].dropna().values.tolist()
        for i in range(clusters.shape[1])
    }


def groupDataByTemperature(
    data: pd.DataFrame, statistic="mean", normalize=True
) -> pd.DataFrame:
    """
    Group data by temperature
    """
    if statistic == "mean":
        data = data.groupby(
            data.columns.str.extract(r"_(\d+)_", expand=False), axis=1
        ).mean()
    elif statistic == "median":
        data = data.groupby(
            data.columns.str.extract(r"_(\d+)_", expand=False), axis=1
        ).median()
    else:
        raise ValueError(f"Unknown statistic: {statistic}")
    if normalize:
        data = data.div(data.max(axis=1), axis=0)
    return data


def computeGeneSilhouettes(clusters, data):
    """
    Compute gene silhouettes for each gene in clusters
    """
    genes_in_clusters, cluster_labels = [], []
    for cluster_id, cluster in clusters.items():
        genes_in_clusters.extend(cluster)
        cluster_labels.extend([cluster_id for _ in range(len(cluster))])

    X = data.loc[genes_in_clusters, :].values
    sil_values = silhouette_samples(X, cluster_labels)

    return {gene_id: sil_values[i] for i, gene_id in enumerate(genes_in_clusters)}


def computeGeneAverageExpression(data: pd.DataFrame) -> dict:
    """Get median expression level across samples

    Args:
        data (pd.DataFrame): _description_
    """
    return data.median(axis=1).to_dict()


def rankGenesWithinClusters(clusters, data, method: str = "silhouette"):
    """Rank genes within each cluster based on their silhouette or average expression level

    Args:
        clusters (_type_): _description_
        data (_type_): _description_
        method (str, optional): _description_. Defaults to "silhouette".
    """
    if "sil" in method:
        return rankGenesWithinClustersUsingSilhouette(clusters, data)
    else:
        return rankGenesWithinClustersAverageExpression(clusters, data)


def rankGenesWithinClustersUsingSilhouette(clusters, data):
    """
    Rank genes within each cluster based on their silhouette
    """

    gene_sil = computeGeneSilhouettes(clusters, data)

    ranked_clusters = {}
    for cluster_id, cluster in clusters.items():
        sil_dict = {gene_id: gene_sil[gene_id] for gene_id in cluster}
        ranked_dict = dict(
            sorted(sil_dict.items(), key=lambda item: item[1], reverse=True)
        )
        ranked_clusters[cluster_id] = ranked_dict

    return ranked_clusters


def rankGenesWithinClustersAverageExpression(clusters, data):
    """
    Rank genes within each cluster based on their average expression level
    """

    gene_expr = computeGeneAverageExpression(data)

    ranked_clusters = {}
    for cluster_id, cluster in clusters.items():
        expr_dict = {gene_id: gene_expr[gene_id] for gene_id in cluster}
        ranked_dict = dict(
            sorted(expr_dict.items(), key=lambda item: item[1], reverse=True)
        )
        ranked_clusters[cluster_id] = ranked_dict

    return ranked_clusters


def writeExcelOfClusterGenes(
    clusters, out_path, gene_info, gene_pathways, gene_systems
):
    """
    Write excel file with gene product and KEGG pathways for genes in each cluster.
    """
    silhouette = type(clusters[list(clusters.keys())[0]]) == dict
    writer = pd.ExcelWriter(out_path, engine="xlsxwriter")

    for cluster_id, cluster in clusters.items():
        cluster_info = {}
        for gene_id in cluster:
            if gene_id in gene_info.keys():
                if silhouette:
                    cluster_info[gene_id] = {
                        "Gene silhouette": cluster[gene_id],
                        "Gene product": gene_info[gene_id],
                        "KEGG system": ", ".join(gene_systems[gene_id]),
                        "KEGG subsystem": ", ".join(gene_pathways[gene_id]),
                    }
                else:
                    cluster_info[gene_id] = {
                        "Gene product": gene_info[gene_id],
                        "KEGG system": ", ".join(gene_systems[gene_id]),
                        "KEGG subsystem": ", ".join(gene_pathways[gene_id]),
                    }

        pd.DataFrame(cluster_info).transpose().to_excel(
            writer, sheet_name=f"Cluster {cluster_id}"
        )
    writer.save()
