from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

from dokdonia.utils import quantile_normalize


def plotClusters(pdata, clusters, outfile: str = None, figsize: tuple = (15, 18)):
    n_rows = int(np.ceil(len(clusters) / 2))
    fig, axes = plt.subplots(nrows=n_rows, ncols=2)
    plt.subplots_adjust(hspace=0.3)
    coords = list(np.ndindex((n_rows, 2)))
    for n, cluster_id in enumerate(clusters):
        i, j = coords[n]
        try:
            axis = axes[i, j]
        except Exception:
            axis = axes[j]
        cluster = clusters[cluster_id]
        pdata[pdata.index.isin(cluster)].transpose().plot(
            legend=False,
            figsize=figsize,
            title=f"{cluster_id}, size={len(cluster)}",
            ax=axis,
            color="#9a9a9a",
            linewidth=0.8,
            marker=".",
            markerfacecolor="#ee9929",
            markersize=12,
        )
    for axis in axes.flatten():
        if not axis.lines:
            axis.set_visible(False)
    if outfile is not None:
        plt.savefig(outfile, dpi=300)
    plt.show()


def plotClusterData(pdata, cluster, ax=None, cluster_id=None):
    pdata[pdata.index.isin(cluster)].transpose().plot(
        legend=False,
        title=f"{cluster_id}, size={len(cluster)}",
        ax=ax,
        color="#9a9a9a",
        linewidth=0.8,
        marker=".",
        markerfacecolor="#ee9929",
        markersize=12,
    )


def plot_density(df, ax=None, show_plot=True, title=None, show_legend=False):
    sns.set(style="whitegrid")
    if ax is None:
        fig, ax = plt.subplots()

    qn_df = quantile_normalize(df)
    qn_sample = qn_df.iloc[:, 0]

    for column in df.columns:
        if pd.api.types.is_numeric_dtype(df[column]):
            sns.kdeplot(df[column], label=column, ax=ax, alpha=1)

    sns.kdeplot(
        qn_sample, label="Quantile Normalized", ax=ax, alpha=0.5, color="black", lw=2.8
    )
    if show_legend:
        ax.legend()
    ax.set_xlabel("log2(expression)")
    ax.set_ylabel("Density")
    if title is not None:
        ax.set_title(title)
    if show_plot:
        plt.show()
    return ax
