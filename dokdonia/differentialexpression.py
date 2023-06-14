from __future__ import annotations
import pandas as pd
import numpy as np
from diffexpr.py_deseq import py_DESeq2
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
from rpy2.robjects import Formula
import logging

rpy2_logger.setLevel(logging.ERROR)


def getAverageStandardRatio(IS_counts, standards_data):
    """
    Compute average standard ratios from standard counts and metadata
    """
    cond_out = "D_25_R1"
    IS = IS_counts["index"].values
    conditions = IS_counts.columns.values.tolist()
    conditions.remove("index")
    conditions.remove(cond_out)

    avg_st_ratios = {}
    for cond_id in conditions:
        st_ratios = []
        for st_id in IS:

            st_copies = standards_data[
                (standards_data["Sample ID"] == cond_id)
                & (standards_data["Standard"] == st_id)
            ]["Standard added (copias)"].values[0]

            st_counts = IS_counts[IS_counts["index"] == st_id][cond_id].values[0]

            st_ratios.append(st_copies / st_counts)

        avg_st_ratios[cond_id] = {
            "average": np.mean(st_ratios),
            "std": np.std(st_ratios),
            "cv": np.std(st_ratios) / np.mean(st_ratios),
        }
    return avg_st_ratios


def getTranscriptsPerCell(counts, avg_st_ratios, abundance_meta):
    """
    Normalize counts by internal standards and cell abundances
    transcripts/cell = (counts * avg_st_ratio) / total_cell_abundance

    """
    cond_out = "D_25_R1"  # no average standard ratio data for this condition
    conditions = abundance_meta["Sample"].values.tolist()
    conditions.remove(cond_out)

    n_counts = counts[counts.columns.intersection(conditions + ["index"])].copy()
    for cond_id in conditions:
        n_cells = abundance_meta[abundance_meta["Sample"] == cond_id][
            "Total_cell_abundance"
        ].values[0]
        avg_ratio = avg_st_ratios[cond_id]["average"]
        n_counts[cond_id] = (n_counts[cond_id] * avg_ratio) / n_cells

    return n_counts


def get_transcript_cell(counts, abundance_meta, cond_outs: list[str] = None):
    """
    Get total cell abundance for a condition
    """
    if cond_outs is None:
        cond_outs = [
            "D_25_R1",
            "L_18_R4",
            "D_18_R4",
        ]  # no average standard ratio data for this condition
    conditions = abundance_meta["Sample"].values.tolist()
    for cond_out in cond_outs:
        conditions.remove(cond_out)

    n_counts = counts[counts.columns.intersection(conditions + ["index"])].copy()
    for cond_id in conditions:
        n_cells = abundance_meta[abundance_meta["Sample"] == cond_id][
            "Total_cell_abundance"
        ].values[0]
        avg_ratio = abundance_meta[abundance_meta["Sample"] == cond_id][
            "Avg_Standard_ratio"
        ].values[0]
        n_counts[cond_id] = (n_counts[cond_id] * avg_ratio) / n_cells

    return n_counts


def get_transcripts_volume(
    counts, abundance_meta, cond_outs: list[str] = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Get total cell abundance for a condition
    """
    if "Sample" not in abundance_meta.columns:
        abundance_meta["Sample"] = (
            abundance_meta["Light/Dark"]
            + "_"
            + abundance_meta["Temperature"].astype(str)
            + "_"
            + abundance_meta["Replicate"]
        )
    if cond_outs is None:
        cond_outs = ["D_25_R1"]  # no average standard ratio data for this condition
    conditions = abundance_meta["Sample"].values.tolist()
    for cond_out in cond_outs:
        conditions.remove(cond_out)

    n_counts = counts[counts.columns.intersection(conditions + ["index"])].copy()
    n_counts_biovol = counts[counts.columns.intersection(conditions + ["index"])].copy()
    for cond_id in conditions:
        n_cells = abundance_meta[abundance_meta["Sample"] == cond_id][
            "Total abundance (cells)"
        ].values[0]
        biovolume = abundance_meta[abundance_meta["Sample"] == cond_id][
            "Total biovolume (um3)"
        ].values[0]
        avg_ratio = abundance_meta[abundance_meta["Sample"] == cond_id][
            "Avg_Standard_ratio"
        ].values[0]
        n_counts[cond_id] = (n_counts[cond_id] * avg_ratio) / n_cells
        n_counts_biovol[cond_id] = (n_counts_biovol[cond_id] * avg_ratio) / biovolume

    return n_counts, n_counts_biovol


def scale_counts_by_internal_standard(counts, abundance_meta):
    """
    Get total cell abundance for a condition
    """
    cond_outs = [
        "D_25_R1",
        "L_18_R4",
        "D_18_R4",
    ]  # no average standard ratio data for this condition
    conditions = abundance_meta["Sample"].values.tolist()
    for cond_out in cond_outs:
        conditions.remove(cond_out)

    n_counts = counts[counts.columns.intersection(conditions + ["index"])].copy()
    for cond_id in conditions:
        avg_ratio = abundance_meta[abundance_meta["Sample"] == cond_id][
            "Avg_Standard_ratio"
        ].values[0]
        n_counts[cond_id] = n_counts[cond_id] * avg_ratio

    return n_counts


def scale_counts_by_cell_abundance(counts, abundance_meta):
    """
    Get total cell abundance for a condition
    """
    cond_outs = [
        "D_25_R1",
        "L_18_R4",
        "D_18_R4",
    ]  # no average standard ratio data for this condition
    conditions = abundance_meta["Sample"].values.tolist()
    for cond_out in cond_outs:
        conditions.remove(cond_out)

    n_counts = counts[counts.columns.intersection(conditions + ["index"])].copy()
    for cond_id in conditions:
        n_cells = abundance_meta[abundance_meta["Sample"] == cond_id][
            "Total_cell_abundance"
        ].values[0]
        n_counts[cond_id] = (n_counts[cond_id]) / n_cells

    return n_counts


def getMetaMatrix(counts):
    return pd.DataFrame(
        {
            "lighting": [s[0] for s in counts.columns[1:]],
            "temperature": [int(s[2:4]) for s in counts.columns[1:]],
            "replicate": [s[5:] for s in counts.columns[1:]],
        },
        index=counts.columns[1:],
    )


def runDEtest(
    counts,
    test="Wald",
    alpha=1e-2,
    formula="~ lighting",
    reduced_formula=None,
    log2fold_cutoff=0,
) -> tuple:
    """
    Runs DeSeq2
    reduced_formula only for LRT test
    log2fold_cutoff: threshold to consider genes as DE when
    pair-wise comparisons with Wald test
    (Schurch et al., 2016 recommends 0.5 for 3 replicates)
    """
    meta = getMetaMatrix(counts)
    dds = py_DESeq2(
        count_matrix=counts,
        design_matrix=meta,
        design_formula=formula,
        gene_column="index",
    )

    if test == "LRT":
        dds.run_deseq(test=test, reduced=Formula(reduced_formula))
        dds.get_deseq_result(alpha=alpha)
    else:
        dds.run_deseq(test=test)
        dds.get_deseq_result(alpha=alpha)  # , lfcThreshold=log2fold_cutoff)
    res = dds.deseq_result
    res = res[res.padj < alpha]
    stats = {
        "DE+": res.log2FoldChange.where(res.log2FoldChange >= log2fold_cutoff)
        .dropna()
        .shape[0],
        "DE-": res.log2FoldChange.where(res.log2FoldChange <= -log2fold_cutoff)
        .dropna()
        .shape[0],
    }
    return (res, stats)


def addProteinNamesToDeseqResult(gbk, deseq_res):
    proteins = [gbk.getGeneInfo(gene_id)["product"][0] for gene_id in deseq_res.index]
    deseq_res["proteins"] = proteins
    return deseq_res


def deseq2Normalize(counts, coldata, gene_column="gene_id", design_formula=f"~ 1"):
    """
    Run deseq2
    """
    dds = py_DESeq2(
        count_matrix=counts,
        design_matrix=coldata,
        design_formula=design_formula,
        gene_column=gene_column,
    )
    dds.run_deseq()
    deseq2_counts = dds.normalized_count()
    return deseq2_counts


def tpmNormalize(counts, patric):
    gene_lengths = {}
    for gene_id in counts["index"]:
        try:
            gene_lengths[gene_id] = patric.loc[patric["RefSeq Locus Tag"] == gene_id][
                "Length"
            ].item()
        except Exception:
            gene_lengths[gene_id] = np.nan

    rpk = counts.iloc[:, 1:].divide(list(gene_lengths.values()), axis=0)
    tpm = rpk.divide(rpk.sum(axis=0).values / 1e6, axis=1)
    tpm.insert(0, "index", counts["index"])
    return tpm
