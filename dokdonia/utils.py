from __future__ import annotations

import pickle
import subprocess
from pathlib import Path

import pandas as pd
import numpy as np


def saveToPickleFile(python_object, path_to_file="object.pkl"):
    """
    Save python object to pickle file
    """
    out_file = open(path_to_file, "wb")
    pickle.dump(python_object, out_file)
    out_file.close()


def readFromPickleFile(path_to_file="object.pkl"):
    """
    Load python object from pickle file.
    Returns python object.
    """
    in_file = open(path_to_file, "rb")
    python_object = pickle.load(in_file)
    return python_object


def terminal_execute(
    command_str: str,
    suppress_shell_output=False,
    work_dir: Path = None,
    return_output=False,
) -> subprocess.STDOUT:
    """Execute given command in terminal through Python.
    Args:
        command_str (str): terminal command to be executed.
        suppress_shell_output (bool, optional): suppress shell output. Defaults to False.
        work_dir (Path, optional): change working directory. Defaults to None.
        return_output (bool, optional): whether to return execution output. Defaults to False.
    Returns:
        subprocess.STDOUT: subprocess output.
    """
    if suppress_shell_output:
        suppress_code = ">/dev/null 2>&1"
        command_str = f"{command_str} {suppress_code}"
    output = subprocess.run(
        command_str, shell=True, cwd=work_dir, capture_output=return_output
    )
    return output

def compute_replicate_averages(
        data: pd.DataFrame,
        method: str = "mean",
        temperatures: list[int] = [10, 18, 25, 34],
        conditions: list[str] = ["L", "D"]
        ) -> pd.DataFrame:
    """
    Compute average expression among replicates and / or light/dark
    """
    data_dict = {}
    for condition in conditions:
        data_dict[condition] = []
        for temp in temperatures:
            colpattern = f"{condition}_{temp}"
            if "mean" in method:
                series = data.loc[:, [col for col in data.columns if colpattern in col]].mean(axis=1)
            else:
                series = data.loc[:, [col for col in data.columns if colpattern in col]].median(axis=1)
            series.name = colpattern
            data_dict[condition].append(
                series
                )
    return (pd.DataFrame(data_dict[c]).transpose() for c in conditions)

def merge_two_conditions(
    data: pd.DataFrame,
    method: str = "mean",
    conditions: list[str] = ["L", "D"]
    ) -> pd.DataFrame:
    """
    Collapse dataframe columns by condition, either computing the mean or the median.
    All columns named equally after condition removed will be merged.
    """
    data_list = []
    unique_left_conditions = np.unique(
        [col.strip(conditions[0]).strip(conditions[1]).strip("_") for col in data.columns]
        )
    for condition in unique_left_conditions:
        colpattern = condition
        if "mean" in method:
            series = data.loc[:, [col for col in data.columns if colpattern in col]].mean(axis=1)
        else:
            series = data.loc[:, [col for col in data.columns if colpattern in col]].median(axis=1)
        series.name = colpattern
        data_list.append(
            series
            )
    return pd.DataFrame(data_list).transpose()
