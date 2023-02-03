import pickle
import subprocess
from pathlib import Path


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
