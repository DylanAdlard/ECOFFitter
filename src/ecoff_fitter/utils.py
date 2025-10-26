import pandas as pd
import yaml
import numpy as np
import os


def read_input(data, sheet_name=None):
    """
    Read MIC input data from a DataFrame or file and validate required columns.

    Args:
        data (str | DataFrame): Input source — a pandas DataFrame or file path
            (.csv, .tsv, .txt, .xlsx, .xls).
        sheet_name (str, optional): Excel sheet name to read if input is an Excel file.

    Returns:
        DataFrame: Cleaned MIC data with columns:
            - "MIC" (str): MIC values.
            - "observations" (int): Observation counts.

    Raises:
        ValueError: If input type is invalid or required columns are missing.
        FileNotFoundError: If the specified file does not exist.
    """

    if isinstance(data, pd.DataFrame):
        df = data.copy()
    elif isinstance(data, str):

        ext = os.path.splitext(data)[-1].lower()

        if ext in [".csv"]:
            df = pd.read_csv(data)
        elif ext in [".tsv", ".txt"]:
            df = pd.read_csv(data, sep="\t")
        elif ext in [".xlsx", ".xls"]:
            df = pd.read_excel(data, sheet_name=sheet_name)
        else:
            raise ValueError(f"Unsupported file type: {ext}")
    else:
        raise ValueError("Input must be a pandas DataFrame or path to file.")

    expected_cols = ["MIC", "observations"]
    missing = [c for c in expected_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df["MIC"] = df["MIC"].astype(str).str.strip()
    df["observations"] = (
        pd.to_numeric(df["observations"], errors="coerce").fillna(0).astype(int)
    )

    df = df.dropna(subset=["MIC"])
    df = df.reset_index(drop=True)

    return df


def read_params(params):
    """
    Read ECOFF model parameters from a file or dictionary.

    Args:
        params (str | dict): File path to a YAML or text parameter file,
            or a dictionary containing configuration values.

    Returns:
        tuple: (dilution_factor, distributions, tail_dilutions)

    Raises:
        FileNotFoundError: If the specified file does not exist.
        ValueError: If file type is unsupported.
        AssertionError: If input is not a valid file path or dictionary.
    """

    if isinstance(params, str):

        if not os.path.exists(params):
            raise FileNotFoundError(f"Parameter file not found: {params}")

        ext = os.path.splitext(params)[-1].lower()

        if ext in [".yaml", "yml"]:
            with open(params, "r") as f:
                params = yaml.safe_load(f)
        elif ext == ".txt":
            params = {}
            with open(params, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    key, val = [x.strip() for x in line.split("=", 1)]
                    if key == "dilution_factor":
                        params[key] = int(val)
                    elif key == "tail_dilutions":
                        if val.lower() == "none":
                            params[key] = None
                        else:
                            params[key] = int(val)
                    elif key == "distributions":
                        params[key] = int(val)
                    else:
                        params[key] = val
        else:
            raise ValueError(f"Unsupported file type: {ext}")

    else:

        assert isinstance(
            params, dict
        ), "params must either be a file path or a dictionary"

    return (
        params["dilution_factor"],
        params["distributions"],
        params["tail_dilutions"],
    )
