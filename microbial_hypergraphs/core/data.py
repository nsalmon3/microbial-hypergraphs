"""
Manages datasets that networks can be built off of
"""

import numpy as np
import pandas as pd

from microbial_hypergraphs.core.env import ROOT_PATH

DATA_PATH = ROOT_PATH / "microbial_hypergraphs" / "data"


def get_taxonomy() -> pd.DataFrame:
    """
    Returns a pandas DataFrame representing the taxonomy information of different OTUs, with some cleaning done from the raw file.
    """
    df = pd.read_csv(DATA_PATH / "raw" / "combined_taxonomy.csv", index_col="OTU_id")

    df["Domain"] = df["Domain"].str.strip().str.removeprefix("k__")
    df["Domain"] = df["Domain"].replace("Eukaryota_kgd_Incertae_sedis", "Eukaryota")
    df["Domain"] = df["Domain"].replace("Unassigned", None)

    df["Phylum"] = df["Phylum"].str.strip().str.removeprefix("p__")
    df["Phylum"] = df["Phylum"].str.removeprefix("[")
    df["Phylum"] = df["Phylum"].str.removesuffix("]")
    df["Phylum"] = df["Phylum"].replace("", None)
    df["Phylum"] = df["Phylum"].replace(np.nan, None)

    df["Class"] = df["Class"].str.strip().str.removeprefix("c__")
    df["Class"] = df["Class"].str.removeprefix("[")
    df["Class"] = df["Class"].str.removesuffix("]")
    df["Class"] = df["Class"].replace("", None)
    df["Class"] = df["Class"].replace(np.nan, None)

    df["Order"] = df["Order"].str.strip().str.removeprefix("o__")
    df["Order"] = df["Order"].str.removeprefix("[")
    df["Order"] = df["Order"].str.removesuffix("]")
    df["Order"] = df["Order"].replace("", None)
    df["Order"] = df["Order"].replace(np.nan, None)

    df["Family"] = df["Family"].str.strip().str.removeprefix("f__")
    df["Family"] = df["Family"].str.removeprefix("[")
    df["Family"] = df["Family"].str.removesuffix("]")
    df["Family"] = df["Family"].replace("", None)
    df["Family"] = df["Family"].replace(np.nan, None)

    df["Genus"] = df["Genus"].str.strip().str.removeprefix("g__")
    df["Genus"] = df["Genus"].str.removeprefix("[")
    df["Genus"] = df["Genus"].str.removesuffix("]")
    df["Genus"] = df["Genus"].replace("", None)
    df["Genus"] = df["Genus"].replace(np.nan, None)

    df["Species"] = df["Species"].str.strip().str.removeprefix("s__")
    df["Species"] = df["Species"].str.removeprefix("[")
    df["Species"] = df["Species"].str.removesuffix("]")
    df["Species"] = df["Species"].replace("", None)
    df["Species"] = df["Species"].replace(np.nan, None)

    return df


def get_sample_info() -> pd.DataFrame:
    """
    Cleans and returns sample information

    :return: Cleaned dataframe
    :rtype: DataFrame
    """

    df = pd.read_csv(DATA_PATH / "raw" / "design_10302021.csv", header=[0, 1])

    df.columns = df.columns.droplevel(0)
    df.set_index("Sample_No", inplace=True)

    return df[
        [
            "Country",
            "Region",
            "Plot_ID",
            "Lat",
            "Long",
            "Treatment_main",
            "Treatment_minor",
        ]
    ].copy()


def get_otu_samples() -> pd.DataFrame:
    """
    Cleans and returns OTU information for each sample spot

    :return: DataFrame of OTU values for each sample
    :rtype: DataFrame
    """

    df = pd.read_csv(DATA_PATH / "raw" / "combined_otu_table.csv", index_col="OTU_id")
    df.drop(columns=["Domain"], inplace=True)
    df.columns = [int(col[1:4]) for col in df.columns]
    df.columns.name = "Sample_No"
    df = df.T

    return df
