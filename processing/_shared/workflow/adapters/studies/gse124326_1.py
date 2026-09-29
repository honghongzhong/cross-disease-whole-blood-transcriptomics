"""Select source expression for gse124326_1 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
import gzip
import numpy as np
import pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse124326_1"
    / "process"
)
RAW = HERE.parent / "raw"


def extract():
    """Select the case and control matrices before effect estimation."""
    meta = pd.read_csv(HERE / "sample_manifest.tsv", sep="\t", dtype=str)
    meta = meta[meta["selected_author_final_analysis"].eq("True")].copy()
    case = meta.loc[meta["group"].isin(["BP1", "BP2"]), "title"].tolist()
    control = meta.loc[meta["group"].eq("Control"), "title"].tolist()
    assert len(case) == 239 and len(control) == 205
    with gzip.open(RAW / "GSE124326_count_matrix.txt.gz", "rt", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
    sample_cols = {c.removesuffix(".counts"): c for c in header[1:]}
    ids = case + control
    assert set(ids).issubset(sample_cols)
    case_cols = [sample_cols[s] for s in case]
    control_cols = [sample_cols[s] for s in control]
    cols = case_cols + control_cols
    x = pd.read_csv(
        RAW / "GSE124326_count_matrix.txt.gz",
        sep="\t",
        compression="gzip",
        usecols=["gene"] + cols,
    )
    x["gene"] = x["gene"].astype(str).str.replace("\\..*$", "", regex=True)
    x[cols] = x[cols].apply(pd.to_numeric, errors="coerce")
    x = x.groupby("gene", sort=True)[cols].sum()
    __raw_counts = x.copy()
    x = np.log2(x.div(x.sum(axis=0), axis=1) * 1000000.0 + 1.0)
    __case_frame = x[case_cols]
    __control_frame = x[control_cols]
    return locals()
