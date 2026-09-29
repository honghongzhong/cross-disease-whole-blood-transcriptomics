"""Select source expression for gse48000_1 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
import numpy as np, pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse48000_1"
    / "process"
)
RAW = HERE.parent / "raw"


def qn(a):
    o = np.argsort(a, axis=0)
    s = np.take_along_axis(a, o, axis=0)
    v = s.mean(1)
    z = np.empty_like(a)
    for j in range(a.shape[1]):
        z[o[:, j], j] = v
    return z


def extract():
    """Select the case and control matrices before effect estimation."""
    m = pd.read_csv(HERE / "sample_manifest.tsv", sep="\t")
    x = pd.read_csv(
        RAW / "GSE48000_non_normalized_set1.txt.gz", sep="\t", compression="gzip"
    ).set_index("ID_REF")
    ids = m.subject_id.tolist()
    x = x[ids].apply(pd.to_numeric, errors="coerce")
    x = pd.DataFrame(
        qn(np.log2(x.clip(lower=1).to_numpy())), index=x.index, columns=ids
    )
    a = (
        pd.read_csv(
            RAW / "GPL10558.annot.gz",
            sep="\t",
            skiprows=28,
            dtype=str,
            encoding="latin1",
            low_memory=False,
        )
        .rename(columns={"ID": "probe_id", "Gene symbol": "gene_id"})[
            ["probe_id", "gene_id"]
        ]
        .dropna()
        .drop_duplicates("probe_id")
    )
    x = (
        x.reset_index(names="probe_id")
        .merge(a, on="probe_id")
        .drop(columns="probe_id")
        .groupby("gene_id")[ids]
        .mean()
    )
    ca = m.loc[m.group.eq("High-risk VTE"), "subject_id"].tolist()
    co = m.loc[m.group.eq("Healthy"), "subject_id"].tolist()
    __case_frame = x[ca]
    A = x[ca].to_numpy()
    __control_frame = x[co]
    return locals()
