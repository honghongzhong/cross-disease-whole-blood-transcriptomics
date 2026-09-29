"""Select source expression for gse111459_1 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
import numpy as np, pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse111459_1"
    / "process"
)
RAW = HERE.parent / "raw"


def extract():
    """Select the case and control matrices before effect estimation."""
    m = pd.read_csv(HERE / "sample_manifest.tsv", sep="\t")
    m = m[m.include_locked.astype(str).str.upper().eq("TRUE")].copy()
    ids = m.sample_key.tolist()
    raw = pd.read_csv(
        RAW / "GSE111459_TBM_Blood_RawCount.txt.gz", sep="\t", compression="gzip"
    )
    if any((x not in raw.columns for x in ids)):
        raise ValueError("sample column missing")
    counts = raw[ids].apply(pd.to_numeric, errors="coerce")
    gene = raw["Gene"].astype(str)
    lib = counts.sum(axis=0)
    __raw_counts = counts.copy()
    x = np.log2(counts.div(lib, axis=1) * 1000000.0 + 1)
    x = x.groupby(gene).mean()
    ca = m.loc[m.disease_group.eq("case"), "sample_key"].tolist()
    co = m.loc[m.disease_group.eq("healthy_control"), "sample_key"].tolist()
    __case_frame = x[ca]
    A = x[ca].to_numpy()
    __control_frame = x[co]
    return locals()
