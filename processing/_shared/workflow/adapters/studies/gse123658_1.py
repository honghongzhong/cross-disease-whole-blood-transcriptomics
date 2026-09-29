"""Select source expression for gse123658_1 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
import numpy as np, pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse123658_1"
    / "process"
)
RAW = HERE.parent / "raw"


def extract():
    """Select the case and control matrices before effect estimation."""
    m = pd.read_csv(HERE / "sample_manifest.tsv", sep="\t")
    m = m[m.include_locked.astype(str).str.upper().eq("TRUE")]
    x = pd.read_csv(
        RAW / "GSE123658_read_counts.gene_level.txt.gz", sep="\t", compression="gzip"
    ).set_index("Samples")
    ids = m.sample_key.tolist()
    x = x[ids].apply(pd.to_numeric, errors="coerce")
    __raw_counts = x.copy()
    x = np.log2(x.div(x.sum(), axis=1) * 1000000.0 + 1)
    ca = m.loc[m.disease_group.eq("case"), "sample_key"].tolist()
    co = m.loc[m.disease_group.eq("healthy_control"), "sample_key"].tolist()
    __case_frame = x[ca]
    A = x[ca].to_numpy()
    __control_frame = x[co]
    return locals()
