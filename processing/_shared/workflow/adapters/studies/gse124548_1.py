"""Select source expression for gse124548_1 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
import numpy as np, pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse124548_1"
    / "process"
)
RAW = HERE.parent / "raw"


def extract():
    """Select the case and control matrices before effect estimation."""
    m = pd.read_csv(HERE / "sample_manifest.tsv", sep="\t")
    m = m[m.gsm.notna()]
    x = pd.read_csv(
        RAW / "GSE124548_raw_baseline_counts.tsv.gz", sep="\t", compression="gzip"
    ).set_index("gene_id")
    ids = m.raw_name.tolist()
    x = x[ids].apply(pd.to_numeric, errors="coerce")
    __raw_counts = x.copy()
    x = np.log2(x.div(x.sum(), axis=1) * 1000000.0 + 1)
    ca = m.loc[m.diagnosis.eq("CF"), "raw_name"].tolist()
    co = m.loc[m.diagnosis.eq("healthy control"), "raw_name"].tolist()
    __case_frame = x[ca]
    A = x[ca].to_numpy()
    __control_frame = x[co]
    return locals()
