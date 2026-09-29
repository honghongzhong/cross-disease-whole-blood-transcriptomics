"""Select source expression for gse165082_1 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
import numpy as np, pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse165082_1"
    / "process"
)
RAW = HERE.parent / "raw"


def extract():
    """Select the case and control matrices before effect estimation."""
    m = pd.read_csv(HERE / "sample_manifest.tsv", sep="\t")
    m["sample_id"] = m.sample_id.astype(str).str.zfill(3)
    m["count_group"] = m.group.replace({"Control": "CC"})
    x = pd.read_csv(
        RAW / "GSE165082_PD-CC.counts.txt.gz", sep="\t", compression="gzip"
    ).set_index("Geneid")
    ids = m.sample_id.str.cat(m.count_group, sep="_").tolist()
    x = x[ids].apply(pd.to_numeric, errors="coerce")
    __raw_counts = x.copy()
    x = np.log2(x.div(x.sum(), axis=1) * 1000000.0 + 1)
    ca = (
        m.loc[m.group.eq("PD"), "sample_id"]
        .str.cat(m.loc[m.group.eq("PD"), "count_group"], sep="_")
        .tolist()
    )
    co = (
        m.loc[m.group.eq("Control"), "sample_id"]
        .str.cat(m.loc[m.group.eq("Control"), "count_group"], sep="_")
        .tolist()
    )
    __case_frame = x[ca]
    A = x[ca].to_numpy()
    __control_frame = x[co]
    return locals()
