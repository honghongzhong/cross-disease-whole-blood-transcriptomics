"""Select source expression for gse112087_1 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
import gzip, io, tarfile
import numpy as np, pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse112087_1"
    / "process"
)
RAW = HERE.parent / "raw"


def extract():
    """Select the case and control matrices before effect estimation."""
    m = pd.read_csv(HERE / "sample_manifest.tsv", sep="\t")
    m = m[m.include_locked.astype(str).str.upper().eq("TRUE")].copy()
    lane_to_donor = {}
    for _, r in m.iterrows():
        for gsm in str(r.lane_gsm_ids).split(";"):
            lane_to_donor[gsm] = r.sample_key
    donor_counts = {}
    with tarfile.open(RAW / "GSE112087_RAW.tar") as t:
        for name in t.getnames():
            gsm = name.split("_", 1)[0]
            if gsm not in lane_to_donor:
                continue
            raw = gzip.decompress(t.extractfile(name).read())
            d = pd.read_csv(
                io.BytesIO(raw), sep="\t", usecols=["EnsemblID", "ReadCount"]
            )
            d["ReadCount"] = pd.to_numeric(d.ReadCount, errors="coerce")
            d = d.groupby("EnsemblID").ReadCount.sum()
            donor_counts[lane_to_donor[gsm]] = (
                d
                if lane_to_donor[gsm] not in donor_counts
                else donor_counts[lane_to_donor[gsm]].add(d, fill_value=0)
            )
    x = pd.DataFrame(donor_counts)
    lib = x.sum()
    __raw_counts = x.copy()
    x = np.log2(x.div(lib, axis=1) * 1000000.0 + 1)
    ca = m.loc[m.disease_group.eq("case"), "sample_key"].tolist()
    co = m.loc[m.disease_group.eq("healthy_control"), "sample_key"].tolist()
    __case_frame = x[ca]
    A = x[ca].to_numpy()
    __control_frame = x[co]
    return locals()
