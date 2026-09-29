"""Select source expression for gse334922_1 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
import gzip
import io
import numpy as np
import pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse334922_1"
    / "process"
)
RAW = HERE.parent / "raw" / "samples"


def read_sample(gsm: str):
    path = RAW / f"{gsm}.xlsx.gz"
    if not path.exists():
        raise FileNotFoundError(path)
    data = gzip.open(path, "rb").read()
    df = pd.read_excel(io.BytesIO(data), sheet_name=0)
    if df.shape[1] != 2 or str(df.columns[0]).lower() != "gene":
        raise RuntimeError(f"Unexpected expression layout for {gsm}: {df.shape}")
    gene = df.iloc[:, 0].astype(str)
    value = pd.to_numeric(df.iloc[:, 1], errors="coerce")
    if gene.duplicated().any():
        df = (
            pd.DataFrame({"gene_id": gene, "value": value})
            .groupby("gene_id", sort=True)["value"]
            .mean()
            .reset_index()
        )
    else:
        df = pd.DataFrame({"gene_id": gene, "value": value})
    return df.set_index("gene_id")["value"]


def extract():
    """Select the case and control matrices before effect estimation."""
    sm = pd.read_csv(HERE / "sample_manifest.tsv", sep="\t", dtype=str)
    included = sm[sm.inclusion == "include"]
    case = included.loc[included.group == "case", "gsm"].tolist()
    control = included.loc[included.group == "healthy_control", "gsm"].tolist()
    if len(case) != 4 or len(control) != 4:
        raise RuntimeError(
            f"Unexpected primary sample counts: cases={len(case)}, controls={len(control)}"
        )
    if set(case) & set(control):
        raise RuntimeError("Case/control overlap")
    series = {gsm: read_sample(gsm) for gsm in case + control}
    x = pd.DataFrame(series, dtype=float)
    x = np.log2(x + 1.0)
    __case_frame = x[case]
    __control_frame = x[control]
    return locals()
