"""Select source expression for gse269497_1 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
import numpy as np
import pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse269497_1"
    / "process"
)
RAW = HERE.parent / "raw"


def extract():
    """Select the case and control matrices before effect estimation."""
    path = RAW / "GSE269497_genes_FPKM_expression.txt.gz"
    x = pd.read_csv(path, sep="\t", compression="gzip", dtype={"Gene ID": str})
    case = [f"AIT-{i} FPKM" for i in range(1, 6)]
    control = [f"Con-{i} FPKM" for i in range(1, 6)]
    genes = x[["Gene Symbol"] + case + control].copy()
    genes["Gene Symbol"] = genes["Gene Symbol"].astype(str).str.strip()
    genes = genes[genes["Gene Symbol"].str.match("^[A-Za-z0-9_.-]+$")]
    for c in case + control:
        genes[c] = pd.to_numeric(genes[c], errors="coerce")
    genes = genes.dropna(subset=case + control)
    genes = genes.groupby("Gene Symbol", sort=True)[case + control].mean()
    genes = np.log2(genes + 1.0)
    __case_frame = genes[case]
    A = genes[case].to_numpy()
    __control_frame = genes[control]
    return locals()
