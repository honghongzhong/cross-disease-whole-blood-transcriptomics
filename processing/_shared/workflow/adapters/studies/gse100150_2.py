"""Select source expression for gse100150_2 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
from io import StringIO
import gzip
import pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse100150_2"
    / "process"
)
RAW = HERE.parent / "raw"


def read_matrix(path):
    lines = []
    with gzip.open(path, "rt", encoding="latin1", errors="replace") as f:
        active = False
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                active = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if active:
                lines.append(line)
    return pd.read_csv(StringIO("".join(lines)), sep="\t", quotechar='"')


def extract():
    """Select the case and control matrices before effect estimation."""
    meta = pd.read_csv(HERE / "sample_manifest.tsv", sep="\t")
    meta = meta[meta.qc_pass.astype(str).str.upper().eq("TRUE")].copy()
    meta["group"] = meta["class"].map({"Case": "Case", "Healthy": "Healthy"})
    matrix = read_matrix(RAW / "GSE100150_series_matrix.txt.gz").rename(
        columns={"ID_REF": "probe_id"}
    )
    gsm = meta.gsm.tolist()
    if any((x not in matrix.columns for x in gsm)):
        raise ValueError("Selected GSM is absent from matrix")
    ann = pd.read_csv(
        RAW / "GPL6884.annot.gz",
        sep="\t",
        skiprows=28,
        dtype=str,
        encoding="latin1",
        low_memory=False,
    )
    ann = (
        ann.rename(columns={"ID": "probe_id", "Gene symbol": "gene_id"})[
            ["probe_id", "gene_id"]
        ]
        .dropna()
        .drop_duplicates("probe_id")
    )
    x = matrix[["probe_id"] + gsm].merge(ann, on="probe_id", how="inner")
    for c in gsm:
        x[c] = pd.to_numeric(x[c], errors="coerce")
    x = (
        x[x.gene_id.str.match("^[A-Za-z0-9_.-]+$", na=False)]
        .groupby("gene_id")[gsm]
        .mean()
    )
    case = meta.loc[meta.group.eq("Case"), "gsm"].tolist()
    ctrl = meta.loc[meta.group.eq("Healthy"), "gsm"].tolist()
    __case_frame = x[case]
    __control_frame = x[ctrl]
    return locals()
