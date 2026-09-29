"""Select source expression for gse48060_1 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
from io import StringIO
import gzip, pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse48060_1"
    / "process"
)
RAW = HERE.parent / "raw"


def extract():
    """Select the case and control matrices before effect estimation."""
    m = pd.read_csv(HERE / "sample_manifest.tsv", sep="\t")
    m = m[m.disease_group.isin(["case", "healthy_control"])]
    lines = []
    on = False
    with gzip.open(
        RAW / "GSE48060_series_matrix.txt.gz", "rt", encoding="latin1", errors="replace"
    ) as f:
        for l in f:
            if l.startswith("!series_matrix_table_begin"):
                on = True
                continue
            if l.startswith("!series_matrix_table_end"):
                break
            if on:
                lines.append(l)
    x = pd.read_csv(StringIO("".join(lines)), sep="\t", quotechar='"').rename(
        columns={"ID_REF": "gene_id"}
    )
    ids = m.sample_id.tolist()
    x = x[["gene_id"] + ids]
    x[ids] = x[ids].apply(pd.to_numeric, errors="coerce")
    x = x.dropna().groupby("gene_id")[ids].mean()
    ca = m.loc[m.disease_group.eq("case"), "sample_id"].tolist()
    co = m.loc[m.disease_group.eq("healthy_control"), "sample_id"].tolist()
    __case_frame = x[ca]
    A = x[ca].to_numpy()
    __control_frame = x[co]
    return locals()
