"""Select source expression for gse51405_2 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
from io import StringIO
import gzip
import pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse51405_2"
    / "process"
)
RAW = HERE.parent / "raw"


def read_matrix(path):
    lines = []
    active = False
    with gzip.open(path, "rt", errors="replace") as f:
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
    assert len(meta) == 83 and meta.gsm.is_unique
    assert meta.source.str.lower().eq("whole blood").all()
    selected = meta[meta.primary_scope.eq("primary")].copy()
    case = selected.loc[selected.illness.str.startswith("CVID"), "gsm"].tolist()
    control = selected.loc[selected.illness.eq("Healthy"), "gsm"].tolist()
    assert len(case) == 59 and len(control) == 24
    x = (
        read_matrix(RAW / "series_matrix.txt.gz")
        .rename(columns={"ID_REF": "gene_id"})
        .set_index("gene_id")
    )
    x.columns = [str(c).strip('"') for c in x.columns]
    x = (
        x[case + control]
        .apply(pd.to_numeric, errors="coerce")
        .dropna()
        .groupby(level=0, sort=True)
        .mean()
    )
    __case_frame = x[case]
    __control_frame = x[control]
    return locals()
