"""Select source expression for gse89594_1 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
import gzip
import pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse89594_1"
    / "process"
)
ROOT = HERE.parent
RAW = ROOT / "raw"


def read_matrix(path):
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                d = pd.read_csv(
                    fh, sep="\t", quotechar='"', comment="!", low_memory=False
                )
                break
        else:
            raise RuntimeError("series matrix table not found")
    d = d.rename(columns={d.columns[0]: "probe_id"})
    d["probe_id"] = d["probe_id"].astype(str)
    return d


def extract():
    """Select the case and control matrices before effect estimation."""
    matrix = read_matrix(RAW / "GSE89594_series_matrix.txt.gz")
    controls = list(matrix.columns[1:31])
    asd = list(matrix.columns[31:63])
    cases = list(matrix.columns[63:95])
    if len(controls) != 30 or len(asd) != 32 or len(cases) != 32:
        raise RuntimeError(
            f"unexpected groups: {len(controls)}, {len(asd)}, {len(cases)}"
        )
    x = matrix.set_index("probe_id")
    __case_frame = x[cases].apply(pd.to_numeric, errors="coerce")
    case = x[cases].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    __control_frame = x[controls].apply(pd.to_numeric, errors="coerce")
    return locals()
