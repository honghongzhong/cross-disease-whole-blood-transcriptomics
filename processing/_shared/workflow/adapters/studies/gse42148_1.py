"""Select source expression for gse42148_1 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
import csv
import gzip
import pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse42148_1"
    / "process"
)
ROOT = HERE.parent
RAW = ROOT / "raw"


def read_matrix(path):
    metadata = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_"):
                parts = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
                metadata.setdefault(parts[0], []).append(parts[1:])
            if line.startswith("!series_matrix_table_begin"):
                d = pd.read_csv(
                    fh, sep="\t", quotechar='"', comment="!", low_memory=False
                )
                break
        else:
            raise RuntimeError("series matrix table not found")
    d = d.rename(columns={d.columns[0]: "probe_id"})
    d["probe_id"] = d["probe_id"].astype(str)
    return (d, metadata)


def extract():
    """Select the case and control matrices before effect estimation."""
    matrix, meta = read_matrix(RAW / "GSE42148_series_matrix.txt.gz")
    accessions = meta["!Sample_geo_accession"][0]
    source = meta["!Sample_source_name_ch1"][0]
    characteristics = meta["!Sample_characteristics_ch1"]
    if not len(accessions) == len(source) == len(matrix.columns) - 1:
        raise RuntimeError("metadata and matrix column counts do not match")

    def value(x):
        return x.split(":", 1)[1].strip().lower() if ":" in x else x.strip().lower()

    tissue = [value(x) for x in characteristics[0]]
    state = [value(x) for x in characteristics[3]]
    if set(tissue) != {"whole blood"}:
        raise RuntimeError("not all selected samples are labelled whole blood")
    cases = [a for a, s in zip(accessions, state) if s == "coronary artery disease"]
    controls = [a for a, s in zip(accessions, state) if s == "control"]
    if len(cases) != 13 or len(controls) != 11:
        raise RuntimeError(
            f"unexpected groups: cases={len(cases)}, controls={len(controls)}"
        )
    if len(set(accessions)) != len(accessions):
        raise RuntimeError("duplicate GEO sample accession")
    x = matrix.set_index("probe_id")
    __case_frame = x[cases].apply(pd.to_numeric, errors="coerce")
    case = x[cases].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    __control_frame = x[controls].apply(pd.to_numeric, errors="coerce")
    return locals()
