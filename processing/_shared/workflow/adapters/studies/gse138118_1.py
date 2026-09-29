"""Select source expression for gse138118_1 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
import csv
import gzip
import pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse138118_1"
    / "process"
)
RAW = HERE.parent / "raw"


def parse_matrix_metadata(path):
    rows = {}
    with gzip.open(path, "rt", errors="replace", newline="") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_"):
                f = next(csv.reader([line], delimiter="\t"))
                rows[f[0]] = [v.strip('"') for v in f[1:]]
    n = len(rows["!Sample_geo_accession"])
    out = []
    for i in range(n):
        title = rows["!Sample_title"][i]
        source = rows["!Sample_source_name_ch1"][i]
        label = f"{title} {source}".lower()
        group = (
            "case"
            if "ucb positive" in label
            else "healthy_control"
            if "health volunteer" in label
            else "excluded"
        )
        out.append(
            {
                "gsm": rows["!Sample_geo_accession"][i],
                "title": title,
                "source": source,
                "group": group,
            }
        )
    return pd.DataFrame(out)


def extract():
    """Select the case and control matrices before effect estimation."""
    matrix = RAW / "GSE138118_series_matrix.txt.gz"
    meta = parse_matrix_metadata(matrix)
    assert len(meta) == 75
    assert meta.title.str.contains("Blood", case=False).all()
    case = meta.loc[meta.group.eq("case"), "gsm"].tolist()
    control = meta.loc[meta.group.eq("healthy_control"), "gsm"].tolist()
    assert len(case) == 46
    assert len(control) == 29
    with gzip.open(matrix, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                header = next(csv.reader([next(fh)], delimiter="\t"))
                break
        data = pd.read_csv(fh, sep="\t", comment="!", header=None, names=header)
    data.columns = [c.strip('"') for c in data.columns]
    data = data.rename(columns={data.columns[0]: "probe_id"}).set_index("probe_id")
    data = data.apply(pd.to_numeric, errors="coerce")
    assert set(case + control).issubset(data.columns)
    __case_frame = data[case]
    __control_frame = data[control]
    return locals()
