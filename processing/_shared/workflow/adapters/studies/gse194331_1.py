"""Select source expression for gse194331_1 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
import csv
import gzip
import numpy as np
import pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse194331_1"
    / "process"
)
RAW = HERE.parent / "raw"


def matrix_metadata(path):
    rows = {}
    with gzip.open(path, "rt", errors="replace", newline="") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            fields = next(csv.reader([line], delimiter="\t"))
            key = fields[0]
            values = [v.strip('"') for v in fields[1:]]
            rows[key] = values
    gsms = rows["!Sample_geo_accession"]
    titles = rows["!Sample_title"]
    source = rows["!Sample_source_name_ch1"]
    characteristics = rows["!Sample_characteristics_ch1"]
    out = []
    for i, gsm in enumerate(gsms):
        text = characteristics[i].lower()
        if "healthy" in text or "control" in text:
            group = "healthy_control"
        elif "ap" in text:
            group = "case"
        else:
            group = "excluded"
        out.append(
            {
                "gsm": gsm,
                "title": titles[i],
                "source": source[i],
                "characteristics": characteristics[i],
                "group": group,
            }
        )
    return pd.DataFrame(out)


def extract():
    """Select the case and control matrices before effect estimation."""
    meta = matrix_metadata(RAW / "GSE194331_series_matrix.txt.gz")
    assert len(meta) == 119, len(meta)
    assert meta.source.str.lower().eq("blood").all()
    case = meta.loc[meta.group.eq("case"), "title"].tolist()
    control = meta.loc[meta.group.eq("healthy_control"), "title"].tolist()
    assert len(case) == 87, len(case)
    assert len(control) == 32, len(control)
    with gzip.open(
        RAW / "GSE194331_HC_PAN_PANSEP_counts.txt.gz", "rt", errors="replace"
    ) as fh:
        counts = pd.read_csv(fh, sep="\t")
    counts = counts.rename(columns={counts.columns[0]: "gene_id"}).set_index("gene_id")
    counts = counts.apply(pd.to_numeric, errors="coerce")
    assert set(case + control).issubset(counts.columns), "count/title mismatch"
    libs = counts.sum(axis=0)
    assert (libs > 0).all()
    __raw_counts = counts.copy()
    logcpm = np.log2(counts.div(libs, axis=1) * 1000000.0 + 1.0)
    __case_frame = logcpm[case]
    A = logcpm[case].to_numpy()
    __control_frame = logcpm[control]
    return locals()
