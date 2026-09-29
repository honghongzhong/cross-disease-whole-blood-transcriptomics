"""Select source expression for gse34404_1 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
from io import StringIO
import gzip
import math
import pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse34404_1"
    / "process"
)
RAW = HERE.parent / "raw"


def extract():
    """Select the case and control matrices before effect estimation."""
    metadata = pd.read_csv(HERE / "sample_manifest.tsv", sep="\t", dtype=str)
    metadata = metadata[metadata["selected_primary"].eq("TRUE")].copy()
    case_subjects = metadata.loc[
        metadata["group"].eq("symptomatic_P_falciparum_malaria"), "subject_id"
    ].tolist()
    control_subjects = metadata.loc[
        metadata["group"].eq("age_matched_healthy_control"), "subject_id"
    ].tolist()
    assert len(case_subjects) == 94 and len(control_subjects) == 61
    with gzip.open(RAW / "GSE34404_non-normalized.txt.gz", "rt", errors="replace") as f:
        header = f.readline().rstrip("\r\n").split("\t")
        signal_columns = [x for x in header[1:] if x.endswith("_AVG_Signal")]
        signal_subjects = [x[: -len("_1_AVG_Signal")] for x in signal_columns]
        assert set(signal_subjects) == set(case_subjects + control_subjects)
        selected_columns = [0] + [
            header.index(s + "_1_AVG_Signal") for s in case_subjects + control_subjects
        ]
        rows = []
        for line in f:
            values = line.rstrip("\r\n").split("\t")
            if len(values) != len(header):
                continue
            try:
                numeric = [float(values[i]) for i in selected_columns[1:]]
            except ValueError:
                continue
            if all((math.isfinite(v) for v in numeric)):
                rows.append([values[0]] + numeric)
    matrix = pd.DataFrame(rows, columns=["probe_id"] + case_subjects + control_subjects)
    with open(RAW / "GPL10558_acc.cgi", encoding="latin1", errors="replace") as f:
        lines = f.read().splitlines()
    begin = lines.index("!platform_table_begin")
    annotation = pd.read_csv(
        StringIO("\n".join(lines[begin + 1 :])), sep="\t", dtype=str, low_memory=False
    )[["ID", "Symbol"]].rename(columns={"ID": "probe_id", "Symbol": "gene_id"})
    annotation["gene_id"] = annotation["gene_id"].fillna("").str.strip()
    annotation = annotation[annotation["gene_id"].str.match("^[A-Za-z0-9_.-]+$")]
    matrix = matrix.merge(annotation, on="probe_id", how="inner")
    expression = matrix.drop(columns="probe_id").groupby("gene_id", sort=True).mean()
    __case_frame = expression[case_subjects]
    A = expression[case_subjects].to_numpy(dtype=float)
    __control_frame = expression[control_subjects]
    return locals()
