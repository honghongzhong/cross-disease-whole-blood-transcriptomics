"""Select source expression for gse205465_1 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
import gzip
import numpy as np
import pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse205465_1"
    / "process"
)
RAW = HERE.parent / "raw"
COUNTS = RAW / "counts.tsv.gz"
SOFT = RAW / "family.soft.gz"


def parse_soft():
    text = gzip.open(SOFT, "rt", errors="replace").read()
    rows = []
    for block in text.split("^SAMPLE = ")[1:]:
        lines = block.splitlines()

        def vals(prefix):
            return [
                x.split(" = ", 1)[1].strip().strip('"')
                for x in lines
                if x.startswith(prefix)
            ]

        chars = vals("!Sample_characteristics_ch1")
        joined = " | ".join(chars)
        diagnosis = next(
            (
                x.split(":", 1)[1].strip()
                for x in chars
                if x.lower().startswith("diagnosis:")
            ),
            "",
        )
        rows.append(
            {
                "gsm": vals("!Sample_geo_accession")[0],
                "title": vals("!Sample_title")[0],
                "source": vals("!Sample_source_name_ch1")[0],
                "characteristics": joined,
                "diagnosis": diagnosis,
                "biospecimen": "whole blood"
                if "whole blood"
                in (joined + " " + " ".join(vals("!Sample_source_name_ch1"))).lower()
                else "",
            }
        )
    sm = pd.DataFrame(rows)
    sm["group"] = np.where(
        sm.diagnosis.str.lower().eq("thrombotic primary aps"), "case", "healthy_control"
    )
    sm["count_column"] = sm.gsm
    return sm


def extract():
    """Select the case and control matrices before effect estimation."""
    sm = parse_soft()
    counts = pd.read_csv(COUNTS, sep="\t", compression="gzip", engine="python")
    counts = counts.rename(columns={counts.columns[0]: "gene_id"})
    available = set(counts.columns[1:])
    sm["inclusion"] = np.where(
        sm.gsm.isin(available), "include", "exclude_missing_count"
    )
    included = sm[sm.inclusion == "include"]
    case = included.loc[included.group == "case", "gsm"].tolist()
    control = included.loc[included.group == "healthy_control", "gsm"].tolist()
    if len(case) != 60 or len(control) != 28:
        raise RuntimeError(
            f"Unexpected available groups: cases={len(case)}, controls={len(control)}"
        )
    if not included.biospecimen.eq("whole blood").all():
        raise RuntimeError("Non-whole-blood sample included")
    x = counts.set_index("gene_id")[case + control].apply(
        pd.to_numeric, errors="coerce"
    )
    libsize = x.sum(axis=0)
    __raw_counts = x.copy()
    expr = np.log2(x.div(libsize, axis=1) * 1000000 + 1.0)
    __case_frame = expr[case]
    __control_frame = expr[control]
    return locals()
