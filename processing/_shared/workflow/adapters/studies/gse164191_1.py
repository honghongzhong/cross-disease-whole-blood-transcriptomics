"""Select source expression for gse164191_1 within the shared workflow."""

"Source study reader staged as an input adapter for reproduction."
from pathlib import Path
import gzip
from io import StringIO
import pandas as pd

HERE = (
    Path(__import__("os").environ["BIOLOGY_WORK_ROOT"])
    / "dataset"
    / "gse164191_1"
    / "process"
)
RAW = HERE.parent / "raw"


def extract():
    """Select the case and control matrices before effect estimation."""
    soft = gzip.open(RAW / "GSE164191_family.soft.gz", "rt", errors="replace").read()
    rows = []
    for block in soft.split("^SAMPLE = ")[1:]:

        def get(prefix):
            for line in block.splitlines():
                if line.startswith(prefix):
                    return line.split(" = ", 1)[1].strip().strip('"')
            return ""

        gsm = get("!Sample_geo_accession")
        title = get("!Sample_title")
        status = (
            get("!Sample_characteristics_ch1 = disease status:")
            .split(":", 1)[-1]
            .strip()
        )
        tissue = get("!Sample_characteristics_ch1 = tissue:").split(":", 1)[-1].strip()
        group = (
            "case"
            if status == "colorectal cancer"
            else "healthy_control"
            if status == "normal"
            else "excluded"
        )
        rows.append((gsm, title, group, tissue))
    meta = pd.DataFrame(rows, columns=["gsm", "title", "group", "tissue"])
    meta = meta[
        (meta.group != "excluded") & meta.tissue.str.lower().eq("peripheral blood")
    ]
    case = meta.loc[meta.group.eq("case"), "gsm"].tolist()
    control = meta.loc[meta.group.eq("healthy_control"), "gsm"].tolist()
    assert len(case) == 59 and len(control) == 62
    with gzip.open(RAW / "GSE164191_series_matrix.txt.gz", "rt", errors="replace") as f:
        lines = []
        in_table = False
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                lines.append(line)
    x = pd.read_csv(StringIO("".join(lines)), sep="\t", quotechar='"').rename(
        columns={"ID_REF": "gene_id"}
    )
    x = x[["gene_id"] + case + control]
    x[case + control] = x[case + control].apply(pd.to_numeric, errors="coerce")
    x = x.dropna(subset=case + control).groupby("gene_id", sort=True).mean()
    __case_frame = x[case]
    __control_frame = x[control]
    return locals()
