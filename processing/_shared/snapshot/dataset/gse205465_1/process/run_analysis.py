from pathlib import Path
import gzip
import hashlib
import re

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RAW = HERE.parent / "raw"
COUNTS = RAW / "counts.tsv.gz"
SOFT = RAW / "family.soft.gz"
OUT = HERE.parent / "result" / "GSE205465_thrPAPS_vs_HC_E_U.tsv"


def parse_soft():
    text = gzip.open(SOFT, "rt", errors="replace").read()
    rows = []
    for block in text.split("^SAMPLE = ")[1:]:
        lines = block.splitlines()
        def vals(prefix):
            return [x.split(" = ", 1)[1].strip().strip('"') for x in lines if x.startswith(prefix)]
        chars = vals("!Sample_characteristics_ch1")
        joined = " | ".join(chars)
        diagnosis = next((x.split(":", 1)[1].strip() for x in chars if x.lower().startswith("diagnosis:")), "")
        rows.append({
            "gsm": vals("!Sample_geo_accession")[0],
            "title": vals("!Sample_title")[0],
            "source": vals("!Sample_source_name_ch1")[0],
            "characteristics": joined,
            "diagnosis": diagnosis,
            "biospecimen": "whole blood" if "whole blood" in (joined + " " + " ".join(vals("!Sample_source_name_ch1"))).lower() else "",
        })
    sm = pd.DataFrame(rows)
    sm["group"] = np.where(sm.diagnosis.str.lower().eq("thrombotic primary aps"), "case", "healthy_control")
    sm["count_column"] = sm.gsm
    return sm


def main():
    sm = parse_soft()
    counts = pd.read_csv(COUNTS, sep="\t", compression="gzip", engine="python")
    counts = counts.rename(columns={counts.columns[0]: "gene_id"})
    available = set(counts.columns[1:])
    sm["inclusion"] = np.where(sm.gsm.isin(available), "include", "exclude_missing_count")
    sm.to_csv(HERE / "sample_manifest.tsv", sep="\t", index=False)
    included = sm[sm.inclusion == "include"]
    case = included.loc[included.group == "case", "gsm"].tolist()
    control = included.loc[included.group == "healthy_control", "gsm"].tolist()
    if len(case) != 60 or len(control) != 28:
        raise RuntimeError(f"Unexpected available groups: cases={len(case)}, controls={len(control)}")
    if not included.biospecimen.eq("whole blood").all():
        raise RuntimeError("Non-whole-blood sample included")
    x = counts.set_index("gene_id")[case + control].apply(pd.to_numeric, errors="coerce")
    libsize = x.sum(axis=0)
    expr = np.log2(x.div(libsize, axis=1) * 1_000_000 + 1.0)
    A, B = expr[case].to_numpy(), expr[control].to_numpy()
    e = A.mean(axis=1) - B.mean(axis=1)
    u = np.sqrt(A.var(axis=1, ddof=1) / len(case) + B.var(axis=1, ddof=1) / len(control))
    out = pd.DataFrame({"gene_id": expr.index.astype(str), "E": e, "U": u})
    out = out[np.isfinite(out.E) & np.isfinite(out.U) & (out.U > 0)]
    out = out[~out.index.duplicated()].sort_values("gene_id")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False, float_format="%.10g")
    print(f"Wrote {len(out)} rows; cases={len(case)}, controls={len(control)}, missing={len(sm)-len(included)}")
    print(f"SHA256={hashlib.sha256(OUT.read_bytes()).hexdigest().upper()}")


if __name__ == "__main__":
    main()
