from pathlib import Path
import csv
import gzip

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RAW = ROOT / "raw"
RESULT = ROOT / "result"


def bh(p):
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    vals = p[ok]
    order = np.argsort(vals)
    ranked = vals[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    tmp = np.empty_like(q)
    tmp[order] = np.minimum(q, 1.0)
    out[ok] = tmp
    return out


def read_matrix(path):
    metadata = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_"):
                parts = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
                metadata.setdefault(parts[0], []).append(parts[1:])
            if line.startswith("!series_matrix_table_begin"):
                d = pd.read_csv(fh, sep="\t", quotechar='"', comment="!", low_memory=False)
                break
        else:
            raise RuntimeError("series matrix table not found")
    d = d.rename(columns={d.columns[0]: "probe_id"})
    d["probe_id"] = d["probe_id"].astype(str)
    return d, metadata


def main():
    matrix, meta = read_matrix(RAW / "GSE42148_series_matrix.txt.gz")
    accessions = meta["!Sample_geo_accession"][0]
    source = meta["!Sample_source_name_ch1"][0]
    characteristics = meta["!Sample_characteristics_ch1"]
    if not (len(accessions) == len(source) == len(matrix.columns) - 1):
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
        raise RuntimeError(f"unexpected groups: cases={len(cases)}, controls={len(controls)}")
    if len(set(accessions)) != len(accessions):
        raise RuntimeError("duplicate GEO sample accession")

    x = matrix.set_index("probe_id")
    case = x[cases].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    ctrl = x[controls].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    n1 = np.isfinite(case).sum(1)
    n0 = np.isfinite(ctrl).sum(1)
    e = np.nanmean(case, 1) - np.nanmean(ctrl, 1)
    u = np.sqrt(np.nanvar(case, 1, ddof=1) / n1 + np.nanvar(ctrl, 1, ddof=1) / n0)
    p = np.full(len(x), np.nan)
    for i in range(len(x)):
        a = case[i, np.isfinite(case[i])]
        b = ctrl[i, np.isfinite(ctrl[i])]
        if len(a) >= 2 and len(b) >= 2:
            p[i] = ttest_ind(a, b, equal_var=False).pvalue
    out = pd.DataFrame({
        "probe_id": x.index,
        "effect": e,
        "uncertainty": u,
        "p_value": p,
        "n_case": n1,
        "n_control": n0,
    }).reset_index(drop=True)
    out["adjusted_p_value"] = bh(out["p_value"].to_numpy())
    out = out[np.isfinite(out["effect"]) & np.isfinite(out["uncertainty"]) & (out["uncertainty"] > 0)].copy()
    out["contrast_id"] = "GSE42148_Coronary_Artery_Disease_vs_HC"
    out["analysis_method"] = "Welch mean difference on GEO processed Agilent whole-blood signal; uncertainty is Welch SE"
    out["source_path"] = "raw/GSE42148_series_matrix.txt.gz"
    out["feature_space"] = "GPL4133_probe"
    cols = ["probe_id", "effect", "uncertainty", "p_value", "adjusted_p_value", "n_case", "n_control",
            "contrast_id", "analysis_method", "source_path", "feature_space"]
    out[cols].to_csv(RESULT / "GSE42148_Coronary_Artery_Disease_vs_HC_E_U.tsv", sep="\t", index=False, float_format="%.10g")
    pd.DataFrame({
        "sample_accession": cases + controls,
        "group": ["case"] * len(cases) + ["control"] * len(controls),
        "tissue": ["whole blood"] * (len(cases) + len(controls)),
        "sex": ["male"] * (len(cases) + len(controls)),
    }).to_csv(ROOT / "sample_manifest.tsv", sep="\t", index=False)
    print(f"Wrote {len(out)} rows; cases={len(cases)}, controls={len(controls)}")


if __name__ == "__main__":
    main()
