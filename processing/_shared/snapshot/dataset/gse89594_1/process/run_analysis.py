from pathlib import Path
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
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                d = pd.read_csv(fh, sep="\t", quotechar='"', comment="!", low_memory=False)
                break
        else:
            raise RuntimeError("series matrix table not found")
    d = d.rename(columns={d.columns[0]: "probe_id"})
    d["probe_id"] = d["probe_id"].astype(str)
    return d


def read_annotation(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                ann = pd.read_csv(fh, sep="\t", comment="!", low_memory=False)
                break
        else:
            raise RuntimeError("platform table not found")
    keep = [c for c in ["ID", "GENE_SYMBOL", "LOCUSLINK_ID", "GENE_NAME"] if c in ann.columns]
    ann = ann[keep].rename(columns={"ID": "probe_id", "GENE_SYMBOL": "canonical_symbol", "LOCUSLINK_ID": "entrez_id"})
    ann["probe_id"] = ann["probe_id"].astype(str)
    return ann.drop_duplicates("probe_id")


def main():
    matrix = read_matrix(RAW / "GSE89594_series_matrix.txt.gz")
    controls = list(matrix.columns[1:31])
    asd = list(matrix.columns[31:63])
    cases = list(matrix.columns[63:95])
    if len(controls) != 30 or len(asd) != 32 or len(cases) != 32:
        raise RuntimeError(f"unexpected groups: {len(controls)}, {len(asd)}, {len(cases)}")
    x = matrix.set_index("probe_id")
    case = x[cases].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    ctrl = x[controls].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    n1 = np.isfinite(case).sum(1)
    n0 = np.isfinite(ctrl).sum(1)
    m1 = np.nanmean(case, 1)
    m0 = np.nanmean(ctrl, 1)
    v1 = np.nanvar(case, 1, ddof=1)
    v0 = np.nanvar(ctrl, 1, ddof=1)
    u = np.sqrt(v1 / n1 + v0 / n0)
    e = m1 - m0
    p = np.full(len(x), np.nan)
    for i in range(len(x)):
        a = case[i, np.isfinite(case[i])]
        b = ctrl[i, np.isfinite(ctrl[i])]
        if len(a) >= 2 and len(b) >= 2:
            p[i] = ttest_ind(a, b, equal_var=False).pvalue
    ann_path = RAW / "GPL16699_platform.txt"
    ann = read_annotation(ann_path) if ann_path.exists() else None
    out = pd.DataFrame({"probe_id": x.index, "effect": e, "uncertainty": u, "p_value": p,
                        "n_case": n1, "n_control": n0})
    if ann is not None:
        out = out.reset_index(drop=True).merge(ann, on="probe_id", how="left")
        out["entrez_id"] = pd.to_numeric(out["entrez_id"], errors="coerce")
        out["canonical_symbol"] = out["canonical_symbol"].fillna("")
        out = out[out["canonical_symbol"].ne("") | out["entrez_id"].notna()].copy()
    out["adjusted_p_value"] = bh(out["p_value"].to_numpy())
    out = out[np.isfinite(out["effect"]) & np.isfinite(out["uncertainty"]) & (out["uncertainty"] > 0)].copy()
    out["contrast_id"] = "GSE89594_Williams_syndrome_vs_healthy_control"
    out["analysis_method"] = "Welch mean difference on GEO series-matrix log2 signal; uncertainty is Welch SE"
    out["source_path"] = "raw/GSE89594_series_matrix.txt.gz"
    out["feature_space"] = "GPL16699_probe"
    cols = ["probe_id", "canonical_symbol", "entrez_id", "effect", "uncertainty", "p_value", "adjusted_p_value", "n_case", "n_control", "contrast_id", "analysis_method", "source_path", "feature_space"]
    for c in cols:
        if c not in out:
            out[c] = ""
    out[cols].to_csv(RESULT / "GSE89594_Williams_vs_HC_E_U.tsv", sep="\t", index=False, float_format="%.10g")
    print(f"Wrote {len(out)} rows; cases={len(cases)}, controls={len(controls)}")


if __name__ == "__main__":
    main()
