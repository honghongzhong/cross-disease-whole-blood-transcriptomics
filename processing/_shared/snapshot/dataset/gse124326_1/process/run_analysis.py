from pathlib import Path
import gzip
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

HERE = Path(__file__).resolve().parent
RAW = HERE.parent / "raw"
OUT = HERE.parent / "result"


def main():
    meta = pd.read_csv(HERE / "sample_manifest.tsv", sep="\t", dtype=str)
    meta = meta[meta["selected_author_final_analysis"].eq("True")].copy()
    case = meta.loc[meta["group"].isin(["BP1", "BP2"]), "title"].tolist()
    control = meta.loc[meta["group"].eq("Control"), "title"].tolist()
    assert len(case) == 239 and len(control) == 205
    with gzip.open(RAW / "GSE124326_count_matrix.txt.gz", "rt", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
    sample_cols = {c.removesuffix(".counts"): c for c in header[1:]}
    ids = case + control
    assert set(ids).issubset(sample_cols)
    case_cols = [sample_cols[s] for s in case]
    control_cols = [sample_cols[s] for s in control]
    cols = case_cols + control_cols
    x = pd.read_csv(RAW / "GSE124326_count_matrix.txt.gz", sep="\t", compression="gzip", usecols=["gene"] + cols)
    x["gene"] = x["gene"].astype(str).str.replace(r"\..*$", "", regex=True)
    x[cols] = x[cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    x = x.groupby("gene", sort=True)[cols].sum()
    x = np.log2(x.div(x.sum(axis=0), axis=1) * 1e6 + 1.0)
    A, B = x[case_cols].to_numpy(), x[control_cols].to_numpy()
    e = A.mean(axis=1) - B.mean(axis=1)
    u = np.sqrt(A.var(axis=1, ddof=1) / len(case) + B.var(axis=1, ddof=1) / len(control))
    p = ttest_ind(A, B, axis=1, equal_var=False, nan_policy="omit").pvalue
    out = pd.DataFrame({"gene_id": x.index, "E": e, "U": u, "p_value": p})
    out = out[np.isfinite(out["E"]) & np.isfinite(out["U"]) & (out["U"] > 0)].sort_values("gene_id")
    pv = out["p_value"].to_numpy(); order = np.argsort(pv); q = np.empty_like(pv)
    q[order] = np.clip(np.minimum.accumulate((pv[order] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1], 0, 1)
    out["q_value"] = q
    OUT.mkdir(exist_ok=True)
    target = OUT / "GSE124326_bipolar_disorder_vs_control_E_U.tsv"
    out.to_csv(target, sep="\t", index=False, float_format="%.10g")
    print(f"Wrote {len(out)} genes to {target}")


if __name__ == "__main__":
    main()
