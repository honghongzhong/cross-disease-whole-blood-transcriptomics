from pathlib import Path
from io import StringIO
import gzip
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

HERE = Path(__file__).resolve().parent
RAW = HERE.parent / "raw"
OUT = HERE.parent / "result"


def read_matrix(path):
    lines = []
    active = False
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                active = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if active:
                lines.append(line)
    return pd.read_csv(StringIO("".join(lines)), sep="\t", quotechar='"')


def main():
    meta = pd.read_csv(HERE / "sample_manifest.tsv", sep="\t")
    assert len(meta) == 83 and meta.gsm.is_unique
    assert meta.source.str.lower().eq("whole blood").all()
    selected = meta[meta.primary_scope.eq("primary")].copy()
    case = selected.loc[selected.illness.str.startswith("CVID"), "gsm"].tolist()
    control = selected.loc[selected.illness.eq("Healthy"), "gsm"].tolist()
    assert len(case) == 59 and len(control) == 24
    x = read_matrix(RAW / "series_matrix.txt.gz").rename(columns={"ID_REF": "gene_id"}).set_index("gene_id")
    x.columns = [str(c).strip('"') for c in x.columns]
    x = x[case + control].apply(pd.to_numeric, errors="coerce").dropna().groupby(level=0, sort=True).mean()
    A, B = x[case].to_numpy(), x[control].to_numpy()
    e = A.mean(axis=1) - B.mean(axis=1)
    u = np.sqrt(A.var(axis=1, ddof=1) / len(case) + B.var(axis=1, ddof=1) / len(control))
    p = ttest_ind(A, B, axis=1, equal_var=False, nan_policy="omit").pvalue
    out = pd.DataFrame({"gene_id": x.index.astype(str), "E": e, "U": u, "p_value": p})
    out = out[np.isfinite(out.E) & np.isfinite(out.U) & (out.U > 0)].sort_values("gene_id")
    pv = out.p_value.to_numpy(); order = np.argsort(pv); q = np.empty_like(pv)
    q[order] = np.clip(np.minimum.accumulate((pv[order] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1], 0, 1)
    out["q_value"] = q
    OUT.mkdir(exist_ok=True)
    target = OUT / "GSE51405_CVID_vs_HC_E_U.tsv"
    out.to_csv(target, sep="\t", index=False, float_format="%.10g")
    print(f"Wrote {len(out)} rows to {target}; case={len(case)} control={len(control)}")


if __name__ == "__main__":
    main()
