from pathlib import Path
import csv
import gzip
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

HERE = Path(__file__).resolve().parent
RAW = HERE.parent / "raw"
OUT = HERE.parent / "result"


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
        group = "case" if "ucb positive" in label else "healthy_control" if "health volunteer" in label else "excluded"
        out.append({
            "gsm": rows["!Sample_geo_accession"][i],
            "title": title,
            "source": source,
            "group": group,
        })
    return pd.DataFrame(out)


def main():
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
    A, B = data[case].to_numpy(), data[control].to_numpy()
    e = A.mean(axis=1) - B.mean(axis=1)
    u = np.sqrt(A.var(axis=1, ddof=1) / len(case) + B.var(axis=1, ddof=1) / len(control))
    p = ttest_ind(A, B, axis=1, equal_var=False, nan_policy="omit").pvalue
    out = pd.DataFrame({"probe_id": data.index.astype(str), "E": e, "U": u, "p_value": p})
    out = out[np.isfinite(out.E) & np.isfinite(out.U) & (out.U > 0)].drop_duplicates("probe_id").sort_values("probe_id")
    pv = out.p_value.to_numpy(); order = np.argsort(pv); q = np.empty_like(pv)
    q[order] = np.clip(np.minimum.accumulate((pv[order] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1], 0, 1)
    out["q_value"] = q
    OUT.mkdir(exist_ok=True)
    target = OUT / "GSE138118_Bladder_Urothelial_Carcinoma_vs_HC_E_U.tsv"
    out.to_csv(target, sep="\t", index=False, float_format="%.10g")
    meta.to_csv(HERE.parent / "sample_manifest.tsv", sep="\t", index=False)
    print(f"Wrote {len(out)} probes to {target}")
    print(f"case={len(case)} control={len(control)} input_probes={len(data)}")


if __name__ == "__main__":
    main()
