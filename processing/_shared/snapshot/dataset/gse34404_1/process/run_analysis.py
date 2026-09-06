from pathlib import Path
from io import StringIO
import gzip
import math
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

HERE = Path(__file__).resolve().parent
RAW = HERE.parent / "raw"
OUT = HERE.parent / "result"


def bh_adjust(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(np.nan_to_num(p, nan=1.0))
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty_like(p)
    q[order] = np.clip(adj, 0, 1)
    return q


def main():
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
        signal_subjects = [x[:-len("_1_AVG_Signal")] for x in signal_columns]
        assert set(signal_subjects) == set(case_subjects + control_subjects)
        selected_columns = [0] + [header.index(s + "_1_AVG_Signal") for s in case_subjects + control_subjects]
        rows = []
        for line in f:
            values = line.rstrip("\r\n").split("\t")
            if len(values) != len(header):
                continue
            try:
                numeric = [float(values[i]) for i in selected_columns[1:]]
            except ValueError:
                continue
            if all(math.isfinite(v) for v in numeric):
                rows.append([values[0]] + numeric)
    matrix = pd.DataFrame(rows, columns=["probe_id"] + case_subjects + control_subjects)

    with open(RAW / "GPL10558_acc.cgi", encoding="latin1", errors="replace") as f:
        lines = f.read().splitlines()
    begin = lines.index("!platform_table_begin")
    annotation = pd.read_csv(
        StringIO("\n".join(lines[begin + 1:])), sep="\t", dtype=str, low_memory=False
    )[["ID", "Symbol"]].rename(columns={"ID": "probe_id", "Symbol": "gene_id"})
    annotation["gene_id"] = annotation["gene_id"].fillna("").str.strip()
    annotation = annotation[annotation["gene_id"].str.match(r"^[A-Za-z0-9_.-]+$")]
    matrix = matrix.merge(annotation, on="probe_id", how="inner")
    expression = matrix.drop(columns="probe_id").groupby("gene_id", sort=True).mean()

    A = expression[case_subjects].to_numpy(dtype=float)
    B = expression[control_subjects].to_numpy(dtype=float)
    e = A.mean(axis=1) - B.mean(axis=1)
    u = np.sqrt(A.var(axis=1, ddof=1) / len(case_subjects) + B.var(axis=1, ddof=1) / len(control_subjects))
    p = ttest_ind(A, B, axis=1, equal_var=False, nan_policy="omit").pvalue
    out = pd.DataFrame({"gene_id": expression.index, "E": e, "U": u, "p_value": p, "q_value": bh_adjust(p)})
    out = out[np.isfinite(out["E"]) & np.isfinite(out["U"]) & (out["U"] > 0)].sort_values("gene_id")
    OUT.mkdir(exist_ok=True)
    path = OUT / "GSE34404_malaria_vs_HC_E_U.tsv"
    out.to_csv(path, sep="\t", index=False, float_format="%.10g")
    print(f"Wrote {len(out)} genes to {path}")
    print(f"case_n={len(case_subjects)} control_n={len(control_subjects)} probe_rows={len(matrix)}")


if __name__ == "__main__":
    main()
