from pathlib import Path
import csv
import gzip
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

HERE = Path(__file__).resolve().parent
RAW = HERE.parent / "raw"
OUT = HERE.parent / "result"


def matrix_metadata(path):
    rows = {}
    with gzip.open(path, "rt", errors="replace", newline="") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            fields = next(csv.reader([line], delimiter="\t"))
            key = fields[0]
            values = [v.strip('"') for v in fields[1:]]
            rows[key] = values
    gsms = rows["!Sample_geo_accession"]
    titles = rows["!Sample_title"]
    source = rows["!Sample_source_name_ch1"]
    characteristics = rows["!Sample_characteristics_ch1"]
    out = []
    for i, gsm in enumerate(gsms):
        text = characteristics[i].lower()
        if "healthy" in text or "control" in text:
            group = "healthy_control"
        elif "ap" in text:
            group = "case"
        else:
            group = "excluded"
        out.append({"gsm": gsm, "title": titles[i], "source": source[i], "characteristics": characteristics[i], "group": group})
    return pd.DataFrame(out)


def main():
    meta = matrix_metadata(RAW / "GSE194331_series_matrix.txt.gz")
    assert len(meta) == 119, len(meta)
    assert meta.source.str.lower().eq("blood").all()
    # The matrix records the whole-blood specimen in a separate characteristic row;
    # the parsed disease row is retained below for the group audit.
    case = meta.loc[meta.group.eq("case"), "title"].tolist()
    control = meta.loc[meta.group.eq("healthy_control"), "title"].tolist()
    assert len(case) == 87, len(case)
    assert len(control) == 32, len(control)

    with gzip.open(RAW / "GSE194331_HC_PAN_PANSEP_counts.txt.gz", "rt", errors="replace") as fh:
        counts = pd.read_csv(fh, sep="\t")
    counts = counts.rename(columns={counts.columns[0]: "gene_id"}).set_index("gene_id")
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    assert set(case + control).issubset(counts.columns), "count/title mismatch"
    libs = counts.sum(axis=0)
    assert (libs > 0).all()
    logcpm = np.log2(counts.div(libs, axis=1) * 1e6 + 1.0)
    A = logcpm[case].to_numpy()
    B = logcpm[control].to_numpy()
    e = A.mean(axis=1) - B.mean(axis=1)
    u = np.sqrt(A.var(axis=1, ddof=1) / len(case) + B.var(axis=1, ddof=1) / len(control))
    p = ttest_ind(A, B, axis=1, equal_var=False, nan_policy="omit").pvalue
    out = pd.DataFrame({"gene_id": logcpm.index.astype(str), "E": e, "U": u, "p_value": p})
    out = out[np.isfinite(out.E) & np.isfinite(out.U) & (out.U > 0)].drop_duplicates("gene_id").sort_values("gene_id")
    pv = out.p_value.to_numpy()
    order = np.argsort(pv)
    q = np.empty_like(pv)
    q[order] = np.clip(np.minimum.accumulate((pv[order] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1], 0, 1)
    out["q_value"] = q
    OUT.mkdir(exist_ok=True)
    target = OUT / "GSE194331_Acute_Pancreatitis_vs_HC_E_U.tsv"
    out.to_csv(target, sep="\t", index=False, float_format="%.10g")
    meta.to_csv(HERE.parent / "sample_manifest.tsv", sep="\t", index=False)
    print(f"Wrote {len(out)} genes to {target}")
    print(f"case={len(case)} control={len(control)} input_genes={len(counts)}")


if __name__ == "__main__":
    main()
