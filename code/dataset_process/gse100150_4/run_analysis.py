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
    with gzip.open(path, "rt", encoding="latin1", errors="replace") as f:
        active = False
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                active = True; continue
            if line.startswith("!series_matrix_table_end"): break
            if active: lines.append(line)
    return pd.read_csv(StringIO("".join(lines)), sep="\t", quotechar='"')

def main():
    meta = pd.read_csv(HERE / "sample_manifest.tsv", sep="\t")
    meta = meta[meta.qc_pass.astype(str).str.upper().eq("TRUE")].copy()
    meta["group"] = meta["class"].map({"Case": "Case", "Healthy": "Healthy"})
    matrix = read_matrix(RAW / "GSE100150_series_matrix.txt.gz").rename(columns={"ID_REF": "probe_id"})
    gsm = meta.gsm.tolist()
    if any(x not in matrix.columns for x in gsm): raise ValueError("Selected GSM is absent from matrix")
    ann = pd.read_csv(RAW / "GPL6884.annot.gz", sep="\t", skiprows=28, dtype=str, encoding="latin1", low_memory=False)
    ann = ann.rename(columns={"ID": "probe_id", "Gene symbol": "gene_id"})[["probe_id", "gene_id"]].dropna().drop_duplicates("probe_id")
    x = matrix[["probe_id"] + gsm].merge(ann, on="probe_id", how="inner")
    for c in gsm: x[c] = pd.to_numeric(x[c], errors="coerce")
    x = x[x.gene_id.str.match(r"^[A-Za-z0-9_.-]+$", na=False)].groupby("gene_id")[gsm].mean()
    case = meta.loc[meta.group.eq("Case"), "gsm"].tolist(); ctrl = meta.loc[meta.group.eq("Healthy"), "gsm"].tolist()
    a, b = x[case].to_numpy(), x[ctrl].to_numpy()
    effect = a.mean(1) - b.mean(1); uncertainty = np.sqrt(a.var(1, ddof=1)/len(case) + b.var(1, ddof=1)/len(ctrl))
    p = ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit").pvalue
    order = np.argsort(np.nan_to_num(p, nan=1.0)); ranked = p[order] * len(p) / np.arange(1, len(p)+1); q = np.empty_like(p); q[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    out = pd.DataFrame({"gene_id": x.index, "E": effect, "U": uncertainty, "p_value": p, "q_value": np.clip(q, 0, 1)})
    out = out[np.isfinite(out.E) & np.isfinite(out.U) & (out.U > 0)].sort_values("gene_id")
    condition = str(meta.loc[meta.group.eq("Case"), "dataset"].iloc[0])
    names = {"Bcell":"B_cell_deficiency", "FLU":"Flu", "JDM":"Juvenile_Dermatomyositis", "MS Patient":"MS", "Transplant":"Liver_Transplant", "HIV":"HIV"}
    condition = names.get(condition, condition)
    OUT.mkdir(exist_ok=True); path = OUT / f"GSE100150_{condition}_vs_HC_E_U.tsv"; out.to_csv(path, sep="\t", index=False, float_format="%.10g")
    print(f"Wrote {len(out)} genes to {path}")

if __name__ == "__main__": main()
