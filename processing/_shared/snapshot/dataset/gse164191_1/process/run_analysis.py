from pathlib import Path
import gzip
from io import StringIO
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

HERE = Path(__file__).resolve().parent
RAW = HERE.parent / "raw"
OUT = HERE.parent / "result"


def main():
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
        status = get("!Sample_characteristics_ch1 = disease status:").split(":", 1)[-1].strip()
        tissue = get("!Sample_characteristics_ch1 = tissue:").split(":", 1)[-1].strip()
        group = "case" if status == "colorectal cancer" else "healthy_control" if status == "normal" else "excluded"
        rows.append((gsm, title, group, tissue))
    meta = pd.DataFrame(rows, columns=["gsm", "title", "group", "tissue"])
    meta = meta[(meta.group != "excluded") & meta.tissue.str.lower().eq("peripheral blood")]
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
    x = pd.read_csv(StringIO("".join(lines)), sep="\t", quotechar='"').rename(columns={"ID_REF": "gene_id"})
    x = x[["gene_id"] + case + control]
    x[case + control] = x[case + control].apply(pd.to_numeric, errors="coerce")
    x = x.dropna(subset=case + control).groupby("gene_id", sort=True).mean()
    A, B = x[case].to_numpy(), x[control].to_numpy()
    e = A.mean(axis=1) - B.mean(axis=1)
    u = np.sqrt(A.var(axis=1, ddof=1) / len(case) + B.var(axis=1, ddof=1) / len(control))
    p = ttest_ind(A, B, axis=1, equal_var=False, nan_policy="omit").pvalue
    out = pd.DataFrame({"gene_id": x.index, "E": e, "U": u, "p_value": p})
    out = out[np.isfinite(out.E) & np.isfinite(out.U) & (out.U > 0)].sort_values("gene_id")
    pv = out.p_value.to_numpy(); order = np.argsort(pv); q = np.empty_like(pv)
    q[order] = np.clip(np.minimum.accumulate((pv[order] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1], 0, 1)
    out["q_value"] = q
    OUT.mkdir(exist_ok=True)
    target = OUT / "GSE164191_colorectal_cancer_vs_HC_E_U.tsv"
    out.to_csv(target, sep="\t", index=False, float_format="%.10g")
    meta.to_csv(HERE / "sample_manifest.tsv", sep="\t", index=False)
    print(f"Wrote {len(out)} genes to {target}")


if __name__ == "__main__":
    main()
