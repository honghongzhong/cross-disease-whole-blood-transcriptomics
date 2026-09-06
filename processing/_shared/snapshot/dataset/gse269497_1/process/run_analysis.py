from pathlib import Path
import gzip
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

HERE = Path(__file__).resolve().parent
RAW = HERE.parent / "raw"
OUT = HERE.parent / "result"


def bh_adjust(p):
    order = np.argsort(np.nan_to_num(p, nan=1.0))
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.empty_like(p)
    q[order] = np.clip(np.minimum.accumulate(ranked[::-1])[::-1], 0, 1)
    return q


def main():
    path = RAW / "GSE269497_genes_FPKM_expression.txt.gz"
    x = pd.read_csv(path, sep="\t", compression="gzip", dtype={"Gene ID": str})
    case = [f"AIT-{i} FPKM" for i in range(1, 6)]
    control = [f"Con-{i} FPKM" for i in range(1, 6)]
    genes = x[["Gene Symbol"] + case + control].copy()
    genes["Gene Symbol"] = genes["Gene Symbol"].astype(str).str.strip()
    genes = genes[genes["Gene Symbol"].str.match(r"^[A-Za-z0-9_.-]+$")]
    for c in case + control:
        genes[c] = pd.to_numeric(genes[c], errors="coerce")
    genes = genes.dropna(subset=case + control)
    genes = genes.groupby("Gene Symbol", sort=True)[case + control].mean()
    # FPKM is converted to log2(FPKM + 1) before effect/uncertainty calculation.
    genes = np.log2(genes + 1.0)
    A = genes[case].to_numpy()
    B = genes[control].to_numpy()
    e = A.mean(axis=1) - B.mean(axis=1)
    u = np.sqrt(A.var(axis=1, ddof=1) / len(case) + B.var(axis=1, ddof=1) / len(control))
    p = ttest_ind(A, B, axis=1, equal_var=False, nan_policy="omit").pvalue
    out = pd.DataFrame({"gene_id": genes.index, "E": e, "U": u, "p_value": p, "q_value": bh_adjust(p)})
    out = out[np.isfinite(out["E"]) & np.isfinite(out["U"]) & (out["U"] > 0)]
    out = out.sort_values("gene_id")
    OUT.mkdir(exist_ok=True)
    target = OUT / "GSE269497_autoimmune_thyroiditis_vs_HC_E_U.tsv"
    out.to_csv(target, sep="\t", index=False, float_format="%.10g")
    print(f"Wrote {len(out)} genes to {target}")


if __name__ == "__main__":
    main()
