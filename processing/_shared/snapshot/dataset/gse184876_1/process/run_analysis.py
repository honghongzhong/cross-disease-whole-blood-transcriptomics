from pathlib import Path
import hashlib
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "raw" / "counts.txt.gz"
OUT = HERE.parent / "result" / "GSE184876_iRHOM2_deficiency_vs_HC_E_U.tsv"


def main():
    sm = pd.read_csv(HERE / "sample_manifest.tsv", sep="\t", dtype=str)
    case = sm.loc[sm.group == "case", "count_column"].tolist()
    control = sm.loc[sm.group == "control", "count_column"].tolist()
    if len(case) != 4 or len(control) != 6:
        raise RuntimeError(f"Unexpected sample counts: cases={len(case)}, controls={len(control)}")
    df = pd.read_csv(DATA, sep="\t", compression="gzip")
    df = df.rename(columns={df.columns[0]: "gene_id"})
    missing = (set(case) | set(control)) - set(df.columns)
    if missing:
        raise RuntimeError(f"Count columns missing: {sorted(missing)}")
    counts = df[case + control].apply(pd.to_numeric, errors="coerce")
    libsize = counts.sum(axis=0)
    expr = np.log2(counts.div(libsize, axis=1) * 1_000_000 + 1)
    X1, X0 = expr[case], expr[control]
    E = X1.mean(axis=1) - X0.mean(axis=1)
    U = np.sqrt(X1.var(axis=1, ddof=1) / len(case) + X0.var(axis=1, ddof=1) / len(control))
    out = pd.DataFrame({"gene_id": df["gene_id"].astype(str), "E": E, "U": U})
    out = out[np.isfinite(out.E) & np.isfinite(out.U) & (out.U > 0)]
    out = out.drop_duplicates("gene_id").sort_values("gene_id")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False, float_format="%.10g")
    print(f"Wrote {len(out)} rows; cases={len(case)}, controls={len(control)}")
    print(f"SHA256={hashlib.sha256(OUT.read_bytes()).hexdigest().upper()}")


if __name__ == "__main__":
    main()
