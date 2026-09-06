from pathlib import Path
import hashlib
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "raw" / "series_matrix.txt.gz"
OUT = HERE.parent / "result" / "GSE201332_MDD_vs_HC_E_U.tsv"


def main():
    sm = pd.read_csv(HERE / "sample_manifest.tsv", sep="\t", dtype=str)
    case = sm.loc[sm["group"] == "case", "gsm"].tolist()
    control = sm.loc[sm["group"] == "control", "gsm"].tolist()
    if len(case) != 20 or len(control) != 20:
        raise RuntimeError(f"Unexpected sample counts: cases={len(case)}, controls={len(control)}")
    # GEO series-matrix metadata lines start with !; the first non-metadata
    # line is the quoted ID_REF/GSM table header.
    df = pd.read_csv(DATA, sep="\t", compression="gzip", comment="!", dtype=str)
    df.columns = [str(c).strip('"') for c in df.columns]
    missing = (set(case) | set(control)) - set(df.columns)
    if missing:
        raise RuntimeError(f"GSMs missing from matrix: {sorted(missing)}")
    X_case = df[case].apply(pd.to_numeric, errors="coerce")
    X_ctrl = df[control].apply(pd.to_numeric, errors="coerce")
    n1 = X_case.notna().sum(axis=1).astype(float)
    n0 = X_ctrl.notna().sum(axis=1).astype(float)
    E = X_case.mean(axis=1) - X_ctrl.mean(axis=1)
    U = np.sqrt(X_case.var(axis=1, ddof=1) / n1 + X_ctrl.var(axis=1, ddof=1) / n0)
    out = pd.DataFrame({"probe_id": df.iloc[:, 0].astype(str), "E": E, "U": U})
    out = out[np.isfinite(out["E"]) & np.isfinite(out["U"]) & (out["U"] > 0)]
    out = out.drop_duplicates("probe_id").sort_values("probe_id")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False, float_format="%.10g")
    print(f"Wrote {len(out)} rows; cases={len(case)}, controls={len(control)}")
    print(f"SHA256={hashlib.sha256(OUT.read_bytes()).hexdigest().upper()}")


if __name__ == "__main__":
    main()
