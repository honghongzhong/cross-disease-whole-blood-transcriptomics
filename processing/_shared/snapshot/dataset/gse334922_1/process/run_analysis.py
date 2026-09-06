from pathlib import Path
import gzip
import hashlib
import io

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RAW = HERE.parent / "raw" / "samples"
OUT = HERE.parent / "result" / "GSE334922_CAMR_vs_HC_E_U.tsv"


def read_sample(gsm: str):
    path = RAW / f"{gsm}.xlsx.gz"
    if not path.exists():
        raise FileNotFoundError(path)
    data = gzip.open(path, "rb").read()
    df = pd.read_excel(io.BytesIO(data), sheet_name=0)
    if df.shape[1] != 2 or str(df.columns[0]).lower() != "gene":
        raise RuntimeError(f"Unexpected expression layout for {gsm}: {df.shape}")
    gene = df.iloc[:, 0].astype(str)
    value = pd.to_numeric(df.iloc[:, 1], errors="coerce")
    if gene.duplicated().any():
        df = pd.DataFrame({"gene_id": gene, "value": value}).groupby("gene_id", sort=True)["value"].mean().reset_index()
    else:
        df = pd.DataFrame({"gene_id": gene, "value": value})
    return df.set_index("gene_id")["value"]


def main():
    sm = pd.read_csv(HERE / "sample_manifest.tsv", sep="\t", dtype=str)
    included = sm[sm.inclusion == "include"]
    case = included.loc[included.group == "case", "gsm"].tolist()
    control = included.loc[included.group == "healthy_control", "gsm"].tolist()
    if len(case) != 4 or len(control) != 4:
        raise RuntimeError(f"Unexpected primary sample counts: cases={len(case)}, controls={len(control)}")
    if set(case) & set(control):
        raise RuntimeError("Case/control overlap")

    series = {gsm: read_sample(gsm) for gsm in case + control}
    x = pd.DataFrame(series, dtype=float)
    x = np.log2(x + 1.0)
    A, B = x[case].to_numpy(), x[control].to_numpy()
    e = A.mean(axis=1) - B.mean(axis=1)
    u = np.sqrt(A.var(axis=1, ddof=1) / len(case) + B.var(axis=1, ddof=1) / len(control))
    out = pd.DataFrame({"gene_id": x.index.astype(str), "E": e, "U": u})
    out = out[np.isfinite(out.E) & np.isfinite(out.U) & (out.U > 0)].sort_values("gene_id")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False, float_format="%.10g")
    digest = hashlib.sha256(OUT.read_bytes()).hexdigest().upper()
    print(f"Wrote {len(out)} rows; cases={len(case)}, controls={len(control)}")
    print(f"SHA256={digest}")


if __name__ == "__main__":
    main()
